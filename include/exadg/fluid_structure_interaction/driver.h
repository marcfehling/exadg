/*  ______________________________________________________________________
 *
 *  ExaDG - High-Order Discontinuous Galerkin for the Exa-Scale
 *
 *  Copyright (C) 2021 by the ExaDG authors
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *  ______________________________________________________________________
 */

#ifndef INCLUDE_EXADG_FLUID_STRUCTURE_INTERACTION_DRIVER_H_
#define INCLUDE_EXADG_FLUID_STRUCTURE_INTERACTION_DRIVER_H_

// application
#include <exadg/fluid_structure_interaction/user_interface/application_base.h>

// utilities
#include <exadg/functions_and_boundary_conditions/interface_coupling.h>
#include <exadg/functions_and_boundary_conditions/verify_boundary_conditions.h>
#include <exadg/matrix_free/matrix_free_data.h>
#include <exadg/utilities/print_general_infos.h>
#include <exadg/utilities/timer_tree.h>

// grid
#include <exadg/grid/grid_motion_elasticity.h>
#include <exadg/grid/grid_motion_poisson.h>
#include <exadg/poisson/spatial_discretization/operator.h>

// IncNS
#include <exadg/incompressible_navier_stokes/postprocessor/postprocessor.h>
#include <exadg/incompressible_navier_stokes/spatial_discretization/operator_coupled.h>
#include <exadg/incompressible_navier_stokes/spatial_discretization/operator_dual_splitting.h>
#include <exadg/incompressible_navier_stokes/spatial_discretization/operator_pressure_correction.h>
#include <exadg/incompressible_navier_stokes/time_integration/time_int_bdf_coupled_solver.h>
#include <exadg/incompressible_navier_stokes/time_integration/time_int_bdf_dual_splitting.h>
#include <exadg/incompressible_navier_stokes/time_integration/time_int_bdf_pressure_correction.h>

// Structure
#include <exadg/structure/spatial_discretization/operator.h>
#include <exadg/structure/time_integration/time_int_gen_alpha.h>

namespace ExaDG
{
namespace FSI
{
/*
 * Own implementation of matrix class.
 */
template<typename Number>
class Matrix
{
public:
  // Constructor.
  Matrix(unsigned int const size) : M(size)
  {
    data.resize(M * M);

    init();
  }

  void
  init()
  {
    for(unsigned int i = 0; i < M; ++i)
      for(unsigned int j = 0; j < M; ++j)
        data[i * M + j] = Number(0.0);
  }

  Number
  get(unsigned int const i, unsigned int const j) const
  {
    AssertThrow(i < M && j < M, dealii::ExcMessage("Index exceeds matrix dimensions."));

    return data[i * M + j];
  }

  void
  set(Number const value, unsigned int const i, unsigned int const j)
  {
    AssertThrow(i < M && j < M, dealii::ExcMessage("Index exceeds matrix dimensions."));

    data[i * M + j] = value;
  }

private:
  // number of rows and columns of matrix
  unsigned int const  M;
  std::vector<Number> data;
};

template<typename VectorType, typename Number>
void
compute_QR_decomposition(std::vector<VectorType> & Q, Matrix<Number> & R, Number const eps = 1.e-2)
{
  for(unsigned int i = 0; i < Q.size(); ++i)
  {
    Number const norm_initial = Number(Q[i].l2_norm());

    // orthogonalize
    for(unsigned int j = 0; j < i; ++j)
    {
      Number r_ji = Q[j] * Q[i];
      R.set(r_ji, j, i);
      Q[i].add(-r_ji, Q[j]);
    }

    // normalize or drop if linear dependent
    Number r_ii = Number(Q[i].l2_norm());
    if(r_ii < eps * norm_initial)
    {
      Q[i] = 0.0;
      for(unsigned int j = 0; j < i; ++j)
        R.set(0.0, j, i);
      R.set(1.0, i, i);
    }
    else
    {
      R.set(r_ii, i, i);
      Q[i] *= 1. / r_ii;
    }
  }
}

/*
 *  Matrix has to be upper triangular with d_ii != 0 for all 0 <= i < n
 */
template<typename Number>
void
backward_substitution(Matrix<Number> const &      matrix,
                      std::vector<Number> &       dst,
                      std::vector<Number> const & rhs)
{
  int const n = dst.size();

  for(int i = n - 1; i >= 0; --i)
  {
    double value = rhs[i];
    for(int j = i + 1; j < n; ++j)
    {
      value -= matrix.get(i, j) * dst[j];
    }

    dst[i] = value / matrix.get(i, i);
  }
}

template<typename Number, typename VectorType>
void
backward_substitution_multiple_rhs(Matrix<Number> const &          matrix,
                                   std::vector<VectorType> &       dst,
                                   std::vector<VectorType> const & rhs)
{
  int const n = dst.size();

  for(int i = n - 1; i >= 0; --i)
  {
    VectorType value = rhs[i];
    for(int j = i + 1; j < n; ++j)
    {
      value.add(-matrix.get(i, j), dst[j]);
    }

    dst[i].equ(1.0 / matrix.get(i, i), value);
  }
}

template<typename VectorType>
void
inv_jacobian_times_residual(VectorType &                                                  b,
                            std::vector<std::shared_ptr<std::vector<VectorType>>> const & D_history,
                            std::vector<std::shared_ptr<std::vector<VectorType>>> const & R_history,
                            std::vector<std::shared_ptr<std::vector<VectorType>>> const & Z_history,
                            VectorType const &                                            residual)
{
  VectorType a = residual;

  // reset
  b = 0.0;

  for(int idx = Z_history.size() - 1; idx >= 0; --idx)
  {
    std::shared_ptr<std::vector<VectorType>> D = D_history[idx];
    std::shared_ptr<std::vector<VectorType>> R = R_history[idx];
    std::shared_ptr<std::vector<VectorType>> Z = Z_history[idx];

    int const           k = Z->size();
    std::vector<double> Z_times_a(k, 0.0);
    for(int i = 0; i < k; ++i)
      Z_times_a[i] = (*Z)[i] * a;

    // add to b
    for(int i = 0; i < k; ++i)
      b.add(Z_times_a[i], (*D)[i]);

    // add to a
    for(int i = 0; i < k; ++i)
      a.add(-Z_times_a[i], (*R)[i]);
  }
}

struct PartitionedData
{
  PartitionedData()
    : method("Aitken"),
      abs_tol(1.e-12),
      rel_tol(1.e-3),
      omega_init(0.1),
      reused_time_steps(0),
      partitioned_iter_max(100),
      geometric_tolerance(1.e-10)
  {
  }

  std::string  method;
  double       abs_tol;
  double       rel_tol;
  double       omega_init;
  unsigned int reused_time_steps;
  unsigned int partitioned_iter_max;

  // tolerance used to locate points at the fluid-structure interface
  double geometric_tolerance;
};

template<int dim, typename Number>
class WrapperStructure
{
public:
  void
  setup(std::shared_ptr<StructureFSI::ApplicationBase<dim, Number>> application,
        MPI_Comm const                                              mpi_comm,
        bool const                                                  is_test);

  // matrix-free
  std::shared_ptr<MatrixFreeData<dim, Number>>     matrix_free_data;
  std::shared_ptr<dealii::MatrixFree<dim, Number>> matrix_free;

  // spatial discretization
  std::shared_ptr<Structure::Operator<dim, Number>> pde_operator;

  // temporal discretization
  std::shared_ptr<Structure::TimeIntGenAlpha<dim, Number>> time_integrator;

  // postprocessor
  std::shared_ptr<Structure::PostProcessor<dim, Number>> postprocessor;
};

template<int dim, typename Number>
class WrapperFluid
{
public:
  WrapperFluid()
  {
    timer_tree = std::make_shared<TimerTree>();
  }

  void
  setup(std::shared_ptr<FluidFSI::ApplicationBase<dim, Number>> application,
        MPI_Comm const                                          mpi_comm,
        bool const                                              is_test);

  void
  solve_ale(std::shared_ptr<FluidFSI::ApplicationBase<dim, Number>> application,
            bool const                                              is_test) const;

  std::shared_ptr<TimerTree>
  get_timings_ale() const
  {
    return timer_tree;
  }

  // matrix-free
  std::shared_ptr<MatrixFreeData<dim, Number>>     matrix_free_data;
  std::shared_ptr<dealii::MatrixFree<dim, Number>> matrix_free;

  // spatial discretization
  std::shared_ptr<IncNS::SpatialOperatorBase<dim, Number>> pde_operator;

  // temporal discretization
  std::shared_ptr<IncNS::TimeIntBDF<dim, Number>> time_integrator;

  // Postprocessor
  std::shared_ptr<IncNS::PostProcessorBase<dim, Number>> postprocessor;

  // moving mapping (ALE)
  std::shared_ptr<GridMotionBase<dim, Number>> ale_grid_motion;

  // use a PDE solver for grid motion
  std::shared_ptr<MatrixFreeData<dim, Number>>     ale_matrix_free_data;
  std::shared_ptr<dealii::MatrixFree<dim, Number>> ale_matrix_free;

  // Poisson-type grid motion
  std::shared_ptr<Poisson::Operator<dim, Number, dim>> ale_poisson_operator;

  // elasticity-type grid motion
  std::shared_ptr<Structure::Operator<dim, Number>> ale_elasticity_operator;

  /*
   * Computation time (wall clock time).
   */
  std::shared_ptr<TimerTree> timer_tree;
};

template<int dim, typename Number>
class Driver
{
private:
  typedef dealii::LinearAlgebra::distributed::Vector<Number> VectorType;

public:
  Driver(std::string const &                           input_file,
         MPI_Comm const &                              comm,
         std::shared_ptr<ApplicationBase<dim, Number>> application,
         bool const                                    is_test);

  static void
  add_parameters(dealii::ParameterHandler & prm, PartitionedData & fsi_data);

  void
  setup();

  void
  solve() const;

  void
  print_performance_results(double const total_time) const;

private:
  void
  setup_interface_coupling();

  void
  set_start_time() const;

  void
  synchronize_time_step_size() const;

  unsigned int
  solve_partitioned_problem() const;

  void
  coupling_structure_to_ale(VectorType const & displacement_structure) const;

  void
  coupling_structure_to_fluid(bool const extrapolate) const;

  void
  coupling_fluid_to_structure(bool const end_of_time_step) const;

  void
  apply_dirichlet_neumann_scheme(VectorType &       d_tilde,
                                 VectorType const & d,
                                 unsigned int       iteration) const;

  bool
  check_convergence(VectorType const & residual) const;

  void
  print_solver_info_header(unsigned int const i) const;

  void
  print_solver_info_converged(unsigned int const i) const;

  void
  print_partitioned_iterations() const;

  // MPI communicator
  MPI_Comm const mpi_comm;

  // output to std::cout
  dealii::ConditionalOStream pcout;

  // do not print wall times if is_test
  bool const is_test;

  // application
  std::shared_ptr<ApplicationBase<dim, Number>> application;

  std::shared_ptr<WrapperStructure<dim, Number>> structure;

  std::shared_ptr<WrapperFluid<dim, Number>> fluid;

  // interface coupling
  std::shared_ptr<InterfaceCoupling<dim, dim, Number>> structure_to_fluid;
  std::shared_ptr<InterfaceCoupling<dim, dim, Number>> structure_to_ale;
  std::shared_ptr<InterfaceCoupling<dim, dim, Number>> fluid_to_structure;

  /*
   * Fixed-point iteration.
   */
  PartitionedData fsi_data;

  // required for quasi-Newton methods
  mutable std::vector<std::shared_ptr<std::vector<VectorType>>> D_history, R_history, Z_history;

  /*
   * Computation time (wall clock time).
   */
  mutable TimerTree timer_tree;

  /*
   * The first number counts the number of time steps, the second number the total number
   * (accumulated over all time steps) of iterations of the partitioned FSI scheme.
   */
  mutable std::pair<unsigned int, unsigned long long> partitioned_iterations;
};

} // namespace FSI
} // namespace ExaDG


#endif /* INCLUDE_EXADG_FLUID_STRUCTURE_INTERACTION_DRIVER_H_ */
