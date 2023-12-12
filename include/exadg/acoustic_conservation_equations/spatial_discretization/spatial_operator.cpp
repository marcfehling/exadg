/*  ______________________________________________________________________
 *
 *  ExaDG - High-Order Discontinuous Galerkin for the Exa-Scale
 *
 *  Copyright (C) 2023 by the ExaDG authors
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

// deal.II
#include <deal.II/numerics/vector_tools.h>

// ExaDG
#include <exadg/acoustic_conservation_equations/spatial_discretization/spatial_operator.h>
#include <exadg/grid/mapping_dof_vector.h>
#include <exadg/operators/finite_element.h>
#include <exadg/operators/quadrature.h>
#include <exadg/utilities/exceptions.h>

namespace ExaDG
{
namespace Acoustics
{
template<int dim, typename Number>
SpatialOperator<dim, Number>::SpatialOperator(
  std::shared_ptr<Grid<dim> const>               grid_in,
  std::shared_ptr<dealii::Mapping<dim> const>    mapping_in,
  std::shared_ptr<BoundaryDescriptor<dim> const> boundary_descriptor_in,
  std::shared_ptr<FieldFunctions<dim> const>     field_functions_in,
  Parameters const &                             parameters_in,
  std::string const &                            field_in,
  MPI_Comm const &                               mpi_comm_in)
  : Interface::SpatialOperator<Number>(),
    grid(grid_in),
    mapping(mapping_in),
    boundary_descriptor(boundary_descriptor_in),
    field_functions(field_functions_in),
    param(parameters_in),
    field(field_in),
    dof_handler_p(*grid_in->triangulation),
    dof_handler_u(*grid_in->triangulation),
    mpi_comm(mpi_comm_in),
    pcout(std::cout, dealii::Utilities::MPI::this_mpi_process(mpi_comm_in) == 0)
{
  pcout << std::endl
        << "Construct acoustic conservation equations operator ..." << std::endl
        << std::flush;

  initialize_dof_handler_and_constraints();

  pcout << std::endl << "... done!" << std::endl << std::flush;
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::fill_matrix_free_data(
  MatrixFreeData<dim, Number> & matrix_free_data) const
{
  // append mapping flags
  matrix_free_data.append_mapping_flags(Operators::Kernel<dim, Number>::get_mapping_flags());

  // dof handler
  matrix_free_data.insert_dof_handler(&dof_handler_p, field + dof_index_p);
  matrix_free_data.insert_dof_handler(&dof_handler_u, field + dof_index_u);

  // constraint
  matrix_free_data.insert_constraint(&constraint_p, field + dof_index_p);
  matrix_free_data.insert_constraint(&constraint_u, field + dof_index_u);

  // quadrature for pressure
  std::shared_ptr<dealii::Quadrature<dim>> quadrature_p =
    create_quadrature<dim>(param.grid.element_type, param.degree_p + 1);
  matrix_free_data.insert_quadrature(*quadrature_p, field + quad_index_p);

  // quadrature for velocity
  std::shared_ptr<dealii::Quadrature<dim>> quadrature_u =
    create_quadrature<dim>(param.grid.element_type, param.degree_u + 1);
  matrix_free_data.insert_quadrature(*quadrature_u, field + quad_index_u);

  // quadrature that works for pressure and velocity
  std::shared_ptr<dealii::Quadrature<dim>> quadrature_p_u =
    create_quadrature<dim>(param.grid.element_type, std::max(param.degree_p, param.degree_u) + 1);
  matrix_free_data.insert_quadrature(*quadrature_p_u, field + quad_index_p_u);
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::setup()
{
  // initialize MatrixFree and MatrixFreeData
  std::shared_ptr<dealii::MatrixFree<dim, Number>> mf =
    std::make_shared<dealii::MatrixFree<dim, Number>>();
  std::shared_ptr<MatrixFreeData<dim, Number>> mf_data =
    std::make_shared<MatrixFreeData<dim, Number>>();

  fill_matrix_free_data(*mf_data);

  mf->reinit(*get_mapping(),
             mf_data->get_dof_handler_vector(),
             mf_data->get_constraint_vector(),
             mf_data->get_quadrature_vector(),
             mf_data->data);

  // Subsequently, call the other setup function with MatrixFree/MatrixFreeData objects as
  // arguments.
  this->setup(mf, mf_data);
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::setup(
  std::shared_ptr<dealii::MatrixFree<dim, Number> const> matrix_free_in,
  std::shared_ptr<MatrixFreeData<dim, Number> const>     matrix_free_data_in)
{
  pcout << std::endl
        << "Setup acoustic conservation equations operator ..." << std::endl
        << std::flush;

  // MatrixFree
  matrix_free      = matrix_free_in;
  matrix_free_data = matrix_free_data_in;

  initialize_operators();

  pcout << std::endl << "... done!" << std::endl << std::flush;
}

template<int dim, typename Number>
dealii::MatrixFree<dim, Number> const &
SpatialOperator<dim, Number>::get_matrix_free() const
{
  return *matrix_free;
}

template<int dim, typename Number>
std::string
SpatialOperator<dim, Number>::get_dof_name_pressure() const
{
  return field + dof_index_p;
}

template<int dim, typename Number>
unsigned int
SpatialOperator<dim, Number>::get_dof_index_pressure() const
{
  return matrix_free_data->get_dof_index(get_dof_name_pressure());
}

template<int dim, typename Number>
std::string
SpatialOperator<dim, Number>::get_dof_name_velocity() const
{
  return field + dof_index_u;
}

template<int dim, typename Number>
unsigned int
SpatialOperator<dim, Number>::get_dof_index_velocity() const
{
  return matrix_free_data->get_dof_index(get_dof_name_velocity());
}

template<int dim, typename Number>
unsigned int
SpatialOperator<dim, Number>::get_quad_index_pressure_velocity() const
{
  return matrix_free_data->get_quad_index(field + quad_index_p_u);
}

template<int dim, typename Number>
unsigned int
SpatialOperator<dim, Number>::get_quad_index_pressure() const
{
  return matrix_free_data->get_quad_index(field + quad_index_p);
}

template<int dim, typename Number>
unsigned int
SpatialOperator<dim, Number>::get_quad_index_velocity() const
{
  return matrix_free_data->get_quad_index(field + quad_index_u);
}

template<int dim, typename Number>
std::shared_ptr<dealii::Mapping<dim> const>
SpatialOperator<dim, Number>::get_mapping() const
{
  return mapping;
}

template<int dim, typename Number>
dealii::FiniteElement<dim> const &
SpatialOperator<dim, Number>::get_fe_p() const
{
  return *fe_p;
}

template<int dim, typename Number>
dealii::FiniteElement<dim> const &
SpatialOperator<dim, Number>::get_fe_u() const
{
  return *fe_u;
}

template<int dim, typename Number>
dealii::DoFHandler<dim> const &
SpatialOperator<dim, Number>::get_dof_handler_p() const
{
  return dof_handler_p;
}

template<int dim, typename Number>
dealii::DoFHandler<dim> const &
SpatialOperator<dim, Number>::get_dof_handler_u() const
{
  return dof_handler_u;
}

template<int dim, typename Number>
dealii::AffineConstraints<Number> const &
SpatialOperator<dim, Number>::get_constraint_p() const
{
  return constraint_p;
}

template<int dim, typename Number>
dealii::AffineConstraints<Number> const &
SpatialOperator<dim, Number>::get_constraint_u() const
{
  return constraint_u;
}

template<int dim, typename Number>
dealii::types::global_dof_index
SpatialOperator<dim, Number>::get_number_of_dofs() const
{
  return dof_handler_u.n_dofs() + dof_handler_p.n_dofs();
}

/*
 * Initialization of vectors.
 */
template<int dim, typename Number>
void
SpatialOperator<dim, Number>::initialize_dof_vector(BlockVectorType & dst) const
{
  dst.reinit(2);

  matrix_free->initialize_dof_vector(dst.block(block_index_pressure), get_dof_index_pressure());
  matrix_free->initialize_dof_vector(dst.block(block_index_velocity), get_dof_index_velocity());

  dst.collect_sizes();
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::prescribe_initial_conditions(BlockVectorType & dst,
                                                           double const      time) const
{
  field_functions->initial_solution_pressure->set_time(time);
  field_functions->initial_solution_velocity->set_time(time);

  // This is necessary if Number == float
  using VectorTypeDouble = dealii::LinearAlgebra::distributed::Vector<double>;

  VectorTypeDouble pressure_double;
  VectorTypeDouble velocity_double;
  pressure_double = dst.block(block_index_pressure);
  velocity_double = dst.block(block_index_velocity);

  dealii::VectorTools::interpolate(*get_mapping(),
                                   dof_handler_p,
                                   *(field_functions->initial_solution_pressure),
                                   pressure_double);

  dealii::VectorTools::interpolate(*get_mapping(),
                                   dof_handler_u,
                                   *(field_functions->initial_solution_velocity),
                                   velocity_double);

  dst.block(block_index_pressure) = pressure_double;
  dst.block(block_index_velocity) = velocity_double;
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::evaluate(BlockVectorType &       dst,
                                       BlockVectorType const & src,
                                       double const            time) const
{
  evaluate_acoustic_operator(dst, src, time);

  // shift to the right-hand side of the equation
  dst *= -1.0;

  apply_inverse_mass_operator(dst, dst);
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::evaluate_acoustic_operator(BlockVectorType &       dst,
                                                         BlockVectorType const & src,
                                                         double const            time) const
{
  acoustic_operator.evaluate(dst, src, time);
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::apply_inverse_mass_operator(BlockVectorType &       dst,
                                                          BlockVectorType const & src) const
{
  inverse_mass_pressure.apply(dst.block(block_index_pressure), src.block(block_index_pressure));
  inverse_mass_velocity.apply(dst.block(block_index_velocity), src.block(block_index_velocity));
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::initialize_dof_handler_and_constraints()
{
  fe_p = create_finite_element<dim>(param.grid.element_type, true, 1, param.degree_p);
  fe_u = create_finite_element<dim>(param.grid.element_type, true, dim, param.degree_u);

  // enumerate degrees of freedom
  dof_handler_p.distribute_dofs(*fe_p);
  dof_handler_u.distribute_dofs(*fe_u);

  // close constraints
  constraint_u.close();
  constraint_p.close();

  // Output DoF information
  pcout << "Pressure:" << std::endl;
  print_parameter(pcout, "degree of 1D polynomials", param.degree_p);
  print_parameter(pcout, "number of dofs per cell", fe_p->n_dofs_per_cell());
  print_parameter(pcout, "number of dofs (total)", dof_handler_p.n_dofs());

  pcout << "Velocity:" << std::endl;
  print_parameter(pcout, "degree of 1D polynomials", param.degree_u);
  print_parameter(pcout, "number of dofs per cell", fe_u->n_dofs_per_cell());
  print_parameter(pcout, "number of dofs (total)", dof_handler_u.n_dofs());

  pcout << "Pressure and velocity:" << std::endl;
  print_parameter(pcout,
                  "number of dofs per cell",
                  fe_p->n_dofs_per_cell() + fe_u->n_dofs_per_cell());
  print_parameter(pcout, "number of dofs (total)", get_number_of_dofs());

  pcout << std::flush;
}

template<int dim, typename Number>
void
SpatialOperator<dim, Number>::initialize_operators()
{
  // inverse mass operator pressure
  {
    InverseMassOperatorData data;
    data.dof_index  = get_dof_index_pressure();
    data.quad_index = get_quad_index_pressure();
    inverse_mass_pressure.initialize(*matrix_free, data);
  }

  // inverse mass operator velocity
  {
    InverseMassOperatorData data;
    data.dof_index  = get_dof_index_velocity();
    data.quad_index = get_quad_index_velocity();
    inverse_mass_velocity.initialize(*matrix_free, data);
  }

  // acoustic operator
  {
    OperatorData<dim> data;
    data.dof_index_pressure   = get_dof_index_pressure();
    data.dof_index_velocity   = get_dof_index_velocity();
    data.quad_index           = get_quad_index_pressure_velocity();
    data.block_index_pressure = block_index_pressure;
    data.block_index_velocity = block_index_velocity;
    data.speed_of_sound       = param.speed_of_sound;
    data.density              = param.density;
    data.formulation          = param.formulation;
    data.bc                   = boundary_descriptor;
    acoustic_operator.initialize(*matrix_free, data);
  }
}

template class SpatialOperator<2, float>;
template class SpatialOperator<3, float>;

template class SpatialOperator<2, double>;
template class SpatialOperator<3, double>;

} // namespace Acoustics
} // namespace ExaDG