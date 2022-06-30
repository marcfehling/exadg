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

// deal.II
#include <deal.II/base/exceptions.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/lac/la_parallel_vector.h>

// ExaDG
#include <exadg/amr/grid_refinement.h>

namespace ExaDG
{
template<int dim, typename VectorType, int spacedim>
GridRefinement<dim, VectorType, spacedim>::GridRefinement(
                AMRData const &amr_data,
                dealii::Mapping<dim, spacedim> const &mapping,
                dealii::DoFHandler<dim, spacedim> const &dof_handler,
                unsigned int const face_quadrature_degree,
                VectorType const &locally_relevant_solution,
                dealii::parallel::distributed::Triangulation<dim, spacedim> & triangulation)
 : amr_data(amr_data)
 , mapping(mapping)
 , dof_handler(dof_handler)
 , face_quadrature(face_quadrature_degree)
 , locally_relevant_solution(locally_relevant_solution)
 , triangulation(triangulation)
{
  // TODO: check dof_handler.get_triangulation() == triangulation

  // AssertThrow(
  //       (dynamic_cast<
  //          typename dealii::parallel::distributed::Triangulation<dim, spacedim> const *>(
  //          &dof_handler.get_triangulation()) != nullptr),
  //       dealii::ExcMessage("Only implemented for distributed Triangulations so far."));
}


template<int dim, typename VectorType, int spacedim>
void
GridRefinement<dim, VectorType, spacedim>::estimate_mark_refine()
{
  estimate_error();
  mark_cells();
  refine_triangulation();
}


template<int dim, typename VectorType, int spacedim>
void
GridRefinement<dim, VectorType, spacedim>::estimate_error()
{
  criteria.reinit(triangulation.n_active_cells());

  // TODO: add neumann bc
  dealii::KellyErrorEstimator<dim, spacedim>::estimate(
    mapping,
    dof_handler,
    face_quadrature,
    /*neumann_bc=*/std::map<dealii::types::boundary_id, dealii::Function<spacedim, typename VectorType::value_type> const *>(),
    locally_relevant_solution,
    criteria);
}


template<int dim, typename VectorType, int spacedim>
void
GridRefinement<dim, VectorType, spacedim>::mark_cells() const
{
  switch (amr_data.refinement_type)
  {
    case RefinementType::None:
      AssertThrow(false, dealii::ExcMessage("Internal error."));
      break;

    case RefinementType::FixedFraction:
      dealii::parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
        triangulation,
        criteria,
        amr_data.refine_fraction,
        amr_data.coarsen_fraction);
      break;

    case RefinementType::FixedNumber:
      dealii::parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
        triangulation,
        criteria,
        amr_data.refine_fraction,
        amr_data.coarsen_fraction);
      break;

    default:
      AssertThrow(false, dealii::ExcMessage("Not implemented."));
      break;
  }
}


template<int dim, typename VectorType, int spacedim>
void
GridRefinement<dim, VectorType, spacedim>::refine_triangulation() const
{
  triangulation.execute_coarsening_and_refinement();
}


// explicit instantiations
template class GridRefinement<2, dealii::LinearAlgebra::distributed::Vector<float>, 2>;
template class GridRefinement<3, dealii::LinearAlgebra::distributed::Vector<float>, 3>;

template class GridRefinement<2, dealii::LinearAlgebra::distributed::Vector<double>, 2>;
template class GridRefinement<3, dealii::LinearAlgebra::distributed::Vector<double>, 3>;

} // namespace ExaDG