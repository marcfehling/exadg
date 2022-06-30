/*  ______________________________________________________________________
 *
 *  ExaDG - High-Order Discontinuous Galerkin for the Exa-Scale
 *
 *  Copyright (C) 2022 by the ExaDG authors
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

#ifndef INCLUDE_EXADG_AMR_GRID_REFINEMENT_H_
#define INCLUDE_EXADG_AMR_GRID_REFINEMENT_H_

// deal.II
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/error_estimator.h>

// ExaDG
#include <exadg/amr/amr_data.h>

namespace ExaDG
{
template<int dim, typename VectorType, int spacedim = dim>
class GridRefinement
{
public:
  // Constructor.
  // TODO: I don't know where to get face quadrature rules from, so we take a quadrature
  //       degree instead.
  GridRefinement(AMRData const &amr_data,
                 dealii::Mapping<dim, spacedim> const & mapping,
                 dealii::DoFHandler<dim, spacedim> const &dof_handler,
                 unsigned int const face_quadrature_degree,
                 VectorType const & locally_relevant_solution,
                 dealii::parallel::distributed::Triangulation<dim, spacedim> &triangulation);

  // Perform error estimation, mark cells, and refinement.
  void estimate_mark_refine();

private:
  // Estimate errors using Kelly estimator.
  void estimate_error();

  // Mark cells for refinement.
  void mark_cells() const;

  // Perform refinement.
  void refine_triangulation() const;

  // TODO: ptr to amrdata? or rather copy?
  AMRData const & amr_data;

  // Reference to Mapping.
  dealii::Mapping<dim, spacedim> const & mapping;

  // Reference to DoFHandler.
  dealii::DoFHandler<dim, spacedim> const & dof_handler;

  // Quadrature rule for error estimator.
  dealii::QGauss<dim - 1> const face_quadrature;

  // Reference to Solution.
  VectorType const & locally_relevant_solution;

  // Triangulation that will be refined. Must be the same one as used in the DoFHandler.
  dealii::parallel::distributed::Triangulation<dim, spacedim> & triangulation;

  // Error estimates on locally owned part of domain.
  dealii::Vector<float> criteria;
};

} // namespace ExaDG

#endif /* INCLUDE_EXADG_AMR_GRID_REFINEMENT_H_ */