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

#ifndef INCLUDE_EXADG_AMR_AMR_DATA_H_
#define INCLUDE_EXADG_AMR_AMR_DATA_H_

// ExaDG
#include <exadg/amr/enum_types.h>

namespace ExaDG
{
struct AMRData
{
  // Type of refinement to be imposed.
  RefinementType refinement_type = RefinementType::None;

  // Fraction of cells to be refined.
  double refine_fraction = 0.;

  // Fraction of cells to be coarsened.
  double coarsen_fraction = 0.;
};

} // namespace ExaDG

#endif /* INCLUDE_EXADG_AMR_AMR_DATA_H_ */