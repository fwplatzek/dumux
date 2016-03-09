// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup ThreePThreeCModel
 */
/*!
 * \file
 *
 * \brief Defines the properties required for the three-phase three-component
 *        fully implicit model.
 */
#ifndef DUMUX_3P3C_PROPERTIES_HH
#define DUMUX_3P3C_PROPERTIES_HH

#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/properties.hh>


namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit three-phase three-component problems
NEW_TYPE_TAG(ThreePThreeC);
NEW_TYPE_TAG(BoxThreePThreeC, INHERITS_FROM(BoxModel, ThreePThreeC));
NEW_TYPE_TAG(CCThreePThreeC, INHERITS_FROM(CCModel, ThreePThreeC));

//! The type tags for the corresponding non-isothermal problems
NEW_TYPE_TAG(ThreePThreeCNI, INHERITS_FROM(ThreePThreeC, NonIsothermal));
NEW_TYPE_TAG(BoxThreePThreeCNI, INHERITS_FROM(BoxModel, ThreePThreeCNI));
NEW_TYPE_TAG(CCThreePThreeCNI, INHERITS_FROM(CCModel, ThreePThreeCNI));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents); //!< Number of fluid components in the system
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(FluidSystem); //!< Type of the multi-component relations
NEW_PROP_TAG(FluidState); //!< Type of the fluid state to be used

NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The parameters of the material law (extracted from the spatial parameters)
NEW_PROP_TAG(EffectiveDiffusivityModel); //!< The employed model for the computation of the effective diffusivity

NEW_PROP_TAG(ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(ImplicitMassUpwindWeight); //!< The value of the upwind parameter for the mobility
NEW_PROP_TAG(ImplicitMobilityUpwindWeight); //!< Weight for the upwind mobility in the velocity calculation
NEW_PROP_TAG(UseConstraintSolver); //!< Determines whether a constraint solver should be used explicitly
NEW_PROP_TAG(SpatialParamsForchCoeff); //!< Property for the forchheimer coefficient
NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
}
}

#endif
