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
 * \ingroup TracerModel
 * \file
 *
 * \brief Defines the properties required for the fully implicit tracer model.
 */

#ifndef DUMUX_TRACER_PROPERTIES_HH
#define DUMUX_TRACER_PROPERTIES_HH

#include <dumux/porousmediumflow/nonisothermal/implicit/properties.hh>

namespace Dumux
{
// \{
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the fully implicit tracer model.
NEW_TYPE_TAG(Tracer);
NEW_TYPE_TAG(BoxTracer, INHERITS_FROM(BoxModel, Tracer));
NEW_TYPE_TAG(CCTracer, INHERITS_FROM(CCTpfaModel, Tracer));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents);   //!< Number of fluid components in the system
NEW_PROP_TAG(Indices); //!< Enumerations for the model
NEW_PROP_TAG(SpatialParams); //!< The type of the spatial parameters
NEW_PROP_TAG(EffectiveDiffusivityModel); //!< The employed model for the computation of the effective diffusivity
NEW_PROP_TAG(FluidSystem); //!< Type of the multi-component relations
NEW_PROP_TAG(UseMoles); //!< Defines whether mole (true) or mass (false) fractions are used
NEW_PROP_TAG(SpatialParamsForchCoeff); //!< Property for the forchheimer coefficient
NEW_PROP_TAG(TauTortuosity); //!< Tortuosity value (tau) used in macroscopic diffusion
NEW_PROP_TAG(ProblemEnableGravity); //!< gravity flag for compatibility

}
// \}
}

#endif