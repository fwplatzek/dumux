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
 * \ingroup RichardsModel
 * \file
 *
 * \brief Contains the default definitions for the properties required
 *        by the Richards fully implicit model.
 */
#ifndef DUMUX_RICHARDS_PROPERTY_DEFAULTS_HH
#define DUMUX_RICHARDS_PROPERTY_DEFAULTS_HH

#include "model.hh"
#include "indices.hh"
#include "volumevariables.hh"
#include "properties.hh"
#include "newtoncontroller.hh"

#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux
{
// \{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties values
//////////////////////////////////////////////////////////////////
//! Number of equations required by the model
SET_INT_PROP(Richards, NumEq, 1);

//! Number of fluid phases considered
//! Although the number of phases is two for the Richards model the
//! non-wetting phase is not balanced and thus we only have one
//! phase with a balance equation
SET_INT_PROP(Richards, NumPhases, 1);

//! Only the water component is balanced for Richards
SET_INT_PROP(Richards, NumComponents, 1);

//! The local residual operator
SET_TYPE_PROP(Richards, LocalResidual, ImmiscibleLocalResidual<TypeTag>);

//! The global model used
SET_TYPE_PROP(Richards, Model, RichardsModel<TypeTag>);

//! The class for the volume averaged quantities
SET_TYPE_PROP(Richards, VolumeVariables, RichardsVolumeVariables<TypeTag>);

//! Smarter newton controller
SET_TYPE_PROP(Richards, NewtonController, RichardsNewtonController<TypeTag>);

//! Enable advection
SET_BOOL_PROP(Richards, EnableAdvection, true);

//! The two-phase model has no molecular diffusion
SET_BOOL_PROP(Richards, EnableMolecularDiffusion, false);

//! Isothermal model by default
SET_BOOL_PROP(Richards, EnableEnergyBalance, false);

//! The class with all index definitions for the model
SET_TYPE_PROP(Richards, Indices, RichardsIndices<TypeTag>);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(Richards, SpatialParams, ImplicitSpatialParams<TypeTag>);

/*!
 * \brief Set type of the parameter objects for the material law
 *
 * By default this is just retrieved from the material law.
 */
SET_TYPE_PROP(Richards, MaterialLawParams, typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the liquid phase fluid system with simple H2O.
 */
SET_PROP(Richards, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar>>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(Richards, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

// enable gravity by default
SET_BOOL_PROP(Richards, ProblemEnableGravity, true);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(Richards, SpatialParamsForchCoeff, 0.55);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(RichardsNI, ThermalConductivityModel)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
public:
    using type = ThermalConductivitySomerton<Scalar, Indices>;
};

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(RichardsNI, IsothermalModel, RichardsModel<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(RichardsNI, IsothermalVolumeVariables, RichardsVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(RichardsNI, IsothermalLocalResidual, ImmiscibleLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(RichardsNI, IsothermalIndices, RichardsIndices<TypeTag>);

//set isothermal NumEq
SET_INT_PROP(RichardsNI, IsothermalNumEq, 1);

// \}
}

} // end namespace

#endif