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

#include "fluxvariables.hh"
#include "model.hh"
#include "problem.hh"
#include "indices.hh"
#include "volumevariables.hh"
#include "properties.hh"
#include "newtoncontroller.hh"
#include "localresidual.hh"

#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/spatialparams/implicit.hh>

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
SET_INT_PROP(Richards, NumPhases, 2);

//! Use pressure [Pa] by default
SET_BOOL_PROP(Richards, UseHead, false);

//! The local residual operator
SET_TYPE_PROP(Richards, LocalResidual, RichardsLocalResidual<TypeTag>);

//! The global model used
SET_TYPE_PROP(Richards, Model, RichardsModel<TypeTag>);

//! The class for the volume averaged quantities
SET_TYPE_PROP(Richards, VolumeVariables, RichardsVolumeVariables<TypeTag>);

//! The class for the quantities required for the flux calculation
SET_TYPE_PROP(Richards, FluxVariables, RichardsFluxVariables<TypeTag>);

//! The class of the newton controller
SET_TYPE_PROP(Richards, NewtonController, RichardsNewtonController<TypeTag>);

//! The upwind weight for the mass conservation equations
SET_SCALAR_PROP(Richards, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(Richards, ImplicitMobilityUpwindWeight, 1.0);

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
SET_PROP(Richards, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

/*!
 * \brief The wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. Please be aware that you
 * should be careful to use the Richards model in conjunction with
 * liquid non-wetting phases. This is only meaningful if the viscosity
 * of the liquid phase is _much_ lower than the viscosity of the
 * wetting phase.
 */
SET_PROP(Richards, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, NullComponent<Scalar> > type;
};

/*!
 * \brief The wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. This doed not need to be
 * specified by the problem for the Richards model to work because the
 * Richards model does not conserve the non-wetting phase.
 */
SET_PROP(Richards, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::GasPhase<Scalar, NullComponent<Scalar> > type;
};

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the immiscible twophase fluid system. The
 * actual fluids used are specified using in the problem definition by
 * the WettingPhase and NonwettingPhase properties. Be aware that
 * using different fluid systems in conjunction with the Richards
 * model only makes very limited sense.
 */
SET_PROP(Richards, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(Richards, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};

// disable velocity output by default

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
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

//! temperature is already written by the isothermal model
SET_BOOL_PROP(RichardsNI, NiOutputLevel, 0);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(RichardsNI, IsothermalModel, RichardsModel<TypeTag>);

// set isothermal FluxVariables
SET_TYPE_PROP(RichardsNI, IsothermalFluxVariables, RichardsFluxVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(RichardsNI, IsothermalVolumeVariables, RichardsVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(RichardsNI, IsothermalLocalResidual, RichardsLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(RichardsNI, IsothermalIndices, RichardsIndices<TypeTag>);

//set isothermal NumEq
SET_INT_PROP(RichardsNI, IsothermalNumEq, 1);

// \}
}

} // end namespace

#endif