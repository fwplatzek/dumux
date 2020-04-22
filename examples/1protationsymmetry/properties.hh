// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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

// ### Header guard
#ifndef DUMUX_ONEP_ROTATION_SYMMETRY_PROPERTIES_HH
#define DUMUX_ONEP_ROTATION_SYMMETRY_PROPERTIES_HH

// This file defines the `TypeTag` used for the single-phase rotation symmetry simulation, for
// which we then define the necessary properties.
//
// ### Include files
// The `TypeTag` defined for this simulation will inherit all properties from the
// `OneP` type tag, a convenience type tag that predefines most of the required
// properties for single-phase flow simulations in DuMuX. The properties that will be
// defined in this file are those that depend on user choices and no meaningful
// default can be set.
#include <dumux/porousmediumflow/1p/model.hh>

// We want to use `YaspGrid`, an implementation of the dune grid interface for structured grids:
#include <dune/grid/yaspgrid.hh>
// In this example, we want to discretize the equations with the cell centered finite volume
// scheme using two-point-flux approximation or with a vertex-centered scheme, the Box scheme:
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
// The fluid properties are specified in the following headers
// (we use the constant component for which its properties can be defined in the input file):
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

// The local residual for incompressible flow is included.
// The one-phase flow model (included above) uses a default implementation of the
// local residual for single-phase flow. However, in this example we are using an
// incompressible fluid phase. Therefore, we are including the specialized local
// residual which contains functionality to analytically compute the entries of
// the Jacobian matrix. We will use this in the main file.
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

// We include the problem and spatial parameters headers used for this simulation.
#include "problem.hh"
#include "spatialparams.hh"

// For rotational symmetric problems we use special geometry traits
#include <dumux/discretization/rotationsymmetricgridgeometrytraits.hh>

// ### Basic property definitions for the 1p problem
// We enter the namespace Dumux::Properties in order to import the entire Dumux namespace for general use:
namespace Dumux:: Properties {

// A `TypeTag` for our simulation is created which inherits from the one-phase flow model
// and the cell centered finite volume scheme with two-point-flux discretization scheme:
namespace TTag {
struct OnePRotSymBox { using InheritsFrom = std::tuple<OneP, BoxModel>; };
}

// We use a structured 1D grid with an offset:
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePRotSymBox> { using type =  Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 1> >; };

// Special grid geometry traits are needed
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::OnePRotSymBox>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using GGTraits = RotationSymmetricGridGeometryTraits<BoxDefaultGridGeometryTraits<GridView>, RotationPolicy::disc>;
public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, GGTraits>;
};

// The problem class specifies initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePRotSymBox> { using type = OnePTestProblem<TypeTag>; };

// We define the spatial parameters for our simulation:
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePRotSymBox>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = OnePTestSpatialParams<GridGeometry, Scalar>;
};

// We use the local residual that contains analytic derivative methods for incompressible flow:
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePRotSymBox> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// In the following we define the fluid system to be used:
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePRotSymBox>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// This enables grid-wide caching of the volume variables.
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePRotSymBox> { static constexpr bool value = false; };
//This enables grid wide caching for the flux variables.
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePRotSymBox> { static constexpr bool value = false; };
// This enables grid-wide caching for the finite volume grid geometry
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePRotSymBox> { static constexpr bool value = false; };
// The caches store values that were already calculated for later usage. This increases the memory demand but makes the simulation faster.

// We leave the namespace Dumux::Properties.
} // end namespace Dumux::Properties
#endif
