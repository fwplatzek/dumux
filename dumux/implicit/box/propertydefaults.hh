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
 * \ingroup BoxModel
 * \file
 * \brief Default properties for box models
 */
#ifndef DUMUX_BOX_PROPERTY_DEFAULTS_HH
#define DUMUX_BOX_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/propertydefaults.hh>
#include <dumux/implicit/fvelementgeometry.hh>
#include <dumux/porousmediumflow/implicit/fluxvariablescache.hh>

#include "fluxvariablescachevector.hh"
#include "volumevariablesvector.hh"
#include "elementboundarytypes.hh"
#include "fvelementgeometryvector.hh"
#include "localresidual.hh"
#include "localjacobian.hh"
#include "assembler.hh"
#include "properties.hh"
#include "stencils.hh"

namespace Dumux {

/*!
 * \brief The box model combines the advantages of the finite-volume (FV) and finite-element (FE) methods on a dual grid
 */
// forward declarations
template<class TypeTag> class BoxLocalResidual;
template<class TypeTag> class BoxElementBoundaryTypes;
template<class TypeTag> class BoxStencilsVector;
template<class TypeTag, bool enableFVElementGeometryCache> class BoxFVElementGeometryVector;

namespace Properties {
//! Set the corresponding discretization method property
SET_INT_PROP(BoxModel, DiscretizationMethod, GET_PROP(TypeTag, DiscretizationMethods)::Box);

//! Set the default for the FVElementGeometry vector
SET_TYPE_PROP(BoxModel, FVElementGeometryVector, BoxFVElementGeometryVector<TypeTag, GET_PROP_VALUE(TypeTag, EnableGlobalFVElementGeometryCache)>);

//! The sub control volume
SET_PROP(BoxModel, SubControlVolume)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using ScvGeometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld>;
    using IndexType = typename GridView::IndexSet::IndexType;
public:
    using type = SubControlVolume<ScvGeometry, IndexType, /*isBox=*/true>;
};

SET_PROP(BoxModel, SubControlVolumeFace)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using ScvfGeometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld>;
    using IndexType = typename GridView::IndexSet::IndexType;
public:
    using type = SubControlVolumeFace<ScvfGeometry, IndexType>;
};

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(BoxModel, ElementBoundaryTypes, BoxElementBoundaryTypes<TypeTag>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(BoxModel, DofMapper, typename GET_PROP_TYPE(TypeTag, VertexMapper));

//! The stencil container
SET_TYPE_PROP(BoxModel, StencilsVector, BoxStencilsVector<TypeTag>);

//! The global current volume variables vector class
SET_TYPE_PROP(BoxModel, CurrentVolumeVariablesVector, Dumux::BoxVolumeVariablesVector<TypeTag, false, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The global previous volume variables vector class
SET_TYPE_PROP(BoxModel, PreviousVolumeVariablesVector, Dumux::BoxVolumeVariablesVector<TypeTag, true, GET_PROP_VALUE(TypeTag, EnableGlobalVolumeVariablesCache)>);

//! The global flux variables cache vector class
SET_TYPE_PROP(BoxModel, FluxVariablesCacheVector, Dumux::BoxFluxVariablesCacheVector<TypeTag>);

//! Set the BaseLocalResidual to BoxLocalResidual
SET_TYPE_PROP(BoxModel, BaseLocalResidual, BoxLocalResidual<TypeTag>);

//! Assembler for the global jacobian matrix
SET_TYPE_PROP(BoxModel, JacobianAssembler, Dumux::BoxAssembler<TypeTag>);

//! The local jacobian operator
SET_TYPE_PROP(BoxModel, LocalJacobian, Dumux::BoxLocalJacobian<TypeTag>);

//! indicate that this is a box discretization
SET_BOOL_PROP(BoxModel, ImplicitIsBox, true);

} // namespace Properties
} // namespace Dumux

#endif
