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
 * \file
 * \brief Class for the sub-control volume face in mpfa schemes
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_SUBCONTROLVOLUMEFACE_HH

#include <dune/common/version.hh>
#include "methods.hh"

namespace Dumux
{

/*!
 * \ingroup Mpfa
 * \brief Default implementation of the class for a sub-control volume face in mpfa methods.
 */
template<class ScvfGeometryTraits>
class CCMpfaDefaultSubControlVolumeFace
{
    using GridIndexType = typename ScvfGeometryTraits::GridIndexType;
    using Scalar = typename ScvfGeometryTraits::Scalar;
    using GlobalPosition = typename ScvfGeometryTraits::GlobalPosition;
    using CornerStorage = typename ScvfGeometryTraits::CornerStorage;
    using Geometry = typename ScvfGeometryTraits::Geometry;

public:
    //! state the traits public and thus export all types
    using Traits = ScvfGeometryTraits;

    /*!
     * \brief Constructor
     *
     * \param geomHelper The mpfa geometry helper
     * \param corners The corners of the scv face
     * \param unitOuterNormal The unit outer normal vector of the scvf
     * \param vIdxGlobal The global vertex index the scvf is connected to
     * \param vIdxLocal The element-local vertex index the scvf is connected to
     * \param scvfIndex The global index of this scv face
     * \param insideScvIdx The inside scv index connected to this face
     * \param outsideScvIndices The outside scv indices connected to this face
     * \param q The parameterization of the quadrature point on the scvf for flux calculation
     * \param boundary Boolean to specify whether or not the scvf is on a boundary
     */
    template<class MpfaHelper>
    explicit CCMpfaDefaultSubControlVolumeFace(const MpfaHelper& helper,
                                               CornerStorage&& corners,
                                               GlobalPosition&& unitOuterNormal,
                                               GridIndexType vIdxGlobal,
                                               unsigned int vIdxLocal,
                                               GridIndexType scvfIndex,
                                               GridIndexType insideScvIdx,
                                               const std::vector<GridIndexType>& outsideScvIndices,
                                               Scalar q,
                                               bool boundary)
             : boundary_(boundary),
               vertexIndex_(vIdxGlobal),
               scvfIndex_(scvfIndex),
               insideScvIdx_(insideScvIdx),
               outsideScvIndices_(outsideScvIndices),
               vIdxInElement_(vIdxLocal),
               corners_(std::move(corners)),
               center_(0.0),
               unitOuterNormal_(std::move(unitOuterNormal))
               {
                     // compute the center of the scvf
                     for (const auto& corner : corners_)
                         center_ += corner;
                     center_ /= corners_.size();

                     // use helper class to obtain area & integration point
                     ipGlobal_ = helper.getScvfIntegrationPoint(corners_, q);
                     area_ = helper.getScvfArea(corners_);
               }

    //! The area of the sub control volume face
    Scalar area() const { return area_; }

    //! returns bolean if the sub control volume face is on the domain boundary
    bool boundary() const { return boundary_; }

    //! The global index of this sub control volume face
    GridIndexType index() const { return scvfIndex_; }

    //! Returns the index of the vertex the scvf is connected to
    GridIndexType vertexIndex() const { return vertexIndex_; }

    //! Returns the element-local vertex index the scvf is connected to
    unsigned int vertexIndexInElement() const { return vIdxInElement_; }

    //! index of the inside sub control volume
    GridIndexType insideScvIdx() const { return insideScvIdx_; }

    //! The number of outside scvs connection via this scv face
    std::size_t numOutsideScvs() const { return outsideScvIndices_.size(); }

    //! index of the outside sub control volume or boundary scv index
    //! returns undefined behaviour if index exceeds numOutsideScvs
    GridIndexType outsideScvIdx(int i = 0) const { return outsideScvIndices_[i]; }

    //! returns the outside scv indices (can be more than one index for dim < dimWorld)
    const std::vector<GridIndexType>& outsideScvIndices() const { return outsideScvIndices_; }

    //! Returns the number of corners
    std::size_t corners() const { return corners_.size(); }

    //! Returns the corner for a given local index
    const GlobalPosition& corner(unsigned int localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

    //! Returns the global position of the vertex the scvf is connected to
    const GlobalPosition& vertexCorner() const { return corners_.back(); }

    //! Returns the global position of the center of the element facet this scvf is embedded in
    const GlobalPosition& facetCorner() const { return corner(0); }

    //! The center of the sub control volume face
    const GlobalPosition& center() const { return center_; }

    //! The integration point for flux evaluations in global coordinates
    const GlobalPosition& ipGlobal() const { return ipGlobal_; }

    //! returns the unit outer normal vector (assumes non-curved geometries)
    const GlobalPosition& unitOuterNormal() const { return unitOuterNormal_; }

    //! The geometry of the sub control volume face
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
    Geometry geometry() const { return Geometry(Dune::GeometryTypes::cube(Geometry::mydimension), corners_); }
#else
    Geometry geometry() const { return Geometry(Dune::GeometryType(Dune::GeometryType::cube, Geometry::mydimension), corners_); }
#endif

private:
    bool boundary_;
    GridIndexType vertexIndex_;
    GridIndexType scvfIndex_;
    GridIndexType insideScvIdx_;
    std::vector<GridIndexType> outsideScvIndices_;
    unsigned int vIdxInElement_;

    CornerStorage corners_;
    GlobalPosition center_;
    GlobalPosition ipGlobal_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;
};

/*!
 * \ingroup Mpfa
 * \brief Class for a sub control volume face in mpfa methods, i.e a part of the boundary
 *        of a control volume we compute fluxes on. Per default, we use the default
 *        implementation of the mpfa scvf class. If a scheme requires a different implementation,
 *        provide a specialization for it.
 *
 * \param M the mpfa method used
 * \param GT the traits class for the geometry type
 */
template<MpfaMethods M, class GT>
using CCMpfaSubControlVolumeFace = CCMpfaDefaultSubControlVolumeFace< GT >;

} // end namespace Dumux

#endif