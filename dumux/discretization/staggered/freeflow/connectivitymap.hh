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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFreeFlowConnectivityMap
 */
#ifndef DUMUX_STAGGERED_FREEFLOW_CONNECTIVITY_MAP_HH
#define DUMUX_STAGGERED_FREEFLOW_CONNECTIVITY_MAP_HH

#include <vector>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Stores the dof indices corresponding to the neighboring cell centers and faces
 *        that contribute to the derivative calculation. Specialization for the staggered free flow model.
 */
template<class FVGridGeometry>
class StaggeredFreeFlowConnectivityMap
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;

    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = std::size_t;

    using CellCenterIdxType = typename FVGridGeometry::DofTypeIndices::CellCenterIdx;
    using FaceIdxType = typename FVGridGeometry::DofTypeIndices::FaceIdx;

    using CellCenterToCellCenterMap = std::vector<std::vector<IndexType>>;
    using CellCenterToFaceMap = std::vector<std::vector<IndexType>>;
    using FaceToCellCenterMap = std::vector<std::vector<IndexType>>;
    using FaceToFaceMap = std::vector<std::vector<IndexType>>;

    using Stencil = std::vector<IndexType>;

public:

    //! Update the map and prepare the stencils
    void update(const FVGridGeometry& fvGridGeometry)
    {
        const auto numDofsCC = fvGridGeometry.gridView().size(0);
        const auto numDofsFace = fvGridGeometry.gridView().size(1);
        const auto numBoundaryFacets = fvGridGeometry.numBoundaryScvf();
        cellCenterToCellCenterMap_.resize(numDofsCC);
        cellCenterToFaceMap_.resize(numDofsCC);
        faceToCellCenterMap_.resize(2*numDofsFace - numBoundaryFacets);
        faceToFaceMap_.resize(2*numDofsFace - numBoundaryFacets);

        std::vector<Stencil> fullFaceToCellCenterStencils;
        fullFaceToCellCenterStencils.resize(numDofsFace);
        std::vector<Stencil> fullfaceToFaceStencils;
        fullfaceToFaceStencils.resize(numDofsFace);

        for(auto&& element: elements(fvGridGeometry.gridView()))
        {
            // restrict the FvGeometry locally and bind to the element
            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bindElement(element);

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto dofIdxCellCenter = fvGridGeometry.elementMapper().index(element);
                computeCellCenterToCellCenterStencil_(cellCenterToCellCenterMap_[dofIdxCellCenter], element, fvGeometry, scvf);
                computeCellCenterToFaceStencil_(cellCenterToFaceMap_[dofIdxCellCenter], element, fvGeometry, scvf);

                const auto scvfIdx = scvf.index();
                computeFaceToCellCenterStencil_(faceToCellCenterMap_[scvfIdx], fvGeometry, scvf);
                computeFaceToFaceStencil_(faceToFaceMap_[scvfIdx], fvGeometry, scvf);
            }
        }
    }

    //! Returns the stencil of a cell center dof w.r.t. other cell center dofs
    const std::vector<IndexType>& operator() (CellCenterIdxType, CellCenterIdxType, const IndexType globalI) const
    {
        return cellCenterToCellCenterMap_[globalI];
    }

    //! Returns the stencil of a cell center dof w.r.t. face dofs
    const std::vector<IndexType>& operator() (CellCenterIdxType, FaceIdxType, const IndexType globalI) const
    {
        return cellCenterToFaceMap_[globalI];
    }

    //! Returns the stencil of a face dof w.r.t. cell center dofs
    const std::vector<IndexType>& operator() (FaceIdxType, CellCenterIdxType, const IndexType globalI) const
    {
        return faceToCellCenterMap_[globalI];
    }

    //! Returns the stencil of a face dof w.r.t. other face dofs
    const std::vector<IndexType>& operator() (FaceIdxType, FaceIdxType, const IndexType globalI) const
    {
        return faceToFaceMap_[globalI];
    }

private:

    /*
     * \brief Computes the stencil for the occupation pattern of the CCToCC block. This is empty.
     */
    void computeCellCenterToCellCenterStencil_(Stencil& stencil,
                                               const Element& element,
                                               const FVElementGeometry& fvGeometry,
                                               const SubControlVolumeFace& scvf)
    {}

    /*
     * \brief Computes the stencil for the occupation pattern of the CCToFace block.
     *        Basically, these are the dof indices of the element's faces.
     */
    void computeCellCenterToFaceStencil_(Stencil& stencil,
                                         const Element& element,
                                         const FVElementGeometry& fvGeometry,
                                         const SubControlVolumeFace& scvf)
    {
        stencil.push_back(scvf.dofIndex());
    }

    /*
     * \brief Computes the stencil for the occupation pattern of the FaceToCC block.
     *        This is the dof index of the element itself.
     */
    void computeFaceToCellCenterStencil_(Stencil& stencil,
                                         const FVElementGeometry& fvGeometry,
                                         const SubControlVolumeFace& scvf)
    {
        stencil.push_back(scvf.insideScvIdx());
    }

    /*
     * \brief Computes the stencil for the occupation pattern of the FaceToFace block.
     */
    void computeFaceToFaceStencil_(Stencil& stencil,
                                   const FVElementGeometry& fvGeometry,
                                   const SubControlVolumeFace& scvf)
    {
        // the first entries are always the face dofIdx itself and the one of the opposing face
        if(stencil.empty())
        {
            stencil.push_back(scvf.dofIndex());
            stencil.push_back(scvf.dofIndexOpposingFace());
        }

        for(const auto& data : scvf.pairData())
        {
            stencil.push_back(data.normalPair.first);
            const auto outerParallelFaceDofIdx = data.outerParallelFaceDofIdx;
            if(outerParallelFaceDofIdx >= 0)
                stencil.push_back(outerParallelFaceDofIdx);
            if(!scvf.boundary())
                stencil.push_back(data.normalPair.second);
        }
    }

    CellCenterToCellCenterMap cellCenterToCellCenterMap_;
    CellCenterToFaceMap cellCenterToFaceMap_;
    FaceToCellCenterMap faceToCellCenterMap_;
    FaceToFaceMap faceToFaceMap_;
};

} // end namespace Dumux

#endif // DUMUX_STAGGERED_FREEFLOW_CONNECTIVITY_MAP_HH
