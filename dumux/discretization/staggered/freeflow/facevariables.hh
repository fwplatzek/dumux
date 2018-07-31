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
 * \copydoc Dumux::StaggeredFaceVariables
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FREEFLOW_FACEVARIABLES_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FREEFLOW_FACEVARIABLES_HH

#include <array>

namespace Dumux
{

/*!
 * \ingroup StaggeredDiscretization
 * \brief The face variables class for free flow staggered grid models.
 *        Contains all relevant velocities for the assembly of the momentum balance.
 */
template<class FacePrimaryVariables, int dim>
class StaggeredFaceVariables
{
    static constexpr int numPairs = (dim == 2) ? 2 : 4;
    using Scalar = typename FacePrimaryVariables::block_type;

public:

    /*!
    * \brief Partial update of the face variables. Only the face itself is considered.
    *
    * \param priVars The face-specific primary variales
    */
    void updateOwnFaceOnly(const FacePrimaryVariables& priVars)
    {
        velocitySelf_ = priVars[0];
    }

    /*!
    * \brief Complete update of the face variables (i.e. velocities for free flow)
    *        for a given face
    *
    * \param faceSol The face-specific solution vector
    * \param problem The problem
    * \param element The element
    * \param fvGeometry The finite-volume geometry
    * \param scvf The sub-control volume face of interest
    */
    template<class SolVector, class Problem, class Element,
             class FVElementGeometry, class SubControlVolumeFace>
    void update(const SolVector& faceSol,
                const Problem& problem,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const SubControlVolumeFace& scvf)
    {
        velocitySelf_ = faceSol[scvf.dofIndex()];
        velocityOpposite_ = faceSol[scvf.dofIndexOpposingFace()];

        // treat the velocity forward of the self face
        if(scvf.hasFrontalNeighbor())
            velocityForward_ = faceSol[scvf.dofIndexForwardFace()];

        // treat the velocity previous to the opposite face
        if(scvf.hasPreviousNeighbor())
            velocityPrevious_ = faceSol[scvf.dofIndexPreviousFace()];

        // handle all sub faces
        for(int i = 0; i < scvf.pairData().size(); ++i)
        {
            const auto& subFaceData = scvf.pairData(i);

            // treat the velocities normal to the face
            velocityNormalInside_[i] = faceSol[subFaceData.normalPair.first];

            if(scvf.hasOuterNormal(i))
                velocityNormalOutside_[i] = faceSol[subFaceData.normalPair.second];

            // treat the velocity parallel to the face
            if(scvf.hasFirstParallelNeighbor(i))
                velocityFirstParallel_[i] = faceSol[subFaceData.firstParallelFaceDofIdx];

            // treat the velocity parallel to the parallel face
            if(scvf.hasSecondParallelNeighbor(i))
                velocitySecondParallel_[i] = faceSol[subFaceData.secondParallelFaceDofIdx];
        }
    }

    /*!
    * \brief Returns the velocity at the face itself
    */
    Scalar velocitySelf() const
    {
        return velocitySelf_;
    }

    /*!
    * \brief Returns the velocity at the opposing face
    */
    Scalar velocityOpposite() const
    {
        return velocityOpposite_;
    }

    /*!
    * \brief Returns the velocity at the previous face
    *
    * \param localSubFaceIdx The local index of the subface
    */
    Scalar velocityPrevious(const int localSubFaceIdx) const
    {
        return velocityPrevious_[localSubFaceIdx];
    }

    /*!
    * \brief Returns the velocity at the forward face
    *
    * \param localSubFaceIdx The local index of the subface
    */
    Scalar velocityForward(const int localSubFaceIdx) const
    {
        return velocityForward_[localSubFaceIdx];
    }

    /*!
    * \brief Returns the velocity at the parallel face
    *
    * \param localSubFaceIdx The local index of the subface
    */
    Scalar velocityFirstParallel(const int localSubFaceIdx) const
    {
        return velocityFirstParallel_[localSubFaceIdx];
    }

    /*!
    * \brief Returns the velocity at the next parallel face
    *
    * \param localSubFaceIdx The local index of the subface
    */
    Scalar velocitySecondParallel(const int localSubFaceIdx) const
    {
        return velocitySecondParallel_[localSubFaceIdx];
    }

    /*!
    * \brief Returns the velocity at the inner normal face
    *
    * \param localSubFaceIdx The local index of the subface
    */
    Scalar velocityNormalInside(const int localSubFaceIdx) const
    {
        return velocityNormalInside_[localSubFaceIdx];
    }

    /*!
    * \brief Returns the velocity at the outer normal face
    *
    * \param localSubFaceIdx The local index of the subface
    */
    Scalar velocityNormalOutside(const int localSubFaceIdx) const
    {
        return velocityNormalOutside_[localSubFaceIdx];
    }

private:

    Scalar velocityForward_;
    Scalar velocitySelf_;
    Scalar velocityOpposite_;
    Scalar velocityPrevious_;
    std::array<Scalar, numPairs> velocityFirstParallel_;
    std::array<Scalar, numPairs> velocitySecondParallel_;
    std::array<Scalar, numPairs> velocityNormalInside_;
    std::array<Scalar, numPairs> velocityNormalOutside_;
};

} // end namespace

#endif
