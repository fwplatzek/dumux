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
 *
 * \brief The primary variable switch base class for compositional models
 */
#ifndef DUMUX_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_PRIMARY_VARIABLE_SWITCH_HH

#include <dune/common/exceptions.hh>

namespace Dumux
{
/*!
 * \ingroup ImplicitModel
 * \brief The primary variable switch controlling the phase presence state variable
 */
template<class TypeTag>
class PrimaryVariableSwitch
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, PrimaryVariableSwitch);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

public:

    void init(const Problem& problem)
    {
        auto numDofs = problem.model().numDofs();
        phasePresence_.resize(numDofs);
        wasSwitched_.resize(numDofs, false);

        for (const auto& element : elements(problem.gridView()))
        {
            const auto& fvGeometry = problem.model().fvGeometries(element);
            for (auto&& scv : fvGeometry.scvs())
            {
                auto dofIdxGlobal = scv.dofIndex();
                phasePresence_[dofIdxGlobal] = problem.initialPhasePresence(scv);
            }
        }

        oldPhasePresence_ = phasePresence_;
    }

    int phasePresence(IndexType dofIdxGlobal) const
    {
        return phasePresence_[dofIdxGlobal];
    }

    void setPhasePresence(IndexType dofIdxGlobal, int phasePresence)
    {
        phasePresence_[dofIdxGlobal] = phasePresence;
    }

    void setOldPhasePresence(IndexType dofIdxGlobal, int phasePresence)
    {
        oldPhasePresence_[dofIdxGlobal] = phasePresence;
    }

    /*!
     * \brief Resets the current phase presence to the old one.
     * This is done after an update failed.
     */
    void resetPhasePresence()
    {
        phasePresence_ = oldPhasePresence_;
        wasSwitched_.assign(phasePresence_.size(), false);
    }

    /*!
     * \brief Sets the old phase presence state to the current one.
     * This is done when advancing to the next time step (problem post processing)
     */
    void updateOldPhasePresence()
    {
        oldPhasePresence_ = phasePresence_;
        wasSwitched_.assign(phasePresence_.size(), false);
    }

    /*!
     * \brief Update the variable switch / phase presence
     * \param problem The problem
     * \param curSol The current solution to be updated / modified
     */
    bool update(const Problem& problem, SolutionVector& curSol)
    {
        bool switched = false;
        visited_.assign(phasePresence_.size(), false);
        for (const auto& element : elements(problem.gridView()))
        {
            const auto& fvGeometry = problem.model().fvGeometries(element);
            for (auto&& scv : fvGeometry.scvs())
            {
                auto dofIdxGlobal = scv.dofIndex();
                if (!visited_[dofIdxGlobal])
                {
                    // Note this implies that volume variables don't differ
                    // in any sub control volume associated with the dof!
                    visited_[dofIdxGlobal] = true;
                    // Compute temporary volVars on which grounds we decide
                    // if we need to switch the primary variables
                    VolumeVariables volVars;
                    volVars.update(curSol[dofIdxGlobal], problem, element, scv);
                    const auto& globalPos = scv.dofPosition();

                    if (asImp_().update_(curSol[dofIdxGlobal], volVars, dofIdxGlobal, globalPos))
                        switched = true;

                }
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag for our partition.
        if (problem.gridView().comm().size() > 1)
            switched = problem.gridView().comm().max(switched);

        return switched;
    }

protected:

    //! return actual implementation (static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    //! return actual implementation (static polymorphism)
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // perform variable switch at a degree of freedom location
    bool update_(PrimaryVariables& priVars,
                 const VolumeVariables& volVars,
                 IndexType dofIdxGlobal,
                 const GlobalPosition& globalPos)
    {
        // evaluate if the primary variable switch would switch
        // to be implemented by the deriving class
        DUNE_THROW(Dune::NotImplemented, "This model seems to use a primary variable switch but none is implemented!");
    }

    std::vector<int> phasePresence_;
    std::vector<int> oldPhasePresence_;
    std::vector<bool> wasSwitched_;
    std::vector<bool> visited_;
};

} // end namespace dumux

#endif
