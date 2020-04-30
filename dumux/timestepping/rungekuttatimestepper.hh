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
/*!
 * \file
 * \brief A time stepper performing a single time step of a transient simulation
 */
#ifndef DUMUX_TIMESTEPPING_RUNGEKUTTA_TIMESTEPPER_HH
#define DUMUX_TIMESTEPPING_RUNGEKUTTA_TIMESTEPPER_HH

#include <iterator>
#include <algorithm>
#include <memory>
#include <dumux/common/timeloop.hh>
#include <dumux/timestepping/rungekuttamethods.hh>

namespace Dumux {

/*!
 * \file
 * \brief Time stepping with a Runge-Kutta method
 * \note We limit ourselves to implicit Runge-Kutta methods where solving
 *       a stage can only depend on the values of the same stage and stages before
 *       not bigger stages (which would require solving larger linear systems)
 */
template<class PDESolver>
class RungeKuttaTimeStepper
{
    using Scalar = typename PDESolver::Scalar;
    using TimeLoop = TimeLoopBase<Scalar>;
    using SolutionVector = typename PDESolver::Assembler::ResidualType;

    class StageSolutions
    {
    public:
        StageSolutions(SolutionVector& oldSol, SolutionVector& newSol, std::size_t numStages)
        : numStages_(numStages)
        {
            sols_.push_back(&oldSol);
            std::generate_n(std::back_inserter(sols_), numStages-1, []{ return new SolutionVector(); });
            sols_.push_back(&newSol);
        }

        SolutionVector& operator[] (std::size_t i)
        { return *sol_; }

        std::size_t numStages() const { return sols_.size()-1; }
    private:
        std::vector<SolutionVector*> sols_;
    };

public:
    RungeKuttaTimeStepper(std::shared_ptr<PDESolver> pdeSolver,
                          std::shared_ptr<const RungeKuttaMethod<Scalar>> rkMethod)
    : pdeSolver_(pdeSolver)
    , rkMethod_(rkMethod)
    {}

    /*!
     * \brief Advance one time step of the given time loop
     */
    void step(SolutionVector& sol, SolutionVector& oldSol, TimeLoop& timeLoop)
    {
        //pdeSolver_->assembler().preStep(*method);

        const auto numStages = method_->numStages();
        auto solutions = StageSolutions(oldSol, sol, numStages);
        for (auto stageIdx = 0U; stageIdx < numStages; ++stageIdx)
        {
            const auto stageTime = timeLoop.time() + rkMethod_->paramD[stageIdx]*timeLoop.timeStepSize();


        }

    }

    /*!
     * \brief Set/change the time step method
     */
    void setMethod(std::shared_ptr<const RungeKuttaMethod<Scalar>> rkMethod)
    { rkMethod_ = rkMethod; }

private:
    std::shared_ptr<PDESolver> pdeSolver_;
    std::shared_ptr<const RungeKuttaMethod<Scalar>> rkMethod_;
};

} // end namespace Dumux

#endif
