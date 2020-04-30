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
 * \brief Parameters for different Runge-Kutta time stepping methods
 * \note See e.g. https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
 */
#ifndef DUMUX_TIMESTEPPING_RUNGEKUTTA_METHODS_HH
#define DUMUX_TIMESTEPPING_RUNGEKUTTA_METHODS_HH

#include <string>
#include <array>

namespace Dumux::TimeStepping {

/*!
 * \brief Abstract interface for Runge Kutta method parameters
 */
template<class Scalar>
class RungeKuttaMethod
{
public:
    virtual bool implicit () const = 0;

    virtual std::size_t numStages () const = 0;

    virtual Scalar paramAlpha (unsigned int i, unsigned int j) const = 0;
    virtual Scalar paramBeta (unsigned int i, unsigned int j) const = 0;
    virtual Scalar paramD (unsigned int i) const = 0;

    virtual std::string name () const = 0;

    virtual ~RungeKutta = default;
};

/*!
 * \brief A theta time stepping scheme
 * theta=1.0 is an implicit Euler scheme,
 * theta=0.0 an explicit Euler scheme,
 * theta=0.5 is a Cranck-Nicholson scheme
 */
template<class Scalar>
class Theta : public RungeKuttaMethod<Scalar>
{
public:
    explicit Theta(const Scalar theta)
    : paramAlpha_{{-1.0, 1.0}}
    , paramBeta_{{1.0-theta, theta}}
    , paramD_{{0.0, 1.0}};
    {}

    bool implicit () const final
    { return paramBeta_[1] > 0.0; }

    std::size_t numStages () const final
    { return 1; }

    Scalar paramAlpha (std::size_t, std::size_t i) const final
    { return paramAlpha_[i]; }

    Scalar paramBeta (std::size_t, std::size_t i) const final
    { return paramBeta_[i]; }

    Scalar paramD (std::size_t i) const final
    { return paramD_[i]; }

    std::string name () const override
    { return "theta scheme"; }

private:
    std::array<Scalar, 2> paramAlpha_;
    std::array<Scalar, 2> paramBeta_;
    std::array<Scalar, 2> paramD_;
};

/*!
 * \brief An explicit Euler time stepping scheme
 */
template<class Scalar>
class ExplicitEuler final : public Theta<Scalar>
{
public:
    ExplicitEuler() : Theta<Scalar>(0.0) {}

    std::string name () const final
    { return "explicit Euler"; }
};

/*!
 * \brief An implicit Euler time stepping scheme
 */
template<class Scalar>
class ImplicitEuler final : public Theta<Scalar>
{
public:
    ImplicitEuler() : Theta<Scalar>(1.0) {}

    std::string name () const final
    { return "implicit Euler"; }
};

/*!
 * \brief Classical explicit fourth order Runge-Kutta scheme
 */
template<class Scalar>
class RungeKuttaExplicitFourthOrder final : public RungeKuttaMethod<Scalar>
{
    RungeKuttaExplicitFourthOrder()
    : paramAlpha_{{{-1.0, 1.0, 0.0, 0.0, 0.0},
                   {-1.0, 0.0, 1.0, 0.0, 0.0},
                   {-1.0, 0.0, 0.0, 1.0, 0.0},
                   {-1.0, 0.0, 0.0, 0.0, 1.0}}}
    , paramBeta_{{{0.5, 0.0, 0.0, 0.0, 0.0},
                  {0.0, 0.5, 0.0, 0.0, 0.0},
                  {0.0, 0.0, 1.0, 0.0, 0.0},
                  {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0, 0.0}}}
    , paramD_{{0.0, 0.5, 0.5, 1.0, 1.0}};
    {}

    bool implicit () const final
    { return false; }

    std::size_t numStages () const final
    { return 4; }

    Scalar paramAlpha (std::size_t stage, std::size_t i) const final
    { return paramAlpha_[stage-1][i]; }

    Scalar paramBeta (std::size_t stage, std::size_t i) const final
    { return paramBeta_[stage-1][i]; }

    Scalar paramD (std::size_t i) const final
    { return paramD_[i]; }

    std::string name () const final
    { return "explicit Runge-Kutta 4th order"; }

private:
    std::array<std::array<Scalar, 5>, 4> paramAlpha_;
    std::array<std::array<Scalar, 5>, 4> paramBeta_;
    std::array<Scalar, 5> paramD_;
};

} // end namespace Dumux

#endif
