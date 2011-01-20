// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Some exceptions thrown in DUMUX
 */
#ifndef DUMUX_EXCEPTIONS_HH
#define DUMUX_EXCEPTIONS_HH

#include <dune/common/exceptions.hh>

#include <string>

namespace Dumux {
/*!
 * \brief Exception thrown if a fixable numerical problem occurs.
 *
 * (e.g. time step too big, etc.)
 */
class NumericalProblem : public Dune::Exception
{
public:
    // copy constructor
    NumericalProblem(const NumericalProblem &v)
        : Dune::Exception(v)
    {}

    // default constructor
    NumericalProblem()
    {}

    // constructor with error message
    NumericalProblem(const std::string &s)
    { this->message(s); }
 };
}

#endif
