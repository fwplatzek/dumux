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
 * \ingroup Linear
 * \brief Dumux sequential linear solver backends
 */
#ifndef DUMUX_SEQ_SOLVER_BACKEND_HH
#define DUMUX_SEQ_SOLVER_BACKEND_HH

#include <type_traits>
#include <tuple>
#include <utility>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/matrixmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/common/version.hh>
#include <dune/common/hybridutilities.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/matrix.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/amgbackend.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A general solver backend allowing arbitrary preconditioners and solvers.
 *
 * This class is used as a base class for specific solver-preconditioner
 * combinations. Several parameters from the group LinearSolver are read to
 * customize the solver and preconditioner:
 *
 * - Verbosity: determines how verbose the linear solver should print output.
 * - MaxIterations: the maximum number of iterations for the linear solver.
 * - ResidualReduction: the threshold for declaration of convergence.
 * - PreconditionerRelaxation: relaxation parameter for the preconditioner.
 * - PreconditionerIterations: usually specifies the number of times the
 *                             preconditioner is applied. In case of ILU(n),
 *                             it specifies the order of the applied ILU.
 */
class IterativePreconditionedSolverImpl
{
public:

    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    static bool solve(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                      const std::string& modelParamGroup = "")
    {
        Preconditioner precond(A, s.precondIter(), s.relaxation());

        // make a linear operator from a matrix
        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter linearOperator(A);

        Solver solver(linearOperator, precond, s.residReduction(), s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }

    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    static bool solveWithGMRes(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                               const std::string& modelParamGroup = "")
    {
        // get the restart threshold
        const int restartGMRes = getParamFromGroup<double>(modelParamGroup, "LinearSolver.GMResRestart");

        Preconditioner precond(A, s.precondIter(), s.relaxation());

        // make a linear operator from a matrix
        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter linearOperator(A);

        Solver solver(linearOperator, precond, s.residReduction(), restartGMRes, s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }

    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    static bool solveWithILU0Prec(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                                  const std::string& modelParamGroup = "")
    {
        Preconditioner precond(A, s.relaxation());

        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter operatorA(A);

        Solver solver(operatorA, precond, s.residReduction(), s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }

    // solve with RestartedGMRes (needs restartGMRes as additional argument)
    template<class Preconditioner, class Solver, class SolverInterface, class Matrix, class Vector>
    static bool solveWithILU0PrecGMRes(const SolverInterface& s, const Matrix& A, Vector& x, const Vector& b,
                                       const std::string& modelParamGroup = "")
    {
        // get the restart threshold
        const int restartGMRes = getParamFromGroup<int>(modelParamGroup, "LinearSolver.GMResRestart");

        Preconditioner precond(A, s.relaxation());

        using MatrixAdapter = Dune::MatrixAdapter<Matrix, Vector, Vector>;
        MatrixAdapter operatorA(A);

        Solver solver(operatorA, precond, s.residReduction(), restartGMRes, s.maxIter(), s.verbosity());

        Vector bTmp(b);

        Dune::InverseOperatorResult result;
        solver.apply(x, bTmp, result);

        return result.converged;
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU(n)-preconditioned BiCSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. It can be applied to
 * nonsymmetric matrices.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: ILU(n) incomplete LU factorization. The order n can be
 * provided by the parameter LinearSolver.PreconditionerIterations and controls
 * the fill-in. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class ILUnBiCGSTABBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqILU<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "ILUn preconditioned BiCGSTAB solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential SOR-preconditioned BiCSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. It can be applied to
 * nonsymmetric matrices.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: SOR successive overrelaxation method. The relaxation is
 * controlled by the parameter LinearSolver.PreconditionerRelaxation. In each
 * preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class SORBiCGSTABBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqSOR<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "SOR preconditioned BiCGSTAB solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential SSOR-preconditioned BiCGSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. While, it can be
 * applied to nonsymmetric matrices, the preconditioner SSOR assumes symmetry.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: SSOR symmetric successive overrelaxation method. The
 * relaxation is controlled by the parameter LinearSolver.PreconditionerRelaxation.
 * In each preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class SSORBiCGSTABBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqSSOR<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "SSOR preconditioned BiCGSTAB solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential GS-preconditioned BiCGSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. While, it can be
 * applied to nonsymmetric matrices, the preconditioner SSOR assumes symmetry.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: GS Gauss-Seidel method. It can be damped by the relaxation
 * parameter LinearSolver.PreconditionerRelaxation. In each preconditioning step,
 * it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class GSBiCGSTABBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqGS<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "SSOR preconditioned BiCGSTAB solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential Jacobi-preconditioned BiCSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. While, it can be
 * applied to nonsymmetric matrices, the preconditioner SSOR assumes symmetry.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: Jacobi method. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation. In each preconditioning step, it is
 * applied as often as given by the parameter LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class JacBiCGSTABBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqJac<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "Jac preconditioned BiCGSTAB solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU(n)-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: ILU(n) incomplete LU factorization. The order n can be
 * provided by the parameter LinearSolver.PreconditionerIterations and controls
 * the fill-in. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class ILUnCGBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqILU<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "ILUn preconditioned CG solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential SOR-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: SOR successive overrelaxation method. The relaxation is
 * controlled by the parameter LinearSolver.PreconditionerRelaxation. In each
 * preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class SORCGBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqSOR<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "SOR preconditioned CG solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential SSOR-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: SSOR symmetric successive overrelaxation method. The
 * relaxation is controlled by the parameter LinearSolver.PreconditionerRelaxation.
 * In each preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class SSORCGBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqSSOR<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "SSOR preconditioned CG solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential GS-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: GS Gauss-Seidel method. It can be damped by the relaxation
 * parameter LinearSolver.PreconditionerRelaxation. In each preconditioning step,
 * it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class GSCGBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqGS<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "GS preconditioned CG solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential Jacobi-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: Jacobi method. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation. In each preconditioning step, it is
 * applied as often as given by the parameter LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class JacCGBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqJac<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solve<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "GS preconditioned CG solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential SSOR-preconditioned GMRes solver.
 *
 * Solver: The GMRes (generalized minimal residual) method is an iterative
 * method for the numerical solution of a nonsymmetric system of linear
 * equations.\n
 * See: Saad, Y., Schultz, M. H. (1986). "GMRES: A generalized minimal residual
 * algorithm for solving nonsymmetric linear systems." SIAM J. Sci. and Stat.
 * Comput. 7: 856–869.
 *
 * Preconditioner: SSOR symmetric successive overrelaxation method. The
 * relaxation is controlled by the parameter LinearSolver.PreconditionerRelaxation.
 * In each preconditioning step, it is applied as often as given by the parameter
 * LinearSolver.PreconditionerIterations.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class SSORRestartedGMResBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqSSOR<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::RestartedGMResSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solveWithGMRes<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "SSOR preconditioned GMRes solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU(0)-preconditioned BiCGSTAB solver.
 *
 * Solver: The BiCGSTAB (stabilized biconjugate gradients method) solver has
 * faster and smoother convergence than the original BiCG. It can be applied to
 * nonsymmetric matrices.\n
 * See: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A Fast and Smoothly Converging
 * Variant of Bi-CG for the Solution of Nonsymmetric Linear Systems".
 * SIAM J. Sci. and Stat. Comput. 13 (2): 631–644. doi:10.1137/0913035.
 *
 * Preconditioner: ILU(0) incomplete LU factorization. The order 0 indicates
 * that no fill-in is allowed. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class ILU0BiCGSTABBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqILU<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::BiCGSTABSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solveWithILU0Prec<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "ILU0 preconditioned BiCGSTAB solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU(0)-preconditioned CG solver.
 *
 * Solver: CG (conjugate gradient) is an iterative method for solving linear
 * systems with a symmetric, positive definite matrix.\n
 * See:  Helfenstein, R., Koko, J. (2010). "Parallel preconditioned conjugate
 * gradient algorithm on GPU", Journal of Computational and Applied Mathematics,
 * Volume 236, Issue 15, Pages 3584–3590, http://dx.doi.org/10.1016/j.cam.2011.04.025.
 *
 * Preconditioner: ILU(0) incomplete LU factorization. The order 0 indicates
 * that no fill-in is allowed. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class ILU0CGBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqILU<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::CGSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solveWithILU0Prec<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "ILU0 preconditioned BiCGSTAB solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU0-preconditioned GMRes solver.
 *
 * Solver: The GMRes (generalized minimal residual) method is an iterative
 * method for the numerical solution of a nonsymmetric system of linear
 * equations.\n
 * See: Saad, Y., Schultz, M. H. (1986). "GMRES: A generalized minimal residual
 * algorithm for solving nonsymmetric linear systems." SIAM J. Sci. and Stat.
 * Comput. 7: 856–869.
 *
 * Preconditioner: ILU(0) incomplete LU factorization. The order 0 indicates
 * that no fill-in is allowed. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class ILU0RestartedGMResBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqILU<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::RestartedGMResSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solveWithILU0PrecGMRes<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }

    std::string name() const
    {
        return "ILU0 preconditioned BiCGSTAB solver";
    }
};

/*!
 * \ingroup Linear
 * \brief Sequential ILU(n)-preconditioned GMRes solver.
 *
 * Solver: The GMRes (generalized minimal residual) method is an iterative
 * method for the numerical solution of a nonsymmetric system of linear
 * equations.\n
 * See: Saad, Y., Schultz, M. H. (1986). "GMRES: A generalized minimal residual
 * algorithm for solving nonsymmetric linear systems." SIAM J. Sci. and Stat.
 * Comput. 7: 856–869.
 *
 * Preconditioner: ILU(n) incomplete LU factorization. The order n can be
 * provided by the parameter LinearSolver.PreconditionerIterations and controls
 * the fill-in. It can be damped by the relaxation parameter
 * LinearSolver.PreconditionerRelaxation.\n
 * See: Golub, G. H., and Van Loan, C. F. (2012). Matrix computations. JHU Press.
 */
class ILUnRestartedGMResBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel: set this to more than one if the matrix to solve is nested multiple times
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using Preconditioner = Dune::SeqILU<Matrix, Vector, Vector, precondBlockLevel>;
        using Solver = Dune::RestartedGMResSolver<Vector>;

        return IterativePreconditionedSolverImpl::template solveWithGMRes<Preconditioner, Solver>(*this, A, x, b, this->paramGroup());
    }
};

#if HAVE_SUPERLU
/*!
 * \ingroup Linear
 * \brief Direct linear solver using the SuperLU library.
 *
 * See: Li, X. S. (2005). "An overview of SuperLU: Algorithms, implementation,
 * and user interface." ACM Transactions on Mathematical Software (TOMS) 31(3): 302-325.
 * http://crd-legacy.lbl.gov/~xiaoye/SuperLU/
 */
class SuperLUBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel is unused and just here for compatibility with iterative solvers
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static_assert(isBCRSMatrix<Matrix>::value, "SuperLU only works with BCRS matrices!");
        using BlockType = typename Matrix::block_type;
        static_assert(BlockType::rows == BlockType::cols, "Matrix block must be quadratic!");
        constexpr auto blockSize = BlockType::rows;

        Dune::SuperLU<Matrix> solver(A, this->verbosity() > 0);

        Vector bTmp(b);
        solver.apply(x, bTmp, result_);

        int size = x.size();
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < blockSize; j++)
            {
                using std::isnan;
                using std::isinf;
                if (isnan(x[i][j]) || isinf(x[i][j]))
                {
                    result_.converged = false;
                    break;
                }
            }
        }

        return result_.converged;
    }

    std::string name() const
    {
        return "SuperLU solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    Dune::InverseOperatorResult result_;
};
#endif // HAVE_SUPERLU

#if HAVE_UMFPACK
/*!
 * \ingroup Linear
 * \brief Direct linear solver using the UMFPack library.
 *
 * See: Davis, Timothy A. (2004). "Algorithm 832". ACM Transactions on
 * Mathematical Software 30 (2): 196–199. doi:10.1145/992200.992206.
 * http://faculty.cse.tamu.edu/davis/suitesparse.html
 */
class UMFPackBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    // precondBlockLevel is unused and just here for compatibility with iterative solvers
    template<int precondBlockLevel = 1, class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        static_assert(isBCRSMatrix<Matrix>::value, "UMFPack only works with BCRS matrices!");
        using BlockType = typename Matrix::block_type;
        static_assert(BlockType::rows == BlockType::cols, "Matrix block must be quadratic!");
        constexpr auto blockSize = BlockType::rows;

        Dune::UMFPack<Matrix> solver(A, this->verbosity() > 0);

        Vector bTmp(b);
        solver.apply(x, bTmp, result_);

        int size = x.size();
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < blockSize; j++)
            {
                using std::isnan;
                using std::isinf;
                if (isnan(x[i][j]) || isinf(x[i][j]))
                {
                    result_.converged = false;
                    break;
                }
            }
        }

        return result_.converged;
    }

    std::string name() const
    {
        return "UMFPack solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    Dune::InverseOperatorResult result_;
};
#endif // HAVE_UMFPACK


/*!
 * \name Solver for MultiTypeBlockMatrix's
 */
// \{

/*
 * A simple ilu0 block diagonal preconditioner
 */
template<class M, class X, class Y, int blockLevel = 2>
class BlockDiagILU0Preconditioner : public Dune::Preconditioner<X, Y>
{
    template<std::size_t i>
    using DiagBlockType = std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;

    template<std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<X>()[Dune::index_constant<i>{}])>;

#if DUNE_VERSION_NEWER(DUNE_ISTL,2,6)
    template<std::size_t i>
    using BlockILU = Dune::SeqILU<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>, blockLevel-1>;
#else
    template<std::size_t i>
    using BlockILU = Dune::SeqILU0<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>, blockLevel-1>;
#endif

    using ILUTuple = typename makeFromIndexedType<std::tuple, BlockILU, std::make_index_sequence<M::size()> >::type;

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = typename std::decay_t<M>;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    /*! \brief Constructor.

       Constructor gets all parameters to operate the prec.
       \param A The (multi type block) matrix to operate on.
       \param w The relaxation factor.
     */
    BlockDiagILU0Preconditioner(const M& m, double w = 1.0)
    : BlockDiagILU0Preconditioner(m, w, std::make_index_sequence<M::size()>{})
    {
        static_assert(blockLevel >= 2, "Only makes sense for MultiTypeBlockMatrix!");
    }

    /*!
       \brief Prepare the preconditioner.

       \copydoc Preconditioner::pre(X&,Y&)
     */
    void pre (X& v, Y& d) final {}

    /*!
     * \brief Apply the preconditoner.
     * \copydoc Preconditioner::apply(X&,const Y&)
     */
    void apply (X& v, const Y& d) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(ilu_)), [&](const auto i)
        {
            std::get<decltype(i)::value>(ilu_).apply(v[i], d[i]);
        });
    }

    /*!
     * \brief Clean up.
     * \copydoc Preconditioner::post(X&)
     */
    void post (X&) final {}

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const final
    {
        return Dune::SolverCategory::sequential;
    }

private:
    template<std::size_t... Is>
    BlockDiagILU0Preconditioner (const M& m, double w, std::index_sequence<Is...> is)
    : ilu_(std::make_tuple(BlockILU<Is>(m[Dune::index_constant<Is>{}][Dune::index_constant<Is>{}], w)...))
    {}

    ILUTuple ilu_;
};


/*
 * \brief A simple ilu0 block diagonal preconditioned BiCGSTABSolver
 * \note expects a system as a multi-type block-matrix
 * | A  B |
 * | C  D |
 */
class BlockDiagILU0BiCGSTABSolver : public LinearSolver
{

public:
    using LinearSolver::LinearSolver;

    // Solve saddle-point problem using a Schur complement based preconditioner
    template<int precondBlockLevel = 2, class Matrix, class Vector>
    bool solve(const Matrix& M, Vector& x, const Vector& b)
    {
        BlockDiagILU0Preconditioner<Matrix, Vector, Vector> preconditioner(M);
        Dune::MatrixAdapter<Matrix, Vector, Vector> op(M);
        Dune::BiCGSTABSolver<Vector> solver(op, preconditioner, this->residReduction(),
                                            this->maxIter(), this->verbosity());
        auto bTmp(b);
        solver.apply(x, bTmp, result_);

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

    std::string name() const
    { return "block-diagonal ILU0 preconditioned BiCGSTAB solver"; }

private:
    Dune::InverseOperatorResult result_;
};

/*
 * A simple ilu0 block diagonal preconditioner
 */
template<class M, class X, class Y, int blockLevel = 2>
class BlockDiagAMGPreconditioner : public Dune::Preconditioner<X, Y>
{
    template<std::size_t i>
    using DiagBlockType = std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;

    template<std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<X>()[Dune::index_constant<i>{}])>;

    template<std::size_t i>
    using Smoother = Dune::SeqSSOR<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>>;

    template<std::size_t i>
    using LinearOperator = Dune::MatrixAdapter<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>>;

    template<std::size_t i>
    using ScalarProduct = Dune::SeqScalarProduct<VecBlockType<i>>;

    template<std::size_t i>
    using BlockAMG = Dune::Amg::AMG<LinearOperator<i>, VecBlockType<i>, Smoother<i>>;

    using AMGTuple = typename makeFromIndexedType<std::tuple, BlockAMG, std::make_index_sequence<M::size()> >::type;

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = typename std::decay_t<M>;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    /*! \brief Constructor.

       Constructor gets all parameters to operate the prec.
       \param A The (multi type block) matrix to operate on.
       \param w The relaxation factor.
     */
    template<class LOP, class Criterion, class SmootherArgs>
    BlockDiagAMGPreconditioner(const LOP& lop, const Criterion& c, const SmootherArgs& sa)
    : BlockDiagAMGPreconditioner(lop, c, sa, std::make_index_sequence<M::size()>{})
    {
        static_assert(blockLevel >= 2, "Only makes sense for MultiTypeBlockMatrix!");
    }

    /*!
       \brief Prepare the preconditioner.

       \copydoc Preconditioner::pre(X&,Y&)
     */
    void pre (X& v, Y& d) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(amg_)), [&](const auto i)
        {
            std::get<decltype(i)::value>(amg_).pre(v[i], d[i]);
        });
    }

    /*!
     * \brief Apply the preconditoner.
     * \copydoc Preconditioner::apply(X&,const Y&)
     */
    void apply (X& v, const Y& d) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(amg_)), [&](const auto i)
        {
            std::get<decltype(i)::value>(amg_).apply(v[i], d[i]);
        });
    }

    /*!
     * \brief Clean up.
     * \copydoc Preconditioner::post(X&)
     */
    void post (X& v) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(amg_)), [&](const auto i)
        {
            std::get<decltype(i)::value>(amg_).post(v[i]);
        });
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const final
    {
        return Dune::SolverCategory::sequential;
    }

private:
    template<class LOP, class Criterion, class SmootherArgs, std::size_t... Is>
    BlockDiagAMGPreconditioner (const LOP& lop, const Criterion& c, const SmootherArgs& sa, std::index_sequence<Is...> is)
    : amg_(std::make_tuple(BlockAMG<Is>(*std::get<Is>(lop), *std::get<Is>(c), *std::get<Is>(sa))...))
    {}

    AMGTuple amg_;
};

/*
 * \brief A simple ilu0 block diagonal preconditioned BiCGSTABSolver
 * \note expects a system as a multi-type block-matrix
 * | A  B |
 * | C  D |
 */
class BlockDiagAMGBiCGSTABSolver : public LinearSolver
{
    template<class M, std::size_t i>
    using DiagBlockType = std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;

    template<class X, std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<X>()[Dune::index_constant<i>{}])>;

    template<class M, class X, std::size_t i>
    using Smoother = Dune::SeqSSOR<DiagBlockType<M, i>, VecBlockType<X, i>, VecBlockType<X, i>>;

    template<class M, class X, std::size_t i>
    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother<M, X, i>>::Arguments;

    template<class M, std::size_t i>
    using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<DiagBlockType<M, i>, Dune::Amg::FirstDiagonal>>;

    template<class M, class X, std::size_t i>
    using LinearOperator = Dune::MatrixAdapter<DiagBlockType<M, i>, VecBlockType<X, i>, VecBlockType<X, i>>;

public:
    using LinearSolver::LinearSolver;

    // Solve saddle-point problem using a Schur complement based preconditioner
    template<int precondBlockLevel = 2, class Matrix, class Vector>
    bool solve(const Matrix& m, Vector& x, const Vector& b)
    {
        //! \todo Check whether the default accumulation mode atOnceAccu is needed.
        //! \todo make parameters changeable at runtime from input file / parameter tree
        Dune::Amg::Parameters params(15, 2000, 1.2, 1.6, Dune::Amg::atOnceAccu);
        params.setDefaultValuesIsotropic(3);
        params.setDebugLevel(this->verbosity());

        auto criterion = makeCriterion_<Criterion, Matrix>(params, std::make_index_sequence<Matrix::size()>{});
        auto smootherArgs = makeSmootherArgs_<SmootherArgs, Matrix, Vector>(std::make_index_sequence<Matrix::size()>{});

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(m)), [&](const auto i)
        {
            auto& args = std::get<decltype(i)::value>(smootherArgs);
            args->iterations = 1;
            args->relaxationFactor = 1;
        });

        auto linearOperator = makeLinearOperator_<LinearOperator, Matrix, Vector>(m, std::make_index_sequence<Matrix::size()>{});

        BlockDiagAMGPreconditioner<Matrix, Vector, Vector> preconditioner(linearOperator, criterion, smootherArgs);

        Dune::MatrixAdapter<Matrix, Vector, Vector> op(m);
        Dune::BiCGSTABSolver<Vector> solver(op, preconditioner, this->residReduction(),
                                            this->maxIter(), this->verbosity());
        auto bTmp(b);
        solver.apply(x, bTmp, result_);

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

    std::string name() const
    { return "block-diagonal ILU0 preconditioned BiCGSTAB solver"; }

private:
    template<template<class M, std::size_t i> class Criterion, class Matrix, class Params, std::size_t... Is>
    auto makeCriterion_(const Params& p, std::index_sequence<Is...>)
    {
        return std::make_tuple(std::make_shared<Criterion<Matrix, Is>>(p)...);
    }

    template<template<class M, class X, std::size_t i> class SmootherArgs, class Matrix, class Vector, std::size_t... Is>
    auto makeSmootherArgs_(std::index_sequence<Is...>)
    {
        return std::make_tuple(std::make_shared<SmootherArgs<Matrix, Vector, Is>>()...);
    }

    template<template<class M, class X, std::size_t i> class LinearOperator, class Matrix, class Vector, std::size_t... Is>
    auto makeLinearOperator_(const Matrix& m, std::index_sequence<Is...>)
    {
        return std::make_tuple(std::make_shared<LinearOperator<Matrix, Vector, Is>>(m[Dune::index_constant<Is>{}][Dune::index_constant<Is>{}])...);
    }

    Dune::InverseOperatorResult result_;
};

/*!
  \brief Schur complement inverse operator.
*/
template<class AType, class BType, class CType, class DType, class InvOpAType, class X, class Y>
class SchurComplement : public Dune::LinearOperator<X,Y>
{
public:
    // export types
    typedef DType matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

    //! constructor: just store a reference to a matrix
    explicit SchurComplement (const AType& A, const BType& B,
                              const CType& C, const DType& D,
                              const std::shared_ptr<InvOpAType>& invA)
    : A_(A), B_(B), C_(C), D_(D), invA_(invA) {}

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply (const X& x, Y& y) const
    {
        // apply B, note x and aTmp1 have different size (Bx)
        X aTmp1(A_.N());
        aTmp1 = 0.0;
        B_.mv(x, aTmp1);

        // apply A^-1 (A^-1Bx)
        auto aTmp2 = aTmp1;

        // if we use the exact Schur complement the eigenvalues are 1 and GMRes converges in 2 steps
        // ... a single V-cycle of AMG instead
        Dune::InverseOperatorResult result;
        invA_->apply(aTmp2, aTmp1, result);
        // invA_.pre(aTmp2, aTmp1);
        // invA_.apply(aTmp2, aTmp1);
        // invA_.post(aTmp2);

        // apply C (CA^-1Bx)
        C_.mv(aTmp2, y);
        // switch signs
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
        auto tmp = y;
        this->apply(x, tmp);
        y.axpy(alpha, tmp);
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
      return Dune::SolverCategory::sequential;
    }

  private:
    const AType& A_;
    const BType& B_;
    const CType& C_;
    const DType& D_;
    std::shared_ptr<InvOpAType> invA_;
};

template <class GridView>
class SchurComplementSolver : public LinearSolver
{
public:
    SchurComplementSolver() = default;

    // expects a system as a multi-type matrix
    // | A  B |
    // | C  D |

    // Solve saddle-point problem using a Schur complement based preconditioner
    template<int precondBlockLevel = 1, class Matrix, class Vector, class FaceVector>
    bool solve(const Matrix& coefficientMatrix, Vector& CCVector, const Vector& b, const FaceVector& sampleFaceVec)
    {
        int verbosity = getParamFromGroup<int>(this->paramGroup(), "LinearSolver2.schurVerbosity");

        using namespace Dune::Indices;
        auto& A = coefficientMatrix[_1][_1];
        auto& B = coefficientMatrix[_1][_0];
        auto& C = coefficientMatrix[_0][_1];
        auto& D = coefficientMatrix[_0][_0];

        using AType = std::decay_t<decltype(A)>;
        using BType = std::decay_t<decltype(B)>;
        using CType = std::decay_t<decltype(C)>;
        using DType = std::decay_t<decltype(D)>;

        using VVector = std::decay_t<decltype(sampleFaceVec)>;
        using PVector = std::decay_t<decltype(CCVector)>;

        // or use an AMG approximation
        static const int precVerbosityA = getParamFromGroup<int>(this->paramGroup(), "LinearSolver2.precVerbosityA");
        static const int precMaxIter = getParamFromGroup<int>(this->paramGroup(), "LinearSolver2.precMaxIter");
        const double precResidReduction = getParamFromGroup<double>(this->paramGroup(), "LinearSolver2.precResidReduction");
        static const int precVerbositySchurAppr = getParamFromGroup<int>(this->paramGroup(), "LinearSolver2.precVerbositySchurAppr");

        Dune::Amg::Parameters params(15, 2000, 1.2, 1.6, Dune::Amg::atOnceAccu);
        params.setDefaultValuesIsotropic(GridView::Traits::Grid::dimension);
        params.setDebugLevel(verbosity);
        params.setGamma(1); // one V-cycle

        // matrix A
        using LinearOperatorA = Dune::MatrixAdapter<AType, VVector, VVector>;
        using ScalarProductA = Dune::SeqScalarProduct<VVector>;
        using SmootherA = Dune::SeqSSOR<AType, VVector, VVector>;
        using CommA = Dune::Amg::SequentialInformation;
        auto Aop = std::make_shared<LinearOperatorA>(coefficientMatrix[_1][_1]);
        auto spA = std::make_shared<ScalarProductA>();
        auto commA = std::make_shared<CommA>();

        using SmootherArgsA = typename Dune::Amg::SmootherTraits<SmootherA>::Arguments;
        using CriterionA = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<AType, Dune::Amg::FirstDiagonal> >;
        CriterionA criterionA(params);
        SmootherArgsA smootherArgsA;
        smootherArgsA.iterations = 1;
        smootherArgsA.relaxationFactor = 1.0;

        using InnerPreconditionerA = Dune::Amg::AMG<LinearOperatorA, VVector, SmootherA, CommA>;
        auto innerPrecA = std::make_shared<InnerPreconditionerA>(*Aop, criterionA, smootherArgsA, *commA);

        using InnerSolverA = Dune::RestartedFlexibleGMResSolver<VVector>;
        auto invOpA = std::make_shared<InnerSolverA>(*Aop, *spA, *innerPrecA, precResidReduction, 100, precMaxIter, precVerbosityA);
        using InvOpVelType = std::decay_t<decltype(*invOpA)>;

        //schurApproximate
        DType schurApproximate = schurApproximate_(coefficientMatrix);

        using LinearOperatorSchurAppr = Dune::MatrixAdapter<DType, PVector, PVector>;
        using ScalarProductSchurAppr = Dune::SeqScalarProduct<PVector>;
        using SmootherSchurAppr = Dune::SeqSSOR<DType, PVector, PVector>;
        using CommSchurAppr = Dune::Amg::SequentialInformation;
        auto schurApprOp = std::make_shared<LinearOperatorSchurAppr>(schurApproximate);
        auto spSchurAppr = std::make_shared<ScalarProductSchurAppr>();
        auto commSchurAppr = std::make_shared<CommSchurAppr>();

        using SmootherArgsSchurAppr = typename Dune::Amg::SmootherTraits<SmootherSchurAppr>::Arguments;
        using CriterionSchurAppr = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<DType, Dune::Amg::FirstDiagonal> >;
        CriterionSchurAppr criterionSchurAppr(params);
        SmootherArgsSchurAppr smootherArgsSchurAppr;
        smootherArgsSchurAppr.iterations = 1;
        smootherArgsSchurAppr.relaxationFactor = 1.0;

        using InnerPreconditionerSchurAppr = Dune::Amg::AMG<LinearOperatorSchurAppr, PVector, SmootherSchurAppr, CommSchurAppr>;
        auto innerPrecASchurAppr = std::make_shared<InnerPreconditionerSchurAppr>(*schurApprOp, criterionSchurAppr, smootherArgsSchurAppr, *commSchurAppr);

        using InnerSolverSchurAppr = Dune::RestartedFlexibleGMResSolver<PVector>;
        auto invOpSchurAppr = std::make_shared<InnerSolverSchurAppr>(*schurApprOp, *spSchurAppr, *innerPrecASchurAppr, precResidReduction, 100, precMaxIter, precVerbositySchurAppr);

        auto schurPrec = std::make_shared<Dune::InverseOperator2Preconditioner<InnerSolverSchurAppr>>(*invOpSchurAppr);

        //schur solver
        using SchurComplementType = SchurComplement<AType, BType, CType, DType, InvOpVelType, PVector, PVector>;
        auto schur = std::make_shared<SchurComplementType>(A, B, C, D, invOpA);

        const int maxIter = getParamFromGroup<double>(this->paramGroup(), "LinearSolver2.schurMaxIter");
        const double residReduction = getParamFromGroup<double>(this->paramGroup(), "LinearSolver2.schurResidualReduction");
        const int restartGMRes = getParamFromGroup<int>(this->paramGroup(), "LinearSolver2.schurGMResRestart");

        Vector bTmp(b);
        using Solver = Dune::RestartedGMResSolver<Vector>;
        Solver solver(*schur, *schurPrec, residReduction, restartGMRes, maxIter, verbosity);
        solver.apply(CCVector, bTmp, result_);

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

private:
    template <class Matrix>
    auto schurApproximate_(const Matrix& coefficientMatrix)
    {
        using namespace Dune::Indices;
        auto& A = coefficientMatrix[_1][_1];
        auto& B = coefficientMatrix[_1][_0];
        auto& C = coefficientMatrix[_0][_1];
        auto& D = coefficientMatrix[_0][_0];

        using AType = std::decay_t<decltype(A)>;
        using BType = std::decay_t<decltype(B)>;
        using DType = std::decay_t<decltype(D)>;

        const std::size_t numDofsFaceReduced = A.N();

        //get a diagonal matrix
        AType invDiagA;
        invDiagA.setBuildMode(AType::random);
        setInvDiagAPattern_(invDiagA, numDofsFaceReduced);
        auto row = A.begin();

        for(; row != A.end(); ++row)
        {
            using size_type = typename AType::size_type;
            size_type rowIdx = row.index();

            //invDigaA = inverse(diagonal(A))
            invDiagA[rowIdx][rowIdx] = 1./(A[rowIdx][rowIdx]);
        }

        BType invAB;
        Dune::matMultMat(invAB, invDiagA, B);

        DType schurApproximate;
        Dune::matMultMat(schurApproximate, C, invAB);

        return schurApproximate;
    }

    /*!
     * \brief Resizes the  matrix invDiagA and sets the matrix' sparsity pattern.
     */
    template <class AType>
    void setInvDiagAPattern_(AType& invDiagA, std::size_t numDofsFaceReduced)
    {
        // set the size of the sub-matrizes
        invDiagA.setSize(numDofsFaceReduced, numDofsFaceReduced);

        // set occupation pattern of the coefficient matrix
        Dune::MatrixIndexSet occupationPatternInvDiagA;
        occupationPatternInvDiagA.resize(numDofsFaceReduced, numDofsFaceReduced);

        // evaluate the acutal pattern
        for (int i = 0; i < invDiagA.N(); ++i)
        {
             occupationPatternInvDiagA.add(i, i);
        }

        occupationPatternInvDiagA.exportIdx(invDiagA);
    }

    Dune::InverseOperatorResult result_;
};


// \}

} // end namespace Dumux

#endif
