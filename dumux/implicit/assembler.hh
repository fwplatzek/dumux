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
 * \brief An assembler for the global Jacobian matrix for fully implicit models.
 */
#ifndef DUMUX_IMPLICIT_ASSEMBLER_HH
#define DUMUX_IMPLICIT_ASSEMBLER_HH

#include "properties.hh"

namespace Dumux {

/*!
 * \ingroup ImplicitModel
 * \brief An assembler for the global Jacobian matrix for fully implicit models.
 */
template<class TypeTag>
class ImplicitAssembler
{
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    enum{ dim = GridView::dimension };
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    // copying the jacobian assembler is not a good idea
    ImplicitAssembler(const ImplicitAssembler &);

public:

    ImplicitAssembler() : problemPtr_(nullptr)
    {}

    /*!
     * \brief Initialize the jacobian assembler.
     *
     * At this point we can assume that all objects in the problem and
     * the model have been allocated. We can not assume that they are
     * fully initialized, though.
     *
     * \param problem The problem object
     */
    void init(Problem& problem)
    {
        problemPtr_ = &problem;

        // initialize the BCRS matrix
        asImp_().createMatrix_();

        // initialize the jacobian matrix with zeros
        *matrix_ = 0;

        // allocate the residual vector
        residual_.resize(problem.model().numDofs());
    }

    /*!
     * \brief Assemble the global Jacobian of the residual and the residual for the current solution.
     *
     * The current state of affairs (esp. the previous and the current
     * solutions) is represented by the model object.
     */
    void assemble()
    {
        bool succeeded;
        try {
            asImp_().assemble_();
            succeeded = true;
            if (gridView_().comm().size() > 1)
                succeeded = gridView_().comm().min(succeeded);
        }
        catch (Dumux::NumericalProblem &e)
        {
            std::cout << "rank " << problem_().gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
            if (gridView_().comm().size() > 1)
                succeeded = gridView_().comm().min(succeeded);
        }

        if (!succeeded) {
            DUNE_THROW(NumericalProblem,
                       "A process did not succeed in linearizing the system");
        }
    }

    /*!
     * \brief Return constant reference to global Jacobian matrix.
     */
    const JacobianMatrix& matrix() const
    { return *matrix_; }
    JacobianMatrix& matrix()
    { return *matrix_; }

    /*!
     * \brief Return constant reference to global residual vector.
     */
    const SolutionVector& residual() const
    { return residual_; }
    SolutionVector& residual()
    { return residual_; }

protected:

    // reset the global linear system of equations. if partial
    // reassemble is enabled, this means that the jacobian matrix must
    // only be erased partially!
    void resetSystem_()
    {
        // reset the right hand side.
        residual_ = 0.0;

        // reset the matrix
        (*matrix_) = 0;
    }

    // linearize the whole system
    void assemble_()
    {
        resetSystem_();

        // reassemble the elements...
        for (const auto& element : elements(gridView_())) {
            if (element.partitionType() == Dune::GhostEntity)
            {
                asImp_().assembleGhostElement_(element);
            }
            else
            {
                asImp_().assembleElement_(element);
            }
        }
    }

    // assemble an interior element
    void assembleElement_(const Element &element)
    {
        model_().localJacobian().assemble(element);

        if (!isBox)
        {
            auto globalI = elementMapper_().index(element);

            // update the right hand side
            residual_[globalI] = model_().localJacobian().residual(0);
            for (int j = 0; j < residual_[globalI].dimension; ++j)
                assert(std::isfinite(residual_[globalI][j]));

            // diagonal entry
            (*matrix_)[globalI][globalI] = model_().localJacobian().mat(0, 0);

            const auto& stencil = model_().stencils(element).neighborStencil();

            unsigned int j = 1;
            for (auto&& globalJ : stencil)
                (*matrix_)[globalI][globalJ] = model_().localJacobian().mat(0, j++);
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented, "Box assembly");
        }
    }

    // "assemble" a ghost element
    void assembleGhostElement_(const Element &element)
    {
        int globalI = elementMapper_().index(element);

        // update the right hand side
        residual_[globalI] = 0.0;

        // update the diagonal entry
        typedef typename JacobianMatrix::block_type BlockType;
        BlockType &J = (*matrix_)[globalI][globalI];
        for (int j = 0; j < BlockType::rows; ++j)
            J[j][j] = 1.0;
    }

protected:
    Problem &problem_()
    { return *problemPtr_; }
    const Problem &problem_() const
    { return *problemPtr_; }
    const Model &model_() const
    { return problem_().model(); }
    Model &model_()
    { return problem_().model(); }
    const GridView &gridView_() const
    { return problem_().gridView(); }
    const VertexMapper &vertexMapper_() const
    { return problem_().vertexMapper(); }
    const ElementMapper &elementMapper_() const
    { return problem_().elementMapper(); }

    Problem *problemPtr_;

    // the jacobian matrix
    std::shared_ptr<JacobianMatrix> matrix_;
    // the right-hand side
    SolutionVector residual_;

private:

    // Construct the BCRS matrix for the global jacobian
    void createMatrix_()
    {
        auto numDofs = problem_().model().numDofs();

        // allocate raw matrix
        matrix_ = std::make_shared<JacobianMatrix>(numDofs, numDofs, JacobianMatrix::random);

        // set the row sizes
        setRowSizes_();

        // set the indices
        addIndices_();
    }

    void setRowSizes_()
    {
        if (!isBox)
        {
            for (const auto& element : elements(gridView_()))
            {
                // the global index of the element at hand
                const auto globalI = elementMapper_().index(element);
                const auto& stencil = model_().stencils(element).elementStencil();

                matrix_->setrowsize(globalI, stencil.size());
            }
            matrix_->endrowsizes();
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented, "Box Assmebly");
        }
    }

    void addIndices_()
    {
        if (!isBox)
        {
            for (const auto& element : elements(gridView_()))
            {
                // the global index of the element at hand
                const auto globalI = elementMapper_().index(element);
                const auto& stencil = model_().stencils(element).elementStencil();


                for (auto&& globalJ : stencil)
                    matrix_->addindex(globalI, globalJ);
            }
            matrix_->endindices();
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented, "Box Assmebly");
        }
    }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

} // namespace Dumux

#endif
