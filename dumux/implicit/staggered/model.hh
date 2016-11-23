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
 * \brief Base class for all models which use the one-phase,
 *        fully implicit model.
 *        Adaption of the fully implicit scheme to the one-phase flow model.
 */

#ifndef DUMUX_STAGGERED_BASEMODEL_HH
#define DUMUX_STAGGERED_BASEMODEL_HH

// #include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include <dumux/implicit/model.hh>
#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup NavierStokesModel
 * \brief A single-phase, isothermal flow model using the fully implicit scheme.
 *
 * Single-phase, isothermal flow model, which uses a standard Darcy approach as the
 * equation for the conservation of momentum:
 * \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 * \f]
 *
 * and solves the mass continuity equation:
 * \f[
 \phi \frac{\partial \varrho}{\partial t} + \text{div} \left\lbrace
 - \varrho \frac{\textbf K}{\mu} \left( \textbf{grad}\, p -\varrho {\textbf g} \right) \right\rbrace = q,
 * \f]
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 * The model supports compressible as well as incompressible fluids.
 */
template<class TypeTag >
class StaggeredBaseModel : public ImplicitModel<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalFVGeometry) GlobalFVGeometry;
//     typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };
    using StencilsVector = typename GET_PROP_TYPE(TypeTag, StencilsVector);
    using Element = typename GridView::template Codim<0>::Entity;
    using Implementation = typename GET_PROP_TYPE(TypeTag, Model);
    using NewtonMethod = typename GET_PROP_TYPE(TypeTag, NewtonMethod);
    using NewtonController = typename GET_PROP_TYPE(TypeTag, NewtonController);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    static_assert(!isBox, "must not be box!");

public:

     /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem &problem)
    {
        this->problemPtr_ = &problem;

        this->updateBoundaryIndices_();

        this->globalFvGeometryPtr_ = std::make_shared<GlobalFVGeometry>(problem.gridView());
        this->globalFvGeometryPtr_->update(problem);

        this->uCur_[cellCenterIdx].resize(asImp_().numCellCenterDofs());
        this->uCur_[faceIdx].resize(asImp_().numFaceDofs());

        // apply initial solution
        // for compositional models initial the phase presence herein
        asImp_().applyInitialSolution_();

        // resize and update the volVars with the initial solution
        this->curGlobalVolVars_.update(problem, this->curSol());

        // update stencils
        this->stencilsVector_.update(problem);

        // update the flux variables caches
        this->globalfluxVarsCache_.update(problem);

        // initialize assembler and create matrix
        this->localJacobian_.init(problem);
        this->jacAsm_ = std::make_shared<JacobianAssembler>();
        this->jacAsm_->init(problem);

        // also set the solution of the "previous" time step to the
        // initial solution.
        this->uPrev_ = this->uCur_;
        this->prevGlobalVolVars_ = this->curGlobalVolVars_;
    }

     /*!
     * \brief Try to progress the model to the next timestep.
     *
     * \param solver The non-linear solver
     * \param controller The controller which specifies the behaviour
     *                   of the non-linear solver
     */
    bool update(NewtonMethod &solver,
                NewtonController &controller)
    {
#if HAVE_VALGRIND
        for (size_t i = 0; i < this->curSol()[cellCenterIdx].size(); ++i)
            Valgrind::CheckDefined(this->curSol()[cellCenterIdx][i]);
        for (size_t i = 0; i < this->curSol()[faceIdx].size(); ++i)
            Valgrind::CheckDefined(this->curSol()[faceIdx][i]);
#endif // HAVE_VALGRIND

        asImp_().updateBegin();

        int converged = solver.execute(controller);

        if (this->gridView_().comm().size() > 1)
        {
            converged = this->gridView_().comm().min(converged);
        }
        if (converged) {
            asImp_().updateSuccessful();
        }
        else
            asImp_().updateFailed();

#if HAVE_VALGRIND
        for (size_t i = 0; i < this->curSol()[cellCenterIdx].size(); ++i)
            Valgrind::CheckDefined(this->curSol()[cellCenterIdx][i]);
        for (size_t i = 0; i < this->curSol()[faceIdx].size(); ++i)
            Valgrind::CheckDefined(this->curSol()[faceIdx][i]);
#endif // HAVE_VALGRIND

        return converged;
    }

     /*!
     * \brief Applies the initial solution for all vertices of the grid.
     *
     * \todo the initial condition needs to be unique for
     *       each vertex. we should think about the API...
     */
    void applyInitialSolution_()
    {
        // first set the whole domain to zero
        this->uCur_ = Scalar(0.0);

        // iterate through leaf grid and evaluate initial
        // condition at the center of each sub control volume
        for (const auto& element : elements(this->gridView_()))
        {
            // deal with the current element, bind FVGeometry to it
            auto fvGeometry = localView(this->globalFvGeometry());
            fvGeometry.bindElement(element);

            // loop over sub control volumes
            for (auto&& scv : scvs(fvGeometry))
            {
                // let the problem do the dirty work of nailing down
                // the initial solution.
                auto initPriVars = this->problem_().initial(scv);

                auto dofIdxGlobal = scv.dofIndex();
                this->uCur_[cellCenterIdx][dofIdxGlobal] += initPriVars;
            }
        }
    }

     /*!
     * \brief Write the current solution for a vertex to a restart
     *        file.
     *
     * \param outstream The stream into which the vertex data should
     *                  be serialized to
     * \param entity The entity which's data should be
     *               serialized, i.e. a vertex for the box method
     *               and an element for the cell-centered method
     */
    template <class Entity>
    void serializeEntity(std::ostream &outstream,
                         const Entity &entity)
    {
        DUNE_THROW(Dune::NotImplemented, "deserializeEntity() not implemented yet");
//         int dofIdxGlobal = dofMapper().index(entity);
//         int dofIdxGlobal = dofMapper().index(entity);
//
//         // write phase state
//         if (!outstream.good()) {
//             DUNE_THROW(Dune::IOError,
//                        "Could not serialize vertex "
//                        << dofIdxGlobal);
//         }
//
//         for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
//             outstream << curSol()[dofIdxGlobal][eqIdx] << " ";
//         }
    }

     /*!
     * \brief Reads the current solution variables for a vertex from a
     *        restart file.
     *
     * \param instream The stream from which the vertex data should
     *                  be deserialized from
     * \param entity The entity which's data should be
     *               serialized, i.e. a vertex for the box method
     *               and an element for the cell-centered method
     */
    template <class Entity>
    void deserializeEntity(std::istream &instream,
                           const Entity &entity)
    {
        DUNE_THROW(Dune::NotImplemented, "deserializeEntity() not implemented yet");
//         int dofIdxGlobal = dofMapper().index(entity);
//
//         for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
//             if (!instream.good())
//                 DUNE_THROW(Dune::IOError,
//                            "Could not deserialize vertex "
//                            << dofIdxGlobal);
//             instream >> curSol()[dofIdxGlobal][eqIdx];
//         }
    }

    /*!
     * \brief \copybrief Dumux::ImplicitModel::addOutputVtkFields
     *
     * Specialization for the NavierStokesModel, adding the pressure and
     * the process rank to the VTK writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        // typedef Dune::BlockVector<Dune::FieldVector<double, dimWorld> > VectorField;

        // create the required scalar fields
        unsigned numDofs = this->numDofs();
        auto *p = writer.allocateManagedBuffer(numDofs);
        // VectorField *velocity = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        // ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        // if (velocityOutput.enableOutput())
        // {
        //     // initialize velocity field
        //     for (unsigned int i = 0; i < numDofs; ++i)
        //     {
        //         (*velocity)[i] = double(0);
        //     }
        // }

        unsigned numElements = this->gridView_().size(0);
        auto *rank = writer.allocateManagedBuffer(numElements);

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
            auto eIdx = this->elementMapper().index(element);
            (*rank)[eIdx] = this->gridView_().comm().rank();

            // get the local fv geometry
            auto fvGeometry = localView(this->globalFvGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(this->curGlobalVolVars());
            elemVolVars.bindElement(element, fvGeometry, this->curSol());

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                auto dofIdxGlobal = scv.dofIndex();

                (*p)[dofIdxGlobal] = volVars.pressure();
            }

            // velocity output
            //velocityOutput.calculateVelocity(*velocity, elemVolVars, fvGeometry, element, /*phaseIdx=*/0);
        }

        writer.attachDofData(*p, "p", isBox);
        // if (velocityOutput.enableOutput())
        // {
        //     writer.attachDofData(*velocity,  "velocity", isBox, dim);
        // }
        writer.attachCellData(*rank, "process rank");
    }

    /*!
     * \brief Add the vector fields for analysing the convergence of
     *        the newton method to the a VTK multi writer.
     *
     * \tparam MultiWriter The type of the VTK multi writer
     *
     * \param writer  The VTK multi writer object on which the fields should be added.
     * \param u       The solution function
     * \param deltaU  The delta of the solution function before and after the Newton update
     */
    template <class MultiWriter>
    void addConvergenceVtkFields(MultiWriter &writer,
                                 const SolutionVector &u,
                                 const SolutionVector &deltaU)
    {
//         typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
//
//         SolutionVector residual(u);
//         asImp_().globalResidual(residual, u);
//
//         // create the required scalar fields
//         unsigned numDofs = asImp_().numDofs();
//
//         // global defect of the two auxiliary equations
//         ScalarField* def[numEq];
//         ScalarField* delta[numEq];
//         ScalarField* x[numEq];
//         for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
//             x[eqIdx] = writer.allocateManagedBuffer(numDofs);
//             delta[eqIdx] = writer.allocateManagedBuffer(numDofs);
//             def[eqIdx] = writer.allocateManagedBuffer(numDofs);
//         }
//
//         for (unsigned int dofIdxGlobal = 0; dofIdxGlobal < u.size(); dofIdxGlobal++)
//         {
//             for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//             {
//                 (*x[eqIdx])[dofIdxGlobal] = u[dofIdxGlobal][eqIdx];
//                 (*delta[eqIdx])[dofIdxGlobal] = - deltaU[dofIdxGlobal][eqIdx];
//                 (*def[eqIdx])[dofIdxGlobal] = residual[dofIdxGlobal][eqIdx];
//             }
//         }
//
//         for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
//             std::ostringstream oss;
//             oss.str(""); oss << "x_" << eqIdx;
//             if (isBox)
//                 writer.attachVertexData(*x[eqIdx], oss.str());
//             else
//                 writer.attachCellData(*x[eqIdx], oss.str());
//             oss.str(""); oss << "delta_" << eqIdx;
//             if (isBox)
//                 writer.attachVertexData(*delta[eqIdx], oss.str());
//             else
//                 writer.attachCellData(*delta[eqIdx], oss.str());
//             oss.str(""); oss << "defect_" << eqIdx;
//             if (isBox)
//                 writer.attachVertexData(*def[eqIdx], oss.str());
//             else
//                 writer.attachCellData(*def[eqIdx], oss.str());
//         }
//
//         asImp_().addOutputVtkFields(u, writer);
    }

     /*!
     * \brief Compute the global residual for an arbitrary solution
     *        vector.
     *
     * \param residual Stores the result
     * \param u The solution for which the residual ought to be calculated
     */
    Scalar globalResidual(SolutionVector &residual,
                          const SolutionVector &u)
    {
        SolutionVector tmp(this->curSol());
        this->curSol() = u;
        Scalar res = globalResidual(residual);
        this->curSol() = tmp;
        return res;
    }

    /*!
     * \brief Compute the global residual for the current solution
     *        vector.
     *
     * \param residual Stores the result
     */
    Scalar globalResidual(SolutionVector &residual)
    {
        residual = 0;

        for (const auto& element : elements(this->gridView_())) {
            this->localResidual().eval(element);


            int globalI = this->elementMapper().index(element);
            residual[cellCenterIdx][globalI] = this->localResidual().residual(0);
        }

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (this->gridView_().comm().size() > 1)
            result2 = this->gridView_().comm().sum(result2);

        return std::sqrt(result2);
    }

     /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    size_t numDofs() const
    {
        return this->gridView_().size(0) + this->gridView_().size(1);
    }

     /*!
     * \brief Returns the number of cell center degrees of freedoms (DOFs)
     */
    size_t numCellCenterDofs() const
    {
        return this->gridView_().size(0);
    }

     /*!
     * \brief Returns the number of cell center degrees of freedoms (DOFs)
     */
    size_t numFaceDofs() const
    {
        return this->gridView_().size(1);
    }

    /*!
    * \brief Returns the size of a complete face dof stencil
    */
    size_t fullFaceToCellCenterStencilSize(const int idx) const
    {
        return this->stencilsVector_.fullFaceToCellCenterStencilSize(idx);
    }

    /*!
    * \brief Returns the size of a complete face dof stencil
    */
    size_t fullfaceToFaceStencilSize(const int idx) const
    {
        return this->stencilsVector_.fullfaceToFaceStencilSize(idx);
    }

    /*!
    * \brief Returns a unique pointer to the complete face dof stencils which is used once for setting up the global matrix and deleted afterwards
    */
    auto getFullFaceToCellCenterStencilsPtr()
    {
        return this->stencilsVector_.getFullFaceToCellCenterStencilsPtr();
    }

    /*!
    * \brief Returns a unique pointer to the complete face dof stencils which is used once for setting up the global matrix and deleted afterwards
    */
    auto getFullfaceToFaceStencilsPtr()
    {
        return this->stencilsVector_.getFullfaceToFaceStencilsPtr();
    }
private:

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

};
}

#include "propertydefaults.hh"

#endif
