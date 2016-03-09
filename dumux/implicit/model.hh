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
 * \brief Base class for fully-implicit models
 */
#ifndef DUMUX_IMPLICIT_MODEL_HH
#define DUMUX_IMPLICIT_MODEL_HH

#include <dune/geometry/type.hh>
#include <dune/istl/bvector.hh>

#include "properties.hh"
#include <dumux/implicit/adaptive/gridadaptproperties.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/parallel/vertexhandles.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief The base class for the vertex centered finite volume
 *        discretization scheme.
 */
template<class TypeTag>
class ImplicitModel
{
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariablesVector) VolumeVariablesVector;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariablesVector) FluxVariablesVector;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension
    };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometryVector) FVElementGeometryVector;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, LocalJacobian) LocalJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename Dune::ReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<CoordScalar, dim> ReferenceElement;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    // copying a model is not a good idea
    ImplicitModel(const ImplicitModel &);

public:
    /*!
     * \brief The constructor.
     */
    ImplicitModel() : problemPtr_(nullptr) {}

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem &problem)
    {
        problemPtr_ = &problem;

        updateBoundaryIndices_();

        fvGeometries_ = std::make_shared<FVElementGeometryVector>(gridView_());
        fvGeometries_.update();

        int numDofs = asImp_().numDofs();
        uCur_.resize(numDofs);
        if (isBox)
            boxVolume_.resize(numDofs);

        localJacobian_.init(problem_());
        jacAsm_ = std::make_shared<JacobianAssembler>();
        jacAsm_->init(problem_());

        asImp_().applyInitialSolution_();

        // update the volVars with the initial solution
        volVarsVector_.update(problem_(), curSol());

        // also set the solution of the "previous" time step to the
        // initial solution.
        uPrev_ = uCur_;
        prevVolVarsVector_ = volVarsVector_;
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
        SolutionVector tmp(curSol());
        curSol() = u;
        Scalar res = globalResidual(residual);
        curSol() = tmp;
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

        for (const auto& element : elements(gridView_())) {
            localResidual().eval(element);

            if (isBox)
            {
                for (int i = 0; i < element.subEntities(dim); ++i)
                {
                    int globalI = vertexMapper().subIndex(element, i, dim);
                    residual[globalI] += localResidual().residual(i);
                }
            }
            else
            {
                int globalI = elementMapper().index(element);
                residual[globalI] = localResidual().residual(0);
            }
        }

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (gridView_().comm().size() > 1)
            result2 = gridView_().comm().sum(result2);

        // add up the residuals on the process borders
        if (isBox && gridView_().comm().size() > 1) {
            VertexHandleSum<PrimaryVariables, SolutionVector, VertexMapper>
                sumHandle(residual, vertexMapper());
            gridView_().communicate(sumHandle,
                                    Dune::InteriorBorder_InteriorBorder_Interface,
                                    Dune::ForwardCommunication);
        }

        return std::sqrt(result2);
    }

    /*!
     * \brief Compute the integral over the domain of the storage
     *        terms of all conservation quantities.
     *
     * \param storage Stores the result
     */
    void globalStorage(PrimaryVariables &storage)
    {
        storage = 0;

        for (const auto& element : elements(gridView_()))
        {
            if(element.partitionType() == Dune::InteriorEntity)
            {
                localResidual().evalStorage(element);

                if (isBox)
                {
                    for (int i = 0; i < element.subEntities(dim); ++i)
                    {
                        storage += localResidual().storageTerm()[i];
                    }
                }
                else
                {
                    storage += localResidual().storageTerm()[0];
                }
            }
        }

        if (gridView_().comm().size() > 1)
            storage = gridView_().comm().sum(storage);
    }

    /*!
     * \brief Returns the volume \f$\mathrm{[m^3]}\f$ of a given control volume.
     *
     * \param vIdxGlobal The global index of the control volume's
     *                  associated vertex
     */
    Scalar boxVolume(const int vIdxGlobal) const
    {
        if (isBox)
        {
            return boxVolume_[vIdxGlobal][0];
        }
        else
        {
            DUNE_THROW(Dune::InvalidStateException,
                       "requested box volume for cell-centered model");
        }
    }

    void adaptVariableSize()
    {
        uCur_.resize(numDofs());
    }

    /*!
     * \brief Reference to the current solution as a block vector.
     */
    const SolutionVector &curSol() const
    { return uCur_; }

    /*!
     * \brief Reference to the current solution as a block vector.
     */
    SolutionVector &curSol()
    { return uCur_; }

    /*!
     * \brief Reference to the previous solution as a block vector.
     */
    const SolutionVector &prevSol() const
    { return uPrev_; }

    /*!
     * \brief Reference to the previous solution as a block vector.
     */
    SolutionVector &prevSol()
    { return uPrev_; }

    /*!
     * \brief Returns the operator assembler for the global jacobian of
     *        the problem.
     */
    JacobianAssembler &jacobianAssembler()
    { return *jacAsm_; }

    /*!
     * \copydoc jacobianAssembler()
     */
    const JacobianAssembler &jacobianAssembler() const
    { return *jacAsm_; }

    /*!
     * \brief Returns the local jacobian which calculates the local
     *        stiffness matrix for an arbitrary element.
     *
     * The local stiffness matrices of the element are used by
     * the jacobian assembler to produce a global linerization of the
     * problem.
     */
    LocalJacobian &localJacobian()
    { return localJacobian_; }
    /*!
     * \copydoc localJacobian()
     */
    const LocalJacobian &localJacobian() const
    { return localJacobian_; }

    /*!
     * \brief Returns the local residual function.
     */
    LocalResidual &localResidual()
    { return localJacobian().localResidual(); }
    /*!
     * \copydoc localResidual()
     */
    const LocalResidual &localResidual() const
    { return localJacobian().localResidual(); }

    /*!
     * \brief Returns the maximum relative shift between two vectors of
     *        primary variables.
     *
     * \param priVars1 The first vector of primary variables
     * \param priVars2 The second vector of primary variables
     */
    Scalar relativeShiftAtDof(const PrimaryVariables &priVars1,
                              const PrimaryVariables &priVars2)
    {
        Scalar result = 0.0;
        for (int j = 0; j < numEq; ++j) {
            Scalar eqErr = std::abs(priVars1[j] - priVars2[j]);
            eqErr /= std::max<Scalar>(1.0, std::abs(priVars1[j] + priVars2[j])/2);

            result = std::max(result, eqErr);
        }
        return result;
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
        for (size_t i = 0; i < curSol().size(); ++i)
            Valgrind::CheckDefined(curSol()[i]);
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
        for (size_t i = 0; i < curSol().size(); ++i) {
            Valgrind::CheckDefined(curSol()[i]);
        }
#endif // HAVE_VALGRIND

        return converged;
    }

    /*!
     * \brief Check the plausibility of the current solution
     *
     *        This has to be done by the actual model, it knows
     *        best, what (ranges of) variables to check.
     *        This is primarily a hook
     *        which the actual model can overload.
     */
    void checkPlausibility() const
    { }

    /*!
     * \brief Called by the update() method before it tries to
     *        apply the newton method. This is primarily a hook
     *        which the actual model can overload.
     */
    void updateBegin()
    {
        if(GET_PROP_VALUE(TypeTag, AdaptiveGrid) && problem_().gridAdapt().wasAdapted())
        {
            uPrev_ = uCur_;
            prevVolVarsVector_ = volVarsVector_;

            updateBoundaryIndices_();

            int numDofs = asImp_().numDofs();

            if (isBox)
                boxVolume_.resize(numDofs);

            jacAsm_->init(problem_());
        }

    }


    /*!
     * \brief Called by the update() method if it was
     *        successful. This is primarily a hook which the actual
     *        model can overload.
     */
    void updateSuccessful()
    { }

    void updateVolVars()
    { volVarsVector_.update(problem_(), curSol()); }

    /*!
     * \brief Called by the update() method if it was
     *        unsuccessful. This is primarily a hook which the actual
     *        model can overload.
     */
    void updateFailed()
    {
        // Reset the current solution to the one of the
        // previous time step so that we can start the next
        // update at a physically meaningful solution.
        uCur_ = uPrev_;
        volVarsVector_ = prevVolVarsVector_;

        jacAsm_->reassembleAll();
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and
     *        the result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        // make the current solution the previous one.
        uPrev_ = uCur_;
        prevVolVarsVector_ = volVarsVector_;
    }

    /*!
     * \brief Serializes the current state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        if (isBox)
            res.template serializeEntities<dim>(asImp_(), this->gridView_());
        else
            res.template serializeEntities<0>(asImp_(), this->gridView_());
    }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        if (isBox)
            res.template deserializeEntities<dim>(asImp_(), this->gridView_());
        else
            res.template deserializeEntities<0>(asImp_(), this->gridView_());

        prevSol() = curSol();
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
        int dofIdxGlobal = dofMapper().index(entity);

        // write phase state
        if (!outstream.good()) {
            DUNE_THROW(Dune::IOError,
                       "Could not serialize vertex "
                       << dofIdxGlobal);
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            outstream << curSol()[dofIdxGlobal][eqIdx] << " ";
        }
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
        int dofIdxGlobal = dofMapper().index(entity);

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good())
                DUNE_THROW(Dune::IOError,
                           "Could not deserialize vertex "
                           << dofIdxGlobal);
            instream >> curSol()[dofIdxGlobal][eqIdx];
        }
    }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    size_t numDofs() const
    {
        if (isBox)
            return gridView_().size(dim);
        else
            return gridView_().size(0);
    }

    /*!
     * \brief Mapper for the entities where degrees of freedoms are
     *        defined to indices.
     *
     * Is the box method is used, this means a mapper
     * for vertices, if the cell centered method is used,
     * this means a mapper for elements.
     */
    template <class T = TypeTag>
    const typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox), VertexMapper>::type &dofMapper() const
    {
        return problem_().vertexMapper();
    }
    template <class T = TypeTag>
    const typename std::enable_if<!GET_PROP_VALUE(T, ImplicitIsBox), ElementMapper>::type &dofMapper() const
    {
        return problem_().elementMapper();
    }

    /*!
     * \brief Mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return problem_().vertexMapper(); }

    /*!
     * \brief Mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return problem_().elementMapper(); }

    /*!
     * \brief Resets the Jacobian matrix assembler, so that the
     *        boundary types can be altered.
     */
    void resetJacobianAssembler ()
    {
        jacAsm_.template reset<JacobianAssembler>(0);
        jacAsm_ = std::make_shared<JacobianAssembler>();
        jacAsm_->init(problem_());
    }

    /*!
     * \brief Update the weights of all primary variables within an
     *        element given the complete set of volume variables
     *
     * \param element The DUNE codim 0 entity
     * \param volVars All volume variables for the element
     */
    void updatePVWeights(const FVElementGeometry& fvGeometry) const
    { }

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
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;

        SolutionVector residual(u);
        asImp_().globalResidual(residual, u);

        // create the required scalar fields
        unsigned numDofs = asImp_().numDofs();

        // global defect of the two auxiliary equations
        ScalarField* def[numEq];
        ScalarField* delta[numEq];
        ScalarField* x[numEq];
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            x[eqIdx] = writer.allocateManagedBuffer(numDofs);
            delta[eqIdx] = writer.allocateManagedBuffer(numDofs);
            def[eqIdx] = writer.allocateManagedBuffer(numDofs);
        }

        for (unsigned int dofIdxGlobal = 0; dofIdxGlobal < u.size(); dofIdxGlobal++)
        {
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                (*x[eqIdx])[dofIdxGlobal] = u[dofIdxGlobal][eqIdx];
                (*delta[eqIdx])[dofIdxGlobal] = - deltaU[dofIdxGlobal][eqIdx];
                (*def[eqIdx])[dofIdxGlobal] = residual[dofIdxGlobal][eqIdx];
            }
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            std::ostringstream oss;
            oss.str(""); oss << "x_" << eqIdx;
            if (isBox)
                writer.attachVertexData(*x[eqIdx], oss.str());
            else
                writer.attachCellData(*x[eqIdx], oss.str());
            oss.str(""); oss << "delta_" << eqIdx;
            if (isBox)
                writer.attachVertexData(*delta[eqIdx], oss.str());
            else
                writer.attachCellData(*delta[eqIdx], oss.str());
            oss.str(""); oss << "defect_" << eqIdx;
            if (isBox)
                writer.attachVertexData(*def[eqIdx], oss.str());
            else
                writer.attachCellData(*def[eqIdx], oss.str());
        }

        asImp_().addOutputVtkFields(u, writer);
    }

    /*!
     * \brief Add the quantities of a time step which ought to be written to disk.
     *
     * This should be overwritten by the actual model if any secondary
     * variables should be written out. Read: This should _always_ be
     * overwritten by well behaved models!
     *
     * \tparam MultiWriter The type of the VTK multi writer
     *
     * \param sol The global vector of primary variable values.
     * \param writer The VTK multi writer where the fields should be added.
     */
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numDofs = asImp_().numDofs();

        // global defect of the two auxiliary equations
        ScalarField* x[numEq];
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            x[eqIdx] = writer.allocateManagedBuffer(numDofs);
        }

        for (int dofIdxGlobal = 0; dofIdxGlobal < sol.size(); dofIdxGlobal++)
        {
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                (*x[eqIdx])[dofIdxGlobal] = sol[dofIdxGlobal][eqIdx];
            }
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            std::ostringstream oss;
            oss << "primaryVar_" << eqIdx;
            if (isBox)
                writer.attachVertexData(*x[eqIdx], oss.str());
            else
                writer.attachCellData(*x[eqIdx], oss.str());
        }
    }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView() const
    { return problem_().gridView(); }

    /*!
     * \brief Returns true if the entity indicated by 'dofIdxGlobal'
     * is located on / touches the grid's boundary.
     *
     * \param dofIdxGlobal The global index of the entity
     */
    bool onBoundary(const int dofIdxGlobal) const
    { return boundaryIndices_[dofIdxGlobal]; }

    /*!
     * \brief Returns true if a vertex is located on the grid's
     *        boundary.
     *
     * \param element A DUNE Codim<0> entity which contains the control
     *             volume's associated vertex.
     * \param vIdx The local vertex index inside element
     */
    bool onBoundary(const SubControlVolume &scv) const
    {
        return asImp_().onBoundary(onBoundary(scv.dofIdxGlobal()));
    }

    /*!
     * \brief Fill the fluid state according to the primary variables.
     *
     * Taking the information from the primary variables,
     * the fluid state is filled with every information that is
     * necessary to evaluate the model's local residual.
     *
     * \param priVars The primary variables of the model.
     * \param problem The problem at hand.
     * \param element The current element.
     * \param fvGeometry The finite volume element geometry.
     * \param scvIdx The index of the subcontrol volume.
     * \param fluidState The fluid state to fill.
     */
    // template <class FluidState>
    // static void completeFluidState(const PrimaryVariables& priVars,
    //                                const Problem& problem,
    //                                const Element& element,
    //                                const FVElementGeometry& fvGeometry,
    //                                const int scvIdx,
    //                                FluidState& fluidState)
    // {
    //     VolumeVariables::completeFluidState(priVars, problem, element,
    //                                         fvGeometry, scvIdx, fluidState);
    // }

    const FluxVariables& fluxVars(const SubControlVolumeFace& scvf) const
    { return fluxVarsVector_[scvf.index()]; }

    const FluxVariables& fluxVars(unsigned int scvfIdx) const
    { return fluxVarsVector_[scvfIdx]; }

    const VolumeVariables& volVars(const SubControlVolume& scv) const
    { return volVarsVector_[scv.index()]; }

    const VolumeVariables& volVars(unsigned int scvIdx) const
    { return volVarsVector_[scvIdx]; }

    const VolumeVariables& prevVolVars(const SubControlVolume& scv) const
    { return prevVolVarsVector_[scv.index()]; }

    const VolumeVariables& prevVolVars(unsigned int scvIdx) const
    { return prevVolVarsVector_[scvIdx]; }

    const FVElementGeometryVector& fvGeometries() const
    { return fvGeometries_; }

    const FVElementGeometry& fvGeometries(const Element& element) const
    { return fvGeometries_.fvGeometry(elementMapper().index(element)); }

    const FVElementGeometry& fvGeometries(unsigned int eIdx) const
    { return fvGeometries_.fvGeometry(eIdx); }

protected:
    // only to be called by friends and deriving classes
    VolumeVariables& volVars_(const SubControlVolume& scv)
    { return volVarsVector_[scv.index()]; }

    VolumeVariables& volVars_(unsigned int scvIdx)
    { return volVarsVector_[scvIdx]; }

    VolumeVariables& prevVolVars_(const SubControlVolume& scv)
    { return prevVolVarsVector_[scv.index()]; }

    VolumeVariables& prevVolVars_(unsigned int scvIdx)
    { return prevVolVarsVector_[scvIdx]; }

    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    Problem &problem_()
    { return *problemPtr_; }
    /*!
     * \copydoc problem_()
     */
    const Problem &problem_() const
    { return *problemPtr_; }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Reference to the local residal object
     */
    LocalResidual &localResidual_()
    { return localJacobian_.localResidual(); }

    /*!
     * \brief Applies the initial solution for all vertices of the grid.
     *
     * \todo the initial condition needs to be unique for
     *       each vertex. we should think about the API...
     */
    void applyInitialSolution_()
    {
        // first set the whole domain to zero
        uCur_ = Scalar(0.0);
        boxVolume_ = Scalar(0.0);

        // iterate through leaf grid and evaluate initial
        // condition at the center of each sub control volume
        for (const auto& element : elements(gridView_()))
        {
            // deal with the current element
            for(auto&& scv : fvGeometries(element).scvs())
            {
                // let the problem do the dirty work of nailing down
                // the initial solution.
                auto initPriVars = problem_().initial(scv);

                auto dofIdxGlobal = scv.dofIndex();
                if (isBox)
                {
                    // add up the initial values of all sub-control
                    // volumes. If the initial values disagree for
                    // different sub control volumes, the initial value
                    // will be the arithmetic mean.
                    initPriVars *= scv.volume();
                    boxVolume_[dofIdxGlobal] += scv.volume();
                }

                uCur_[dofIdxGlobal] += initPriVars;
            }
        }

        // add up the primary variables and the volumes of the boxes
        // which cross process borders
        if (isBox && gridView_().comm().size() > 1) {
            VertexHandleSum<Dune::FieldVector<Scalar, 1>,
                Dune::BlockVector<Dune::FieldVector<Scalar, 1> >,
                VertexMapper> sumVolumeHandle(boxVolume_, vertexMapper());
            gridView_().communicate(sumVolumeHandle,
                                    Dune::InteriorBorder_InteriorBorder_Interface,
                                    Dune::ForwardCommunication);

            VertexHandleSum<PrimaryVariables, SolutionVector, VertexMapper>
                sumPVHandle(uCur_, vertexMapper());
            gridView_().communicate(sumPVHandle,
                                    Dune::InteriorBorder_InteriorBorder_Interface,
                                    Dune::ForwardCommunication);
        }

        if (isBox)
        {
            // divide all primary variables by the volume of their boxes
            for (unsigned int i = 0; i < uCur_.size(); ++i) {
                uCur_[i] /= boxVolume(i);
            }
        }
    }

    /*!
     * \brief Find all indices of boundary vertices (box) / elements (cell centered).
     */
    void updateBoundaryIndices_()
    {
        boundaryIndices_.resize(numDofs());
        std::fill(boundaryIndices_.begin(), boundaryIndices_.end(), false);

        for (const auto& element : elements(gridView_())) {
            Dune::GeometryType geomType = element.geometry().type();
            const ReferenceElement &refElement = ReferenceElements::general(geomType);

            for (const auto& intersection : intersections(gridView_(), element)) {
                if (intersection.boundary()) {
                    if (isBox)
                    {
                        // add all vertices on the intersection to the set of
                        // boundary vertices
                        int fIdx = intersection.indexInInside();
                        int numFaceVerts = refElement.size(fIdx, 1, dim);
                        for (int faceVertexIdx = 0;
                             faceVertexIdx < numFaceVerts;
                             ++faceVertexIdx)
                        {
                            int vIdx = refElement.subEntity(fIdx,
                                                            1,
                                                            faceVertexIdx,
                                                            dim);
                            int vIdxGlobal = vertexMapper().subIndex(element, vIdx, dim);
                            boundaryIndices_[vIdxGlobal] = true;
                        }
                    }
                    else
                    {
                        int eIdxGlobal = elementMapper().index(element);
                        boundaryIndices_[eIdxGlobal] = true;
                    }
                }
            }
        }
    }

    // the problem we want to solve. defines the constitutive
    // relations, matxerial laws, etc.
    Problem *problemPtr_;

    // calculates the local jacobian matrix for a given element
    LocalJacobian localJacobian_;
    // Linearizes the problem at the current time step using the
    // local jacobian
    std::shared_ptr<JacobianAssembler> jacAsm_;

    // the set of all indices of vertices on the boundary
    std::vector<bool> boundaryIndices_;

    // cur is the current iterative solution, prev the converged
    // solution of the previous time step
    SolutionVector uCur_;
    SolutionVector uPrev_;

    // the current and previous variables (primary and secondary variables)
    VolumeVariablesVector volVarsVector_;
    VolumeVariablesVector prevVolVarsVector_;

    // the flux variables vector
    FluxVariablesVector fluxVarsVector_;

    // the finite volume element geometries
    std::shared_ptr<FVElementGeometryVector> fvGeometries_;

    Dune::BlockVector<Dune::FieldVector<Scalar, 1> > boxVolume_;

private:
    /*!
     * \brief Returns whether messages should be printed
     */
    bool verbose_() const
    { return gridView_().comm().rank() == 0; }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

};

} // end namespace Dumux

#include "propertydefaults.hh"

#endif
