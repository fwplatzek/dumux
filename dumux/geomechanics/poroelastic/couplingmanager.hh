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
 * \ingroup MultiDomain
 * \ingroup PoroElastic
 * \ingroup PorousMediumFlow
 * \brief \copydoc Dumux::PMFlowFlowPoroMechanicsCouplingManager
 */

#ifndef DUMUX_POROMECHANICS_COUPLING_MANAGER_HH
#define DUMUX_POROMECHANICS_COUPLING_MANAGER_HH

#include <algorithm>
#include <type_traits>

#include <dumux/common/properties.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/common/timeloop.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup PoroElastic
 * \ingroup PorousMediumFlow
 * \brief Coupling manager for porous medium flow problems coupled to a poro-mechanical
 *        problem both defined on the same grid. Coupling occurs via the change of porosity
 *        and permeability due to mechanical deformations and the influence of the pore
 *        pressure on the effecive stresses acting on the porous medium.
 *
 * \tparam PMFlowId The porous medium flow domain id
 * \tparam PoroMechId The poro-mechanical domain id
 */
template< class MDTraits,
          std::size_t PMFlowId = 0,
          std::size_t PoroMechId = PMFlowId+1 >
class PoroMechanicsCouplingManager : public CouplingManager< MDTraits >
{
    using ParentType = CouplingManager< MDTraits >;

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;

    // further types specific to the sub-problems
    template<std::size_t id> using Scalar = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Scalar);
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);
    template<std::size_t id> using LocalResidual = typename GET_PROP_TYPE(SubDomainTypeTag<id>, LocalResidual);
    template<std::size_t id> using GridVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVariables);
    template<std::size_t id> using PrimaryVariables = typename GridVariables<id>::PrimaryVariables;
    template<std::size_t id> using GridVolumeVariables = typename GridVariables<id>::GridVolumeVariables;
    template<std::size_t id> using ElementVolumeVariables = typename GridVolumeVariables<id>::LocalView;
    template<std::size_t id> using FVGridGeometry = typename GridVariables<id>::GridGeometry;
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using IndexType = typename GridView<id>::IndexSet::IndexType;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;
    template<std::size_t id> using GlobalPosition = typename Element<id>::Geometry::GlobalCoordinate;

    //! we assume that the two sub-problem operate on the same grid
    static_assert(std::is_same< GridView<PMFlowId>, GridView<PoroMechId> >::value,
                  "The grid types of the two sub-problems have to be equal!");

    //! this coupling manager is for cc - box only
    static_assert(FVGridGeometry<PoroMechId>::discMethod == DiscretizationMethod::box,
                  "Poro-mechanical problem must be discretized with the box scheme for this coupling manager!");

    static_assert(FVGridGeometry<PMFlowId>::discMethod == DiscretizationMethod::cctpfa ||
                  FVGridGeometry<PMFlowId>::discMethod == DiscretizationMethod::ccmpfa,
                  "Porous medium flow problem must be discretized with a cell-centered scheme for this coupling manager!");

    //! this does not work for enabled grid volume variables caching (update of local view in context has no effect)
    static_assert(!GET_PROP_VALUE(SubDomainTypeTag<PMFlowId>, EnableGridVolumeVariablesCache),
                  "Poromechanics framework does not yet work for enabled grid volume variables caching");

    //! Types used for coupling stencils
    template<std::size_t id>
    using CouplingIndexType = typename std::conditional< id == PMFlowId,
                                                         IndexType<PoroMechId>,
                                                         IndexType<PMFlowId> >::type;

    /*!
     * \brief Porous medium flow domain data required for the residual calculation of an
     *        element of the poro-mechanidal domain. We store the data required to do an
     *        update of all primary/secondary variables at the dof of the element.
     */
    struct PoroMechanicsCouplingContext
    {
        // We need unique ptrs because the local views have no default constructor
        Element<PMFlowId> pmFlowElement;
        std::unique_ptr< FVElementGeometry<PMFlowId> > pmFlowFvGeometry;
        std::unique_ptr< ElementVolumeVariables<PMFlowId> > pmFlowElemVolVars;
    };

public:

    /*!
     * \brief The types used for coupling stencils. An element of the poro-mechanical
     *        domain always only couples to the single dof (because we use a cell-centered
     *        scheme in the porous medium flow domain) of this same element.
     */
    template<std::size_t i, std::size_t j = (i == PMFlowId) ? PoroMechId : PMFlowId>
    using CouplingStencilType = typename std::conditional< i == PMFlowId,
                                                           std::vector< CouplingIndexType<i> >,
                                                           std::array< CouplingIndexType<i>, 1> >::type;

    //! the type of the solution vector
    using SolutionVector = typename MDTraits::SolutionVector;

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param pmFlowProblem The porous medium flow problem
     * \param poroMechanicalProblem The poro-mechanical problem
     * \param curSol The current solution
     */
    void init(std::shared_ptr< Problem<PMFlowId> > pmFlowProblem,
              std::shared_ptr< Problem<PoroMechId> > poroMechanicalProblem,
              const SolutionVector& curSol)
    {
        // set up tuple containing the sub-problems
        problemTuple_ = std::make_tuple(pmFlowProblem, poroMechanicalProblem);
        // copy the solution vector
        ParentType::updateSolution(curSol);
        // set up the coupling map pmfow -> poromechanics
        initializeCouplingMap_();
    }

    /*!
     * \brief Return the coupling stencil for a given porous medium flow domain element
     */
    const CouplingStencilType<PMFlowId>& couplingStencil(Dune::index_constant<PMFlowId> pmFlowDomainId,
                                                         const Element<PMFlowId>& element,
                                                         Dune::index_constant<PoroMechId> poroMechDomainId) const
    {
        return pmFlowCouplingMap_[ problem<PMFlowId>().fvGridGeometry().elementMapper().index(element) ];
    }

    /*!
     * \brief Return the coupling element stencil for a given poro-mechanical domain element
     */
    const CouplingStencilType<PoroMechId> couplingStencil(Dune::index_constant<PoroMechId> poroMechDomainId,
                                                          const Element<PoroMechId>& element,
                                                          Dune::index_constant<PMFlowId> pmFlowDomainId) const
    {
        const auto eIdx = problem<PMFlowId>().fvGridGeometry().elementMapper().index(element);
        return CouplingStencilType<PoroMechId>{ {eIdx} };
    }

    //! Pull up the base class' default implementations
    using ParentType::bindCouplingContext;
    using ParentType::updateCoupledVariables;

    /*!
     * \brief For the assembly of the element residual of an element of the poro-mechanics domain,
     *        we have to prepare the element variables of the porous medium flow domain.
     */
    template< class Assembler >
    void bindCouplingContext(Dune::index_constant<PoroMechId> poroMechDomainId,
                             const Element<PoroMechId>& element,
                             const Assembler& assembler)
    {
        // first reset the context
        poroMechCouplingContext_.pmFlowFvGeometry.reset(nullptr);
        poroMechCouplingContext_.pmFlowElemVolVars.reset(nullptr);

        // prepare the fvGeometry and the element volume variables
        // these quantities will be used later to obtain the effective pressure
        auto fvGeometry = localView( problem<PMFlowId>().fvGridGeometry() );
        auto elemVolVars = localView( assembler.gridVariables(Dune::index_constant<PMFlowId>()).curGridVolVars() );

        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, this->curSol()[Dune::index_constant<PMFlowId>()]);
        poroMechCouplingContext_.pmFlowElement = element;
        poroMechCouplingContext_.pmFlowFvGeometry = std::make_unique< FVElementGeometry<PMFlowId> >(fvGeometry);
        poroMechCouplingContext_.pmFlowElemVolVars = std::make_unique< ElementVolumeVariables<PMFlowId> >(elemVolVars);
    }

    /*!
     * \brief After deflection of the solution in the porous medium flow domain during
     *        element residual assembly in the poromechanics domain, we have to update
     *        the porous medium flow element variables of the coupling context
     */
    template< class PoroMechLocalAssembler >
    void updateCouplingContext(Dune::index_constant<PoroMechId> poroMechDomainId,
                               const PoroMechLocalAssembler& poroMechLocalAssembler,
                               Dune::index_constant<PMFlowId> pmFlowDomainId,
                               IndexType<PMFlowId> dofIdxGlobalJ,
                               const PrimaryVariables<PMFlowId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate the deflected pm flow domain primary variable
        ParentType::updateCouplingContext(poroMechDomainId, poroMechLocalAssembler, pmFlowDomainId, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // now, update the coupling context (i.e. elemVolVars)
        const auto& element = poroMechCouplingContext_.pmFlowElement;
        const auto& fvGeometry = *poroMechCouplingContext_.pmFlowFvGeometry;
        poroMechCouplingContext_.pmFlowElemVolVars->bindElement(element, fvGeometry, this->curSol()[pmFlowDomainId]);
    }

    /*!
     * \brief After deflection of the solution in the poromechanics domain during
     *        element residual assembly in that same domain, we have to update the
     *        porous medium flow element variables of the coupling context because
     *        the porosity/permeability might depend on the mechanical deformation
     */
    template< class PoroMechLocalAssembler >
    void updateCouplingContext(Dune::index_constant<PoroMechId> poroMechDomainIdI,
                               const PoroMechLocalAssembler& poroMechLocalAssembler,
                               Dune::index_constant<PoroMechId> poroMechDomainIdJ,
                               IndexType<PoroMechId> dofIdxGlobalJ,
                               const PrimaryVariables<PoroMechId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate the deflected displacement
        ParentType::updateCouplingContext(poroMechDomainIdI, poroMechLocalAssembler, poroMechDomainIdJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);

        // now, update the coupling context (i.e. elemVolVars)
        (*poroMechCouplingContext_.pmFlowElemVolVars).bindElement(poroMechCouplingContext_.pmFlowElement,
                                                                  *poroMechCouplingContext_.pmFlowFvGeometry,
                                                                  this->curSol()[Dune::index_constant<PMFlowId>()]);
    }

    /*!
     * \brief We need this overload to avoid ambiguity. However, for the porous medium flow
     *        domain weonly have to update the solution, which is done in the parent class.
     */
    template< std::size_t j, class PMFlowLocalAssembler >
    void updateCouplingContext(Dune::index_constant<PMFlowId> pmFlowDomainId,
                               const PMFlowLocalAssembler& pmFlowLocalAssembler,
                               Dune::index_constant<j> domainIdJ,
                               IndexType<j> dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        // communicate the deflected displacement
        ParentType::updateCouplingContext(pmFlowDomainId, pmFlowLocalAssembler, domainIdJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief Update the porous medium flow domain volume variables and flux variables cache
     *        after the coupling context has been updated. This has to be done because the
     *        mechanical deformation enters the porosity/permeability relationships.
     */
    template< class PMFlowLocalAssembler, class UpdatableElementVolVars, class UpdatableFluxVarCache >
    void updateCoupledVariables(Dune::index_constant<PMFlowId> pmFlowDomainId,
                                const PMFlowLocalAssembler& pmFlowLocalAssembler,
                                ElementVolumeVariables<PMFlowId>& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        // update the element volume variables to obtain the updated porosity/permeability
        elemVolVars.bind(pmFlowLocalAssembler.element(),
                         pmFlowLocalAssembler.fvGeometry(),
                         this->curSol()[pmFlowDomainId]);

        // update the transmissibilities subject to the new permeabilities
        elemFluxVarsCache.update(pmFlowLocalAssembler.element(),
                                 pmFlowLocalAssembler.fvGeometry(),
                                 elemVolVars);
    }

    /*!
     * \brief Update the poro-mechanics volume variables after the coupling context has been updated.
     *        This is necessary because the fluid density is stored in them and which potentially is
     *        solution-dependent.
     */
    template< class PoroMechLocalAssembler, class UpdatableElementVolVars, class UpdatableFluxVarCache >
    void updateCoupledVariables(Dune::index_constant<PoroMechId> poroMechDomainId,
                                const PoroMechLocalAssembler& poroMechLocalAssembler,
                                ElementVolumeVariables<PoroMechId>& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        elemVolVars.bind(poroMechLocalAssembler.element(),
                         poroMechLocalAssembler.fvGeometry(),
                         this->curSol()[poroMechDomainId]);
    }

    /*!
     * \brief Compute the divergence of u for a given element and position. This is possibly called
     *        from the user-defined spatial parameters interface for the porosity and/or permeability.
     */
    Scalar<PoroMechId> computeDivU(const Element<PoroMechId>& element, const GlobalPosition<PoroMechId>& globalPos) const
    {
        const auto& poroMechGridGeom = problem<PoroMechId>().fvGridGeometry();
        const auto poroMechElemSol = elementSolution(element, this->curSol()[Dune::index_constant<PoroMechId>()], poroMechGridGeom);

        Scalar<PoroMechId> divU = 0.0;
        const auto gradU = evalGradients(element, element.geometry(), poroMechGridGeom, poroMechElemSol, globalPos);
        for (int i = 0; i < GridView<PoroMechId>::dimension; ++i)
            divU += gradU[i][i];
        return divU;
    }

    /*!
     * \brief Evaluates the coupling element residual of the porous medium flow domain with
     *        respect to the poro-mechanical domain. The deformation might has an effect on
     *        both the permeability as well as the porosity. Thus, we need to compute fluxes
     *        and the storage term here.
     */
    template< class PMFlowLocalAssembler >
    typename LocalResidual<PMFlowId>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<PMFlowId> pmFlowDomainId,
                         const PMFlowLocalAssembler& pmFlowLocalAssembler,
                         Dune::index_constant<PoroMechId> poroMechDomainId,
                         IndexType<PoroMechId> dofIdxGlobalJ)
    {
        auto res = pmFlowLocalAssembler.localResidual().evalFluxAndSource(pmFlowLocalAssembler.element(),
                                                                          pmFlowLocalAssembler.fvGeometry(),
                                                                          pmFlowLocalAssembler.curElemVolVars(),
                                                                          pmFlowLocalAssembler.elemFluxVarsCache(),
                                                                          pmFlowLocalAssembler.elemBcTypes());

        // If the residual instationary, evaluate storage
        if (!pmFlowLocalAssembler.localResidual().isStationary())
            res += pmFlowLocalAssembler.localResidual().evalStorage(pmFlowLocalAssembler.element(),
                                                                    pmFlowLocalAssembler.fvGeometry(),
                                                                    pmFlowLocalAssembler.prevElemVolVars(),
                                                                    pmFlowLocalAssembler.curElemVolVars());

        return res;
    }

    /*!
     * \brief Evaluates the coupling element residual of the poromechanical domain with
     *        respect to the porous medium flow domain. The pressure has an effect on the
     *        mechanical stresses as well as the body forces. Thus, we have to compute
     *        the fluxes as well as the source term here.
     */
    template< class PoroMechLocalAssembler >
    typename LocalResidual<PoroMechId>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<PoroMechId> poroMechDomainId,
                         const PoroMechLocalAssembler& pmFlowLocalAssembler,
                         Dune::index_constant<PMFlowId> pmFlowDomainId,
                         IndexType<PMFlowId> dofIdxGlobalJ)
    {
        return pmFlowLocalAssembler.localResidual().evalFluxAndSource(pmFlowLocalAssembler.element(),
                                                                      pmFlowLocalAssembler.fvGeometry(),
                                                                      pmFlowLocalAssembler.curElemVolVars(),
                                                                      pmFlowLocalAssembler.elemFluxVarsCache(),
                                                                      pmFlowLocalAssembler.elemBcTypes());
    }

    //! Return a const reference to one of the problems
    template<std::size_t id> const Problem<id>& problem() const
    { return *std::get<id>(problemTuple_); }

    //! Return reference to one of the problems
    template<std::size_t id> Problem<id>& problem()
    { return *std::get<id>(problemTuple_); }

    //! Return the coupling context (used in mechanical sub-problem to compute effective pressure)
    const PoroMechanicsCouplingContext& poroMechanicsCouplingContext() const
    { return poroMechCouplingContext_; }

private:
    /*!
     * \brief Initializes the pm flow domain coupling map. Since the elements
     *        of the poro-mechanical domain only couple to the same elements, we
     *        don't set up the maps here but return copy of the "stencil" always.
     */
    void initializeCouplingMap_()
    {
        // some references for convenience
        const auto& pmFlowGridGeom = problem<PMFlowId>().fvGridGeometry();
        const auto& poroMechGridGeom = problem<PoroMechId>().fvGridGeometry();

        // make sure the two grids are really the same. Note that if the two grids
        // happen to have equal number of elements by chance, we don't detect this source of error.
        if (pmFlowGridGeom.gridView().size(0) != poroMechGridGeom.gridView().size(0))
            DUNE_THROW(Dune::InvalidStateException, "The two sub-problems are assumed to operate on the same mesh!");

        pmFlowCouplingMap_.resize(pmFlowGridGeom.gridView().size(0));
        static constexpr int dim = GridView<PMFlowId>::dimension;
        for (const auto& element : elements(pmFlowGridGeom.gridView()))
        {
            const auto eIdx = pmFlowGridGeom.elementMapper().index(element);

            // firstly, the element couples to the nodal dofs in itself
            for (int i = 0; i < element.geometry().corners(); ++i)
                pmFlowCouplingMap_[eIdx].push_back( pmFlowGridGeom.vertexMapper().subIndex(element, i , dim) );

            // the pm flow problem couples to the same elements as in its own stencil
            // due to the dependency of the residual on all permeabilities in its stencil,
            // which in turn depend on the mechanical deformations.
            const auto& inverseConnectivity = pmFlowGridGeom.connectivityMap()[eIdx];
            for (const auto& dataJ : inverseConnectivity)
            {
                const auto elemJ = pmFlowGridGeom.element(dataJ.globalJ);
                for (int i = 0; i < elemJ.geometry().corners(); ++i)
                    pmFlowCouplingMap_[dataJ.globalJ].push_back( pmFlowGridGeom.vertexMapper().subIndex(elemJ, i , dim) );
            }
        }

        // make stencils unique
        for (auto& stencil : pmFlowCouplingMap_)
        {
            std::sort(stencil.begin(), stencil.end());
            stencil.erase(std::unique(stencil.begin(), stencil.end()), stencil.end());
        }
    }

    // tuple for storing pointers to the sub-problems
    using PMFlowProblemPtr = std::shared_ptr< Problem<PMFlowId> >;
    using PoroMechanicalProblemPtr = std::shared_ptr< Problem<PoroMechId> >;
    std::tuple<PMFlowProblemPtr, PoroMechanicalProblemPtr> problemTuple_;

    // Container for storing the coupling element stencils for the pm flow domain
    std::vector< CouplingStencilType<PMFlowId> > pmFlowCouplingMap_;

    // the coupling context of the poromechanics domain
    PoroMechanicsCouplingContext poroMechCouplingContext_;
};

} //end namespace Dumux

#endif