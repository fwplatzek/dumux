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
 * \ingroup EmbeddedCoupling
 * \brief Coupling manager for low-dimensional domains embedded in the bulk
 *        domain. Point sources on each integration point are computed by an AABB tree.
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGERBASE_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_COUPLINGMANAGERBASE_HH

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <unordered_map>

#include <dune/common/timer.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/geometry/intersectingentities.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/embedded/mixeddimensionglue.hh>
#include <dumux/multidomain/embedded/pointsourcedata.hh>
#include <dumux/multidomain/embedded/integrationpointsource.hh>

namespace Dumux {

//! the default point source traits
template<class MDTraits>
struct DefaultPointSourceTraits
{
private:
    template<std::size_t i> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<i>;
    template<std::size_t i> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<i>, FVGridGeometry);
    template<std::size_t i> using NumEqVector = typename GET_PROP_TYPE(SubDomainTypeTag<i>, NumEqVector);
public:
    //! export the point source type for domain i
    template<std::size_t i>
    using PointSource = IntegrationPointSource<typename FVGridGeometry<i>::GlobalCoordinate, NumEqVector<i>>;

    //! export the point source helper type  for domain i
    template<std::size_t i>
    using PointSourceHelper = IntegrationPointSourceHelper;

    //! export the point source data type
    using PointSourceData = Dumux::PointSourceData<MDTraits>;
};

/*!
 * \ingroup MultiDomain
 * \ingroup EmbeddedCoupling
 * \brief Manages the coupling between bulk elements and lower dimensional elements
 *        Point sources on each integration point are computed by an AABB tree.
 */
template<class MDTraits, class Implementation, class PSTraits = DefaultPointSourceTraits<MDTraits>>
class EmbeddedCouplingManagerBase
: public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto bulkIdx = typename MDTraits::template DomainIdx<0>();
    static constexpr auto lowDimIdx = typename MDTraits::template DomainIdx<1>();
    using SolutionVector = typename MDTraits::SolutionVector;
    using PointSourceData = typename PSTraits::PointSourceData;

    // the sub domain type tags
    template<std::size_t id> using PointSource = typename PSTraits::template PointSource<id>;
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;
    template<std::size_t id> using GridView = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridView);
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);
    template<std::size_t id> using PrimaryVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, PrimaryVariables);
    template<std::size_t id> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry);
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    enum {
        bulkDim = GridView<bulkIdx>::dimension,
        lowDimDim = GridView<lowDimIdx>::dimension,
        dimWorld = GridView<bulkIdx>::dimensionworld
    };

    template<std::size_t id>
    static constexpr bool isBox()
    { return FVGridGeometry<id>::discMethod == DiscretizationMethod::box; }

    using CouplingStencil = std::vector<std::size_t>;
    using GlobalPosition = typename Element<bulkIdx>::Geometry::GlobalCoordinate;
    using GlueType = MixedDimensionGlue<GridView<bulkIdx>, GridView<lowDimIdx>>;

public:
    //! export the point source traits
    using PointSourceTraits = PSTraits;
    //! export stencil types
    using CouplingStencils = std::unordered_map<std::size_t, CouplingStencil>;

    /*!
    * \brief call this after grid adaption
    */
    void updateGrid(std::shared_ptr<const FVGridGeometry<bulkIdx>> bulkFvGridGeometry,
                    std::shared_ptr<const FVGridGeometry<lowDimIdx>> lowDimFvGridGeometry)
    {
        glue_ = std::make_shared<GlueType>();
    }

    /*!
     * \brief Constructor
     */
    EmbeddedCouplingManagerBase(std::shared_ptr<const FVGridGeometry<bulkIdx>> bulkFvGridGeometry,
                                std::shared_ptr<const FVGridGeometry<lowDimIdx>> lowDimFvGridGeometry)
    {
        updateGrid(bulkFvGridGeometry, lowDimFvGridGeometry);
    }

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    void init(std::shared_ptr<Problem<bulkIdx>> bulkProblem,
              std::shared_ptr<Problem<lowDimIdx>> lowDimProblem,
              const SolutionVector& curSol)
    {
        this->updateSolution(curSol);
        problemTuple_ = std::make_tuple(bulkProblem, lowDimProblem);

        integrationOrder_ = getParam<int>("MixedDimension.IntegrationOrder", 1);
        asImp_().computePointSourceData(integrationOrder_);
    }

    // \}

    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param elementI the coupled element of domain í
     * \param domainJ the domain index of domain j
     *
     * \note  The element residual definition depends on the discretization scheme of domain i
     *        box: a container of the residuals of all sub control volumes
     *        cc : the residual of the (sub) control volume
     *        fem: the residual of the element
     * \note  This function has to be implemented by all coupling managers for all combinations of i and j
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const Element<i>& element,
                                           Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");

        const auto eIdx = problem(domainI).fvGridGeometry().elementMapper().index(element);
        if (couplingStencils(domainI).count(eIdx))
            return couplingStencils(domainI).at(eIdx);
        else
            return emptyStencil_;
    }

    /*!
     * \brief evaluates the element residual of a coupled element of domain i which depends on the variables
     *        at the degree of freedom with index dofIdxGlobalJ of domain j
     *
     * \param domainI the domain index of domain i
     * \param localAssemblerI the local assembler assembling the element residual of an element of domain i
     * \param domainJ the domain index of domain j
     * \param dofIdxGlobalJ the index of the degree of freedom of domain j which has an influence on the element residual of domain i
     *
     * \note  we only need to evaluate the source contribution to the residual here as the coupling term is the source
     * \return the element residual
     */
    template<std::size_t i, std::size_t j, class LocalAssemblerI>
    decltype(auto) evalCouplingResidual(Dune::index_constant<i> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ)
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");

        typename LocalAssemblerI::LocalResidual::ElementResidualVector residual;

        const auto& element = localAssemblerI.element();
        const auto& fvGeometry = localAssemblerI.fvGeometry();
        const auto& curElemVolVars = localAssemblerI.curElemVolVars();

        residual.resize(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
        {
            auto couplingSource = problem(domainI).scvPointSources(element, fvGeometry, curElemVolVars, scv);
            couplingSource *= -scv.volume()*curElemVolVars[scv].extrusionFactor();
            residual[scv.indexInElement()] = couplingSource;
        }
        return residual;
    }

    // \}

    /* \brief Compute integration point point sources and associated data
     *
     * This method uses grid glue to intersect the given grids. Over each intersection
     * we later need to integrate a source term. This method places point sources
     * at each quadrature point and provides the point source with the necessary
     * information to compute integrals (quadrature weight and integration element)
     * \param order The order of the quadrature rule for integration of sources over an intersection
     * \param verbose If the point source computation is verbose
     */
    void computePointSourceData(std::size_t order = 1, bool verbose = false)
    {
        static_assert(FVGridGeometry<bulkIdx>::discMethod == DiscretizationMethod::cctpfa, "Currently only works for tpfa discretization!");

        // initilize the maps
        // do some logging and profiling
        Dune::Timer watch;
        std::cout << "Initializing the point sources..." << std::endl;

        // clear all internal members like pointsource vectors and stencils
        // initializes the point source id counter
        clear();

        const auto& bulkFvGridGeometry = problem(bulkIdx).fvGridGeometry();
        const auto& lowDimFvGridGeometry = problem(lowDimIdx).fvGridGeometry();

        // intersect the bounding box trees
        glue_->build(bulkFvGridGeometry.boundingBoxTree(),
                     lowDimFvGridGeometry.boundingBoxTree());

        pointSourceData_.reserve(glue_->size());
        averageDistanceToBulkCC_.reserve(glue_->size());
        for (const auto& is : intersections(*glue_))
        {
            // all inside elements are identical...
            const auto& inside = is.inside(0);
            // get the intersection geometry for integrating over it
            const auto intersectionGeometry = is.geometry();

            // get the Gaussian quadrature rule for the local intersection
            const auto& quad = Dune::QuadratureRules<Scalar, lowDimDim>::rule(intersectionGeometry.type(), order);
            const std::size_t lowDimElementIdx = lowDimFvGridGeometry.elementMapper().index(inside);

            // iterate over all quadrature points
            for (auto&& qp : quad)
            {
                // compute the coupling stencils
                for (std::size_t outsideIdx = 0; outsideIdx < is.neighbor(0); ++outsideIdx)
                {
                    const auto& outside = is.outside(outsideIdx);
                    const std::size_t bulkElementIdx = bulkFvGridGeometry.elementMapper().index(outside);

                    // each quadrature point will be a point source for the sub problem
                    const auto globalPos = intersectionGeometry.global(qp.position());
                    const auto id = idCounter_++;
                    const auto qpweight = qp.weight();
                    const auto ie = intersectionGeometry.integrationElement(qp.position());
                    pointSources(bulkIdx).emplace_back(globalPos, id, qpweight, ie, std::vector<std::size_t>({bulkElementIdx}));
                    pointSources(bulkIdx).back().setEmbeddings(is.neighbor(0));
                    pointSources(lowDimIdx).emplace_back(globalPos, id, qpweight, ie, std::vector<std::size_t>({lowDimElementIdx}));
                    pointSources(lowDimIdx).back().setEmbeddings(is.neighbor(0));

                    // pre compute additional data used for the evaluation of
                    // the actual solution dependent source term
                    PointSourceData psData;
                    psData.addLowDimInterpolation(lowDimElementIdx);
                    psData.addBulkInterpolation(bulkElementIdx);

                    // publish point source data in the global vector
                    pointSourceData_.emplace_back(std::move(psData));
                    averageDistanceToBulkCC_.push_back(computeDistance_(outside.geometry(), globalPos));

                    // compute the coupling stencils
                    bulkCouplingStencils_[bulkElementIdx].push_back(lowDimElementIdx);
                    // add this bulk element to the low dim coupling stencil
                    lowDimCouplingStencils_[lowDimElementIdx].push_back(bulkElementIdx);
                }
            }
        }

        // make stencils unique
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(problemTuple_)), [&](const auto domainIdx)
        {
            for (auto&& stencil : this->couplingStencils(domainIdx))
            {
                std::sort(stencil.second.begin(), stencil.second.end());
                stencil.second.erase(std::unique(stencil.second.begin(), stencil.second.end()), stencil.second.end());
            }
        });

        std::cout << "took " << watch.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Methods to be accessed by the subproblems
     */
    // \{

    //! Return a reference to the pointSource data
    const PointSourceData& pointSourceData(std::size_t id) const
    { return pointSourceData_[id]; }

    //! Return a reference to the bulk problem
    template<std::size_t id>
    const Problem<id>& problem(Dune::index_constant<id> domainIdx) const
    { return *std::get<id>(problemTuple_); }

    //! Return a reference to the bulk problem
    template<std::size_t id>
    const GridView<id>& gridView(Dune::index_constant<id> domainIdx) const
    { return problem(domainIdx).fvGridGeometry().gridView(); }

    //! Return data for a bulk point source with the identifier id
    PrimaryVariables<bulkIdx> bulkPriVars(std::size_t id) const
    {
        auto& data = pointSourceData_[id];
        return data.interpolateBulk(this->curSol()[bulkIdx]);
    }

    //! Return data for a low dim point source with the identifier id
    PrimaryVariables<lowDimIdx> lowDimPriVars(std::size_t id) const
    {
        auto& data = pointSourceData_[id];
        return data.interpolateLowDim(this->curSol()[lowDimIdx]);
    }

    //! return the average distance to the coupled bulk cell center
    Scalar averageDistance(std::size_t id) const
    { return averageDistanceToBulkCC_[id]; }

    //! Return reference to bulk point sources
    const std::vector<PointSource<bulkIdx>>& bulkPointSources() const
    { return std::get<bulkIdx>(pointSources_); }

    //! Return reference to low dim point sources
    const std::vector<PointSource<lowDimIdx>>& lowDimPointSources() const
    { return std::get<lowDimIdx>(pointSources_); }

    //! Return the point source if domain i
    template<std::size_t i>
    const std::vector<PointSource<i>>& pointSources(Dune::index_constant<i> dom) const
    { return std::get<i>(pointSources_); }

    //! Return reference to bulk coupling stencil member of domain i
    template<std::size_t i>
    const CouplingStencils& couplingStencils(Dune::index_constant<i> dom) const
    { return (i == bulkIdx) ? bulkCouplingStencils_ : lowDimCouplingStencils_; }

    //! Return reference to point source data vector member
    const std::vector<PointSourceData>& pointSourceData() const
    { return pointSourceData_; }

protected:

    //! computes the bulk vertex indices per element for the box method
    void preComputeBulkVertexIndices()
    {
        // fill helper structure for box discretization
        if (isBox<bulkIdx>())
        {
            bulkVertexIndices_.resize(gridView(bulkIdx).size(0));
            for (const auto& element : elements(gridView(bulkIdx)))
            {
                const auto eIdx = problem(bulkIdx).fvGridGeometry().elementMapper().index(element);
                bulkVertexIndices_[eIdx].resize(element.subEntities(bulkDim));
                for (int i = 0; i < element.subEntities(bulkDim); ++i)
                    bulkVertexIndices_[eIdx][i] = problem(bulkIdx).fvGridGeometry().vertexMapper().subIndex(element, i, bulkDim);
            }
        }
    }

    //! Clear all internal data members
    void clear()
    {
        pointSources(bulkIdx).clear();
        pointSources(lowDimIdx).clear();
        pointSourceData_.clear();
        bulkCouplingStencils_.clear();
        lowDimCouplingStencils_.clear();

        idCounter_ = 0;
    }

    //! Return a reference to the bulk problem
    template<std::size_t id>
    Problem<id>& problem(Dune::index_constant<id> domainIdx)
    { return *std::get<id>(problemTuple_); }

    //! Return reference to point source data vector member
    std::vector<PointSourceData>& pointSourceData()
    { return pointSourceData_; }

    //! Return the point source if domain i
    template<std::size_t i>
    std::vector<PointSource<i>>& pointSources(Dune::index_constant<i> dom)
    { return std::get<i>(pointSources_); }

    //! Return reference to bulk coupling stencil member of domain i
    template<std::size_t i>
    CouplingStencils& couplingStencils(Dune::index_constant<i> dom)
    { return (i == 0) ? bulkCouplingStencils_ : lowDimCouplingStencils_; }

    //! Return a reference to the bulk vertex indices
    const std::vector<std::size_t>& bulkVertexIndices(std::size_t eIdx)
    { return bulkVertexIndices_[eIdx]; }

    //! Return a reference to an empty stencil
    const CouplingStencil& emptyStencil() const
    { return emptyStencil_; }

    const GlueType& glue() const
    { return *glue_; }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    //! id generator for point sources
    std::size_t idCounter_ = 0;

private:

    template<class Geometry, class GlobalPosition>
    Scalar computeDistance_(const Geometry& geometry, const GlobalPosition& p)
    {
        Scalar avgDist = 0.0;
        const auto& quad = Dune::QuadratureRules<Scalar, bulkDim>::rule(geometry.type(), 5);
        for (auto&& qp : quad)
        {
            const auto globalPos = geometry.global(qp.position());
            avgDist += (globalPos-p).two_norm()*qp.weight();
        }
        return avgDist;
    }

    std::tuple<std::shared_ptr<Problem<0>>, std::shared_ptr<Problem<1>>> problemTuple_;

    //! the point source in both domains
    std::tuple<std::vector<PointSource<bulkIdx>>, std::vector<PointSource<lowDimIdx>>> pointSources_;
    mutable std::vector<PointSourceData> pointSourceData_;
    std::vector<Scalar> averageDistanceToBulkCC_;

    //! Stencil data
    std::vector<std::vector<std::size_t>> bulkVertexIndices_;
    CouplingStencils bulkCouplingStencils_;
    CouplingStencils lowDimCouplingStencils_;
    CouplingStencil emptyStencil_;

    //! The glue object
    std::shared_ptr<GlueType> glue_;

    //! integration order for coupling source
    int integrationOrder_ = 1;
};

} // end namespace Dumux

#endif