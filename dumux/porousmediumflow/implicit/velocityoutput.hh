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
 * \brief Velocity output for implicit (porous media) models
 */
#ifndef DUMUX_IMPLICIT_VELOCITYOUTPUT_HH
#define DUMUX_IMPLICIT_VELOCITYOUTPUT_HH

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

//At the moment this property is defined in the individual models -> should be changed
namespace Properties
{
    NEW_PROP_TAG(VtkAddVelocity); //!< Returns whether velocity vectors are written into the vtk output
}

/*!
 * \brief Velocity output for implicit (porous media) models
 */
template<class TypeTag>
class ImplicitVelocityOutput
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    static const bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using CoordScalar = typename GridView::ctype;
    using Stencil = std::vector<IndexType>;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using ReferenceElements = Dune::ReferenceElements<CoordScalar, dim>;

public:
    /*!
     * \brief Constructor initializes the static data with the initial solution.
     *
     * \param problem The problem to be solved
     */
    ImplicitVelocityOutput(const Problem& problem)
    : problem_(problem)
    {
        // check, if velocity output can be used (works only for cubes so far)
        velocityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddVelocity);
        if (velocityOutput_)
        {
            // set the number of scvs the vertices are connected to
            if (isBox)
            {
                // resize to the number of vertices of the grid
                cellNum_.assign(problem.gridView().size(dim), 0);

                for (const auto& vertex : vertices(problem.gridView()))
                    cellNum_[problem.vertexMapper().index(vertex)] = getStencil(vertex).size();
            }
        }
    }

    bool enableOutput()
    { return velocityOutput_; }

    // The following SFINAE enable_if usage allows compilation, even if only a
    //
    // boundaryTypes(const Element&, const scv&)
    //
    // is provided in the problem file. In that case, the compiler cannot detect
    // (without additional measures like "using...") the signature
    //
    // boundaryTypes(const Element&, const scvf&)
    //
    // in the problem base class. Therefore, calls to this method trigger a
    // compiler error. However, that call is needed for calculating velocities
    // if the cell-centered discretization is used. By proceeding as in the
    // following lines, that call will only be compiled if cell-centered
    // actually is used. For the same reason we also provide a isBox-specific
    // implementation of the getStencil method below.
    template <class T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, ImplicitIsBox), BoundaryTypes>::type
    problemBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    { return problem_.boundaryTypes(element, scvf); }

    //! we should never call this method for box models
    template <class T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox), BoundaryTypes>::type
    problemBoundaryTypes(const Element& element, const SubControlVolume& scv) const
    { return BoundaryTypes(); }

    //! returns the elements connected to a vertex
    template<class T = TypeTag>
    const typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox), Stencil>::type&
    getStencil(const Vertex& vertex) const
    { return problem_.model().stencils(vertex).elementIndices(); }

    //! we should never call this method for cc models
    template<class T = TypeTag>
    const typename std::enable_if<!GET_PROP_VALUE(T, ImplicitIsBox), Stencil>::type
    getStencil(const Vertex& vertex) const
    { return Stencil(0); }

    //! Calculate the velocities for the scvs in the element
    //! We assume the local containers to be bound to the complete stencil
    template<class VelocityVector>
    void calculateVelocity(VelocityVector& velocity,
                           const ElementVolumeVariables& elemVolVars,
                           const FVElementGeometry& fvGeometry,
                           const Element& element,
                           int phaseIdx)
    {
        if (!velocityOutput_) return;

        const auto geometry = element.geometry();
        Dune::GeometryType geomType = geometry.type();

        // get the transposed Jacobian of the element mapping
        const auto& referenceElement = ReferenceElements::general(geomType);
        const auto& localPos = referenceElement.position(0, 0);
        const auto jacobianT2 = geometry.jacobianTransposed(localPos);

        // bind the element flux variables cache
        auto elemFluxVarsCache = localView(problem_.model().globalFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

        // the upwind term to be used for the volume flux evaluation
        auto upwindTerm = [phaseIdx](const VolumeVariables& volVars) { return volVars.mobility(phaseIdx); };

        if (isBox)
        {
            using ScvVelocities = Dune::BlockVector<Dune::FieldVector<Scalar, dimWorld> >;
            ScvVelocities scvVelocities(fvGeometry.numScv());
            scvVelocities = 0;

            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                    continue;

                // local position of integration point
                const auto localPosIP = geometry.local(scvf.ipGlobal());

                // Transformation of the global normal vector to normal vector in the reference element
                const auto jacobianT1 = geometry.jacobianTransposed(localPosIP);
                const auto globalNormal = scvf.unitOuterNormal();
                GlobalPosition localNormal(0);
                jacobianT1.mv(globalNormal, localNormal);
                localNormal /= localNormal.two_norm();

                // insantiate the flux variables
                FluxVariables fluxVars;
                fluxVars.initAndComputeFluxes(problem_,
                                              element,
                                              fvGeometry,
                                              elemVolVars,
                                              scvf,
                                              elemFluxVarsCache);

                // get the volume flux divided by the area of the
                // subcontrolvolume face in the reference element
                // TODO: Divide by extrusion factor!!?
                Scalar localArea = scvfReferenceArea_(geomType, scvf.index());
                Scalar flux = fluxVars.advectiveFlux(phaseIdx, upwindTerm) / localArea;

                // transform the volume flux into a velocity vector
                GlobalPosition tmpVelocity = localNormal;
                tmpVelocity *= flux;

                scvVelocities[scvf.insideScvIdx()] += tmpVelocity;
                scvVelocities[scvf.outsideScvIdx()] += tmpVelocity;
            }

            // transform vertex velocities from local to global coordinates
            for (auto&& scv : scvs(fvGeometry))
            {
                int vIdxGlobal = scv.dofIndex();

                // calculate the subcontrolvolume velocity by the Piola transformation
                Dune::FieldVector<CoordScalar, dimWorld> scvVelocity(0);

                jacobianT2.mtv(scvVelocities[scv.index()], scvVelocity);
                scvVelocity /= geometry.integrationElement(localPos)*cellNum_[vIdxGlobal];
                // add up the wetting phase subcontrolvolume velocities for each vertex
                velocity[vIdxGlobal] += scvVelocity;
            }
        }
        else
        {
            if (fvGeometry.numScvf() > element.subEntities(1))
                DUNE_THROW(Dune::NotImplemented, "Velocity output for non-conforming grids");

            if (!geomType.isCube() && !geomType.isSimplex())
                DUNE_THROW(Dune::NotImplemented, "Velocity output for other geometry types than cube and simplex");

            // first we extract the corner indices for each scv for the CIV method
            // for network grids there might be multiple intersection with the same geometryInInside
            // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
            // here we keep track of them
            std::vector<bool> handledScvf;
            if (dim < dimWorld) handledScvf.resize(element.subEntities(1), false);

            // find the local face indices of the scvfs (for conforming meshes)
            std::vector<unsigned int> scvfIndexInInside(element.subEntities(1));
            int localScvfIdx = 0;
            for (const auto& intersection : intersections(problem_.gridView(), element))
            {
                if (dim < dimWorld) if (handledScvf[intersection.indexInInside()]) continue;

                if (intersection.neighbor() || intersection.boundary())
                {
                    scvfIndexInInside[localScvfIdx++] = intersection.indexInInside();
                    // for surface and network grids mark that we handled this face
                    if (dim < dimWorld) handledScvf[intersection.indexInInside()] = true;
                }
            }

            std::vector<Scalar> scvfFluxes(element.subEntities(1), 0.0);
            localScvfIdx = 0;
            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (!scvf.boundary())
                {
                    FluxVariables fluxVars;
                    fluxVars.initAndComputeFluxes(problem_, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
                    scvfFluxes[scvfIndexInInside[localScvfIdx]] = fluxVars.advectiveFlux(phaseIdx, upwindTerm);
                }
                else
                {
                    auto bcTypes = problem_.boundaryTypes(element, scvf);
                    if (bcTypes.hasOnlyDirichlet())
                    {
                        FluxVariables fluxVars;
                        fluxVars.initAndComputeFluxes(problem_, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
                        scvfFluxes[scvfIndexInInside[localScvfIdx]] = fluxVars.advectiveFlux(phaseIdx, upwindTerm);
                    }
                }

                // increment scvf counter
                localScvfIdx++;
            }

            // Correct boundary fluxes in case of Neumann conditions.
            // In this general setting, it would be very difficult to
            // calculate correct phase, i.e., volume, fluxes from arbitrary
            // Neumann conditions. We approximate the Neumann flux by the
            // flux on the opposite face. For extremely distorted grids this can
            // lead to unexpected results (but then TPFA also leads to unexpected results).
            localScvfIdx = 0;
            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                {
                    auto bcTypes = problem_.boundaryTypes(element, scvf);
                    if (bcTypes.hasNeumann())
                    {
                        // cubes
                        if (dim == 1 || geomType.isCube())
                        {
                            const auto fIdx = scvfIndexInInside[localScvfIdx];
                            const auto fIdxOpposite = fIdx%2 ? fIdx-1 : fIdx+1;
                            scvfFluxes[fIdx] = -scvfFluxes[fIdxOpposite];
                        }
                        // simplices
                        else if (geomType.isSimplex())
                            scvfFluxes[scvfIndexInInside[localScvfIdx]] = 0;
                    }
                }

                // increment scvf counter
                localScvfIdx++;
            }

            Dune::FieldVector <Scalar, dim> refVelocity;
            // cubes: On the reference element simply average over opposite fluxes
            // note that this is equal to a corner velocity interpolation method
            if (dim == 1 || geomType.isCube())
            {
                for (int i = 0; i < dim; i++)
                    refVelocity[i] = 0.5 * (scvfFluxes[2*i + 1] - scvfFluxes[2*i]);
            }
            // simplices: Raviart-Thomas-0 interpolation evaluated at the cell center
            else if (geomType.isSimplex())
            {
                for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                {
                    refVelocity[dimIdx] = -scvfFluxes[dim - 1 - dimIdx];
                    for (int fIdx = 0; fIdx < dim + 1; fIdx++)
                        refVelocity[dimIdx] += scvfFluxes[fIdx]/(dim + 1);
                }
            }
            // 3D prism and pyramids
            else
                DUNE_THROW(Dune::NotImplemented, "velocity output for cell-centered and prism/pyramid");

            Dune::FieldVector<Scalar, dimWorld> scvVelocity(0);
            jacobianT2.mtv(refVelocity, scvVelocity);

            scvVelocity /= geometry.integrationElement(localPos);

            int eIdxGlobal = problem_.elementMapper().index(element);

            velocity[eIdxGlobal] = scvVelocity;

        } // cell-centered
    }

private:
    // The area of a subcontrolvolume face in a reference element.
    // The 3d non-cube values have been calculated with quadrilateralArea3D
    // of boxfvelementgeometry.hh.
    static Scalar scvfReferenceArea_(Dune::GeometryType geomType, int fIdx)
    {
        if (dim == 1 || geomType.isCube())
        {
            return 1.0/(1 << (dim-1));
        }
        else if (geomType.isTriangle())
        {
            static const Scalar faceToArea[] = {0.372677996249965,
                                                0.372677996249965,
                                                0.235702260395516};
            return faceToArea[fIdx];
        }
        else if (geomType.isTetrahedron())
        {
            static const Scalar faceToArea[] = {0.102062072615966,
                                                0.102062072615966,
                                                0.058925565098879,
                                                0.102062072615966,
                                                0.058925565098879,
                                                0.058925565098879};
            return faceToArea[fIdx];
        }
        else if (geomType.isPyramid())
        {
            static const Scalar faceToArea[] = {0.130437298687488,
                                                0.130437298687488,
                                                0.130437298687488,
                                                0.130437298687488,
                                                0.150923085635624,
                                                0.1092906420717,
                                                0.1092906420717,
                                                0.0781735959970571};
            return faceToArea[fIdx];
        }
        else if (geomType.isPrism())
        {
            static const Scalar faceToArea[] = {0.166666666666667,
                                                0.166666666666667,
                                                0.166666666666667,
                                                0.186338998124982,
                                                0.186338998124982,
                                                0.117851130197758,
                                                0.186338998124982,
                                                0.186338998124982,
                                                0.117851130197758};
            return faceToArea[fIdx];
        }
        else {
            DUNE_THROW(Dune::NotImplemented, "scvf area for unknown GeometryType");
        }
    }

private:
    const Problem& problem_;
    bool velocityOutput_;
    std::vector<int> cellNum_;
};

} // end namespace Dumux

#endif