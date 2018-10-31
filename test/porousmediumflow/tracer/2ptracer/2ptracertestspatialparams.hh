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
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem
 */
#ifndef DUMUX_TWOP_TRACER_TEST_SPATIAL_PARAMS_HH
#define DUMUX_TWOP_TRACER_TEST_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <iostream>

namespace Dumux {

/*!
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem
 */
template<class FVGridGeometry, class Scalar>
class TwoPTracerTestSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         TwoPTracerTestSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar,
                                       TwoPTracerTestSpatialParams<FVGridGeometry, Scalar>>;

    //static const int numPhases = ModelTraits::numPhases();
    //static const int numScvf = FVGridGeometry::numScvf();
    //using FieldMatrix = Dune::FieldMatrix <Scalar, numPhases, numScvf>;
    //using PhaseVector = std::vector<Scalar>;
    //using FieldVector = Dune::FieldVector<PhaseVector, 2>;

    //static const int numPhases = ModelTraits::numPhases();
    //using ScvfVector  = std::vector<Scalar>;
    //using FieldVector = Dune::FieldVector<ScvfVector, numPhases>;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;


public:

    TwoPTracerTestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry) {}

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief Define the dispersivity.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     * \param elemSol The solution for all dofs of the element
     */
    template<class ElementSolution>
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    { return 0; }

    //! Fluid properties that are spatial params in the tracer model
    //! They can possible vary with space but are usually constants

    //! fluid density

    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return density_[this->fvGridGeometry().elementMapper().index(element)];; }

    void setDensity(const std::vector<Scalar>& d)
    { density_ = d; }

    //! fluid molar mass
    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    { return 0.018; }

    Scalar fluidMolarMass(const GlobalPosition &globalPos) const
    { return 0.018; }

    //! velocity field
    template<class ElementVolumeVariables>
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
//         std::cout<<"SpatialParams: volumeFlux_[scvf.index()="<<scvf.index()<<"]="<<volumeFlux_[scvf.index()]<<std::endl;
            return volumeFlux_[scvf.index()];
    }

    void setVolumeFlux(const std::vector<Scalar>& f)
    { volumeFlux_ = f; }

    //! saturation from twoPProblem
    Scalar saturation(const Element &element,
                      const SubControlVolume& scv) const
    { return saturation_[scv.dofIndex()]; }

    void setSaturation(const std::vector<Scalar>& s)
    { saturation_ = s; }
private:

//     std::vector<Scalar> volumeFlux_(fvGridGeometry->numScvf(), 0.0);
    std::vector<Scalar> volumeFlux_;
    std::vector<Scalar> density_;
    std::vector<Scalar> saturation_;

    //FieldVector volumeFlux_;
    //std::vector<Scalar> volumeFlux_;
};

} // end namespace Dumux

#endif
