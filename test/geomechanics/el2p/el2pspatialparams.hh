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
 * \brief The spatial parameters for the El2P_TestProblem which uses the
 *        linear elastic two-phase model
 */
#ifndef DUMUX_ELTWOPSPARAMETERS_HH
#define DUMUX_ELTWOPSPARAMETERS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
// #include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/geomechanics/el2p/model.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class El2PSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(El2PSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(El2PSpatialParams, SpatialParams, El2PSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(El2PSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedVanGenuchten<Scalar> EffectiveLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> type;
};
}
/*!
 * \ingroup ElTwoPBoxModel
 * \brief The spatial parameters for the El2P_TestProblem which uses the
 *        linear elastic two-phase model
 */
template<class TypeTag>
class El2PSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;


    El2PSpatialParams(const GridView &gridView)
    : ParentType(gridView)
    {
        // episode index
        episode_ = 0;
        // intrinsic permeabilities [m^2]
        Kinit_ = Scalar(0.0); // init permeability
//         K_ = Scalar(0.0); // permeability
        k_ = Scalar(0.0); // conductivity.khodam
        ks_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.ks);//khodam [m/s] we assume this is for standard temp and presure
        for (int i = 0; i < dim; i++){
            Kinit_[i][i] = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.Kinit);
            K_ = (ks_*1.E-3)/(1000*9.81);//khodam [m/s] dynamic viscosity and density are defined from simpleh2o.hh in dumux/material/component , and didnot make them pressure and temp dependent, becasue Ks is constant therefore I consider the standard tem,press for it. k=ks*Krw
        }

        // porosities [-]
        phi_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.Phi);

        // rock density [kg/m^3]
        rockDensity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, FailureParameters.rockDensity);

        // Young's modulus [Pa]
        E_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.E);
        // Poisson's ratio [-]
        nu_ = GET_RUNTIME_PARAM(TypeTag, Scalar,ElasticParameters.nu);
        // Lame parameters [Pa]
        lambda_ = (E_ * nu_) / ((1 + nu_)*(1 - 2 * nu_));
        mu_ = E_ / (2 * (1 + nu_));


        // given Van Genuchten m
        n_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.n);

        m_ = 1.0 - (1.0 / n_);

        alpha_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.alpha);
        Scalar swr_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.swr); //residual water content
        Scalar snr_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.snr);

        // residual saturations
        MaterialParams_.setSwr(swr_);
        MaterialParams_.setSnr(snr_);

        MaterialParams_.setVgAlpha(alpha_);
        MaterialParams_.setVgn(n_);


     }

    ~El2PSpatialParams()
    {}

    /*!
     * \brief This function sets the private variable episode_ to the current episode index
     * which is checked in the hydraulic parameter functions to identify if we are still in the
     * initialization run (episode_ == 1)
     *
     * \param episode The episode index
     */
    void setEpisode(const int& episode)
    {
        episode_ = episode;
        std::cout<< "episode set to: "<< episode_<<std::endl;
    }

    /*!
     * \brief Apply the intrinsic permeability tensor \f$[m^2]\f$ to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     *
     * During the initialization period the intrinsic permeability can be set to a larger
     * value in order to accelerate the calculation of the hydrostatic pressure distribution.
     */
    const DimMatrix intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       int scvIdx) const
    {
        if(episode_ <= 1)
            return Kinit_; // intrinsic permeability applied during initialization
        else
             return K_; // intrinsic permeability
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the soil
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    double porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {
            return phi_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the soil
     *
     * \param globalPos The global position of the vertex
     */
    double porosity(const GlobalPosition& globalPos) const
    {
            return phi_;
    }

    /*!
     * \brief Define the density \f$[kg/m^3]\f$ of the rock
     *
     * \param element The finite element
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    const Scalar rockDensity(const Element &element,
                                        int scvIdx) const
    {
        return rockDensity_;
    }

    /*!
     * \brief Define the density \f$[kg/m^3]\f$ of the rock
     *
     * \param globalPos The global position of the vertex
     */
    const Scalar rockDensity(const GlobalPosition &globalPos) const
    {
        return rockDensity_;
    }

    /*!
     * \brief Define the Lame parameters \f$[Pa]\f$ linear elastic rock
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    const Dune::FieldVector<Scalar,2> lameParams(const Element &element,
                                           const FVElementGeometry &fvGeometry,
                                           int scvIdx) const
    {
        // Lame parameters
        Dune::FieldVector<Scalar, 2> param;

        param[0] = lambda_;
        param[1] = mu_;

        return param;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     *
     * \param element The current element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume.
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvGeometry,
                                               int scvIdx) const
    {
        return MaterialParams_;
    }
    const MaterialLawParams& materialLawParams(const GlobalPosition& globalPos) const//I added this then I can read the materialLawParams via spatialparams in problem to calculate sw just based on globalPos and no need for fvGeometry, element and scvIdx. then I dont need to read via model, as in localoperatopr.hh, and similar to porosirty I just read it via spatialparams.
    {
        return MaterialParams_;
    }

private:
    Dune::FieldMatrix<Scalar,dim,dim> K_, Kinit_;
    Scalar layerBottom_;
    Scalar rockDensity_;
    Scalar phi_, phiInit_;
    Scalar lambda_;
    Scalar mu_;
    Scalar E_;
    Scalar nu_;
    Scalar BrooksCoreyLambda_, m_, alpha_, n_;
    Scalar k_, ks_;//khodam
    MaterialLawParams MaterialParams_;
    static constexpr Scalar eps_ = 3e-6;
    int episode_;

};
}
#endif