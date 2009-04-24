// $Id$

#ifndef DUNE_FVTRANSPORT_DEPRECATED_HH
#define DUNE_FVTRANSPORT_DEPRECATED_HH

#warning this file is deprecated, please use fvtransport.hh instead

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/istl/bvector.hh>
#include "dumux/transport/transport_deprecated.hh"
#include "dumux/transport/fv/numericalflux.hh"
#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/transport/fv/parfvdatahandle.hh"

namespace Dune {
//! \ingroup transport
//! The finite volume model for the solution of the transport equation
template<class G, class RT, class VC> class DeprecatedFVTransport :
        public DeprecatedTransport< G, RT, VC> {
    template<int dim> struct ElementLayout {
        bool contains(Dune::GeometryType gt) {
            return gt.dim() == dim;
        }
    };

    enum {dim = G::dimension};
    enum {dimworld = G::dimensionworld};
    typedef BlockVector< Dune::FieldVector<RT,1> > PressType;
    typedef BlockVector< FieldVector<FieldVector<RT, G::dimension>, 2*G::dimension> >
    VelType;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::LevelGridView GV;
    typedef typename GV::IndexSet IS;
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV,ElementLayout> EM;
    typedef typename G::template Codim<0>::EntityPointer EntityPointer;
    typedef typename G::LevelGridView::IntersectionIterator
    IntersectionIterator;
    typedef typename G::ctype ct;
    typedef BlockVector< Dune::FieldVector<RT,dim> > SlopeType;

public:
    typedef BlockVector< Dune::FieldVector<RT,1> > RepresentationType;
    /*!
     *  \param t time
     *  \param dt time step size to estimate
     *  \param update vector to be filled with the update
     *
     *  This method calculates the update vector, i.e., the FV discretization
     *  of \f$\text{div}\, (f_\text{w} \boldsymbol{v}_t)\f$. The total velocity
     *  \f$\boldsymbol{v}_t)\f$ has to be given by the block vector \a velocity,
     *  containing values at each cell face. The fractional flow function \f$f_\text{w}\f$
     *  is realized by the numerical flux function \a numericalFlux such that the \f$i\f$-th entry
     *  of the vector update is obtained by calculating
     *  \f[ \sum_{j \in \mathcal{N}(i)} v_{ij}g(S_i, S_j) - v_{ji}g(S_j, S_i), \f]
     *  where \f$\mathcal{N}(i)\f$ denotes the set of neighbors of the cell \f$i\f$ and
     *  \f$g\f$ stands for \a numericalFlux. The normal velocities \f$v_{ij}\f$ and \f$v_{ji}\f$
     *  are given by
     *  \f{align*} v_{ij} = \int_{\Gamma_{ij}} \max(\boldsymbol{v}_\text{t}{\cdot}\boldsymbol{n}_{ij}, \, 0), \qquad
     *  v_{ji} = \int_{\Gamma_{ij}} \min(\boldsymbol{v}_\text{t}{\cdot}\boldsymbol{n}_{ij}, \, 0), \f}
     *  where \f$\boldsymbol{n}_{ij}\f$ denotes the unit normal vector from cell \f$i\f$ towards cell \f$j\f$.
     *
     *  Additionally to the \a update vector, the recommended time step size \a dt is calculated
     *  employing the usual CFL condition.
     */
    int update(const RT t, RT& dt, RepresentationType& updateVec, RT& cFLFac);

    void initialTransport();


    /*! @brief constructor
     *
     * @param g a DUNE grid object
     * @param prob an object of class DeprecatedTransportProblem or derived
     * @param lev the grid level on which the DeprecatedTransport equation is to be solved.
     * @param diffPart an object of class DiffusivePart or derived. This determines the diffusive flux incorporated in the transport.
     * @param rec flag to switch on linear reconstruction (second order TVD)
     * @param amax alphamax parameter for slope limiter in TVD
     * @param numFl an object of class Numerical Flux or derived
     */
    DeprecatedFVTransport(G& g, DeprecatedTransportProblem<G, RT, VC>& prob, int lev = 0,
                          DiffusivePart<G,RT>& diffPart = *(new DiffusivePart<G, RT>), bool rec = false,
                          double amax = 0.8, const NumericalFlux<RT>& numFl = *(new Upwind<RT>)) :
        DeprecatedTransport<G, RT, VC>(g, prob, lev),
        elementmapper( g.levelView(lev)), gridview(g.levelView(lev)),
        indexset(gridview.indexSet()), reconstruct(rec),
        numFlux(numFl), diffusivePart(diffPart), alphamax(amax)
    {}

private:

    void CalculateSlopes(SlopeType& slope, RT t, RT& cFLFactor);

private:
    EM elementmapper;
    const GV& gridview;
    const IS& indexset;
    bool reconstruct;
    const NumericalFlux<RT>& numFlux;
    const DiffusivePart<G, RT>& diffusivePart;
    double alphamax;
};

template<class G, class RT, class VC> int DeprecatedFVTransport<G, RT, VC>::update(const RT t, RT& dt,
                                                                                   RepresentationType& updateVec, RT& cFLFac = 1) {
    // initialize dt very large
    dt = 1E100;

    // TODO: remove DEBUG--->
    double maxAd = 0;
    double maxDiff = 0;
    // <---DEBUG

    // set update vector to zero
    updateVec = 0;

    SlopeType slope(elementmapper.size());
    CalculateSlopes(slope, t,cFLFac);

    // compute update vector
    Iterator eendit = this->grid().template lend<0>(this->level());
    for (Iterator it = this->grid().template lbegin<0>(this->level()); it != eendit; ++it) {
        // cell geometry type
        Dune::GeometryType gt = it->geometry().type();

        // cell center in reference element
        const Dune::FieldVector<ct,dim>
            &local = Dune::ReferenceElements<ct,dim>::general(gt).position(0, 0);

        // yufei temmp
        Dune::FieldVector<ct,dim> global = it->geometry().global(local);

        // cell volume, assume linear map here
        double volume = it->geometry().integrationElement(local)
            *Dune::ReferenceElements<ct,dim>::general(gt).volume();

        // cell index
        int indexi = elementmapper.map(*it);

        // for time step calculation
        double sumfactor = 0;
        double sumfactor2 = 0;
        double sumDiff = 0;
        double sumDiff2 = 0;

        // run through all intersections with neighbors and boundary
        IntersectionIterator
            endit = it->ilevelend();
        for (IntersectionIterator
                 is = it->ilevelbegin(); is
                 !=endit; ++is) {
            // local number of facet
            int numberInSelf = is->indexInInside();

            // get geometry type of face
            Dune::GeometryType gtf = is->geometryInInside().type();

            // center in face's reference element
            const Dune::FieldVector<ct,dim-1>&
                facelocal = Dune::ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

            // center of face inside volume reference element
            const Dune::FieldVector<ct,dim>&
                facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(is->indexInInside(),1);

            // get normal vector scaled with volume
            Dune::FieldVector<ct,dimworld> integrationOuterNormal = is->integrationOuterNormal(facelocal);
            integrationOuterNormal *= Dune::ReferenceElements<ct,dim-1>::general(gtf).volume();

            // compute factor occuring in flux formula
            double velocityIJ = std::max(this->transproblem.variables.vTotal(*it, numberInSelf)*integrationOuterNormal/(volume), 0.0);

            double factor=0, diffFactor=0, totfactor=0;

            // handle interior face
            if (is->neighbor())
            {
                // access neighbor
                EntityPointer outside = is->outside();
                int indexj = elementmapper.map(*outside);

                // compute flux from one side only
                // this should become easier with the new IntersectionIterator functionality!
                if ( it->level()>=outside->level() )
                {
                    // compute factor in neighbor
                    Dune::GeometryType nbgt = outside->geometry().type();
                    const Dune::FieldVector<ct,dim>&
                        nblocal = Dune::ReferenceElements<ct,dim>::general(nbgt).position(0,0);

                    double velocityJI = std::max(-(this->transproblem.variables.vTotal(*it, numberInSelf)*integrationOuterNormal/volume), 0.0);

                    // cell center in global coordinates
                    Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);

                    // neighbor cell center in global coordinates
                    Dune::FieldVector<ct,dimworld> nbglobal = outside->geometry().global(nblocal);

                    // distance vector between barycenters
                    Dune::FieldVector<ct,dimworld> distVec = global - nbglobal;

                    // compute distance between cell centers
                    double dist = distVec.two_norm();

                    // get saturation value at cell center
                    double satI = this->transproblem.variables.saturation[indexi];

                    // get saturation value at neighbor cell center
                    double satJ = this->transproblem.variables.saturation[indexj];

                    // calculate the saturation gradient
                    Dune::FieldVector<ct,dim> satGradient = distVec;
                    satGradient *= (satJ - satI)/(dist*dist);

                    // the arithmetic average
                    double satAvg = 0.5*(satI + satJ);

                    // get the diffusive part
                    double diffPart = this->diffusivePart(*it, numberInSelf, satAvg, satGradient, t, satI, satJ)*integrationOuterNormal;

                    // CAREFUL: works only for axisymmetric grids
                    if (reconstruct) {
                        for (int k = 0; k < dim; k++)
                            if (fabs(distVec[k]) > 0.5*dist)
                            {
                                satI -= fabs(distVec[k])/distVec[k]*0.5*dist*slope[indexi][k];
                                satJ += fabs(distVec[k])/distVec[k]*0.5*dist*slope[indexj][k];
                            }
                    }

                    double fI = this->transproblem.materialLaw.fractionalW(satI);
                    double fJ = this->transproblem.materialLaw.fractionalW(satJ);

                    diffFactor = diffPart / volume;
                    factor = diffFactor
                        + velocityJI*numFlux(satJ, satI, fJ, fI)
                        - velocityIJ*numFlux(satI, satJ, fI, fJ);
                    totfactor = velocityJI - velocityIJ;
                    totfactor *= (fI-fJ)/(satI-satJ);
                }
            }

            // handle boundary face
            if (is->boundary())
            {
                // center of face in global coordinates
                Dune::FieldVector<ct,dimworld> faceglobal = is->geometry().global(facelocal);

                //get boundary type
                BoundaryConditions::Flags bctype = this->transproblem.bctype(faceglobal, *it, facelocalDim);

                if (bctype == BoundaryConditions::dirichlet)
                {
                    // get saturation value at cell center
                    double satI = this->transproblem.variables.saturation[indexi];

                    double velocityJI = std::max(-(this->transproblem.variables.vTotal(*it, numberInSelf)*integrationOuterNormal/volume), 0.0);

                    double satBound = this->transproblem.g(faceglobal, *it, facelocalDim);

                    // cell center in global coordinates
                    Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);

                    // distance vector between barycenters
                    Dune::FieldVector<ct,dimworld> distVec = global - faceglobal;

                    // compute distance between cell centers
                    double dist = distVec.two_norm();

                    // calculate the saturation gradient
                    Dune::FieldVector<ct,dim> satGradient = distVec;
                    satGradient *= (satBound - satI)/(dist*dist);

                    // the arithmetic average
                    double satAvg = 0.5*(satI + satBound);

                    // get the diffusive part
                    double diffPart = this->diffusivePart(*it, numberInSelf, satAvg, satGradient, t, satI, satBound)*integrationOuterNormal;

                    // CAREFUL: works only for axisymmetric grids
                    if (reconstruct) {
                        for (int k = 0; k < dim; k++)
                            if (fabs(distVec[k]) > 0.5*dist)
                            {
                                //TODO remove DEBUG--->
                                //double gabagabahey = slope[indexi][k];
                                //<---DEBUG
                                satI -= fabs(distVec[k])/distVec[k]*dist*slope[indexi][k];
                            }
                    }

                    double fI = this->transproblem.materialLaw.fractionalW(satI);
                    double fBound = this->transproblem.materialLaw.fractionalW(satBound);
                    diffFactor = diffPart / volume;
                    factor = diffFactor
                        + velocityJI*numFlux(satBound, satI, fBound, fI)
                        - velocityIJ*numFlux(satI, satBound, fI, fBound);
                    totfactor = velocityJI - velocityIJ;
                    totfactor *= (fI-fBound)/(satI-satBound);
                }
                else
                {
                    //double J = this->transproblem.J(faceglobal, *it, facelocalDim);
                    //factor = J*faceVol;
                    factor = 0;
                    totfactor = 0;
                }
            }

            // add to update vector
            updateVec[indexi] += factor;

            // for time step calculation
            if (totfactor>=0)
                sumfactor += totfactor;
            else
                sumfactor2 += (-totfactor);
            if (diffFactor>=0)
                sumDiff += diffFactor;
            else
                sumDiff += (-diffFactor);
        }
        // end all intersections
        // compute dt restriction
        sumfactor = std::max(sumfactor, sumfactor2);
        sumDiff = std::max(sumDiff, sumDiff2);
        sumfactor = std::max(sumfactor, 10*sumDiff);
        dt = std::min(dt, 1.0/sumfactor);

        // TODO: remove DEBUG--->
        maxAd = std::max(maxAd, sumfactor);
        maxDiff = std::max(maxDiff, sumDiff);
        // <---DEBUG

    } // end grid traversal

    //Correct maximal available volume in the CFL-Criterium
    RT ResSaturationFactor = 1-this->transproblem.materialLaw.wettingPhase.Sr()-this->transproblem.materialLaw.nonwettingPhase.Sr();
    dt = dt*this->transproblem.porosity()*ResSaturationFactor;

    // global min over all partitions
    dt = this->grid().comm().min(dt);

    updateVec /= this->transproblem.porosity();

    //TODO remove DEBUG--->
    //std::cout<<"maxAd "<< maxAd << "\t maxDiff "<< maxDiff<<std::endl;
    // <---DEBUG

    // exchange update
    typedef VectorExchange<EM,RepresentationType> VecExchange;
    VecExchange dataHandle(elementmapper,updateVec);
    this->grid().template communicate<VecExchange>(dataHandle, InteriorBorder_All_Interface,
                                                   ForwardCommunication);

    return 0;
}

template<class G, class RT, class VC> void DeprecatedFVTransport<G, RT, VC>::initialTransport() {
    //    std::cout<<"initsat = "<<&this->transproblem.variables.saturation<<std::endl;
    // iterate through leaf grid an evaluate c0 at cell center
    Iterator eendit = this->grid().template lend<0>(this->level());
    for (Iterator it = this->grid().template lbegin<0>(this->level()); it != eendit; ++it) {
        // get geometry type
        Dune::GeometryType gt = it->geometry().type();

        // get cell center in reference element
        const Dune::FieldVector<ct,dim>
            &local = Dune::ReferenceElements<ct,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);

        // initialize cell concentration
        this->transproblem.variables.saturation[elementmapper.map(*it)] = this->transproblem.initSat(global, *it, local);
    }
    return;
}

template<class G, class RT, class VC> void DeprecatedFVTransport<G, RT, VC>::CalculateSlopes(
                                                                                             SlopeType& slope, RT t, RT& cFLFactor) {

    double stabilityFactor = 1.0 - cFLFactor*sqrt(cFLFactor);

    Iterator endit = this->grid().template lend<0>(this->level());
    for (Iterator it = this->grid().template lbegin<0>(this->level()); it!=endit; ++it) {
        // get some cell properties
        Dune::GeometryType gt = it->geometry().type();
        const Dune::FieldVector<ct,dim>
            &local = Dune::ReferenceElements<ct,dim>::general(gt).position(0, 0);
        Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);
        int indexi = elementmapper.map(*it);

        // vector containing the distances to the neighboring cells
        Dune::FieldVector<double, 2*dim> dist;

        // vector containing the saturations of the neighboring cells
        Dune::FieldVector<double, 2*dim> saturation;

        // location[k], k = 0,...,2dim-1, contains the local index w.r.t. IntersectionIterator,
        // i.e. the numberInSelf, which is east, west, north, south, top, bottom to the cell center
        Dune::FieldVector<int, 2*dim> location(0);

        // run through all intersections with neighbors and boundary
        IntersectionIterator
            isend = it->ilevelend();
        for (IntersectionIterator
                 is = it->ilevelbegin(); is
                 !=isend; ++is) {
            // local number of facet
            int numberInSelf = is->indexInInside();

            // handle interior face
            if (is->neighbor()) {
                // access neighbor
                EntityPointer outside = is->outside();
                int indexj = elementmapper.map(*outside);

                // get saturation value
                saturation[numberInSelf] = this->transproblem.variables.saturation[indexj];

                // compute factor in neighbor
                Dune::GeometryType nbgt = outside->geometry().type();
                const Dune::FieldVector<ct,dim>
                    &nblocal = Dune::ReferenceElements<ct,dim>::general(nbgt).position(0, 0);

                // neighboring cell center in global coordinates
                Dune::FieldVector<ct,dimworld>nbglobal = outside->geometry().global(nblocal);

                // distance vector between barycenters
                Dune::FieldVector<ct,dimworld> distVec = global - nbglobal;

                // compute distance between cell centers
                dist[numberInSelf] = distVec.two_norm();

                // CAREFUL: works only for axiparallel grids
                for (int k = 0; k < dim; k++)
                    if (nbglobal[k] - global[k]> 0.5*dist[numberInSelf]) {
                        location[2*k] = numberInSelf;
                    } else if (nbglobal[k] - global[k]< -0.5*dist[numberInSelf]) {
                        location[2*k + 1] = numberInSelf;
                    }
            }

            // handle boundary face
            if (is->boundary()) {
                // get geometry type of face
                Dune::GeometryType gtf = is->geometryInInside().type();

                // center in face's reference element
                const Dune::FieldVector<ct,dim-1>&
                    facelocal = Dune::ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

                // center of face in global coordinates
                Dune::FieldVector<ct,dimworld>
                    faceglobal = is->geometry().global(facelocal);

                // get saturation value
                saturation[numberInSelf] = this->transproblem.variables.saturation[indexi];//this->transproblem.g(faceglobal, *it, facelocalDim);

                // distance vector between barycenters
                Dune::FieldVector<ct,dimworld> distVec = global - faceglobal;

                // compute distance
                dist[numberInSelf] = distVec.two_norm();

                // CAREFUL: works only for axiparallel grids
                for (int k = 0; k < dim; k++)
                    if (faceglobal[k] - global[k] > 0.5*dist[numberInSelf])
                    {
                        location[2*k] = numberInSelf;
                    }
                    else if (faceglobal[k] - global[k] < -0.5*dist[numberInSelf])
                    {
                        location[2*k + 1] = numberInSelf;
                    }
            }
        } // end all intersections

        for (int k = 0; k < dim; k++) {
            double slopeIK = (saturation[location[2*k]]
                              - saturation[location[2*k + 1]])/(dist[location[2*k]]
                                                                + dist[location[2*k + 1]]);

            double alphaIK = 1.0;
            if (fabs(slopeIK) > 1e-8*dist[location[2*k]]) {
                double satI = this->transproblem.variables.saturation[indexi];
                double alphaRight = stabilityFactor*fabs(2.0
                                                         /(dist[location[2*k]]*slopeIK)
                                                         *(saturation[location[2*k]] - satI));
                double alphaLeft = stabilityFactor*fabs(2.0
                                                        /(dist[location[2*k + 1]]*slopeIK)*(satI
                                                                                            - saturation[location[2*k + 1]]));
                alphaIK = std::min(std::min(std::max(alphaRight, 0.0),
                                            std::max(alphaLeft, 0.0)), alphamax);
            }

            slope[indexi][k] = alphaIK*slopeIK;

        }
    } // end grid traversal

    return;
}

}
#endif
