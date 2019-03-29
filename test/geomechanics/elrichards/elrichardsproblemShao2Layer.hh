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
/**
 * \file
 * \brief Definition of a problem, for the Richards flow linear elasticity problem:
 * Problem definition for the deformation of an elastic solid.
 */
#ifndef DUMUX_ELRICHARDS_TESTPROBLEM_HH
#define DUMUX_ELRICHARDS_TESTPROBLEM_HH

#include <dune/pdelab/finiteelementmap/qkfem.hh>

#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/geomechanics/elrichards/model.hh>//already included
#include <dumux/geomechanics/elrichards/amgbackend.hh>


//#include "elrichardsco2tables.hh"
#include "elrichardsspatialparamsShao2Layer.hh"

namespace Dumux
{
template<class TypeTag>
class ElRichards_TestProblem;


// initial conditions for momentum balance equation
template<class TypeTag, int dim>
class InitialDisplacement;

// initial conditions for mass balance equations
template<class TypeTag>
class InitialPressSat;

namespace Properties {
NEW_TYPE_TAG(ElRichards_TestProblem, INHERITS_FROM(BoxElasticRichards, ElRichardsSpatialParams));
NEW_PROP_TAG(InitialDisplacement); //!< The initial displacement function
NEW_PROP_TAG(InitialPressSat); //!< The initial pressure and saturation function

// Set the grid type
SET_TYPE_PROP(ElRichards_TestProblem, Grid, Dune :: UGGrid <3>);


SET_PROP(ElRichards_TestProblem, PressureFEM)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    typedef Dune::PDELab::QkLocalFiniteElementMap<GridView,Scalar,Scalar,1>  type;
};

SET_PROP(ElRichards_TestProblem, DisplacementFEM)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    typedef Dune::PDELab::QkLocalFiniteElementMap<GridView,Scalar,Scalar,1>  type;
};

// Set the problem property
SET_TYPE_PROP(ElRichards_TestProblem, Problem, ElRichards_TestProblem<TypeTag>);

// Set fluid configuration
SET_PROP(ElRichards_TestProblem, FluidSystem)
{
    typedef H2OAirFluidSystem<TypeTag> type;
};

// Set the CO2 table to be used; in this case not the the default table
//SET_TYPE_PROP(ElRichards_TestProblem, CO2Table, ElRichards::CO2Tables);
// Set the salinity mass fraction of the brine in the reservoir
//SET_SCALAR_PROP(ElRichards_TestProblem, ProblemSalinity, 1e-1);

// Set the soil properties
SET_TYPE_PROP(ElRichards_TestProblem, SpatialParams, ElRichardsSpatialParams<TypeTag>);

// Set the initial displacement function
SET_PROP(ElRichards_TestProblem, InitialDisplacement)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
public:
    typedef InitialDisplacement<TypeTag, dim> type;
};

// Set the initial pressure and saturation function
SET_PROP(ElRichards_TestProblem, InitialPressSat)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    typedef InitialPressSat<TypeTag> type;
};

SET_SCALAR_PROP(ElRichards_TestProblem, NewtonMaxRelativeShift, 1e-5);

// use the algebraic multigrid
//SET_TYPE_PROP(ElRichards_TestProblem, LinearSolver, ElRichardsAMGBackend<TypeTag>);
SET_TYPE_PROP(ElRichards_TestProblem, LinearSolver, SuperLUBackend<TypeTag>);

// central differences to calculate the jacobian by default
SET_INT_PROP(ElRichards_TestProblem, ImplicitNumericDifferenceMethod, 0);

// write the stress and displacement output according to rock mechanics
// sign convention (compressive stresses > 0)
SET_BOOL_PROP(ElRichards_TestProblem, VtkRockMechanicsSignConvention, true);
}

/*!
 * \ingroup ElRichardsBoxProblems
 *
 * \brief Problem definition for a Richards flow process
 * in an elastic deformable matrix.
 *
 * This problem simulates an injection of CO2 into the center of a cube with 1000 m x 1000 m x 1000 m
 * dimension. The bottom boundary of this cube is in 2000 m depth. The initialization period is 1e6 s,
 * the real injection period is 1e6 s, the initial timestep is 10 s.
 * Apart from the pressure and the saturation distribution this problems solves for the changes in
 * solid displacement (ux, uy, uz [m]) due to injection. Based on the solid displacement vector
 * the injection-induced changes in the strain and stress tensors are evaluated.
 * Further the porosity and permeability are functions of the solid displacement.
 *
 * During an initialization period of length tInit [s] the pressure field is initialized.
 *
 * After the initialization the real simulation starts and the pressure field from the initialization
 * period is applied as initial condition and for the definition of the lateral Dirichlet
 * boundary conditions. The solid  displacement field is set to zero and the CO2 injection is started.
 */
template<class TypeTag = TTAG(ElRichards_TestProblem)>
class ElRichards_TestProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;//khodam
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;//khodam
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;//khodam baraie pc
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;//khodam baraie pc

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    enum {
        // indices of the primary variables
            uxIdx = Indices::uxIdx,
            uyIdx = Indices::uyIdx,
            uzIdx = Indices::uzIdx,
            pwIdx = Indices::pwIdx,
            hIdx = Indices::hIdx,

            // equation indices (in Boundary condition we define that which equation should be used. then in case of using momentum balance, it will use displacement u for dirichlet and take stress for neuman. and in case of having conti equation, takes pressure for dirichlet and flux for neumann!!!!!!!!!!!!)
            contiEqIdx = Indices::contiEqIdx,
            momentumXEqIdx = Indices::momentumXEqIdx,
            momentumYEqIdx = Indices::momentumYEqIdx,
            momentumZEqIdx = Indices::momentumZEqIdx,

            // phase indices
            wPhaseIdx = Indices::wPhaseIdx //last line doesnt need ,

    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::BlockVector<GlobalPosition> InitialStressField;

    typedef typename GET_PROP_TYPE(TypeTag, LocalFEMSpace) LocalFEMSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

//    typedef typename GET_PROP_TYPE(TypeTag, CO2Table) CO2Table;
//    typedef Dumux::CO2<Scalar, CO2Table> CO2;
//    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;


public:
    /*!
     * \brief The constructor //we cant add a function to the constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param tInitEnd End of initialization period
     */
    ElRichards_TestProblem(TimeManager &timeManager,
                    const GridView &gridView)
        : ParentType(timeManager, gridView),
        gridView_(gridView)
    {
        std::cout << "ElRichards_TestProblem: Initializing the fluid system for the elrichards model\n";


         //eps_ = 3e-6;
        pnRef_ = 1e5;


        // initialize the tables of the fluid system
//         FluidSystem::init(); //this one or the next init function. in el2P this part was commented, mazbe because it was defined by default
         FluidSystem::init(/*Tmin=*/273,
                           /*Tmax=*/300,
                           /*nT=*/5,
                           /*pmin=*/0,
                           /*pmax=*/1e6,
                           /*np=*/10);
        temperature_ = 283.15;
        // resize the pressure field vector with the number of vertices
         pInit_.resize(gridView_.size(dimWorld));
        // fill the pressure field vector with zeros
        std::fill( pInit_.begin(), pInit_.end(), 0.0 );

        // variable which determines if output should be written (initially set to false)
        output_ = false;
        // define if current run is initialization run
        // (initially set to true, will be set to false if initialization is over)
        initializationRun_ = true;
        // defines if feedback from geomechanics on flow is taken into account or not
        // (usually the coupling is switched off for the initialization run)
        coupled_ = false;
        // set initial episode length equal to length of initialization period
        Scalar tInitEnd = GET_RUNTIME_PARAM(TypeTag, Scalar,TimeManager.TInitEnd);
        this->timeManager().startNextEpisode(tInitEnd);
        // transfer the episode index to spatial parameters
        // (during intialization episode hydraulic different parameters might be applied)
        this->spatialParams().setEpisode(this->timeManager().episodeIndex());

//         depthBOR_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Injection.DepthBOR);//later I need to define a function for depth
//          depthBOR(globalPos) = yMax(globalPos)-globalPos[1];
        episodeLength_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TimeManager.EpisodeLength);

        dt_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TimeManager.DtInitial);
     }

    void init()
    {
        if (this->timeManager().time() < 1e-8)
        {
            // set the initial approximated hydrostatic pressure distribution
            // based on an averaged brine density
            // or based on a pressure polynomial
            this->initializePressure();
            // output is written
            this->setOutput(true);
        }

        ParentType::init();
    }


     //instead of bBoxMax which was constant in the whole domain, we define the max elevation or the ground surface at each point
    Scalar zMax(const GlobalPosition &globalPos) const
    {
      if (globalPos[0]<15+eps_)
          return 0.0;
      else if ((globalPos[0] > 15-eps_) && (globalPos[0] < 30 +eps_) )
          return (-0.4*(globalPos[0]-15) - eps_);
      else if (globalPos[0]>30-eps_)
          return -6.0;
    }

        Scalar zMin(const GlobalPosition &globalPos) const
    {
          return -25.0;
    }

    Scalar depth(const GlobalPosition &globalPos) const
    {
        return  zMax(globalPos) - globalPos[2];//distance from the surface
    }

//     Scalar waterTable = -15;

    Scalar wTdepth(const GlobalPosition &globalPos) const
    {
        return (-9.0 - globalPos[2]);
    }

    // note: pInit is < 0 (just due to geomechanics sign convention applied here)
    // initialize the pressure field for initialization run
    // first an approximate hydrostatic pressure field is calculated based on an
    // averaged density. Then the model runs for the initialization period and
    // calculates the real hydrostatic pressure distribution based on the real
    // density distribution. The calculated pressure field is than applied for
    // initialization of the actual model run and for the pressure Dirichlet boundary values.

    void initializePressure()
    {
         pInit_.resize(gridView_.size(dimWorld));
        for(const auto& vertex : vertices(gridView_))
        {
            int vIdxGlobal = this->vertexMapper().index(vertex);
            GlobalPosition globalPos = vertex.geometry().corner(0);

            // initial approximate pressure distribution at start of initialization run
//             pInit_[vIdxGlobal] = -(1.013e5 + (depthBOR_ - globalPos[dimWorld-1]) * waterDensity_ * 9.81);

                pInit_[vIdxGlobal] = -pInit(globalPos); //initial 'water'pressure
        }
    }

    Scalar pInit(const GlobalPosition globalPos)//then we just can cal it by globalPos
    {
        return pnRef_+waterDensity_*9.81*wTdepth(globalPos); //initial 'water'pressure
    }

    // allows to change the coupled_ variable which defines if geomechanical feedback on flow is taken
    // into account
    void setCoupled(bool coupled)
    {
        coupled_ = coupled;
    }

    // returns the coupled_ variable which defines if geomechanical feedback on flow is taken
    // into account
    const bool coupled() const
    {
        return coupled_;
    }

    // allows to change the output_ variable which defines if output is written
    void setOutput(bool output)
    {
        output_ = output;
    }

    // note: pInit is < 0 (just due to geomechanics sign convention applied here)
    // function which is called after the initialization run in
    // order to fill the pressure field vector pInit_ with the
    // pressure result of the initialization
    void setPressure()
    {
        initializationRun_ = false; // initialization run is now finished

        this->setInitializationRun(initializationRun_);
        std::cout<<"ElRichards_TestProblem: initialized pressure field copied to pInit_"<<std::endl;
        for(const auto& vertex : vertices(gridView_))
        {
            int vIdxGlobal = this->vertexMapper().index(vertex);
            pInit_[vIdxGlobal] = -this->model().curSol().base()[vIdxGlobal*2][0];
        }
    }

    // returns the initializationRun_ variable which defines if this is an initialization run
    bool initializationRun()
    {
        return initializationRun_;
    }



    // allows to set the initializationRun_ variable which defines if this is an initialization run
    void setInitializationRun(bool initializationRun)
    {
        initializationRun_ = initializationRun;
    }

    // function which returns an in-situ stress field that needs to be provided
    // for the principal stress calculation
    GlobalPosition initialStress(const GlobalPosition globalPos, const int dofIdxGlobal) const
    {
      GlobalPosition stress;
      Scalar porosity, rockDensity, gravity;
      gravity = -this->gravity()[dimWorld-1];
      porosity = this->spatialParams().porosity(globalPos);
      rockDensity = this->spatialParams().rockDensity(globalPos);
//here we have considered that the system is fully saturated. else the water weight should be multiplied by the saturation. and is the same for hydrostatic pressure= rho*g*h*Se
      // initial total stress field here assumed to be isotropic, lithostatic
//       stress[0] = waterDensity_ * porosity * gravity * (- globalPos[dim-1])
//                   + (1 - porosity) * rockDensity * gravity * ( - globalPos[dim-1]);
//       if(dim >=2)
//       stress[1] = waterDensity_ * porosity * gravity * ( - globalPos[dim-1])
//                   + (1 - porosity) * rockDensity * gravity * ( - globalPos[dim-1]);
//       if(dim == 3)
//       stress[2] = waterDensity_ * porosity * gravity * ( - globalPos[dim-1])
//                   + (1 - porosity) * rockDensity * gravity * ( - globalPos[dim-1]);

     const MaterialLawParams& materialParams = this->spatialParams().materialLawParams(globalPos);

     Scalar i_ = 3; //number of elements along the depth, or number of devisions that we wanna measure sw along the depth, therefore we have i+1 measurements for each globalPos
//      Scalar z_ = (zMax(globalPos)-zMin(globalPos))/i_;//distance of two points that we measure in them or size of the steps .
     Scalar z_ = (zMax(globalPos)-globalPos[2])/i_;//distance of two points that we measure in them or size of the steps .
//      std::cout<< "z= " << z_ <<std::endl;
     int n_ = floor((zMax(globalPos)-globalPos[2])/z_);//number of sections into which the lenght is diveded.  depth dependent. floor0( -7.7 ) --> -7 floor0( 7.7 ) --> 7
//      std::cout<< "n("<< globalPos << ") is " << n_ <<std::endl;
//      std::cout<< "distance "<<  zMax(globalPos)-globalPos[2] <<std::endl;


     Scalar pc=pnRef_-(pnRef_+waterDensity_*9.81*wTdepth(globalPos));
//      Scalar swnew = MaterialLaw::sw(materialParams, pc);
     Scalar swnew =0.;
     GlobalPosition currentGlobalPos = globalPos;

    // loop over poins, number of points is number of sections + 1
     for (int j=0; j<n_ + 1; ++j)
       {

//               std::cout<< "j("<< currentGlobalPos << ") is " << j <<std::endl;
              currentGlobalPos[2] = globalPos[2]+j*z_;//each time we add one step size to the location for the next measurement location
//               std::cout<< "cgp2="<< currentGlobalPos[2]  <<std::endl;
              Scalar pInit= +pnRef_+waterDensity_*9.81*wTdepth(currentGlobalPos);

               Scalar pc=pnRef_-pInit;
               Scalar sw = MaterialLaw::sw(materialParams, pc);
               if (sw > 1.+eps_){ //this has not been defined in the regularizedVanGenuchten.hh because of causing some numerical problem. therefore I define it here.
                   sw = 1.0;
               }
//               std::cout<< "sw("<< currentGlobalPos[2] << ") is " << sw <<std::endl;
               swnew =  (swnew + sw);
        }

           swnew =  (swnew)/(n_ + 1.);// averaged along the considered height


//          std::cout<< "swnew("<< globalPos << ") is " << swnew <<std::endl;

         stress[0] =((1 - porosity) * rockDensity * gravity * ( zMax(globalPos)- globalPos[dim-1])) + (waterDensity_*9.81*porosity* swnew* ( zMax(globalPos)- globalPos[dim-1]));
         if(dim >=2)
          stress[1] = ((1 - porosity) * rockDensity * gravity * ( zMax(globalPos)- globalPos[dim-1])) + (waterDensity_*9.81*porosity* swnew* ( zMax(globalPos)- globalPos[dim-1]));
         if(dim == 3)
          stress[2] = ((1 - porosity) * rockDensity * gravity * ( zMax(globalPos)- globalPos[dim-1])) + (waterDensity_*9.81*porosity* swnew*( zMax(globalPos)- globalPos[dim-1]));


      return stress;
    }
    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    {
        return "elrichards";
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius at the ground surface
     * and a geothermal gradient of 0.03 K/m.
     */
//     Scalar temperatureAtPos(const GlobalPosition &globalPos) const
//     {
//         Scalar T;
//         T = 283.15 + (depthBOR_ - globalPos[dimWorld-1]) * 0.03;
//
//         return T;
//     };
    Scalar temperature() const
    { return temperature_; }; // -> 10°C

          /*!
     * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     *
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The sub control volume index inside the finite
     *               volume geometry
     */
    Scalar referencePressure(const Element &element,
                             const FVElementGeometry &fvGeometry,
                             const int scvIdx) const
    { return pnRef_; };



    // note: pInit is < 0 (just due to geomechanics sign convention applied here)
    // function which returns the initialized pressure at an arbitrary location within the element
    // called from finite element method (elrichardslocaloperator.hh) and evaluated at Gauss points
    Scalar pInit(const GlobalPosition& globalPos, const GlobalPosition& localPos, const Element& element) const
    {
        Scalar pValue = 0.0;

        typename ElRichards_TestProblem<TypeTag>::LocalFEMSpace feMap(this->gridView());
        const typename LocalFEMSpace::Traits::FiniteElementType
        &localFiniteElement = feMap.find(element.geometry().type());
        typedef Dune::FieldVector<CoordScalar, 1> ShapeValue;
        std::vector<ShapeValue> shapeVal;
        localFiniteElement.localBasis().evaluateFunction(localPos, shapeVal);

        for (int i = 0; i < element.subEntities(dim); i++)
        {
            int vIdxGlobal = this->vertexMapper().subIndex(element, i, dim);
            pValue += pInit_[vIdxGlobal] * shapeVal[i];
        }

        return pValue;
    }

    // note: pInit is < 0
    // function which returns initial pressure distribution
    std::vector<Scalar> pInit()
    {
        return pInit_;
    }

    // returns true if the current solution should be written to
    // disk (i.e. as a VTK file)
    // during initialization no output is written
    // during actual simulation output is written initially and
    // at episode/simulation end
    bool shouldWriteOutput()
    {
        return output_;
    }

    // returns true if the current solution should be written to
    // disk (i.e. as a drs file)
    bool shouldWriteRestartFile() const
    {
        return output_;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The center of the finite volume which ought to be set.
     *
     *               This function is called directly from dumux/geomechanics/el2p/localoperator.hh
     *               If it is renamed to boundaryTypesAtPos it should be adjusted there as well.
     */
    void boundaryTypesAtPos(BoundaryTypes &values, const GlobalPosition& globalPos) const
    {
        values.setAllNeumann();// means when we define flux, zero or not

        // The solid displacement normal to the lateral boundaries is fixed.
        if(onInlet_(globalPos))
        {
//             values.setDirichlet(contiEqIdx);
            if (initializationRun_ == true)
            {
//                 values.setDirichlet(contiEqIdx);
                values.setDirichlet(momentumXEqIdx);
                values.setDirichlet(momentumYEqIdx);
                  values.setDirichlet(momentumZEqIdx);
            }
//             else {
//                 values.setDirichlet(contiEqIdx);
//             }
        }

        if (leftBoundaryO_(globalPos)|| rightBoundaryO_(globalPos)|| frontBoundaryO_(globalPos)|| backBoundaryO_(globalPos)){
            values.setDirichlet(momentumXEqIdx);//this part is aplied all the time
            values.setDirichlet(momentumYEqIdx);

            if (initializationRun_ == true)
            {
//                 values.setDirichlet(contiEqIdx);
                values.setDirichlet(momentumZEqIdx);
            }
        }

        if (leftBoundaryU_(globalPos)|| frontBoundaryU_(globalPos)|| backBoundaryU_(globalPos)){
            values.setDirichlet(momentumXEqIdx);//this part is aplied all the time
            values.setDirichlet(momentumYEqIdx);

            if (initializationRun_ == true)
            {
                values.setDirichlet(contiEqIdx);
                values.setDirichlet(momentumZEqIdx);
            }
        }

        if ( rightBoundaryU_(globalPos)){
            values.setDirichlet(contiEqIdx);
            values.setDirichlet(momentumXEqIdx);//this part is aplied all the time
            values.setDirichlet(momentumYEqIdx);
            if (initializationRun_ == true)
            {
                values.setDirichlet(momentumZEqIdx);
            }
         }

        if (onLowerBoundary_(globalPos)){

             values.setDirichlet(momentumZEqIdx);

            if(initializationRun_ == true)
            {
                values.setDirichlet(momentumXEqIdx);
                values.setDirichlet(momentumYEqIdx);
           }
        }

        if (uGBoundary_(globalPos)){

            if(initializationRun_ == true)
            {
             values.setDirichlet(contiEqIdx);
           }
        }

    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex representing the "half volume on the boundary"
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const//based on vertex and element
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        dirichletAtPos(values, globalPos);
        values[contiEqIdx] = -pInit_[this->vertexMapper().index(vertex)];
     }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * This function is called directly from dumux/geomechanics/el2p/localoperator.hh
     * If it is renamed to dirichletAtPos it should be adjusted there as well.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const//based on globalpos or center of the element
    {
        values = 0.0;

          values[contiEqIdx] =pnRef_+waterDensity_*9.81*wTdepth(globalPos);
          values[momentumXEqIdx] = 0;
          values[momentumYEqIdx] = 0;
          values[momentumZEqIdx] = 0;

//         if(onInlet_(globalPos)){
//              values[contiEqIdx] = pnRef_;//pnRef_; We assume that the pressure at the surface is in equilibrium with the air pressure (from paper:Prediction of Evaporation Losses from Shallow Water Table using a numerical Model,C P Kumar, equation 12)
//         }
//         else {
//             values[contiEqIdx] =pnRef_+waterDensity_*9.81*wTdepth(globalPos);
//         }
//          if (rightBoundary_(globalPos)){
//              values[contiEqIdx] =pnRef_+waterDensity_*9.81*depth;
//         }

    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * This function is called directly from dumux/geomechanics/el2p/localoperator.hh
     * If it is renamed to neumannAtPos it should be adjusted there as well.
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
//     void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
     void solDependentNeumann(PrimaryVariables &values,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const Intersection &intersection,
                      const int scvIdx,
                      const int boundaryFaceIdx,
                      const ElementVolumeVariables &elemVolVars) const
    {
        values = 0.0;
        GlobalPosition globalPos = fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;
          Scalar ks_ = GET_RUNTIME_PARAM(TypeTag, Scalar, HydraulicParameters.ksOver);
//             Scalar ks_ = 1.39e-6;
            values[contiEqIdx] = 0;//no flow
            values[momentumXEqIdx] = 0;
            values[momentumYEqIdx] = 0;
            values[momentumZEqIdx] = 0;

//              if (rightBoundary_(globalPos)){
//              values[contiEqIdx] =ks_*waterDensity_;
//         }
        if (onInlet_(globalPos))
        {
                 values[contiEqIdx] = -ks_*waterDensity_; //  kg / (m2 * s) in 2D is kg/(m*s) inflow should be negative and outflow positive
//                 values[momentumXEqIdx] = 0;//means stress
//                 values[momentumYEqIdx] = 0;
//                 values[momentumZEqIdx] = 0;
        }
        else if (rightBoundaryO_(globalPos))
        {
             if(initializationRun_ == true)
                 values[contiEqIdx] = 0;
             else
             {
                Scalar pressure = elemVolVars[scvIdx].pressure(wPhaseIdx);
                if (pressure <0.0+eps_)
                    values[contiEqIdx] =0.0;
                else
                     values[contiEqIdx] = +ks_*waterDensity_; //  kg / (m2 * s) in 2D is kg/(m*s) inflow
//                      Scalar pn = 1e5; //khodam
//                      Scalar pc = 1e5-pwIdx; //pc [pa]
//                      Scalar sw = MaterialLaw::sw(this->spatialParams().materialLawParams(globalPos),pc);
//                      Scalar kh_=ks_*pow (sw,0.5)*pow(1-pow(1-pow(sw,1/0.27),0.27),2);//Moalem hydraulic conductivity
//                      values[contiEqIdx] = +kh_*waterDensity_;
              }
         }
    }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param values Storage for all primary variables of the initial condition
     * \param globalPos The position for which the boundary type is set
     */
//     void initialAtPos(PrimaryVariables &values,
//                  const GlobalPosition &globalPos) const
//     { initial_(values, globalPos); }; //from richards

    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
            const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        const GlobalPosition globalPos = element.geometry().corner(scvIdx);

        values = 0.0;
        //source(values, globalPos);
    }


    /*!
     * \brief Transfer episode index to spatial parameters in order to apply different
     * hydraulic parameters during pressure initialization
     */
    void preTimeStep()
    {
        this->spatialParams().setEpisode(this->timeManager().episodeIndex());
    }

    /*!
     * \brief Write mass balance information for both fluid phases
     */
    void postTimeStep()
    {
        PrimaryVariables mass;
        this->model().globalStorage(mass);
        double time = this->timeManager().time()+this->timeManager().timeStepSize();

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "TIME, MASS NPhase (kg), MASS WPhase (kg): \n"
            <<"mass: "
            <<time<< " , "
            <<mass[1] << " , "
            <<mass[0]
            <<"\n"
            <<"***************************************" <<std::endl;
        }


    }

    /*!
     * \brief Define length of next episode
     */
    void episodeEnd()
    {
        this->timeManager().startNextEpisode(episodeLength_);
        // At the end of the initializationRun
        if (this->timeManager().time() == GET_RUNTIME_PARAM(TypeTag, Scalar,TimeManager.TInitEnd))
        {
            this->timeManager().setTimeStepSize(dt_);

            this->setCoupled(true);
            // pressure field resulting from the initialization period is applied for the initial
            // and the Dirichlet boundary conditions
            this->setPressure();
            // output is written
            this->setOutput(true);
        }
    }


private:
    void initial_(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        Scalar depth = zMax(globalPos) - globalPos[2];
/*
        if (useHead)
            values[hIdx] = (-pc)/1000./ g *(100.);
        else
            values[pwIdx] = pnRef_ - pc;*/
    }
//     bool topLeftBoundary_(const GlobalPosition &globalPos) const{
//         return ((globalPos[0]<15+eps_) && (globalPos[2]>0-eps_));
//     }
//
//     bool topRightBoundary_(const GlobalPosition &globalPos) const{
//         return (globalPos[0]>30-eps_ && (globalPos[2]>-6-eps_));
//     }
//
//     bool onSlopeBoundary_(const GlobalPosition &globalPos) const{
//         if ((globalPos[0] > 15-eps_) && (globalPos[0] < 30 +eps_) && (globalPos[2] > -0.4*(globalPos[0]-15) - eps_))
//         return 1;
//         else return 0;
//     }

    bool onInlet_(const GlobalPosition &globalPos) const //is that small part of the upper boundary
    {
        return (globalPos[2]>zMax(globalPos)-eps_);
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
         return  (globalPos[2] < zMin(globalPos) + eps_);
    }

    bool leftBoundaryO_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] < 0. + eps_ && globalPos[2] > -9. - eps_);
    }

    bool leftBoundaryU_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] < 0. + eps_ && globalPos[2] < -9. + eps_);
    }

    bool rightBoundaryO_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] >45. - eps_ && globalPos[2] > -9. - eps_);
    }

    bool rightBoundaryU_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] > 45. - eps_ && globalPos[2] < -9. + eps_);
    }


    bool frontBoundaryO_(const GlobalPosition &globalPos) const
    {
        return (globalPos[1] > 20. - eps_ && globalPos[2] > -9. - eps_);
    }

    bool frontBoundaryU_(const GlobalPosition &globalPos) const
    {
        return (globalPos[1] > 20. - eps_ && globalPos[2] < -9. + eps_);
    }

    bool backBoundaryO_(const GlobalPosition &globalPos) const
    {
        return (globalPos[1] < 0. + eps_ && globalPos[2] > -9. - eps_);
    }

    bool backBoundaryU_(const GlobalPosition &globalPos) const
    {
        return (globalPos[1] < 0. + eps_ && globalPos[2] < -9. + eps_);
    }

    bool uGBoundary_(const GlobalPosition &globalPos) const
    {
        return (globalPos[2] < -9. + eps_ );
    }

//     }
    static constexpr Scalar eps_ = 1e-6;
    Scalar pnRef_;
//     Scalar depthBOR_;
    static constexpr Scalar waterDensity_ = 1000;
    Scalar episodeLength_;// = GET_RUNTIME_PARAM(TypeTag, Scalar, TimeManager.EpisodeLength);

    std::vector<Scalar> pInit_;
    GridView gridView_;
    Scalar dt_;
    Scalar temperature_;



public:
    bool initializationRun_, coupled_, output_;
    InitialStressField initialStressField_;
};


/*!
 * \ingroup ElRichardsBoxProblems
 *
 * \brief Initial conditions for momentum balance equation.
 *
 * Set initial conditions for solution of momentum balance equation
 * i.e. initialize solid displacement
 * This function is called from dumux/geomechanics/elrichards/model.hh
 *  * This function is called from dumux/geomechanics/elrichards/model.hh.
 *
 * The primary variables are initialized two times:
 * 1. before the initialization run.
 * 2. at the start of the actual simulation the solid displacement values which have
 *    changed during initialization of the pressure field are set to zero again.
 *
 */
template<class TypeTag, int dim>
class InitialDisplacement :
public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<typename GET_PROP_TYPE(TypeTag, GridView),typename GET_PROP_TYPE(TypeTag, Scalar),dim>,
    InitialDisplacement<TypeTag,dim> >
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::PDELab::AnalyticGridFunctionTraits<typename GET_PROP_TYPE(TypeTag, GridView),typename GET_PROP_TYPE(TypeTag, Scalar),dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InitialDisplacement<TypeTag,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    InitialDisplacement(const GridView & gridView) : BaseT(gridView) {}

    /*!
     * \brief Evaluate initial conditions for the momentum balance equation.
     *
     * \param position The position of the vertex
     * \param values The initial solid displacement vector at the vertex
     */
    inline void evaluateGlobal(const DomainType & position, RangeType & values) const
    {
        values = 0;
    }
};

/*!
 * \ingroup ElRichardsBoxProblems
 *
 * \brief Initial conditions for mass balance equations.
 *
 * Set initial conditions for solution of the mass balance equations
 * i.e. initialize wetting phase pressure and nonwetting phase saturation
 *
 * This function is called from dumux/geomechanics/elrichards/model.hh.
 * The primary variables are initialized two times:
 * 1. before the initialization run.
 * 2. at the start of the actual simulation applying pressure field
 *    calculated during initialization
 *
 */
template<class TypeTag>
class InitialPressSat :
public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<typename GET_PROP_TYPE(TypeTag, GridView),typename GET_PROP_TYPE(TypeTag, Scalar),1>,
    InitialPressSat<TypeTag> >
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GridView,Scalar,1> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InitialPressSat<TypeTag>> BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        // indices of the primary variables
        pwIdx = Indices::pwIdx,
       // saturationIdx = Indices::snIdx,

        dimWorld = GridView::dimensionworld
    };

    typedef typename Dune::MultipleCodimMultipleGeomTypeMapper<GridView,
                        Dune::MCMGVertexLayout> VertexMapper;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    InitialPressSat(const GridView & gridView) : BaseT(gridView) , gridView_(gridView), vertexMapper_(gridView)
    {
        // resize the pressure field vector with the number of vertices
        pInit_.resize(gridView.size(GridView::dimension));
        // fill the pressure field vector with zeros
        std::fill(pInit_.begin(), pInit_.end(), 0.0);
    }

    /*!
     * \brief Evaluate initial conditions for the mass balance equations.
     *
     * \param position The position of the vertex
     * \param values The initial pressure and saturation values at the vertex
     *
     * This function applies the pressure field pInit_ which is defined
     * in the problem.
     */
    inline void evaluateGlobal(const DomainType & position, RangeType & values) const
    {
            bool valueSet;
            valueSet = false;

            // loop over all vertices
            for (const auto& vertex : vertices(gridView_))
            {
                // get global index of current vertex
                int vIdxGlobal = vertexMapper_.index(vertex);
                Dune::FieldVector<double, dimWorld> globalPos =
                                (vertex).geometry().corner(0);

                // compare coordinates of current vertex with position coordinates
                if (globalPos[0] >= position[0] - eps_ && globalPos[0] <= position[0] + eps_
                                && globalPos[1] >= position[1] - eps_ && globalPos[1]
                                <= position[1] + eps_ && globalPos[dimWorld-1] >= position[dimWorld-1] - eps_
                                && globalPos[dimWorld-1] <= position[dimWorld-1] + eps_)
                {
                    // if coordinates are identical write the pressure value for this
                    // vertex (with index vIdxGlobal) into the values vector
                    values[pwIdx] = pInit_[vIdxGlobal];
                    // the value of this vertex is set
                    valueSet = true;
                }
            }

            // check if the pressure value for this vertex has been initialized
            if (valueSet == false)
            {
                std::cout << " pressure value not initialized correctly "
                                << std::endl;
            }

            // initialize saturation values
            //values[saturationIdx] = 0;
        }

    /*!
     * \brief Fill the vector pInit_ for initialization
     *
     * \param pInit The pressure field vector defined in the problem class
     *
     * This function is called from dumux/geomechanics/elrichards/model.hh.
     */
        void setPressure(std::vector<Scalar> pInit)
        {
            std::cout << "InitialPressSat: setPressure function called" << std::endl;
            for (const auto& vertex : vertices(gridView_))
            {
                int vIdxGlobal = vertexMapper_.index(vertex);
                pInit_[vIdxGlobal] = -pInit[vIdxGlobal];
            }
        }

private:
    static constexpr Scalar eps_ = 3e-6;
    Scalar pnRef_;
//     Scalar depthBOR_;
    std::vector<Scalar> pInit_;
    GridView gridView_;
    VertexMapper vertexMapper_;
};

} //end namespace

#endif