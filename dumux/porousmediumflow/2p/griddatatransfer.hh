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
#ifndef DUMUX_TWOP_GRIDDATA_TRANSFER_HH
#define DUMUX_TWOP_GRIDDATA_TRANSFER_HH

#include <dune/grid/utility/persistentcontainer.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/adaptive/griddatatransfer.hh>

namespace Dumux {

/*!
 * \brief Class performing the transfer of data on a grid from before to after adaptation.
 */
template<class TypeTag>
class TwoPGridDataTransfer : public GridDataTransfer
{
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    struct AdaptedValues
    {
        AdaptedValues() : associatedMass(0.0) {}
        ElementSolution u;
        int count = 0;
        PrimaryVariables associatedMass;
        bool wasLeaf = false;
    };

    using PersistentContainer = Dune::PersistentContainer<Grid, AdaptedValues>;

    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;
    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, DiscretizationMethod) == DiscretizationMethods::Box;

    //! export some indices
    enum {
        // index of saturation in primary variables
        saturationIdx = Indices::saturationIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        // formulations
        pwsn = Indices::pwsn,
        pnsw = Indices::pnsw,

        // the formulation that is actually used
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    //! This won't work (mass conservative) for compressible fluids
    static_assert(!FluidSystem::isCompressible(wPhaseIdx)
                  && !FluidSystem::isCompressible(nPhaseIdx),
                  "This adaption helper is only mass conservative for incompressible fluids!");

    //! check if the used formulation is implemented here
    static_assert(formulation == pwsn || formulation == pnsw, "Chosen formulation not known to the TwoPGridDataTransfer");

public:
    /*! \brief Constructor
     *
     *  \param problem The DuMuX problem to be solved
     *  \param fvGridGeometry The finite volume grid geometry
     *  \param gridVariables The secondary variables on the grid
     *  \param sol The solution (primary variables) on the grid
     */
    TwoPGridDataTransfer(std::shared_ptr<const Problem> problem,
                         std::shared_ptr<FVGridGeometry> fvGridGeometry,
                         std::shared_ptr<const GridVariables> gridVariables,
                         SolutionVector& sol)
    : GridDataTransfer()
    , problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , gridVariables_(gridVariables)
    , sol_(sol)
    , adaptionMap_(fvGridGeometry->gridView().grid(), 0)
    {}

    /*! \brief Stores primary variables and additional data
     *
     *  To reconstruct the solution in father elements, problem properties might
     *  need to be accessed. From upper level on downwards, the old solution is stored
     *  into a container object, before the grid is adapted. Father elements hold averaged
     *  information from the son cells for the case of the sons being coarsened.
     */
    void store() override
    {
        adaptionMap_.resize();

        const auto& grid = fvGridGeometry_->gridView().grid();
        for (auto level = grid.maxLevel(); level >= 0; level--)
        {
            for (const auto& element : elements(grid.levelGridView(level)))
            {
                //! get map entry
                auto& adaptedValues = adaptionMap_[element];

                //! put values in the map for leaf elements
                if (element.isLeaf())
                {
                    auto fvGeometry = localView(*fvGridGeometry_);
                    fvGeometry.bindElement(element);

                    //! store current element solution
                    adaptedValues.u = ElementSolution(element, sol_, *fvGridGeometry_);

                    //! compute mass in the scvs
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        VolumeVariables volVars;
                        volVars.update(adaptedValues.u, *problem_, element, scv);

                        const auto poreVolume = scv.volume()*volVars.porosity();
                        adaptedValues.associatedMass[nPhaseIdx] += poreVolume * volVars.density(nPhaseIdx) * volVars.saturation(nPhaseIdx);
                        adaptedValues.associatedMass[wPhaseIdx] += poreVolume * volVars.density(wPhaseIdx) * volVars.saturation(wPhaseIdx);
                    }

                    //! leaf elements always start with count = 1
                    adaptedValues.count = 1;
                    adaptedValues.wasLeaf = true;
                }
                //! Average in father elements
                if (element.level() > 0)
                {
                    auto& adaptedValuesFather = adaptionMap_[element.father()];
                    // For some grids the father element is identical to the son element.
                    // In that case averaging is not necessary.
                    if(&adaptedValues != &adaptedValuesFather)
                        storeAdaptionValues(adaptedValues, adaptedValuesFather);
                }

                //! The vertices of the non-leaf elements exist on the leaf as well
                //! This element solution constructor uses the vertex mapper to obtain
                //! the privars at the vertices, thus, this works for non-leaf elements!
                if(isBox && !element.isLeaf())
                    adaptedValues.u = ElementSolution(element, sol_, *fvGridGeometry_);
            }
        }
    }

    /*! \brief Reconstruct missing primary variables (where elements are created/deleted)
     *
     *  To reconstruct the solution in father elements, problem properties might
     *  need to be accessed.
     *  Starting from the lowest level, the old solution is mapped on the new grid:
     *  Where coarsened, new cells get information from old father element.
     *  Where refined, a new solution is reconstructed from the old father cell,
     *  and then a new son is created. That is then stored into the general data
     *  structure (AdaptedValues).
     */
    void reconstruct() override
    {
        //! resize stuff (grid might have changed)
        adaptionMap_.resize();
        fvGridGeometry_->update();
        sol_.resize(fvGridGeometry_->numDofs());

        //! vectors storing the mass associated with each vertex, when using the box method
        std::vector<Scalar> massCoeff;
        std::vector<Scalar> associatedMass;

        if(isBox)
        {
            massCoeff.resize(fvGridGeometry_->numDofs(), 0.0);
            associatedMass.resize(fvGridGeometry_->numDofs(), 0.0);
        }

        //! iterate over leaf and reconstruct the solution
        for (const auto& element : elements(fvGridGeometry_->gridView().grid().leafGridView(), Dune::Partitions::interior))
        {
            if (!element.isNew())
            {
                const auto& adaptedValues = adaptionMap_[element];

                auto fvGeometry = localView(*fvGridGeometry_);
                fvGeometry.bindElement(element);

                //! obtain element solution from map (divide by count!)
                auto elemSol = adaptedValues.u;
                if (!isBox)
                    elemSol[0] /= adaptedValues.count;

                const auto elementVolume = element.geometry().volume();
                for (const auto& scv : scvs(fvGeometry))
                {
                    VolumeVariables volVars;
                    volVars.update(elemSol, *problem_, element, scv);

                    //! write solution at dof in current solution vector
                    sol_[scv.dofIndex()] = elemSol[scv.indexInElement()];

                    const auto dofIdxGlobal = scv.dofIndex();
                    //! For cc schemes, overwrite the saturation by a mass conservative one here
                    if (!isBox)
                    {
                        //! only recalculate the saturations if element hasn't been leaf before adaptation
                        if (!adaptedValues.wasLeaf)
                        {
                            if (formulation == pwsn)
                            {
                                sol_[dofIdxGlobal][saturationIdx] = adaptedValues.associatedMass[nPhaseIdx];
                                sol_[dofIdxGlobal][saturationIdx] /= elementVolume * volVars.density(nPhaseIdx) * volVars.porosity();
                            }
                            else if (formulation == pnsw)
                            {
                                sol_[dofIdxGlobal][saturationIdx] = adaptedValues.associatedMass[wPhaseIdx];
                                sol_[dofIdxGlobal][saturationIdx] /= elementVolume * volVars.density(wPhaseIdx) * volVars.porosity();
                            }
                        }
                    }

                    //! For the box scheme, add mass & mass coefficient to container (saturations are recalculated at the end)
                    else
                    {
                        const auto scvVolume = scv.volume();
                        if (formulation == pwsn)
                        {
                            massCoeff[dofIdxGlobal] += scvVolume * volVars.density(nPhaseIdx) * volVars.porosity();
                            associatedMass[dofIdxGlobal] += scvVolume / elementVolume * adaptedValues.associatedMass[nPhaseIdx];
                        }
                        else if (formulation == pnsw)
                        {
                            massCoeff[dofIdxGlobal] += scvVolume * volVars.density(wPhaseIdx) * volVars.porosity();
                            associatedMass[dofIdxGlobal] += scvVolume / elementVolume * adaptedValues.associatedMass[wPhaseIdx];
                        }
                    }
                }
            }
            else
            {
                //! value is not in map, interpolate from father element
                assert(element.hasFather() && "new element does not have a father element!");

                //! find the ancestor element that existed on the old grid already
                auto fatherElement = element.father();
                while(fatherElement.isNew() && fatherElement.level() > 0)
                    fatherElement = fatherElement.father();

                if(!isBox)
                {
                    const auto& adaptedValuesFather = adaptionMap_[fatherElement];

                    //! obtain the mass contained in father
                    Scalar massFather = 0.0;
                    if (formulation == pwsn)
                        massFather = adaptedValuesFather.associatedMass[nPhaseIdx];
                    else if (formulation == pnsw)
                        massFather = adaptedValuesFather.associatedMass[wPhaseIdx];

                    // obtain the element solution through the father
                    auto elemSolSon = adaptedValuesFather.u;
                    elemSolSon[0] /= adaptedValuesFather.count;

                    auto fvGeometry = localView(*fvGridGeometry_);
                    fvGeometry.bindElement(element);

                    for (const auto& scv : scvs(fvGeometry))
                    {
                        VolumeVariables volVars;
                        volVars.update(elemSolSon, *problem_, element, scv);

                        //! store constructed values of son in the current solution
                        sol_[scv.dofIndex()] = elemSolSon[0];

                        //! overwrite the saturation by a mass conservative one here
                        Scalar massCoeffSon = 0.0;
                        if (formulation == pwsn)
                            massCoeffSon = scv.volume() * volVars.density(nPhaseIdx) * volVars.porosity();
                        else if (formulation == pnsw)
                            massCoeffSon = scv.volume() * volVars.density(wPhaseIdx) * volVars.porosity();
                        sol_[scv.dofIndex()][saturationIdx] = (scv.volume() / fatherElement.geometry().volume() * massFather)/massCoeffSon;
                    }
                }
                else
                {
                    auto& adaptedValuesFather = adaptionMap_[fatherElement];

                    auto fvGeometry = localView(*fvGridGeometry_);
                    fvGeometry.bindElement(element);

                    //! interpolate solution in the father to the vertices of the new son
                    ElementSolution elemSolSon(element, sol_, *fvGridGeometry_);
                    const auto fatherGeometry = fatherElement.geometry();
                    for (const auto& scv : scvs(fvGeometry))
                        elemSolSon[scv.indexInElement()] = evalSolution(fatherElement,
                                                                        fatherGeometry,
                                                                        adaptedValuesFather.u,
                                                                        scv.dofPosition());

                    //! compute mass & mass coeffients for the scvs (saturations are recalculated at the end)
                    const auto fatherElementVolume = fatherGeometry.volume();
                    for (const auto& scv : scvs(fvGeometry))
                    {
                        VolumeVariables volVars;
                        volVars.update(elemSolSon, *problem_, element, scv);

                        const auto dofIdxGlobal = scv.dofIndex();
                        const auto scvVolume = scv.volume();
                        if (int(formulation) == pwsn)
                        {
                            massCoeff[dofIdxGlobal] += scvVolume * volVars.density(nPhaseIdx) * volVars.porosity();
                            associatedMass[dofIdxGlobal] += scvVolume / fatherElementVolume * adaptedValuesFather.associatedMass[nPhaseIdx];
                        }
                        else if (int(formulation) == pnsw)
                        {
                            massCoeff[dofIdxGlobal] += scvVolume * volVars.density(wPhaseIdx) * volVars.porosity();
                            associatedMass[dofIdxGlobal] += scvVolume / fatherElementVolume * adaptedValuesFather.associatedMass[wPhaseIdx];
                        }

                        //! store constructed (pressure) values of son in the current solution (saturation comes later)
                        sol_[dofIdxGlobal] = elemSolSon[scv.indexInElement()];
                    }
                }
            }
        }

        if(isBox)
        {
            for(std::size_t dofIdxGlobal = 0; dofIdxGlobal < fvGridGeometry_->numDofs(); dofIdxGlobal++)
                sol_[dofIdxGlobal][saturationIdx] = associatedMass[dofIdxGlobal] / massCoeff[dofIdxGlobal];
        }

        //! reset entries in adaptation map
        adaptionMap_.resize( typename PersistentContainer::Value() );
        adaptionMap_.shrinkToFit();
        adaptionMap_.fill( typename PersistentContainer::Value() );

//! TODO: fix adaptive simulations in parallel
//#if HAVE_MPI
//        // communicate ghost data
//        using SolutionTypes = typename GET_PROP(TypeTag, SolutionTypes);
//        using ElementMapper = typename SolutionTypes::ElementMapper;
//        using DataHandle = VectorExchange<ElementMapper, std::vector<CellData> >;
//        DataHandle dataHandle(problem.elementMapper(), this->cellDataGlobal());
//        problem.gridView().template communicate<DataHandle>(dataHandle,
//                                                            Dune::InteriorBorder_All_Interface,
//                                                            Dune::ForwardCommunication);
//#endif
    }

  private:

    /*! \brief Stores sons entries into father element for averaging
     *
     *  Sum up the adaptedValues (sons values) into father element. We store from leaf
     *  upwards, so sons are stored first, then cells on the next leaf (=fathers)
     *  can be averaged.
     *
     *  \param adaptedValues Container for model-specific values to be adapted
     *  \param adaptedValuesFather Values to be adapted of father cell
     */
    static void storeAdaptionValues(AdaptedValues& adaptedValues,
                                    AdaptedValues& adaptedValuesFather)
    {
        //! Add associated mass of the child to the one of the father
        adaptedValuesFather.associatedMass += adaptedValues.associatedMass;

        if(!isBox)
        {
            //! add the child's primary variables to the ones of father
            //! we have to divide the child's ones in case it was composed
            //! of several children as well!
            auto values = adaptedValues.u[0];
            values /= adaptedValues.count;
            adaptedValuesFather.u[0] += values;

            //! keep track of the number of children that composed this father
            adaptedValuesFather.count += 1;

            //! A father element is never leaf
            adaptedValuesFather.wasLeaf = false;
        }
        else
        {
            //! For the box scheme, scaling of primary variables by count is obsolete
            //! Thus, we always want count = 1
            adaptedValuesFather.count = 1;

            //! A father element is never leaf
            adaptedValuesFather.wasLeaf = false;
        }
    }

    std::shared_ptr<const Problem> problem_;
    std::shared_ptr<FVGridGeometry> fvGridGeometry_;
    std::shared_ptr<const GridVariables> gridVariables_;
    SolutionVector& sol_;
    PersistentContainer adaptionMap_;
};

} // end namespace Dumux

#endif /* DUMUX_TWOP_GRIDDATA_TRANSFER_HH */