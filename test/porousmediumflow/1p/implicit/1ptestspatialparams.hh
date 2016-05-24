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
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>
#include <dumux/material/spatialparams/gstatrandomfield.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class OnePTestSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(OnePTestSpatialParams);

// Set properties of the porous medium
NEW_PROP_TAG(SpatialParamsRandomField);
SET_BOOL_PROP(OnePTestSpatialParams, SpatialParamsRandomField, false);
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
template<class TypeTag>
class OnePTestSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    typedef ImplicitSpatialParamsOneP<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef std::vector<Scalar> ScalarVector;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

public:
    OnePTestSpatialParams(const GridView& gridView)
        : ParentType(gridView),
          randomPermeability_(gridView.size(dim), 0.0),
          indexSet_(gridView.indexSet())
    {
        randomField_ = GET_PARAM_FROM_GROUP(TypeTag, bool, SpatialParams, RandomField);
        permeability_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Permeability);
        if(!randomField_)
            permeabilityLens_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PermeabilityLens);
        else
            initRandomField(gridView);

        lensLowerLeft_ = GET_RUNTIME_PARAM(TypeTag, GlobalPosition, SpatialParams.LensLowerLeft);
        lensUpperRight_ = GET_RUNTIME_PARAM(TypeTag, GlobalPosition, SpatialParams.LensUpperRight);
    }

    /*!
     * \brief Return the intrinsic permeability for the current sub-control volume in [m^2].
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index sub-control volume face where the
     *                      intrinsic velocity ought to be calculated.
     */
    Scalar intrinsicPermeability(const Element &element,
                                 const FVElementGeometry &fvGeometry,
                                 const int scvIdx) const
    {
        const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;

        if (isInLens_(globalPos))
        {
            if(randomField_)
                return randomPermeability_[indexSet_.index(element.template subEntity<dim> (0))];
            else
                return permeabilityLens_;
        }
        else
            return permeability_;
    }

    /*! \brief Define the porosity in [-].
   *
   * \param element The finite element
   * \param fvGeometry The finite volume geometry
   * \param scvIdx The local index of the sub-control volume where
   */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx) const
    { return 0.4; }

    /*!
     * \brief This method allows the generation of a statistical field using gstat
     *
     * \param gridView The GridView used by the problem
     */
    void initRandomField(const GridView& gridView)
    {
        std::string gStatControlFile = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Gstat, ControlFile);
        std::string gStatInputFile = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Gstat, InputFile);
        std::string outputFilePrefix = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Gstat, OutputFilePrefix);

        // create random permeability object
        GstatRandomField<GridView, Scalar> randomPermeabilityField(gridView,
                                                                   gStatControlFile,
                                                                   gStatInputFile,
                                                                   outputFilePrefix + ".dat",
                                                                   true);
        randomPermeabilityField.createPowField();
        randomPermeability_.resize(gridView.size(dim), 0.0);

        ElementIterator eItEnd = gridView.template end<0> ();
        for (ElementIterator eIt = gridView.template begin<0> (); eIt != eItEnd; ++eIt)
        {
            randomPermeability_[indexSet_.index(eIt->template subEntity<dim> (0/*scvIdx*/))]
                = randomPermeabilityField.data(*eIt);;
        }

        randomPermeabilityField.writeVtk(outputFilePrefix, "absolute permeability");
    }

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i) {
            if (globalPos[i] < lensLowerLeft_[i] || globalPos[i] > lensUpperRight_[i])
                return false;
        }
        return true;
    }

    bool randomField_;
    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar permeability_, permeabilityLens_;
    ScalarVector randomPermeability_;

    const IndexSet& indexSet_;
};
} // end namespace
#endif

