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
 * \brief A point source class,
 *        i.e. sources located at a single point in space
 */

#ifndef DUMUX_POINTSOURCE_HH
#define DUMUX_POINTSOURCE_HH

#include <dumux/common/boundingboxtree.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/propertysystem.hh>

namespace Dumux
{

namespace Properties
{
// Property forward declarations
NEW_PROP_TAG(ElementVolumeVariables);
NEW_PROP_TAG(FVElementGeometry);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(ImplicitIsBox);
NEW_PROP_TAG(PointSource);
NEW_PROP_TAG(PrimaryVariables);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(TimeManager);
NEW_PROP_TAG(SubControlVolume);
} // end namespace Properties

// forward declarations
template<class TypeTag>
class PointSourceHelper;

/*!
 * \ingroup Common
 * \brief A point source base class
 */
template<class TypeTag>
class PointSource
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GridView::template Codim<0>::Entity Element;

    static const int dimworld = GridView::dimensionworld;
    typedef typename Dune::FieldVector<Scalar, dimworld> GlobalPosition;

    friend class Dumux::PointSourceHelper<TypeTag>;

public:
    //! Constructor for constant point sources
    PointSource(GlobalPosition pos, PrimaryVariables values)
      : values_(values), pos_(pos), embeddings_(1) {}

    //! Constructor for sol dependent point sources, when there is no
    // value known at the time of initialization
    PointSource(GlobalPosition pos)
      : values_(0), pos_(pos), embeddings_(1) {}

    //! Convenience += operator overload modifying only the values
    PointSource& operator+= (Scalar s)
    {
        values_ += s;
        return *this;
    }

    //! Convenience -= operator overload modifying only the values
    PointSource& operator-= (Scalar s)
    {
        values_ -= s;
        return *this;
    }

    //! Convenience *= operator overload modifying only the values
    PointSource& operator*= (Scalar s)
    {
        values_ *= s;
        return *this;
    }

    //! Convenience /= operator overload modifying only the values
    PointSource& operator/= (Scalar s)
    {
        values_ /= s;
        return *this;
    }

    //! Convenience = operator overload modifying only the values
    PointSource& operator= (const PrimaryVariables& values)
    {
        values_ = values;
        return *this;
    }

    //! Convenience = operator overload modifying only the values
    PointSource& operator= (Scalar s)
    {
        values_ = s;
        return *this;
    }

    //! return the source values
    // don't forget to call this when it's overloaded
    // in the derived class
    PrimaryVariables values() const
    {
        auto values = PrimaryVariables(values_);
        values /= embeddings_;
        return values;
    }

    //! return the source position
    const GlobalPosition& position() const
    { return pos_; }

    //! an update function called before adding the value
    // to the local residual in the problem in scvPointSources
    // to be overloaded by derived classes
    void update(const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {}

protected:
    PrimaryVariables values_; //! value of the point source for each equation
private:
    GlobalPosition pos_; //! position of the point source
    std::size_t embeddings_; //! how many SCVs the point source is associated with
};

/*!
 * \ingroup Common
 * \brief A point source class with an identifier to attach data
 */
template<class TypeTag, typename IdType>
class IdPointSource : public Dumux::PointSource<TypeTag>
{
    typedef typename Dumux::PointSource<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    static const int dimworld = GridView::dimensionworld;
    typedef typename Dune::FieldVector<Scalar, dimworld> GlobalPosition;

public:
    //! Constructor for constant point sources
    IdPointSource(GlobalPosition pos, PrimaryVariables values, IdType id)
      :  ParentType(pos, values), id_(id) {}

    //! Constructor for sol dependent point sources, when there is no
    // value known at the time of initialization
    IdPointSource(GlobalPosition pos, IdType id)
      : ParentType(pos, 0), id_(id) {}

    //! return the sources identifier
    IdType id() const
    { return id_; }

private:
    IdType id_;
};

/*!
 * \ingroup Common
 * \brief A point source class for time dependent point sources
 */
template<class TypeTag>
class TimeDependentPointSource : public Dumux::PointSource<TypeTag>
{
    typedef typename Dumux::PointSource<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GridView::template Codim<0>::Entity Element;

    static const int dimworld = GridView::dimensionworld;
    typedef typename Dune::FieldVector<Scalar, dimworld> GlobalPosition;
    // a function that takes a TimeManager and a GlobalPosition
    // and returns the PointSource values as PrimaryVariables
    typedef typename std::function<PrimaryVariables(const TimeManager&, const GlobalPosition&)> ValueFunction;

public:
    //! Constructor for constant point sources
    TimeDependentPointSource(GlobalPosition pos, PrimaryVariables values,
                             ValueFunction valueFunction)
      : ParentType(pos, values), valueFunction_(valueFunction) {}

    //! Constructor for sol dependent point sources, when there is no
    // value known at the time of initialization
    TimeDependentPointSource(GlobalPosition pos,
                             ValueFunction valueFunction)
      : ParentType(pos, 0), valueFunction_(valueFunction) {}

    //! an update function called before adding the value
    // to the local residual in the problem in scvPointSources
    // to be overloaded by derived classes
    void update(const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    { this->values_ = valueFunction_(problem.timeManager(), this->position()); }

private:
    ValueFunction valueFunction_;
};

/*!
 * \ingroup Common
 * \brief A helper class calculating a DOF-index to point source map
 */
template<class TypeTag>
class PointSourceHelper
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PointSource) PointSource;

    static const int dim = GridView::dimension;
    static const int dimworld = GridView::dimensionworld;

    typedef Dumux::BoundingBoxTree<GridView> BoundingBoxTree;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    //! calculate a DOF index to point source map from given vector of point sources
    static void computePointSourceMap(const Problem& problem,
                                      const BoundingBoxTree& boundingBoxTree,
                                      std::vector<PointSource>& sources,
                                      std::map<std::pair<unsigned int, unsigned int>, std::vector<PointSource> >& pointSourceMap)
    {
        for (auto&& source : sources)
        {
            // compute in which elements the point source falls
            std::vector<unsigned int> entities = boundingBoxTree.computeEntityCollisions(source.position());
            // split the source values equally among all concerned entities
            source.embeddings_ *= entities.size();
            // loop over all concernes elements
            for (unsigned int eIdx : entities)
            {
                if(isBox)
                {
                    // check in which subcontrolvolume(s) we are
                    // TODO mapper/problem in bboxtree would allow to make this much better
                    const auto element = boundingBoxTree.entity(eIdx);
                    auto fvGeometry = problem.model().fvGeometries(element);
                    const auto globalPos = source.position();
                    // loop over all sub control volumes and check if the point source is inside
                    std::vector<unsigned int> scvs;
                    for (auto&& scv : fvGeometry.scvs())
                    {
                        if (BoundingBoxTreeHelper<dimworld>::pointInGeometry(scv.geometry(), globalPos))
                            scvs.push_back(scv.indexInElement());
                    }
                    // for all scvs that where tested positiv add the point sources
                    // to the element/scv to point source map
                    for (auto scvIdx : scvs)
                    {
                        const auto key = std::make_pair(eIdx, scvIdx);
                        if (pointSourceMap.count(key))
                            pointSourceMap.at(key).push_back(source);
                        else
                            pointSourceMap.insert({key, {source}});
                        // split equally on the number of matched scvs
                        pointSourceMap.at(key).back().embeddings_ *= scvs.size();
                    }
                }
                else
                {
                    // add the pointsource to the DOF map
                    const auto key = std::make_pair(eIdx, /*scvIdx=*/ 0);
                    if (pointSourceMap.count(key))
                        pointSourceMap.at(key).push_back(source);
                    else
                        pointSourceMap.insert({key, {source}});
                }
            }
        }
    }
};

} // end namespace Dumux

#endif
