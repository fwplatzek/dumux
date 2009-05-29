// $Id$

#ifndef FLUXFUNCTION_HH
#define FLUXFUNCTION_HH

//! \ingroup transport
//! \defgroup flux function transport
/**
 * @file
 * @brief  Base class for defining the flux function of an advection-diffusion equation
 * @author Yufei Cao
 */

namespace Dune
{
//! \ingroup transport
//! Compute the flux function of the transport equation
template<class GridView, class Scalar>
class FluxFunction {
public:
    /*! \brief Realizes the numerical flux function.
     *
     *  \param element cell I
     *  \param saturationW the saturation of the wetting phase
     *  \param T temperature
     *  \param p pressure
     */
private:
    enum{dim = GridView::dimension, dimWorld = GridView::dimensionworld};
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    virtual Scalar operator() (Scalar saturationW, const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, Scalar T=283.15, Scalar p=1e5) const = 0;

    virtual ~FluxFunction()
    { }
};
}

#endif
