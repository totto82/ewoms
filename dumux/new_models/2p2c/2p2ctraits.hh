/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/

#ifndef DUMUX_2P2CTRAITS_HH
#define DUMUX_2P2CTRAITS_HH

#include <dumux/new_models/boxscheme/p1boxtraits.hh>

#include "2p2cnewtoncontroller.hh"

namespace Dune
{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class TwoPTwoCBoxModel;

template<class TypeTag>
class TwoPTwoCBoxJacobian;

template <class Scalar>
class TwoPTwoCPnSwTraits;

template <class TwoPTwoCTraits, class Problem>
class TwoPTwoCVertexData;

template <class TwoPTwoCTraits, class ProblemT>
class TwoPTwoCElementData;

template <class TwoPTwoCTraits, class ProblemT, class VertexData>
class TwoPTwoCFluxData;

////////////////////////////////
// properties
////////////////////////////////

namespace Properties
{
//! set the number of equations to 2
SET_PROP_INT(BoxTwoPTwoC, NumEq, 2);

//! Use the pw-Sn formulation by default
SET_PROP(BoxTwoPTwoC, TwoPTwoCTraits)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
  
public:
    typedef TwoPTwoCPnSwTraits<Scalar> type;
};

//! Use the 2p2c local jacobian operator for the 2p2c model
SET_PROP(BoxTwoPTwoC, LocalJacobian)
{
    typedef TwoPTwoCBoxJacobian<TypeTag> type;
};

//! Use the 2p2c specific newton controller for the 2p2c model
SET_PROP(BoxTwoPTwoC, NewtonController)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod))  NewtonMethod;
  
public:
    typedef TwoPTwoCNewtonController<NewtonMethod> type;
};

//! the Model property
SET_PROP(BoxTwoPTwoC, Model)
{
    typedef Dune::TwoPTwoCBoxModel<TypeTag> type;
};

//! the VertexData property
SET_PROP(BoxTwoPTwoC, VertexData)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))        Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCTraits)) Traits;
    
public:
    typedef TwoPTwoCVertexData<Traits, Problem> type;
};

//! the ElementData property
SET_PROP(BoxTwoPTwoC, ElementData)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))        Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCTraits)) Traits;
    
public:
    typedef TwoPTwoCElementData<Traits, Problem> type;
};

//! the FluxData property
SET_PROP(BoxTwoPTwoC, FluxData)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))    Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCTraits)) Traits;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData)) VertexData;
    
public:
    typedef TwoPTwoCFluxData<Traits, Problem, VertexData> type;
};

}

/*!
 * \brief Generic 2P-2C traits.
 *
 * This class specifies the exact behaviour of the two-phase
 * two-component model. By using a different traits class for the
 * model, the model can change its behaviour considerably.
 */
template <class Scalar>
class TwoPTwoCBaseTraits
{
public:
    static const int numEq = 2;         //!< Number of primary variables / equations
    static const int numPhases = 2;     //!< Number of fluid phases
    static const int numComponents = 2; //!< Number of fluid components within a phase
    
    // Primary variable indices
    static const int pressureIdx = 0;     //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx   = 1;     //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    // present phases
    static const int nPhaseOnly = 0; //!< Only the non-wetting phase is present
    static const int wPhaseOnly = 1; //!< Only the wetting phase is present
    static const int bothPhases = 2;  //!< Both phases are present
    
    // formulation
    static const int pWsN = 0; //!< Pw and Sn as primary variables
    static const int pNsW = 1;  //!< Pn and Sw as primary variables

    // Upwind parameter
    static const Scalar upwindAlpha = 1.0; //!< Upwind parameter. 1.0 means fully the upstream.
};

/*!
 * \brief The traits for the pw-Sn formulation of the 2p2c model.
 */
template <class Scalar>
class TwoPTwoCPwSnTraits : public TwoPTwoCBaseTraits<Scalar>
{
    typedef TwoPTwoCBaseTraits<Scalar>     ParentT;

public:
    static const int formulation = ParentT::pWsN; //!< Formulation to use
    static const int wPhase      = 0;    //!< Index of the wetting phase in a phase vector
    static const int nPhase      = 1;    //!< Index of the non-wetting phase in a phase vector

    static const int wComp       = 0; //!< Index of the wetting component in a component vector
    static const int nComp       = 1; //!< Index of the non-wetting component in a compent vector
};

/*!
 * \brief The traits for the pn-Sw formulation of the 2p2c model.
 */
template <class Scalar>
class TwoPTwoCPnSwTraits : public TwoPTwoCBaseTraits<Scalar>
{
    typedef TwoPTwoCBaseTraits<Scalar>     ParentT;

public:
    static const int formulation = ParentT::pNsW; //!< Formulation to use
    static const int wPhase      = 1;    //!< Index of the wetting phase in a phase vector
    static const int nPhase      = 0;    //!< Index of the non-wetting phase in a phase vector

    static const int wComp       = 1; //!< Index of the wetting component in a solution vector
    static const int nComp       = 0; //!< Index of the non-wetting component in a solution vector
};

/*!
 * \brief The traits for the non-isothermal 2p2c model formulation.
 */
template <class Scalar,
          class BaseTraits = TwoPTwoCPwSnTraits<Scalar> >
class TwoPTwoCNITraits : public BaseTraits
{
public:
    static const int numEq = 3;  //!< Override the mumber of primary variables: We also have temperature
    // Primary variable indices
    static const int temperatureIdx = 2; //! The index for temperature in solution vectors.
};

}

#endif
