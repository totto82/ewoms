// $Id: richardslocalresidual.hh 3738 2010-06-15 14:01:09Z lauser $
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_RICHARDS_LOCAL_RESIDUAL_HH
#define DUMUX_RICHARDS_LOCAL_RESIDUAL_HH

#include <dumux/boxmodels/common/boxlocalresidual.hh>

#include "richardsvolumevariables.hh"

#include "richardsfluxvariables.hh"

namespace Dumux
{
/*!
 * \ingroup RichardsModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the Richards box model.
 */
template<class TypeTag>
class RichardsLocalResidual : public BoxLocalResidual<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;
    enum {
        dimWorld = GridView::dimensionworld,

        contiEqIdx = Indices::contiEqIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
    };

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    static const Scalar mobilityUpwindAlpha = GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));
    
public:
    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element for the Richards
     *        model.
     *
     * This function should not include the source and sink terms.
     */
    void computeStorage(PrimaryVariables &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VolumeVariables &volVars = 
            usePrevSol ? 
            this->prevVolVars_(scvIdx) :
            this->curVolVars_(scvIdx);

        // partial time derivative of the wetting phase mass
        result[contiEqIdx] =
            volVars.density(wPhaseIdx)
            * volVars.saturation(wPhaseIdx)
            * volVars.porosity();
    }


    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(PrimaryVariables &flux, int faceId) const
    {
        FluxVariables fluxVars(this->problem_(),
                               this->elem_(),
                               this->fvElemGeom_(),
                               faceId,
                               this->curVolVars_());

        // calculate the flux in the normal direction of the
        // current sub control volume face
        Vector tmpVec;
        fluxVars.intrinsicPermeability().mv(fluxVars.potentialGradW(),
                                            tmpVec);
        Scalar normalFlux = - (tmpVec*fluxVars.face().normal);
        
        // data attached to upstream and the downstream vertices
        // of the current phase
        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx(normalFlux));
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx(normalFlux));
        
        flux[contiEqIdx] =
            normalFlux
            *
            ((    mobilityUpwindAlpha)*up.density(wPhaseIdx)*up.mobility(wPhaseIdx) 
             +
             (1 - mobilityUpwindAlpha)*dn.density(wPhaseIdx)*dn.mobility(wPhaseIdx));

        // we need the flux from i to j instead of the other way
        flux *= -1;
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(PrimaryVariables &q, int localVertexIdx)
    {
        this->problem_().source(q,
                                this->elem_(),
                                this->fvElemGeom_(),
                                localVertexIdx);
    }
};

};

#endif
