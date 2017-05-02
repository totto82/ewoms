// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief This file contains the flux module which is used for ECL problems
 *
 * This approach to fluxes is very specific to two-point flux approximation and applies
 * what the Eclipse Technical Description calls the "NEWTRAN" tramsmissibilty approach.
 */
#ifndef EWOMS_ECL_FLUX_MODULE_HH
#define EWOMS_ECL_FLUX_MODULE_HH

#include <ewoms/disc/common/fvbaseproperties.hh>
#include <ewoms/common/signum.hh>

#include <opm/common/Valgrind.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(MaterialLaw);
}

template <class TypeTag>
class EclTransIntensiveQuantities;

template <class TypeTag>
class EclTransExtensiveQuantities;

template <class TypeTag>
class EclTransBaseProblem;

/*!
 * \ingroup EclBlackOilSimulator
 * \brief Specifies a flux module which uses ECL transmissibilities.
 */
template <class TypeTag>
struct EclTransFluxModule
{
    typedef EclTransIntensiveQuantities<TypeTag> FluxIntensiveQuantities;
    typedef EclTransExtensiveQuantities<TypeTag> FluxExtensiveQuantities;
    typedef EclTransBaseProblem<TypeTag> FluxBaseProblem;

    /*!
     * \brief Register all run-time parameters for the flux module.
     */
    static void registerParameters()
    { }
};

/*!
 * \ingroup EclBlackOilSimulator
 * \brief Provides the defaults for the parameters required by the
 *        transmissibility based volume flux calculation.
 */
template <class TypeTag>
class EclTransBaseProblem
{ };

/*!
 * \ingroup EclBlackOilSimulator
 * \brief Provides the intensive quantities for the ECL flux module
 */
template <class TypeTag>
class EclTransIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
protected:
    void update_(const ElementContext& elemCtx OPM_UNUSED, unsigned dofIdx OPM_UNUSED, unsigned timeIdx OPM_UNUSED)
    { }
};

/*!
 * \ingroup EclBlackOilSimulator
 * \brief Provides the ECL flux module
 */
template <class TypeTag>
class EclTransExtensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

    enum { dimWorld = GridView::dimensionworld };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { numPhases = FluidSystem::numPhases };

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> EvalDimVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \brief Returns transmissibility for a given sub-control volume face.
     */
    Scalar transmissibility() const
    { return trans_; }

    /*!
     * \brief Return the intrinsic permeability tensor at a face [m^2]
     */
    const DimMatrix& intrinsicPermeability() const
    {
        OPM_THROW(Opm::NotImplemented,
                  "The ECL transmissibility module does not provide an explicit intrinsic permeability");
    }

    /*!
     * \brief Return the pressure potential gradient of a fluid phase at the
     *        face's integration point [Pa/m]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const EvalDimVector& potentialGrad(unsigned phaseIdx OPM_UNUSED) const
    {
        OPM_THROW(Opm::NotImplemented,
                  "The ECL transmissibility module does not provide explicit potential gradients");
    }

    /*!
     * \brief Return the gravity corrected pressure difference between the interior and
     *        the exterior of a face.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Evaluation& pressureDifference(unsigned phaseIdx) const
    { return pressureDifference_[phaseIdx]; }

    /*!
     * \brief Return the filter velocity of a fluid phase at the face's integration point
     *        [m/s]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const EvalDimVector& filterVelocity(unsigned phaseIdx OPM_UNUSED) const
    {
        OPM_THROW(Opm::NotImplemented,
                  "The ECL transmissibility module does not provide explicit filter velocities");
    }

    /*!
     * \brief Return the volume flux of a fluid phase at the face's integration point
     *        \f$[m^3/s / m^2]\f$
     *
     * This is the fluid volume of a phase per second and per square meter of face
     * area.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Evaluation& volumeFlux(unsigned phaseIdx) const
    { return volumeFlux_[phaseIdx]; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    /*!
     * \brief Returns the local index of the degree of freedom in which is
     *        in upstream direction.
     *
     * i.e., the DOF which exhibits a higher effective pressure for
     * the given phase.
     */
    unsigned upstreamIndex_(unsigned phaseIdx) const
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return upIdx_[phaseIdx];
    }

    /*!
     * \brief Returns the local index of the degree of freedom in which is
     *        in downstream direction.
     *
     * i.e., the DOF which exhibits a lower effective pressure for the
     * given phase.
     */
    unsigned downstreamIndex_(unsigned phaseIdx) const
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return dnIdx_[phaseIdx];
    }

    /*!
     * \brief Update the required gradients for interior faces
     */
    void calculateGradients_(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    {
        Opm::Valgrind::SetUndefined(*this);

        const auto& problem = elemCtx.problem();
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        interiorDofIdx_ = scvf.interiorIndex();
        exteriorDofIdx_ = scvf.exteriorIndex();
        assert(interiorDofIdx_ != exteriorDofIdx_);

        unsigned I = stencil.globalSpaceIndex(interiorDofIdx_);
        unsigned J = stencil.globalSpaceIndex(exteriorDofIdx_);
        trans_ = problem.transmissibility(elemCtx, interiorDofIdx_, exteriorDofIdx_);
        faceArea_ = scvf.area();
        thpres_ = problem.thresholdPressure(I, J);

        // estimate the gravity correction: for performance reasons we use a simplified
        // approach for this flux module that assumes that gravity is constant and always
        // acts into the downwards direction. (i.e., no centrifuge experiments, sorry.)
        Scalar g = elemCtx.problem().gravity()[dimWorld - 1];

        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx_, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx_, timeIdx);

        // this is quite hacky because the dune grid interface does not provide a
        // cellCenterDepth() method (so we ask the problem to provide it). The "good"
        // solution would be to take the Z coordinate of the element centroids, but since
        // ECL seems to like to be inconsistent on that front, it needs to be done like
        // here...
        Scalar zIn = problem.dofCenterDepth(elemCtx, interiorDofIdx_, timeIdx);
        Scalar zEx = problem.dofCenterDepth(elemCtx, exteriorDofIdx_, timeIdx);

        // the distances from the DOF's depths. (i.e., the additional depth of the
        // exterior DOF)
        Scalar distZ = zIn - zEx;

        for (unsigned phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            // check shortcut: if the mobility of the phase is zero in the interior as
            // well as the exterior DOF, we can skip looking at the phase.
#warning hack
            if (false && intQuantsIn.mobility(phaseIdx) < 1e-18 && intQuantsEx.mobility(phaseIdx) < 1e-18) {
                upIdx_[phaseIdx] = interiorDofIdx_;
                dnIdx_[phaseIdx] = exteriorDofIdx_;
                pressureDifference_[phaseIdx] = 0.0;
                volumeFlux_[phaseIdx] = 0.0;
                continue;
            }

            // do the gravity correction: compute the hydrostatic pressure for the
            // external at the depth of the internal one
            const Evaluation& rhoIn = intQuantsIn.fluidState().density(phaseIdx);
            Scalar rhoEx = Toolbox::value(intQuantsEx.fluidState().density(phaseIdx));
            Evaluation rhoAvg = (rhoIn + rhoEx)/2;

            const Evaluation& pressureInterior = intQuantsIn.fluidState().pressure(phaseIdx);
            Evaluation pressureExterior = Toolbox::value(intQuantsEx.fluidState().pressure(phaseIdx));
            pressureExterior += rhoAvg*(distZ*g);

            pressureDifference_[phaseIdx] = pressureExterior - pressureInterior;

            // decide the upstream index for the phase. for this we make sure that the
            // degree of freedom which is regarded upstream if both pressures are equal
            // is always the same: if the pressure is equal, the DOF with the lower
            // global index is regarded to be the upstream one.
            if (pressureDifference_[phaseIdx] > 0.0) {
                upIdx_[phaseIdx] = exteriorDofIdx_;
                dnIdx_[phaseIdx] = interiorDofIdx_;
            }
            else if (pressureDifference_[phaseIdx] < 0.0) {
                upIdx_[phaseIdx] = interiorDofIdx_;
                dnIdx_[phaseIdx] = exteriorDofIdx_;
            }
            else {
                // if the pressure difference is zero, we chose the DOF which has the
                // larger volume associated to it as upstream DOF
                Scalar Vin = elemCtx.dofVolume(interiorDofIdx_, /*timeIdx=*/0);
                Scalar Vex = elemCtx.dofVolume(exteriorDofIdx_, /*timeIdx=*/0);
                if (Vin > Vex) {
                    upIdx_[phaseIdx] = interiorDofIdx_;
                    dnIdx_[phaseIdx] = exteriorDofIdx_;
                }
                else if (Vin < Vex) {
                    upIdx_[phaseIdx] = exteriorDofIdx_;
                    dnIdx_[phaseIdx] = interiorDofIdx_;
                }
                else {
                    assert(Vin == Vex);
                    // if the volumes are also equal, we pick the DOF which exhibits the
                    // smaller global index
                    if (I < J) {
                        upIdx_[phaseIdx] = interiorDofIdx_;
                        dnIdx_[phaseIdx] = exteriorDofIdx_;
                    }
                    else {
                        upIdx_[phaseIdx] = exteriorDofIdx_;
                        dnIdx_[phaseIdx] = interiorDofIdx_;
                    }
                }
            }

            // apply the threshold pressure for the intersection. note that the concept
            // of threshold pressure is a quite big hack that only makes sense for ECL
            // datasets. (and even there, its physical justification is quite
            // questionable IMO.)
            if (std::abs(Toolbox::value(pressureDifference_[phaseIdx])) > thpres_) {
                if (pressureDifference_[phaseIdx] < 0.0)
                    pressureDifference_[phaseIdx] += thpres_;
                else
                    pressureDifference_[phaseIdx] -= thpres_;
            }
            else {
                pressureDifference_[phaseIdx] = 0.0;
                volumeFlux_[phaseIdx] = 0.0;
                continue;
            }

            // this is slightly hacky because in the automatic differentiation case, it
            // only works for the element centered finite volume method. for ebos this
            // does not matter, though.
            unsigned upstreamIdx = upstreamIndex_(phaseIdx);
            const auto& up = elemCtx.intensiveQuantities(upstreamIdx, timeIdx);
            if (upstreamIdx == interiorDofIdx_)
                volumeFlux_[phaseIdx] =
                    pressureDifference_[phaseIdx]*up.mobility(phaseIdx)*(-trans_/faceArea_);
            else
                volumeFlux_[phaseIdx] =
                    pressureDifference_[phaseIdx]*(Toolbox::value(up.mobility(phaseIdx))*(-trans_/faceArea_));
        }

        asImp_().updateTrans(pressureDifference_[gasPhaseIdx],
                             trans_/faceArea_,
                             elemCtx,
                             scvfIdx,
                             timeIdx);
    }

    /*!
     * \brief Update the volumetric fluxes for all fluid phases on the interior faces of the context
     */
    void calculateFluxes_(const ElementContext& elemCtx OPM_UNUSED, unsigned scvfIdx OPM_UNUSED, unsigned timeIdx OPM_UNUSED)
    { }

    // transmissibility [m^3 s]
    Scalar trans_;

    // the area of the face between the DOFs [m^2]
    Scalar faceArea_;

    // threshold pressure [Pa]
    Scalar thpres_;

    // the volumetric flux of all phases [m^3/s]
    Evaluation volumeFlux_[numPhases];

    // the difference in effective pressure between the exterior and the interior degree
    // of freedom [Pa]
    Evaluation pressureDifference_[numPhases];

    // the local indices of the interior and exterior degrees of freedom
    unsigned short interiorDofIdx_;
    unsigned short exteriorDofIdx_;
    unsigned short upIdx_[numPhases];
    unsigned short dnIdx_[numPhases];
};

} // namespace Ewoms

#endif
