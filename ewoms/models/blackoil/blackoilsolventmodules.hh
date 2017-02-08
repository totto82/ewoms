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
 * \brief Contains the classes required to extend the black-oil model by solvents.
 */
#ifndef EWOMS_BLACK_OIL_SOLVENT_MODULE_HH
#define EWOMS_BLACK_OIL_SOLVENT_MODULE_HH

#include "blackoilproperties.hh"
#include <ewoms/io/vtkblackoilsolventmodule.hh>

#include <opm/material/fluidsystems/blackoilpvt/SolventPvt.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

#include <opm/common/Valgrind.hpp>
#include <opm/common/Unused.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

#include <string>

namespace Ewoms {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by solvents.
 */
template <class TypeTag, bool enableSolventV = GET_PROP_VALUE(TypeTag, EnableSolvent)>
class BlackOilSolventModule
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Opm::SolventPvt<Scalar> SolventPvt;

    static constexpr unsigned solventSaturationIdx = Indices::solventSaturationIdx;
    static constexpr unsigned contiSolventEqIdx = Indices::contiSolventEqIdx;
    static constexpr unsigned enableSolvent = enableSolventV;
    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:
#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize all internal data structures needed by the solvent module
     */
    static void initFromDeck(const Opm::Deck& deck, const Opm::EclipseState& eclState)
    {
        // some sanity checks: if solvents are enabled, the SOLVENT keyword must be
        // present, if solvents are disabled the keyword must not be present.
        if (enableSolvent && !deck.hasKeyword("SOLVENT")) {
            OPM_THROW(std::runtime_error,
                      "Non-trivial solvent treatment requested at compile time, but "
                      "the deck does not contain the SOLVENT keyword");
        }
        else if (!enableSolvent && deck.hasKeyword("SOLVENT")) {
            OPM_THROW(std::runtime_error,
                      "Solvent treatment disabled at compile time, but the deck "
                      "contains the SOLVENT keyword");
        }

        if (!deck.hasKeyword("SOLVENT"))
            return; // solvent treatment is supposed to be disabled

        solventPvt_.initFromDeck(deck, eclState);
    }

    // TODO(?): full init without opm-parser
#endif

    /*!
     * \brief Register all run-time parameters for the black-oil solvent module.
     */
    static void registerParameters()
    {
        if (!enableSolvent)
            // solvents have disabled at compile time
            return;

        VtkBlackOilSolventModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all solvent specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enableSolvent)
            // solvents have disabled at compile time
            return;

        model.addOutputModule(new Ewoms::VtkBlackOilSolventModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if (!enableSolvent)
            // solvents have disabled at compile time
            return false;

        return pvIdx == solventSaturationIdx;
    }

    static std::string primaryVarName(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        return "saturation_solvent";
    }

    static Scalar primaryVarWeight(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableSolvent)
            return false;

        return eqIdx == contiSolventEqIdx;
    }

    static std::string eqName(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        return "conti^solvent";
    }

    static Scalar eqWeight(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enableSolvent)
            return;

        storage[contiSolventEqIdx] +=
            Toolbox::template decay<LhsEval>(intQuants.porosity())
            * Toolbox::template decay<LhsEval>(intQuants.solventSaturation())
            * Toolbox::template decay<LhsEval>(intQuants.solventDensity());
    }

    template <class UpEval>
    static void addFlux(RateVector& flux,
                        const ExtensiveQuantities& extQuants,
                        const IntensiveQuantities& upQuants)
    {
        if (!enableSolvent)
            return;

        const auto& volFlux = extQuants.solventVolumeFlux();

        const UpEval& upRhoSol =
            Toolbox::template decay<UpEval>(upQuants.solventDensity());

        flux[contiSolventEqIdx] += volFlux*upRhoSol;
    }

    /*!
     * \brief Assign the solvent specific primary variables to a PrimaryVariables object
     */
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  Scalar solventSaturation)
    {
        if (!enableSolvent)
            return;

        priVars[solventSaturationIdx] = solventSaturation;
    }

    /*!
     * \brief Do a Newton-Raphson update the primary variables of the solvents.
     */
    static void updatePrimaryVars(PrimaryVariables& newPv,
                                  const PrimaryVariables& oldPv,
                                  const EqVector& delta)
    {
        if (!enableSolvent)
            return;

        // do a plain unchopped Newton update
        newPv[solventSaturationIdx] = oldPv[solventSaturationIdx] - delta[solventSaturationIdx];
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv OPM_UNUSED,
                                     const EqVector& delta OPM_UNUSED)
    {
        // do not consider consider the cange of solvent primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    /*!
     * \brief Return how much a residual is considered an error
     */
    static Scalar computeResidualError(const EqVector& resid)
    {
        // do not weight the residual of solvents when it comes to convergence
        return std::abs(Toolbox::scalarValue(resid[contiSolventEqIdx]));
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if (!enableSolvent)
            return;

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
        unsigned dofIdx = model.dofMapper().index(dof);
#else
        unsigned dofIdx = model.dofMapper().map(dof);
#endif

        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
        outstream << priVars[solventSaturationIdx];
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if (!enableSolvent)
            return;

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
        unsigned dofIdx = model.dofMapper().index(dof);
#else
        unsigned dofIdx = model.dofMapper().map(dof);
#endif

        PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        instream >> priVars0[solventSaturationIdx];

        // set the primary variables for the beginning of the current time step.
        priVars1 = priVars0[solventSaturationIdx];
    }

    static const SolventPvt& solventPvt()
    { return solventPvt_; }

private:
    static SolventPvt solventPvt_;
};

template <class TypeTag, bool enableSolventV>
typename BlackOilSolventModule<TypeTag, enableSolventV>::SolventPvt
BlackOilSolventModule<TypeTag, enableSolventV>::solventPvt_;

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilSolventIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        solvents extension of the black-oil model.
 */
template <class TypeTag, bool enableSolventV = GET_PROP_VALUE(TypeTag, EnableSolvent)>
class BlackOilSolventIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef BlackOilSolventModule<TypeTag> SolventModule;

    static constexpr int solventSaturationIdx = Indices::solventSaturationIdx;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;

public:
    /*!
     * \brief Called before the saturation functions are doing their magic
     *
     * At this point, the saturations of the fluid state correspond to those if the phases
     * were pure hydrocarbons.
     */
    void preSatFuncUpdate_(const ElementContext& elemCtx,
                           unsigned scvIdx,
                           unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(scvIdx, timeIdx);
        auto& fs = asImp_().fluidState_;
        solventSaturation_ = priVars.makeEvaluation(solventSaturationIdx, timeIdx);
        hydrocarbonSaturation_ = fs.saturation(gasPhaseIdx);

        // make the saturation of the gas phase which is used by the saturation functions
        // the sum of the solvent "saturation" and the saturation the hydrocarbon gas.
        fs.setSaturation(gasPhaseIdx, hydrocarbonSaturation_ + solventSaturation_);
    }

    /*!
     * \brief Called after the saturation functions have been doing their magic
     *
     * At this point, the pressures of the fluid state are final. After this function,
     * all saturations and relative permeabilities must be final. (i.e., the "hydrocarbon
     * saturations".)
     */
    void postSatFuncUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                            unsigned scvIdx OPM_UNUSED,
                            unsigned timeIdx OPM_UNUSED)
    {
        // revert the gas "saturation" of the fluid state back to the saturation of the
        // hydrocarbon gas.
        auto& fs = asImp_().fluidState_;
        fs.setSaturation(gasPhaseIdx, hydrocarbonSaturation_);

        // compute the mobility of the solvent "phase". this only covers the "immiscible"
        // case.
        solventMobility_ = 0.0;
        Evaluation Stot = hydrocarbonSaturation_ + solventSaturation_;
        if (Stot > 0.0) {
            Evaluation Fhydgas = hydrocarbonSaturation_/Stot;
            Evaluation Fsolgas = solventSaturation_/Stot;

#warning TODO: implement SSFN
            Evaluation& krg = asImp_().mobility_[gasPhaseIdx];
            solventMobility_ = krg * Fsolgas;
            krg *= Fhydgas;
        }
    }

    /*!
     * \brief Update the intensive PVT properties needed to handle solvents from the
     *        primary variables.
     *
     * At this point the pressures and saturations of the fluid state are correct.
     */
    void pvtUpdate_(const ElementContext& elemCtx,
                    unsigned scvIdx,
                    unsigned timeIdx)
    {
        const auto& iq = elemCtx.intensiveQuantities(scvIdx, timeIdx);
        const auto& fs = iq.fluidState();
        const auto& solventPvt = SolventModule::solventPvt();

        unsigned pvtRegionIdx = iq.pvtRegionIndex();
        Scalar rhoRef = solventPvt.referenceDensity(pvtRegionIdx);
        const Evaluation& T = fs.temperature(gasPhaseIdx);
        const Evaluation& p = fs.pressure(gasPhaseIdx);
        const Evaluation& b = solventPvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p);

#warning TODO: implement "miscibility"
        solventDensity_ = b*rhoRef;
        solventMobility_ /= solventPvt.viscosity(pvtRegionIdx, T, p);
    }

    const Evaluation& solventSaturation() const
    { return solventSaturation_; }

    const Evaluation& solventDensity() const
    { return solventDensity_; }

    const Evaluation& solventMobility() const
    { return solventMobility_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation hydrocarbonSaturation_;
    Evaluation solventSaturation_;
    Evaluation solventDensity_;
    Evaluation solventMobility_;
};

template <class TypeTag>
class BlackOilSolventIntensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    void preSatFuncUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                           unsigned scvIdx OPM_UNUSED,
                           unsigned timeIdx OPM_UNUSED)
    { }

    void postSatFuncUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                            unsigned scvIdx OPM_UNUSED,
                            unsigned timeIdx OPM_UNUSED)
    { }

    void pvtUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                    unsigned scvIdx OPM_UNUSED,
                    unsigned timeIdx OPM_UNUSED)
    { }

    const Evaluation& solventSaturation() const
    { OPM_THROW(std::runtime_error, "solventSaturation() called but solvents are disabled"); }

    const Evaluation& solventDensity() const
    { OPM_THROW(std::runtime_error, "solventDensity() called but solvents are disabled"); }

    const Evaluation& solventViscosity() const
    { OPM_THROW(std::runtime_error, "solventViscosity() called but solvents are disabled"); }
};

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilSolventExtensiveQuantities
 *
 * \brief Provides the solvent specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableSolventV = GET_PROP_VALUE(TypeTag, EnableSolvent)>
class BlackOilSolventExtensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;

public:
    /*!
     * \brief Generic update method which uses the pressure potential gradients
     */
    void update(const ElementContext& elemCtx,
                unsigned scvfIdx,
                unsigned timeIdx)
    {
#warning HACK
#if 0
        const auto& scvf = elemCtx.stencil(timeIdx).interiorFace(scvfIdx);
        const ExtensiveQuantities& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned upIdx = extQuants.upstreamIndex(gasPhaseIdx);
        const IntensiveQuantities& up = elemCtx.intensiveQuantities(upIdx, timeIdx);

        Evaluation pGradNormal = 0.0;
        const auto& pGrad = asImp_().potentialGrad(gasPhaseIdx);
        const auto& normal = scvf.normal();
        for (unsigned dimIdx = 0; dimIdx < normal.size(); ++dimIdx)
            pGradNormal += pGrad[dimIdx]*normal[dimIdx];

        solventVolumeFlux_ = pGradNormal*up.solventMobility();
#endif
    }

    /*!
     * \brief Update method which assumes that TPFA, pressure differences and
     *        transmissibilities are used.
     */
    void updateTrans(const Evaluation& pDiffGas,
                     Scalar trans,
                     const ElementContext& elemCtx,
                     unsigned scvfIdx,
                     unsigned timeIdx)
    {
        const ExtensiveQuantities& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned upIdx = extQuants.upstreamIndex(gasPhaseIdx);
        const IntensiveQuantities& up = elemCtx.intensiveQuantities(upIdx, timeIdx);

        solventVolumeFlux_ = pDiffGas*up.solventMobility()*(-trans);
    }

    const Evaluation& solventVolumeFlux() const
    { return solventVolumeFlux_; }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation solventVolumeFlux_;
};

template <class TypeTag>
class BlackOilSolventExtensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;

public:
    void update(const ElementContext& context OPM_UNUSED,
                unsigned scvfIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED)
    { }

    const Evaluation& solventVolumeFlux() const
    { OPM_THROW(std::runtime_error, "solventVolumeFlux() called but solvents are disabled"); }
};

} // namespace Ewoms

#endif
