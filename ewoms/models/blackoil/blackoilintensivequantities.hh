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
 * \copydoc Ewoms::BlackOilIntensiveQuantities
 */
#ifndef EWOMS_BLACK_OIL_INTENSIVE_QUANTITIES_HH
#define EWOMS_BLACK_OIL_INTENSIVE_QUANTITIES_HH

#include "blackoilproperties.hh"
#include "blackoilfluidstate.hh"
#include "blackoilsolventmodules.hh"
#include "blackoilpolymermodules.hh"


#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/common/Valgrind.hpp>

#include <dune/common/fmatrix.hh>

#include <cstring>
#include <utility>

namespace Ewoms {
/*!
 * \ingroup BlackOilModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the black-oil model.
 */
template <class TypeTag>
class BlackOilIntensiveQuantities
    : public GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities)
    , public GET_PROP_TYPE(TypeTag, FluxModule)::FluxIntensiveQuantities
    , public BlackOilSolventIntensiveQuantities<TypeTag>
    , public BlackOilPolymerIntensiveQuantities<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluxModule) FluxModule;

    enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };
    enum { enablePolymer = GET_PROP_VALUE(TypeTag, EnablePolymer) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { dimWorld = GridView::dimensionworld };

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef typename FluxModule::FluxIntensiveQuantities FluxIntensiveQuantities;
    typedef Ewoms::BlackOilFluidState<TypeTag> FluidState;

public:
    BlackOilIntensiveQuantities()
    {
        fluidState_.setRs(0.0);
        fluidState_.setRv(0.0);
    }

    BlackOilIntensiveQuantities(const BlackOilIntensiveQuantities& other)
        : ParentType()
    { std::memcpy(this, &other, sizeof(other)); }

    BlackOilIntensiveQuantities& operator=(const BlackOilIntensiveQuantities& other)
    { std::memcpy(this, &other, sizeof(other)); return *this; }

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);

        //fluidState_.setTemperature(elemCtx.problem().temperature(elemCtx, dofIdx, timeIdx));

        const auto& problem = elemCtx.problem();
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        unsigned globalSpaceIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
        unsigned pvtRegionIdx = priVars.pvtRegionIndex();
        fluidState_.setPvtRegionIndex(pvtRegionIdx);

        // extract the water and the gas saturations for convenience
        Evaluation Sw = priVars.makeEvaluation(Indices::waterSaturationIdx, timeIdx);

        Evaluation Sg;
        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg)
            // -> threephase case
            Sg = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
        else if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv) {
            // -> gas-water case
            Sg = 1 - Sw;

            // deal with solvent
            if (enableSolvent)
                Sg -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx);
        }
        else {
            assert(priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Rs);
            // -> oil-water case
            Sg = 0.0;
        }

        Opm::Valgrind::CheckDefined(Sg);
        Opm::Valgrind::CheckDefined(Sw);

        Evaluation So = 1.0 - Sw - Sg;

        // deal with solvent
        if (enableSolvent)
            So -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx);

        fluidState_.setSaturation(waterPhaseIdx, Sw);
        fluidState_.setSaturation(gasPhaseIdx, Sg);
        fluidState_.setSaturation(oilPhaseIdx, So);

        asImp_().solventPreSatFuncUpdate_(elemCtx, dofIdx, timeIdx);
        asImp_().polymerPreSatFuncUpdate_(elemCtx, dofIdx, timeIdx);


        // now we compute all phase pressures
        Evaluation pC[numPhases];
        const auto& materialParams = problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);

        //oil is the reference phase for pressure
        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv) {
            const Evaluation& pg = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fluidState_.setPressure(phaseIdx, pg + (pC[phaseIdx] - pC[gasPhaseIdx]));
        }

        else {
            const Evaluation& po = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fluidState_.setPressure(phaseIdx, po + (pC[phaseIdx] - pC[oilPhaseIdx]));
        }

        // calculate relative permeabilities. note that we store the result into the
        // mobility_ class attribute. the division by the phase viscosity happens later.
        MaterialLaw::relativePermeabilities(mobility_, materialParams, fluidState_);
        Opm::Valgrind::CheckDefined(mobility_);

        asImp_().solventPostSatFuncUpdate_(elemCtx, dofIdx, timeIdx);
        asImp_().polymerPostSatFuncUpdate_(elemCtx, dofIdx, timeIdx);


        Scalar SoMax = elemCtx.model().maxOilSaturation(globalSpaceIdx);

        // take the meaning of the switiching primary variable into account for the gas
        // and oil phase compositions
        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
            // in the threephase case, gas and oil phases are potentially present, i.e.,
            // we use the compositions of the gas-saturated oil and oil-saturated gas.
            if (FluidSystem::enableDissolvedGas()) {
                const Evaluation& RsSat =
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            oilPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);
                fluidState_.setRs(RsSat);
            }
            else
                fluidState_.setRs(0.0);

            if (FluidSystem::enableVaporizedOil()) {
                const Evaluation& RvSat =
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);
                fluidState_.setRv(RvSat);
            }
            else
                fluidState_.setRv(0.0);
        }
        else if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Rs) {
            // if the switching variable is the mole fraction of the gas component in the
            // oil phase, we can directly set the composition of the oil phase
            const auto& Rs = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            fluidState_.setRs(Rs);

            if (FluidSystem::enableVaporizedOil()) {
                // the gas phase is not present, but we need to compute its "composition"
                // for the gravity correction anyway
                const auto& RvSat =
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);

                fluidState_.setRv(RvSat);
            }
            else
                fluidState_.setRv(0.0);
        }
        else {
            assert(priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv);

            const auto& Rv = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            fluidState_.setRv(Rv);

            if (FluidSystem::enableDissolvedGas()) {
                // the oil phase is not present, but we need to compute its "composition" for
                // the gravity correction anyway
                const auto& RsSat =
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            oilPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);

                fluidState_.setRs(RsSat);
            }
            else
                fluidState_.setRs(0.0);
        }

        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.setRegionIndex(pvtRegionIdx);
        paramCache.setMaxOilSat(SoMax);
        paramCache.updateAll(fluidState_);

        // compute the phase densities and transform the phase permeabilities into mobilities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            const auto& b = FluidSystem::inverseFormationVolumeFactor(fluidState_, phaseIdx, pvtRegionIdx);
            fluidState_.setInvB(phaseIdx, b);

            const auto& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            mobility_[phaseIdx] /= mu;
        }
        Opm::Valgrind::CheckDefined(mobility_);

        // calculate the phase densities
        Evaluation rho;
        if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
            rho = fluidState_.invB(waterPhaseIdx);
            rho *= FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
            fluidState_.setDensity(waterPhaseIdx, rho);
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            rho = fluidState_.invB(gasPhaseIdx);
            rho *= FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableVaporizedOil()) {
                rho +=
                    fluidState_.invB(gasPhaseIdx) *
                    fluidState_.Rv() *
                    FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            }
            fluidState_.setDensity(gasPhaseIdx, rho);
        }

        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            rho = fluidState_.invB(oilPhaseIdx);
            rho *= FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableDissolvedGas()) {
                rho +=
                    fluidState_.invB(oilPhaseIdx) *
                    fluidState_.Rs() *
                    FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            }
            fluidState_.setDensity(oilPhaseIdx, rho);
        }

        // retrieve the porosity from the problem
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);

        // the porosity must be modified by the compressibility of the
        // rock...
        Scalar rockCompressibility = problem.rockCompressibility(elemCtx, dofIdx, timeIdx);
        if (rockCompressibility > 0.0) {
            Scalar rockRefPressure = problem.rockReferencePressure(elemCtx, dofIdx, timeIdx);
            Evaluation x = rockCompressibility*(fluidState_.pressure(oilPhaseIdx) - rockRefPressure);
            porosity_ *= 1.0 + x + 0.5*x*x;
        }

        asImp_().solventPvtUpdate_(elemCtx, dofIdx, timeIdx);
        asImp_().polymerPropertiesUpdate_(elemCtx, dofIdx, timeIdx);


        // update the quantities which are required by the chosen
        // velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

#ifndef NDEBUG
        // some safety checks in debug mode
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            assert(std::isfinite(Toolbox::value(fluidState_.density(phaseIdx))));
            assert(std::isfinite(Toolbox::value(fluidState_.saturation(phaseIdx))));
            assert(std::isfinite(Toolbox::value(fluidState_.temperature(phaseIdx))));
            assert(std::isfinite(Toolbox::value(fluidState_.pressure(phaseIdx))));
            assert(std::isfinite(Toolbox::value(fluidState_.invB(phaseIdx))));
        }
        assert(std::isfinite(Toolbox::value(fluidState_.Rs())));
        assert(std::isfinite(Toolbox::value(fluidState_.Rv())));
#endif
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::mobility
     */
    const Evaluation& mobility(unsigned phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the index of the PVT region used to calculate the thermodynamic
     *        quantities.
     *
     * This allows to specify different Pressure-Volume-Temperature (PVT) relations in
     * different parts of the spatial domain. Note that this concept should be seen as a
     * work-around of the fact that the black-oil model does not capture the
     * thermodynamics well enough. (Because there is, err, only a single real world with
     * in which all substances follow the same physical laws and hence the same
     * thermodynamics.) Anyway: Since the ECL file format uses multiple PVT regions, we
     * support it as well in our black-oil model. (Note that, if it is not explicitly
     * specified, the PVT region index is 0.)
     */
    auto pvtRegionIndex() const
        -> decltype(std::declval<FluidState>().pvtRegionIndex())
    { return fluidState_.pvtRegionIndex(); }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
     */
    Evaluation relativePermeability(unsigned phaseIdx) const
    {
        // warning: slow
        return fluidState_.viscosity(phaseIdx)*mobility(phaseIdx);
    }

private:
    friend BlackOilSolventIntensiveQuantities<TypeTag>;
    friend BlackOilPolymerIntensiveQuantities<TypeTag>;


    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    FluidState fluidState_;
    Evaluation porosity_;
    Evaluation mobility_[numPhases];
};

} // namespace Ewoms

#endif
