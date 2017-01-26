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
    static void initFromDeck(const Opm::Deck& deck, const Opm::EclipseState& eclState OPM_UNUSED)
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

    static std::string primaryVarName(unsigned pvIdx)
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

    static std::string eqName(unsigned eqIdx)
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

        const auto& fs = intQuants.fluidState();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            unsigned eqIdx = contiSolventEqIdx;

            storage[eqIdx] +=
                Toolbox::template decay<LhsEval>(intQuants.porosity())
                * Toolbox::template decay<LhsEval>(fs.saturation(phaseIdx))
                * Toolbox::template decay<LhsEval>(intQuants.solventDensity(phaseIdx));
        }
    }

    template <class UpEval>
    static void addPhaseFlux(RateVector& flux,
                             const ExtensiveQuantities& extQuants,
                             const IntensiveQuantities& upQuants,
                             unsigned phaseIdx)
    {
        if (!enableSolvent)
            return;

        const auto& volFlux = extQuants.volumeFlux(phaseIdx);

        const UpEval& upRhoSol =
            Toolbox::template decay<UpEval>(upQuants.solventDensity(phaseIdx));

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
        // do not consider solvents for convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    /*!
     * \brief Return how much a residual is considered an error
     */
    static Scalar computeResidualError(const EqVector& resid OPM_UNUSED)
    {
        // do not consider solvents for convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
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
};

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilSolventIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        solvents extension of the black-oil model.
 */
template <class TypeTag>
class BlackOilSolventIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    static constexpr int solventSaturationIdx = Indices::solventSaturationIdx;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;

public:
    /*!
     * \brief Update the intensive quantities needed to handle solvents from the primary variables.
     */
    template <class Context>
    void update_(const Context& context,
                 unsigned spaceIdx,
                 unsigned timeIdx)
    {
        const PrimaryVariables& priVars = context.primaryVars(spaceIdx, timeIdx);

        solventSaturation_ = priVars.makeEvaluation(solventSaturationIdx, timeIdx);
    }

    const Evaluation& solventSaturation() const
    {
        return solventSaturation_;
    }

    const Evaluation& solventDensity(unsigned phaseIdx OPM_UNUSED) const
    {
        return solventSaturation();
    }

protected:
    Evaluation solventSaturation_;
    Evaluation solventDensity_;
};

} // namespace Ewoms

#endif
