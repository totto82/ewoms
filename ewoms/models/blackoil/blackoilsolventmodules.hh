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
template <class TypeTag, unsigned numSolventsV = GET_PROP_VALUE(TypeTag, NumSolvents)>
class BlackOilSolventModule
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef Opm::MathToolbox<Evaluation> Toolbox;

    static constexpr unsigned solvent0PrimaryVarIdx = Indices::solvent0PrimaryVarIdx;
    static constexpr unsigned contiSolvent0EqIdx = Indices::contiSolvent0EqIdx;
    static constexpr unsigned numSolvents = numSolventsV;
    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:
    /*!
     * \brief Register all run-time parameters for the black-oil solvent module.
     */
    static void registerParameters()
    {
        // TODO: (if any)
    }

    /*!
     * \brief Register all solvent specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model OPM_UNUSED)
    {
        // TODO: (if any)
    }

    static bool primaryVarApplies(unsigned pvIdx)
    { return pvIdx <= solvent0PrimaryVarIdx && solvent0PrimaryVarIdx + numSolvents < pvIdx; }

    static std::string primaryVarName(unsigned pvIdx)
    {
        assert(primaryVarApplies(pvIdx));

        unsigned solCompIdx = pvIdx - solvent0PrimaryVarIdx;
        // TODO: What is chosen as primary variables for the solvents? This could be
        // e.g. the mole fractions in a given phase, total mass fractions, or
        // concentrations...
        return "c^sol," + std::to_string(solCompIdx);
    }

    static Scalar primaryVarWeight(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    { return eqIdx <= contiSolvent0EqIdx && contiSolvent0EqIdx + numSolvents < eqIdx; }

    static std::string eqName(unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        unsigned solCompIdx = eqIdx - contiSolvent0EqIdx;
        return "conti^sol," + std::to_string(solCompIdx);
    }

    static Scalar eqWeight(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    template <class LhsEval>
    static void addPhaseStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                                const IntensiveQuantities& intQuants,
                                unsigned phaseIdx)
    {
        for (unsigned solCompIdx = 0; solCompIdx < numSolvents; ++ solCompIdx) {
            unsigned eqIdx = contiSolvent0EqIdx + solCompIdx;

            storage[eqIdx] +=
                Toolbox::template decay<LhsEval>(intQuants.porosity())
                * Toolbox::template decay<LhsEval>(intQuants.fluidState().saturation(phaseIdx))
                * Toolbox::template decay<LhsEval>(intQuants.solventConcentration(phaseIdx, solCompIdx));
        }
    }

    static void addPhaseFlux(RateVector& flux,
                             const ExtensiveQuantities& extQuants,
                             const IntensiveQuantities& upQuants,
                             unsigned phaseIdx)
    {
        const auto& volFlux = extQuants.volumeFlux(phaseIdx);

        for (unsigned solCompIdx = 0; solCompIdx < numSolvents; ++ solCompIdx) {
            int eqIdx = contiSolvent0EqIdx + solCompIdx;
            const auto& upC = upQuants.solventConcentration(phaseIdx, solCompIdx);

            flux[eqIdx] += volFlux*upC;
        }
    }

    /*!
     * \brief Assign the solvent specific primary variables to a PrimaryVariables object
     */
    template <class SolventValues>
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  const SolventValues& solventPvValues)
    {
        for (unsigned solCompIdx = 0; solCompIdx < numSolvents; ++ solCompIdx) {
            unsigned pvIdx = solvent0PrimaryVarIdx + solCompIdx;
            priVars[pvIdx] = solventPvValues[solCompIdx];
        }
    }

    /*!
     * \brief Do a Newton-Raphson update the primary variables of the solvents.
     */
    static void updatePrimaryVars(PrimaryVariables& newPv,
                                  const PrimaryVariables& oldPv,
                                  const EqVector& delta)
    {
        // do a plain unchopped Newton update
        for (unsigned solCompIdx = 0; solCompIdx < numSolvents; ++ solCompIdx) {
            unsigned pvIdx = solvent0PrimaryVarIdx + solCompIdx;

            newPv[pvIdx] = oldPv[pvIdx] - delta[pvIdx];
        }
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
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
        unsigned dofIdx = model.dofMapper().index(dof);
#else
        unsigned dofIdx = model.dofMapper().map(dof);
#endif

        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];

        for (unsigned solCompIdx = 0; solCompIdx < numSolvents; ++ solCompIdx) {
            unsigned pvIdx = solvent0PrimaryVarIdx + solCompIdx;
            outstream << priVars[pvIdx];
        }
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
        unsigned dofIdx = model.dofMapper().index(dof);
#else
        unsigned dofIdx = model.dofMapper().map(dof);
#endif

        PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        for (unsigned solCompIdx = 0; solCompIdx < numSolvents; ++ solCompIdx) {
            unsigned pvIdx = solvent0PrimaryVarIdx + solCompIdx;
            instream >> priVars0[pvIdx];

            // set the primary variables for the beginning of the current time step.
            priVars1 = priVars0[pvIdx];
        }
    }
};

// the specialization of the solvent module for the trivial case (i.e., no solvents)
template <class TypeTag>
class BlackOilSolventModule<TypeTag, /*numSolvents=*/0>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);

public:
    static void registerParameters()
    {}

    static void registerOutputModules(Model& model OPM_UNUSED)
    { }

    static bool primaryVarApplies(unsigned pvIdx OPM_UNUSED)
    { return false; }

    static std::string primaryVarName(unsigned pvIdx OPM_UNUSED)
    { OPM_THROW(std::logic_error, "Solvents are disabled!"); }

    static Scalar primaryVarWeight(unsigned pvIdx OPM_UNUSED)
    { OPM_THROW(std::logic_error, "Solvents are disabled!"); }

    static bool eqApplies(unsigned eqIdx OPM_UNUSED)
    { return false; }

    static std::string eqName(unsigned eqIdx OPM_UNUSED)
    { OPM_THROW(std::logic_error, "Solvents are disabled!"); }

    static Scalar eqWeight(unsigned eqIdx OPM_UNUSED)
    { OPM_THROW(std::logic_error, "Solvents are disabled!"); }

    template <class LhsEval>
    static void addPhaseStorage(Dune::FieldVector<LhsEval, numEq>& storage OPM_UNUSED,
                                const IntensiveQuantities& intQuants OPM_UNUSED,
                                unsigned phaseIdx OPM_UNUSED)
    { }

    static void addPhaseFlux(RateVector& flux OPM_UNUSED,
                             const ExtensiveQuantities& extQuants OPM_UNUSED,
                             const IntensiveQuantities& upQuants OPM_UNUSED,
                             unsigned phaseIdx OPM_UNUSED)
    { }

    template <class SolventValues>
    static void assignPrimaryVars(PrimaryVariables& priVars OPM_UNUSED,
                                  const SolventValues& solventPvValues OPM_UNUSED)
    { }

    static void updatePrimaryVars(PrimaryVariables& newPv OPM_UNUSED,
                                  const PrimaryVariables& oldPv OPM_UNUSED,
                                  const EqVector& delta OPM_UNUSED)
    { }

    static Scalar computeUpdateError(const PrimaryVariables& oldPv OPM_UNUSED,
                                     const EqVector& delta OPM_UNUSED)
    { return static_cast<Scalar>(0.0); }

    static Scalar computeResidualError(const EqVector& resid OPM_UNUSED)
    { return static_cast<Scalar>(0.0); }

    template <class DofEntity>
    static void serializeEntity(const Model& model OPM_UNUSED,
                                std::ostream& outstream OPM_UNUSED,
                                const DofEntity& dof OPM_UNUSED)
    { }

    template <class DofEntity>
    static void deserializeEntity(Model& model OPM_UNUSED,
                                  std::istream& instream OPM_UNUSED,
                                  const DofEntity& dof OPM_UNUSED)
    { }
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
/*
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
*/
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    static constexpr int numSolvents = GET_PROP_VALUE(TypeTag, NumSolvents);
    static constexpr int solvent0PrimaryVarIdx = Indices::solvent0PrimaryVarIdx;
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

        for (unsigned solCompIdx = 0; solCompIdx < numSolvents; ++solCompIdx) {
            solventConcentration_[solCompIdx] =
                priVars.makeEvaluation(solvent0PrimaryVarIdx + solCompIdx, timeIdx);
        }
    }

    const Evaluation& solventConcentration(unsigned /*phaseIdx*/, unsigned solCompIdx) const
    {
#warning "TODO: when is which phase active? (In each case a non-zero value must be returned for some phase!)"
        //static const Evaluation zero(0.0);
        //if (phaseIdx != oilPhaseIdx)
        //    return zero;

        return solventConcentration_[solCompIdx];
    }

protected:
    // TODO: solvent specific intensive quantities
    std::array<Evaluation, numSolvents> solventConcentration_;
};

} // namespace Ewoms

#endif
