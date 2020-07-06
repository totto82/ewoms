#ifndef TWOPHASEFLUIDSYSTEM_HH
#define TWOPHASEFLUIDSYSTEM_HH

#include <iostream>
#include <cassert>
#include <stdexcept>  // invalid_argument
#include <sstream>
#include <iostream>
#include <string>
#include <random>    // mt19937, normal_distribution
#include <limits>    // epsilon
#include <boost/format.hpp>  // boost::format

#include <opm/common/Exceptions.hpp>
#include <opm/material/IdealGas.hpp>

#include <opm/material/components/Component.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/components/CO2.hpp>
#include <opm/material/components/Brine.hpp>
#include <opm/material/eos/PengRobinsonMixture.hpp>
#include <opm/material/eos/PengRobinsonParamsMixture.hpp>
#include "ChiParameterCache.hpp"

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/fluidsystems/NullParameterCache.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>

#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>

namespace Opm {
template <class Scalar>
class ChiwomsCO2 : public Opm::SimpleCO2<Scalar>
{
public:
    /// Acentric factor
    static Scalar acentricFactor() { return 0.225; }
};

template <class Scalar>
class Octane : public Opm::Component<Scalar, Octane<Scalar> >
{
public:
        /// Chemical name
        static const char* name() { return "C8"; }

        /// Molar mass in \f$\mathrm{[kg/mol]}\f$
        static Scalar molarMass() { return 0.11423; }

        /// Critical temperature in \f$\mathrm[K]}\f$
        static Scalar criticalTemperature() { return 568.7; }

        /// Critical pressure in \f$\mathrm[Pa]}\f$
        static Scalar criticalPressure() { return 2.49e6; }

        /// Acentric factor
        static Scalar acentricFactor() { return 0.398; }
};
/*!
 * \ingroup Fluidsystems
 *
 * \brief A two-phase fluid system with brine and octane as the main components
 * in each their phase, and CO2 as solvent in both.
 */
template <class Scalar>
class TwoPhaseCo2OctaneFluidSystem
        : public Opm::BaseFluidSystem<Scalar, TwoPhaseCo2OctaneFluidSystem<Scalar> >
{
    typedef TwoPhaseCo2OctaneFluidSystem<Scalar> ThisType;
    typedef Opm::BaseFluidSystem<Scalar, ThisType> Base;
    typedef typename Opm::PengRobinson<Scalar> PengRobinson;
    typedef typename Opm::PengRobinsonMixture<Scalar, ThisType> PengRobinsonMixture;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    template <class Evaluation>
    using ParameterCache = Opm::ChiParameterCache<Evaluation, ThisType>;

    /****************************************
     * Fluid phase related static parameters
    ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 2;

    //! Index of the liquid phase
    static const int oilPhaseIdx = 0;
    static const int gasPhaseIdx = 1;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        static const char* name[] = {"o",  // oleic phase
                                    "g"};  // gas phase

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned phaseIdx)
    {
        if (phaseIdx == oilPhaseIdx)
            return false;

        // CO2 have associative effects with octane
        return true;
    }


    /****************************************
    * Component related static parameters
    ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 2;  // octane, co2

    //! The component index of the oil; octane
    static const int OctaneIdx = 0;

    //! The component index of the solvent; co2
    static const int CO2Idx = 1;

    //! The component for pure oil
    typedef Opm::Octane<Scalar> Octane;

    //! The component for pure solvent
    typedef Opm::ChiwomsCO2<Scalar> CO2;


    static void init(Scalar minT = 273.15,
                     Scalar maxT = 373.15,
                     Scalar minP = 1e4,
                     Scalar maxP = 100e6)
    {
        Opm::PengRobinsonParamsMixture<Scalar, ThisType, oilPhaseIdx, /*useSpe5=*/true> prParams;

        // find envelopes of the 'a' and 'b' parameters for the range
        // minT <= T <= maxT and minP <= p <= maxP. For
        // this we take advantage of the fact that 'a' and 'b' for
        // mixtures is just a convex combination of the attractive and
        // repulsive parameters of the pure components

        Scalar minA = 1e30, maxA = -1e30;
        Scalar minB = 1e30, maxB = -1e30;

        prParams.updatePure(minT, minP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(maxT, minP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(minT, maxP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };

        prParams.updatePure(maxT, maxP);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            minA = std::min(prParams.pureParams(compIdx).a(), minA);
            maxA = std::max(prParams.pureParams(compIdx).a(), maxA);
            minB = std::min(prParams.pureParams(compIdx).b(), minB);
            maxB = std::max(prParams.pureParams(compIdx).b(), maxB);
        };
        PengRobinson::init(/*aMin=*/minA, /*aMax=*/maxA, /*na=*/100,
                           /*bMin=*/minB, /*bMax=*/maxB, /*nb=*/200);
    }

    //! \copydoc BaseFluidSystem::componentName
    static const char* componentName(unsigned compIdx)
    {
            static const char* name[] = {
                    Octane::name(),
                    CO2::name()
            };
            assert(0 <= compIdx && compIdx < numComponents);
            return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned compIdx)
    {
        return (compIdx == OctaneIdx)
            ? Octane::molarMass()
            : (compIdx == CO2Idx)
            ? CO2::molarMass()
            : throw std::invalid_argument("Molar mass component index");
    }

    /*!
     * \brief Critical temperature of a component [K].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar criticalTemperature(unsigned compIdx)
    {
            return (compIdx == OctaneIdx)
                    ? Octane::criticalTemperature()
                    : (compIdx == CO2Idx)
                    ? CO2::criticalTemperature()
                    : throw std::invalid_argument("Critical temperature component index");
    }

    /*!
     * \brief Critical pressure of a component [Pa].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar criticalPressure(unsigned compIdx)
    {
            return (compIdx == OctaneIdx)
                    ? Octane::criticalPressure()
                    : (compIdx == CO2Idx)
                    ? CO2::criticalPressure()
                    : throw std::invalid_argument("Critical pressure component index");
    }

    /*!
     * \brief The acentric factor of a component [].
     *
     * \copydetails Doxygen::compIdxParam
     */
    static Scalar acentricFactor(unsigned compIdx)
    {
            return (compIdx == OctaneIdx)
                    ? Octane::acentricFactor()
                    : (compIdx == CO2Idx)
                    ? CO2::acentricFactor()
                    : throw std::invalid_argument("Molar mass component index");
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \copydoc BaseFluidSystem::density
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& paramCache,
                           unsigned phaseIdx)
    {
        if (false)// set to true if you want constant density
        {
            if(phaseIdx == oilPhaseIdx) {
                return 670; 
            } else {
                return 1.7; 
            }
        } else {
            assert(0 <= phaseIdx && phaseIdx < numPhases);
            return fluidState.averageMolarMass(phaseIdx)/paramCache.molarVolume(phaseIdx);
        }
    }

        //! \copydoc BaseFluidSystem::viscosity
        template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
        static LhsEval viscosity(const FluidState& fluidState,
                                 const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                 unsigned phaseIdx)
        {
            assert(0 <= phaseIdx && phaseIdx < numPhases);

            const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
            const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
            const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));

            #warning We use constant viscosity. These needs to be checked
            if(phaseIdx == oilPhaseIdx) {
                return 5e-4; //EOS::oleic_viscosity(T, p, x);
            } else {
                return 1e-5; //EOS::aqueous_viscosity(T, p, x);
            }
        }

        //! \copydoc BaseFluidSystem::enthalpy
        template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
        static LhsEval enthalpy(const FluidState& fluidState,
                                const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                unsigned phaseIdx)
        {
            const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
            const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
            const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));
            throw std::runtime_error("We don't use the enthalpy for non-isothermal runs");
        }

        // according to <https://srdata.nist.gov/solubility/sol_detail.aspx?sysID=38_103>
        // the solubility of octane in aqueous phase is 2.0g/100g sln. since octane
        // (C8H18) has a molecular weight of 114.23 g/mol and water (H2O) of 18.01528 g/mol
        // we have a total of 2.0g/114.23 g/mol ~= 0.0175 mol of octane and
        // (100g-2.0g)/18.01528 g/mol ~= 5.44 mol of water, for a total of 5.45 mol,
        // of which the mole fraction of octane is 0.0175/5.45 ~= 3.2e-3
        //constexpr static Scalar sol_aqueous_oil = 3.208e-3;  // solution of octane in brine

        // the solubility of water in the oleic (octane) phase is according the same
        // reference above, 7.3g/100g sln, giving 7.3g/18.01528 g/mol ~= 0.41 mol of
        // water and (100g-7.3g)/114.23 g/mol ~= 0.81 mol of octane, in a 100 g solution,
        // yielding a mole fraction of 0.41/(0.41 + 0.81) = 0.33 for water.
        //constexpr static Scalar sol_oleic_water = sol_aqueous_oil; // 3.330e-1;

        // partition coefficients when both oleic and aqueous phases are saturated with
        // maximum dissoluted amount of the other component. these coefficients should
        // give fugacities with maximum volatility of the component in its non-native phase
        //constexpr static Scalar k_aqueous_oil = (1 - sol_oleic_water) / sol_aqueous_oil;  // ~200
        //constexpr static Scalar k_oleic_water = (1 - sol_aqueous_oil) / sol_oleic_water;  // ~3


        //! \copydoc BaseFluidSystem::fugacityCoefficient
        template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
        static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& paramCache,
                                           unsigned phaseIdx,
                                           unsigned compIdx)
        {
            assert(0 <= phaseIdx && phaseIdx < numPhases);
            assert(0 <= compIdx && compIdx < numComponents);

            if (phaseIdx == oilPhaseIdx) {
#if 1
#warning HACK We use henry's law
                if (compIdx == OctaneIdx)
                    return 40e3/fluidState.pressure(oilPhaseIdx);
                else
                return 500e3/fluidState.pressure(oilPhaseIdx);
#else
                if (compIdx == CO2Idx)
                    return 500e3/fluidState.pressure(oilPhaseIdx);
                else {
                    return PengRobinsonMixture::computeFugacityCoefficient(fluidState,
								     paramCache,
								     phaseIdx,
								     compIdx);
                }
#endif
            } else if (phaseIdx == gasPhaseIdx) {
                return 1.0;
            } else {
                throw std::invalid_argument("expects oil or gas phase!");
            }
        }

        //! \copydoc BaseFluidSystem::diffusionCoefficient
        template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
        static LhsEval diffusionCoefficient(const FluidState& /*fluidState*/,
                                            const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                            unsigned /*phaseIdx*/,
                                            unsigned /*compIdx*/)
        {
            // The diffusionCoefficient. Needs to be checked
            return 1e-9; 
        }

    /*!
     * \brief Returns the interaction coefficient for two components.
     *
     * The values are from Ivar
     */
    static Scalar interactionCoefficient(unsigned comp1Idx, unsigned comp2Idx)
    {
        unsigned i = std::min(comp1Idx, comp2Idx);
        unsigned j = std::max(comp1Idx, comp2Idx);
        if (i == OctaneIdx && j == CO2Idx)
            return 0.1089;

        return 0;
    }

};

};//namespace opm

#endif // TWOPHASEFLUIDSYSTEM_HH
