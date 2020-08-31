#ifndef ONEPHASEFLUIDSYSTEM_HH
#define ONEPHASEFLUIDSYSTEM_HH

//#include "components.hh"
//#include "chiwoms.h"

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
//#include "ChiParameterCache.hpp"

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/fluidsystems/NullParameterCache.hpp>
#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidsystems/BrineCO2FluidSystem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>

#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>

#include <opm/material/binarycoefficients/H2O_CO2.hpp>
#include <opm/material/binarycoefficients/Brine_CO2.hpp>

namespace Opm {
/*!
 * \ingroup Fluidsystems
 *
 * \brief A one-phase fluid system with co2 and octane as components.
 */
template <class Scalar>
class OnePhaseCo2OctaneFluidSystem
    : public Opm::BaseFluidSystem<Scalar, OnePhaseCo2OctaneFluidSystem<Scalar> >
{
    typedef OnePhaseCo2OctaneFluidSystem<Scalar> ThisType;
    typedef Opm::BaseFluidSystem<Scalar, ThisType> Base;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    template <class Evaluation>
    using ParameterCache = Opm::NullParameterCache<Evaluation>;

    /****************************************
     * Fluid phase related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 1;

    //! Index of the liquid phase
    static const int onlyPhase = 0;

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        static const char* name[] = {"o"};  // oleic phase

        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned /*phaseIdx*/)
    {
        // both CO2 and Octane have associative effects, separately and together
        return false;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 2;

    //! The component index of the oil; octane
    static const int OctaneIdx = 0;

    //! The component index of the solvent; co2
    static const int CO2Idx = 1;

    //! The component for pure oil
    typedef Opm::Octane<Scalar> Octane;

    //! The component for pure solvent
    typedef Opm::CO2<Scalar, Scalar> CO2;

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
            : throw std::invalid_argument("Acentric factor component index");
    }

    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \copydoc BaseFluidSystem::density
     */
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& /*paramCache*/,
                           unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        const auto& x = Opm::decay<LhsEval>(fluidState.moleFraction(phaseIdx, CO2Idx));

        return EOS::oleic_density(T, p, x);
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
        assert(T == (TEMPERATURE + 273.15));

        return EOS::oleic_viscosity(T, p, x);
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

        return EOS::oleic_enthalpy(T, p, x);
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        // thermodynamic function represented by density
        const LhsEval rho = (compIdx == OctaneIdx ?
                             EOS::oleic_density(T, p, LhsEval(0.)) :  // x_CO2 == 0.
                             EOS::oleic_density(T, p, LhsEval(1.)));  // x_CO2 == 1.

        // molar mass
        const Scalar M = molarMass(compIdx);

        // fugacity coefficient calculated
        const LhsEval phi = rho * Opm::IdealGas<Scalar>::R * T / (M * p);

        return phi;
    }

    //! \copydoc BaseFluidSystem::diffusionCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval diffusionCoefficient(const FluidState& /*fluidState*/,
                                        const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                        unsigned /*phaseIdx*/,
                                        unsigned /*compIdx*/)
    {
        return DIFFUSIVITY;
    }
};

};



#endif // ONEPHASEFLUIDSYSTEM_HH
