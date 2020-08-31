#ifndef COMPONENTS_HH
#define COMPONENTS_HH

#include "chiwoms.h"

#include <opm/material/IdealGas.hpp>
#include <opm/material/components/Component.hpp>
#include <opm/material/components/SimpleCO2.hpp>
#include <opm/material/components/CO2.hpp>
#include <opm/material/components/Brine.hpp>
#include <opm/material/components/H2O.hpp>

#include "padua.h"
extern struct poly2d_t co2_c8_dens;
extern struct poly2d_t co2_c8_visc;
extern struct poly2d_t co2_c8_enth;
extern struct poly2d_t co2_h2o_dens;
extern struct poly2d_t co2_h2o_visc;
extern struct poly2d_t co2_h2o_enth;

namespace Opm {
/*!
 * \ingroup Components
 *
 * \brief A simple representation of linear octane
 *
 * \tparam Scalar The type used for scalar values
 */
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

template <class Scalar>
class ChiwomsCO2 : public Opm::SimpleCO2<Scalar>
{
public:
    /// Acentric factor
    static Scalar acentricFactor() { return 0.225; }
};

template <class Scalar>
class ChiwomsBrine : public Opm::H2O<Scalar>
{
public:
    /// Acentric factor
    static Scalar acentricFactor() { return 0.344; }
};

struct EOS
{
        template<typename LhsEval>
        static LhsEval eval(const struct poly2d_t & poly,
                            LhsEval T, LhsEval p, LhsEval x) {
                LhsEval val;
                // if the linear solver overshoots and attempts to pick a solution
                // outside of the domain, just pretend that there is no change from the
                // boundary of the domain; hopefully the solver will regain its senses
                // and bring us back into the domain if we don't trip it up with
                // outlandish values.
                LhsEval fenced_p = Opm::min(Opm::max(p, poly.lo_x), poly.up_x);
                LhsEval fenced_x = Opm::min(Opm::max(x, poly.lo_y), poly.up_y);
                padua_eval(&poly, 1, &fenced_p, &fenced_x, &val);
                return val;
        }

        template<typename LhsEval>
        static LhsEval oleic_density(LhsEval T, LhsEval p, LhsEval x) {
                assert(T == (TEMPERATURE + 273.15));
                return 650;
        }

        template<typename LhsEval>
        static LhsEval aqueous_density(LhsEval T, LhsEval p, LhsEval x) {
                assert(T == (TEMPERATURE + 273.15));
                return 1000;
        }

        template<typename LhsEval>
        static LhsEval oleic_viscosity(LhsEval T, LhsEval p, LhsEval x) {
                assert(T == (TEMPERATURE + 273.15));
                return 5e-3;
        }

        template<typename LhsEval>
        static LhsEval aqueous_viscosity(LhsEval T, LhsEval p, LhsEval x) {
                assert(T == (TEMPERATURE + 273.15));
                return 1e-3;
        }

        template<typename LhsEval>
        static LhsEval oleic_enthalpy(LhsEval T, LhsEval p, LhsEval x) {
                assert(T == (TEMPERATURE + 273.15));
                return 0;
        }

        template<typename LhsEval>
        static LhsEval aqueous_enthalpy(LhsEval T, LhsEval p, LhsEval x) {
                assert(T == (TEMPERATURE + 273.15));
                return 0;
        }
};

} // namespace opm

#endif // COMPONENTS_HH
