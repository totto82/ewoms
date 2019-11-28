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
 * \copydoc Opm::Co2InjectionProblem
 */
#ifndef EWOMS_CO2_INJECTION_PROBLEM_HH
#define EWOMS_CO2_INJECTION_PROBLEM_HH

#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/simulators/linalg/parallelamgbackend.hh>

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
#include <opm/material/binarycoefficients/Brine_CO2.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/Unused.hpp>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <iostream>
#include <string>

namespace Opm {
//! \cond SKIP_THIS
template <class TypeTag>
class Co2InjectionProblem;

namespace Co2Injection {
//#include <opm/material/components/co2tables.inc>
#include <opm/material/components/fineCo2TablesPureWater.inc>
}
//! \endcond
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(Co2InjectionBaseProblem);

// declare the CO2 injection problem specific property tags
NEW_PROP_TAG(FluidSystemPressureLow);
NEW_PROP_TAG(FluidSystemPressureHigh);
NEW_PROP_TAG(FluidSystemNumPressure);
NEW_PROP_TAG(FluidSystemTemperatureLow);
NEW_PROP_TAG(FluidSystemTemperatureHigh);
NEW_PROP_TAG(FluidSystemNumTemperature);

NEW_PROP_TAG(MaxDepth);
NEW_PROP_TAG(Temperature);
NEW_PROP_TAG(SimulationName);

// Set the grid type
//SET_TYPE_PROP(Co2InjectionBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the grid type
#if HAVE_DUNE_ALUGRID
// use dune-alugrid if available
#warning "using dune-alugrid. adaptive grid refinement will be available, but parallelism won't"
SET_TYPE_PROP(Co2InjectionBaseProblem,
              Grid,
              Dune::ALUGrid</*dim=*/2,
                            /*dimWorld=*/2,
                            Dune::cube,
                            Dune::nonconforming>);
#else
SET_TYPE_PROP(Co2InjectionBaseProblem, Grid, Dune::YaspGrid<2>);
#endif


// Set the problem property
SET_TYPE_PROP(Co2InjectionBaseProblem, Problem,
              Opm::Co2InjectionProblem<TypeTag>);

// Set fluid configuration
SET_PROP(Co2InjectionBaseProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::Co2Injection::CO2Tables CO2Tables;

public:
    typedef Opm::BrineCO2FluidSystem<Scalar, CO2Tables> type;
    //typedef Opm::H2ON2FluidSystem<Scalar, /*useComplexRelations=*/false> type;
};

// Set the material Law
SET_PROP(Co2InjectionBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::liquidPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

    // define the material law which is parameterized by effective
    // saturations
    typedef Opm::RegularizedBrooksCorey<Traits> EffMaterialLaw;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::EffToAbsLaw<EffMaterialLaw> type;
};

// Set the thermal conduction law
SET_PROP(Co2InjectionBaseProblem, ThermalConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::SomertonThermalConductionLaw<FluidSystem, Scalar> type;
};

// set the energy storage law for the solid phase
SET_TYPE_PROP(Co2InjectionBaseProblem, SolidEnergyLaw,
              Opm::ConstantSolidHeatCapLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Use the algebraic multi-grid linear solver for this problem
SET_TAG_PROP(Co2InjectionBaseProblem, LinearSolverSplice, ParallelAmgLinearSolver);

// Write the Newton convergence behavior to disk?
SET_BOOL_PROP(Co2InjectionBaseProblem, NewtonWriteConvergence, false);

// Enable gravity
SET_BOOL_PROP(Co2InjectionBaseProblem, EnableGravity, true);

SET_BOOL_PROP(Co2InjectionBaseProblem, EnableDiffusion, true);

// set the defaults for the problem specific properties
SET_SCALAR_PROP(Co2InjectionBaseProblem, FluidSystemPressureLow, 0.8e7);
SET_SCALAR_PROP(Co2InjectionBaseProblem, FluidSystemPressureHigh, 1.2e7);
SET_INT_PROP(Co2InjectionBaseProblem, FluidSystemNumPressure, 100);
SET_SCALAR_PROP(Co2InjectionBaseProblem, FluidSystemTemperatureLow, 290);
SET_SCALAR_PROP(Co2InjectionBaseProblem, FluidSystemTemperatureHigh, 350);
SET_INT_PROP(Co2InjectionBaseProblem, FluidSystemNumTemperature, 100);

SET_SCALAR_PROP(Co2InjectionBaseProblem, MaxDepth, 1);
SET_SCALAR_PROP(Co2InjectionBaseProblem, Temperature, 323.15);
SET_STRING_PROP(Co2InjectionBaseProblem, SimulationName, "uncover");

// The default for the end time of the simulation is 2 hours
SET_SCALAR_PROP(Co2InjectionBaseProblem, EndTime, 2*3600);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(Co2InjectionBaseProblem, InitialTimeStepSize, 1e-3);

// The default DGF file to load
SET_STRING_PROP(Co2InjectionBaseProblem, GridFile, "data/uncover.dgf");

END_PROPERTIES

namespace Opm {
/*!
 * \ingroup TestProblems
 *
 * \brief Problem where \f$CO_2\f$ is injected under a low permeable
 *        layer at a depth of 2700m.
 *
 * The domain is sized 0.2m times 0.2m and consists of a homogenious porous media
 *
 * \f$CO_2\f$ gets injected by means of a forced-flow boundary
 * condition into water-filled aquifer, which is situated 2700m below
 * sea level, at the lower-right boundary (\f$5m<y<15m\f$) and
 * migrates upwards due to buoyancy. It accumulates and eventually
 * enters the lower permeable aquitard.
 *
 * The boundary conditions applied by this problem are no-flow
 * conditions on the top bottom and right boundaries and a free-flow
 * boundary condition on the left. For the free-flow condition,
 * hydrostatic pressure is assumed.
 */
template <class TypeTag>
class Co2InjectionProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { numPhases = FluidSystem::numPhases };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { CO2Idx = FluidSystem::CO2Idx };
    enum { BrineIdx = FluidSystem::BrineIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiCO2EqIdx = conti0EqIdx + CO2Idx };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductionLaw) ThermalConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, SolidEnergyLawParams) SolidEnergyLawParams;
    typedef typename ThermalConductionLaw::Params ThermalConductionLawParams;

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    Co2InjectionProblem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        unsigned refinement = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);
        eps_ = 0.2 /2 / std::pow(2,refinement);

        pressure_ = 100 *1e5;

        temperatureLow_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemTemperatureLow);
        temperatureHigh_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemTemperatureHigh);
        nTemperature_ = EWOMS_GET_PARAM(TypeTag, unsigned, FluidSystemNumTemperature);

        pressureLow_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureLow);
        pressureHigh_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureHigh);
        nPressure_ = EWOMS_GET_PARAM(TypeTag, unsigned, FluidSystemNumPressure);

        maxDepth_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxDepth);
        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);

        // initialize the tables of the fluid system
        // FluidSystem::init();
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);


        // intrinsic permeabilities
        K_ = this->toDimMatrix_(76 * 9.8692 * 1e-13);

        // porosities
        size_t numDof = this->model().numGridDof();
        porosity_.resize(numDof,0.4);
        //some noise to generate fingers
        for (size_t i = 0; i < numDof; ++i) {
            double noise = this->norm_dist(this->rand_gen);
            porosity_[i] += 1 * noise;
        }

        molEps_.resize(numDof, 0.0);
        for (size_t i = 0; i < numDof; ++i) {
            double noise = this->norm_dist(this->rand_gen);
            molEps_[i] += noise;
        }


        // residual saturations
        materialParams_.setResidualSaturation(liquidPhaseIdx, 0.2);
        materialParams_.setResidualSaturation(gasPhaseIdx, 0.0);

        // parameters for the Brooks-Corey law
        materialParams_.setEntryPressure(5e3);
        materialParams_.setLambda(0.5);


        materialParams_.finalize();

        // assume constant heat capacity and granite
        solidEnergyLawParams_.setSolidHeatCapacity(790.0 // specific heat capacity of granite [J / (kg K)]
                                                   * 2700.0); // density of granite [kg/m^3]
        solidEnergyLawParams_.finalize();


        enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
        std::cout << "enableDiffusion "<< enableDiffusion << std::endl;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureLow,
                             "The lower temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureHigh,
                             "The upper temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, FluidSystemNumTemperature,
                             "The number of intervals between the lower and "
                             "upper temperature");

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureLow,
                             "The lower pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureHigh,
                             "The upper pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, FluidSystemNumPressure,
                             "The number of intervals between the lower and "
                             "upper pressure");

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, MaxDepth,
                             "The maximum depth [m] of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SimulationName,
                             "The name of the simulation used for the output "
                             "files");
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << EWOMS_GET_PARAM(TypeTag, std::string, SimulationName)
            << "_" << Model::name();
        if (GET_PROP_VALUE(TypeTag, EnableEnergy))
            oss << "_ni";
        oss << "_" << Model::discretizationName();
        return oss.str();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        Scalar tol = this->model().newtonMethod().tolerance()*1e5;
        this->model().checkConservativeness(tol);

        // Calculate storage terms
        PrimaryVariables storageL, storageG;
        this->model().globalPhaseStorage(storageL, /*phaseIdx=*/0);
        this->model().globalPhaseStorage(storageG, /*phaseIdx=*/1);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: liquid=[" << storageL << "]"
                      << " gas=[" << storageG << "]\n" << std::flush;
        }
#endif // NDEBUG
    }

    /// Constant temperature
    template <class Context>
    Scalar temperature(const Context& context OPM_UNUSED,
                       unsigned spaceIdx OPM_UNUSED,
                       unsigned timeIdx OPM_UNUSED) const
    { return temperature_; }

    /// Constant permeability
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context OPM_UNUSED,
                                           unsigned spaceIdx OPM_UNUSED,
                                           unsigned timeIdx OPM_UNUSED) const
    { return K_; }

    /// Constant porosity
    template <class Context>
    Scalar porosity(const Context& context ,
                    unsigned spaceIdx,
                    unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return porosity_[globalSpaceIdx];
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        return materialParams_;
    }

    /*!
     * \brief Return the parameters for the heat storage law of the rock
     *
     * In this case, we assume the rock-matrix to be granite.
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyLawParams(const Context& context OPM_UNUSED,
                         unsigned spaceIdx OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED) const
    { return solidEnergyLawParams_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams &
    thermalConductionLawParams(const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        return thermalCondParams_;
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector& values, const Context& context,
                  unsigned spaceIdx, unsigned timeIdx) const
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (onTop_(pos)) {
            Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
            const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

            //////
            // set temperature
            //////
            fs.setTemperature(temperature(context, spaceIdx, timeIdx));

            //////
            // set saturations
            //////
            fs.setSaturation(FluidSystem::liquidPhaseIdx, 1.0);
            fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);

            //////
            // set pressures
            //////
            /// \brief pl
            ///
            ///

            //Scalar densityL = FluidSystem::CO2::gasDensity(temperature_, Scalar(pressure_));
            Scalar densityL = FluidSystem::Brine::gasDensity(temperature_, Scalar(pressure_));

            Scalar depth = pos[dim - 1];
            Scalar pl = pressure_ + densityL * this->gravity()[dim - 1] * depth;

            Scalar pC[numPhases];
            const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
            MaterialLaw::capillaryPressures(pC, matParams, fs);

            fs.setPressure(liquidPhaseIdx, pl + (pC[liquidPhaseIdx] - pC[liquidPhaseIdx]));
            fs.setPressure(gasPhaseIdx, pl + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));

            //////
            // set composition of the liquid phase
            //////
            unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

            //Scalar tmp = 0.9 + 1*molEps_[globalSpaceIdx];
            Scalar tmp = std::min(0.01921, 0.01921 + 0.1*molEps_[globalSpaceIdx]);
            fs.setMoleFraction(liquidPhaseIdx, CO2Idx, tmp);
            fs.setMoleFraction(liquidPhaseIdx, BrineIdx,
                               1.0 - fs.moleFraction(liquidPhaseIdx, CO2Idx));


//            fs.setMoleFraction(gasPhaseIdx, CO2Idx, 1.0);
//            fs.setMoleFraction(gasPhaseIdx, BrineIdx,
//                               1.0 - fs.moleFraction(gasPhaseIdx, CO2Idx));

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
            CFRP::solve(fs, paramCache,
                        /*refPhaseIdx=*/liquidPhaseIdx,
                        /*setViscosity=*/true,
                        /*setEnthalpy=*/true);

            fs.checkDefined();

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
            //std::cout << globalSpaceIdx << " " << pos[0] << " " << pos[1] << std::endl;
            //for (unsigned compIdx = 0; compIdx < 2; ++compIdx) {
            //    std::cout  << compIdx << " " << values[compIdx] << " " <<fs.moleFraction(liquidPhaseIdx, compIdx) << std::endl;
            //}
        }
        else
            // no flow on top and bottom
            values.setNoFlow();
    }

    // \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        // const auto& matParams = this->materialLawParams(context, spaceIdx,
        // timeIdx);
        // values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
        values.assignNaive(fs);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    { rate = Scalar(0.0); }

    //! \}

private:
    template <class Context, class FluidState>
    void initialFluidState_(FluidState& fs,
                            const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        //////
        // set temperature
        //////
        fs.setTemperature(temperature(context, spaceIdx, timeIdx));

        //////
        // set saturations
        //////

        // This is used to test the whole cirle simulation
        if (true && onTopCell_(pos)) {
            fs.setSaturation(FluidSystem::liquidPhaseIdx, 1.0);
            fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);

            //////
            // set pressures
            //////
            Scalar pC[numPhases];
            const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
            MaterialLaw::capillaryPressures(pC, matParams, fs);

            Scalar densityL = FluidSystem::Brine::gasDensity(temperature_, Scalar(pressure_));

            Scalar depth = pos[dim - 1];
            // add a slight underpressure to start the fingers
            Scalar pl = pressure_ + densityL * this->gravity()[dim - 1] * depth;

            fs.setPressure(liquidPhaseIdx, pl + (pC[liquidPhaseIdx] - pC[liquidPhaseIdx]));
            fs.setPressure(gasPhaseIdx, pl + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));

            //////
            // set composition of the liquid phase
            //////
            double noise = this->norm_dist(this->rand_gen);
            double co2fraction = std::min(0.01921, 0.01921 * (1.0 + 1 * noise));
            std::cout << co2fraction << std::endl;
            fs.setMoleFraction(liquidPhaseIdx, CO2Idx, co2fraction);
            fs.setMoleFraction(liquidPhaseIdx, BrineIdx,
                               1.0 - fs.moleFraction(liquidPhaseIdx, CO2Idx));

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
            CFRP::solve(fs, paramCache,
                        /*refPhaseIdx=*/liquidPhaseIdx,
                        /*setViscosity=*/true,
                        /*setEnthalpy=*/true);
        } else {

            fs.setSaturation(FluidSystem::liquidPhaseIdx, 1.0);
            fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);

            //////
            // set pressures
            //////
            Scalar pC[numPhases];
            const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
            MaterialLaw::capillaryPressures(pC, matParams, fs);

            Scalar densityL = FluidSystem::Brine::gasDensity(temperature_, Scalar(pressure_));

            Scalar depth = pos[dim - 1];
            // add a slight underpressure to start the fingers
            Scalar pl = pressure_ + densityL * this->gravity()[dim - 1] * depth;

            fs.setPressure(liquidPhaseIdx, pl + (pC[liquidPhaseIdx] - pC[liquidPhaseIdx]));
            fs.setPressure(gasPhaseIdx, pl + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));

            //////
            // set composition of the liquid phase
            //////

            fs.setMoleFraction(liquidPhaseIdx, CO2Idx, 0.00000);
            fs.setMoleFraction(liquidPhaseIdx, BrineIdx,
                               1.0 - fs.moleFraction(liquidPhaseIdx, CO2Idx));

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
            CFRP::solve(fs, paramCache,
                        /*refPhaseIdx=*/liquidPhaseIdx,
                        /*setViscosity=*/true,
                        /*setEnthalpy=*/true);

        }
    }

    bool onTopCell_(const GlobalPosition& pos) const
    { return pos[1] > -eps_; }

    bool onTop_(const GlobalPosition& pos) const
    { return pos[1] > -eps_*0.01; }

    void computeThermalCondParams_(ThermalConductionLawParams& params, Scalar poro)
    {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet = std::pow(lambdaGranite, (1 - poro))
                           * std::pow(lambdaWater, poro);
        Scalar lambdaDry = std::pow(lambdaGranite, (1 - poro));

        params.setFullySaturatedLambda(gasPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(liquidPhaseIdx, lambdaWet);
        params.setVacuumLambda(lambdaDry);
    }

    //bool isFineMaterial_(const GlobalPosition& pos) const
    //{ return pos[dim - 1] > fineLayerBottom_; }

    DimMatrix K_;
    std::vector<Scalar> porosity_;
    std::vector<Scalar> molEps_;

    MaterialLawParams materialParams_;

    ThermalConductionLawParams thermalCondParams_;
    SolidEnergyLawParams solidEnergyLawParams_;

    Scalar temperature_;
    Scalar maxDepth_;
    Scalar eps_;
    Scalar pressure_;

    unsigned nTemperature_;
    unsigned nPressure_;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;

    mutable std::mt19937 rand_gen;
    mutable std::normal_distribution<double> norm_dist{0., 1e-4 / 3.};
};
} // namespace Opm

#endif
