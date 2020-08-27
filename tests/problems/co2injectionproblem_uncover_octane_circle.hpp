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

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/Unused.hpp>
#include "twophasefluidsystem.hh"

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
#include <fstream>


namespace Opm {
//! \cond SKIP_THIS
template <class TypeTag>
class Co2InjectionProblem;

namespace Co2Injection {
#include <opm/material/components/co2tables.inc>
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
NEW_PROP_TAG(EpisodeLength);

NEW_PROP_TAG(FluxOutputFileName);

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

public:
    typedef Opm::TwoPhaseCo2OctaneFluidSystem<Scalar> type;
};

// Set the material Law
SET_PROP(Co2InjectionBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
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

// Newton convergence tolerance --newton-tolerance=
SET_SCALAR_PROP(Co2InjectionBaseProblem, NewtonTolerance, 1e-6);

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
SET_STRING_PROP(Co2InjectionBaseProblem, SimulationName, "uncover_circle");

// The default for the end time of the simulation is 2 hours
SET_SCALAR_PROP(Co2InjectionBaseProblem, EndTime, 2*3600);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(Co2InjectionBaseProblem, InitialTimeStepSize, 1e-3);

// The default DGF file to load
SET_STRING_PROP(Co2InjectionBaseProblem, GridFile, "data/co2injectionBox.dgf");

// write restart for every hour
SET_SCALAR_PROP(Co2InjectionBaseProblem, EpisodeLength, 60. * 60.);

SET_STRING_PROP(Co2InjectionBaseProblem, FluxOutputFileName, "uncover_swelling.txt");

END_PROPERTIES

namespace Opm {

/*!
 * \ingroup TestProblems
 *
 * \brief Problem where \f$CO_2\f$ is mixing in octane 
 * The domain is a sized 0.2m times 0.1m and consists of a permeable homogenious porous media
 * within the half-circle and inpermable on the outside of the half-circle
 * \f$CO_2\f$ gets injected by means of a diffusion at the top boundary 
 *
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
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { CO2Idx = FluidSystem::CO2Idx };
    enum { OctaneIdx = FluidSystem::OctaneIdx };
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
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
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
    {     	
        Scalar epi_len = EWOMS_GET_PARAM(TypeTag, Scalar, EpisodeLength);
	    simulator.setEpisodeLength(epi_len);
    }

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
        outputFileName_ = EWOMS_GET_PARAM(TypeTag, std::string, FluxOutputFileName);
        

        // initialize the tables of the fluid system
        FluidSystem::init();

        // intrinsic permeabilities
        //Scalar perm = 76*10; //in darcy
        Scalar perm = 76; //in darcy
        K_ = this->toDimMatrix_(perm * 9.8692 * 1e-13);
		
        //K_ = this->toDimMatrix_(2.11e3 * 9.8692 * 1e-13);
        // impermeable at the outside of the domain
        KK_ = this->toDimMatrix_(perm * 9.8692 * 1e-13 * 1e-6);


        // porosities
        size_t numDof = this->model().numGridDof();

        //porosity_.resize(numDof,0.4);
        porosity_.resize(numDof,0.4);

        //some noise to generate fingers
        for (size_t i = 0; i < numDof; ++i) {
            double noise = this->norm_dist(this->rand_gen);
            porosity_[i] += 10 * noise;
            porosity_[i] = std::min(porosity_[i], 1.0);
        }

        molEps_.resize(numDof, 0.0);
        for (size_t i = 0; i < numDof; ++i) {
            double noise = this->norm_dist(this->rand_gen);
            molEps_[i] += noise;
        }

        // residual saturations
        materialParams_.setResidualSaturation(oilPhaseIdx, 0.2);
        materialParams_.setResidualSaturation(gasPhaseIdx, 0.0);

        // parameters for the Brooks-Corey law
        materialParams_.setEntryPressure(5e3);
        materialParams_.setLambda(0.5);

        materialParams_.finalize();

        // assume constant heat capacity and granite (Not in use) 
        solidEnergyLawParams_.setSolidHeatCapacity(790.0 // specific heat capacity of granite [J / (kg K)]
                                                   * 2700.0); // density of granite [kg/m^3]
        solidEnergyLawParams_.finalize();

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
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EpisodeLength,
                             "Time interval [s] for episode length");
                             
        EWOMS_REGISTER_PARAM(TypeTag, std::string, FluxOutputFileName,
                             "Name of the output file for where the out flux is stored");
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
    
    // This method must be overridden for the simulator to continue with
    // a new episode. We just start a new episode with the same length as
    // the old one.
    void endEpisode() {
	    Scalar epi_len = EWOMS_GET_PARAM(TypeTag, Scalar, EpisodeLength);
	    this->simulator().startNextEpisode(epi_len);
    }

    // only write output when episodes change, aka. report steps, and
    // include the initial timestep too
    bool shouldWriteOutput() {
        
        // set to false to only write at end of episodes
        if (true)
            return true; // write all 
	    else 
            return this->simulator().episodeWillBeOver()
		    || (this->simulator().timeStepIndex() == -1);
    }

    // we don't care about doing restarts from every fifth timestep, it
    // will just slow us down
    bool shouldWriteRestartFile() {
	    return false;
    }
	void beginTimeStep()
	{ 	
	//oilFlux_ = 0.0;
	//gasFlux_ = 0.0;
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
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (!inCircle_(pos)) {
            return KK_;
        } else {
            return K_;
        }
    }

    /// Constant porosity
    template <class Context>
    Scalar porosity(const Context& context ,
                    unsigned spaceIdx,
                    unsigned timeIdx) const
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);
        if (!inCircle_(pos)) {
            return 1e-6;
        } else {
            unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
            return porosity_[globalSpaceIdx];
        }
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
            fs.setSaturation(FluidSystem::oilPhaseIdx, 1.0);
            fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);

            Scalar pC[numPhases];
            const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
            MaterialLaw::capillaryPressures(pC, matParams, fs);

            // The density of octane at (p = 150 bar, T=50C)
            Scalar densityL = 670.25;
            Scalar densityG = 384.4;

            Scalar depth = pos[dim - 1] - 0.1; //depth in gas zone
            Scalar pl = pressure_ + this->gravity()[dim - 1] * ( densityG * depth + densityL * 0.1);
            fs.setPressure(oilPhaseIdx, pl + (pC[oilPhaseIdx] - pC[oilPhaseIdx]));
            fs.setPressure(gasPhaseIdx, pl + (pC[gasPhaseIdx] - pC[oilPhaseIdx]));

            //////
            // set composition of the liquid phase
            //////
            fs.setMoleFraction(oilPhaseIdx, CO2Idx, 1.0);
            fs.setMoleFraction(oilPhaseIdx, OctaneIdx,
                               1.0 - fs.moleFraction(oilPhaseIdx, CO2Idx));

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
            CFRP::solve(fs, paramCache,
                        /*refPhaseIdx=*/oilPhaseIdx,
                        /*setViscosity=*/true,
                        /*setEnthalpy=*/false);

            fs.checkDefined();

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else
            // no flow on top and bottom
            values.setNoFlow();
    }

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        const auto& matParams = this->materialLawParams(context, spaceIdx,
        timeIdx);
        values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
        //values.assignNaive(fs);
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

private:
    template <class Context, class FluidState>
    void initialFluidState_(FluidState& fs,
                            const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);

        if (onTopCircle_(pos)) {
            //////
            // set temperature
            //////
            fs.setTemperature(temperature(context, spaceIdx, timeIdx));

            //////
            // set saturations
            //////
            fs.setSaturation(FluidSystem::oilPhaseIdx, 1.0);
            fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);

            //////
            // set pressures
            //////
            Scalar pC[numPhases];
            const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
            MaterialLaw::capillaryPressures(pC, matParams, fs);

            // The density of octane at (p = 150 bar, T=50C)
            Scalar densityL = 670.25;
            Scalar densityG = 384.4;

            Scalar depth = pos[dim - 1] - 0.1; //depth in gas zone
            Scalar pl = pressure_ + this->gravity()[dim - 1] * ( densityG * depth + densityL * 0.1);
            fs.setPressure(oilPhaseIdx, pl + (pC[oilPhaseIdx] - pC[oilPhaseIdx]));
            fs.setPressure(gasPhaseIdx, pl + (pC[gasPhaseIdx] - pC[oilPhaseIdx]));

            //////
            // set composition of the liquid phase
            //////
            //Scalar co2molfrac = 0.99999;
            //if (depth < 0.05)
            //    co2molfrac = 0.00001;
            Scalar co2molfrac = 1.0 - 0.001;
            fs.setMoleFraction(oilPhaseIdx, CO2Idx, co2molfrac);
            fs.setMoleFraction(oilPhaseIdx, OctaneIdx,
                               1.0 - fs.moleFraction(oilPhaseIdx, CO2Idx));

            typename FluidSystem::template ParameterCache<Scalar> paramCache;
            typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
            CFRP::solve(fs, paramCache,
                        /*refPhaseIdx=*/oilPhaseIdx,
                        /*setViscosity=*/true,
                        /*setEnthalpy=*/false);
            return;
        }

        //////
        // set temperature
        //////
        fs.setTemperature(temperature(context, spaceIdx, timeIdx));

        //////
        // set saturations
        //////
        fs.setSaturation(FluidSystem::oilPhaseIdx, 1.0);
        fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);

        //////
        // set pressures
        //////
        Scalar pC[numPhases];
        const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        // The density of octane at (p = 150 bar, T=50C)
        Scalar densityL = 670.25;

        Scalar depth = pos[dim - 1];
        Scalar pl = pressure_ + densityL * this->gravity()[dim - 1] * depth;
        fs.setPressure(oilPhaseIdx, pl + (pC[oilPhaseIdx] - pC[oilPhaseIdx]));
        fs.setPressure(gasPhaseIdx, pl + (pC[gasPhaseIdx] - pC[oilPhaseIdx]));

        //////
        // set composition of the liquid phase
        //////
        //Scalar co2molfrac = 0.99999;
        //if (depth < 0.05)
        //    co2molfrac = 0.00001;
        Scalar co2molfrac = 0.001;
        fs.setMoleFraction(oilPhaseIdx, CO2Idx, co2molfrac);
        fs.setMoleFraction(oilPhaseIdx, OctaneIdx,
                           1.0 - fs.moleFraction(oilPhaseIdx, CO2Idx));

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem> CFRP;
        CFRP::solve(fs, paramCache,
                    /*refPhaseIdx=*/oilPhaseIdx,
                    /*setViscosity=*/true,
                    /*setEnthalpy=*/false);
    }
    bool onTopCircle_(const GlobalPosition& pos) const
    { return pos[1] > 0.1; }

    bool onTop_(const GlobalPosition& pos) const
    { return pos[1] > 0.2 - eps_*0.01; }

    bool inCircle_(const GlobalPosition& pos) const
    {
        return true; // change this to always return true if you want square
        //const std::vector<Scalar> origo = {0.1,0.1};
        //return ((pos[0]-origo[0])*(pos[0]-origo[0]) + (pos[1]-origo[1])*(pos[1]-origo[1]) ) < 0.01;
    }

    void computeThermalCondParams_(ThermalConductionLawParams& params, Scalar poro)
    {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet = std::pow(lambdaGranite, (1 - poro))
                           * std::pow(lambdaWater, poro);
        Scalar lambdaDry = std::pow(lambdaGranite, (1 - poro));

        params.setFullySaturatedLambda(gasPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(oilPhaseIdx, lambdaWet);
        params.setVacuumLambda(lambdaDry);
    }

    DimMatrix K_;
    DimMatrix KK_;
    std::vector<Scalar> porosity_;
    std::vector<Scalar> molEps_;

    mutable std::vector<Scalar> fluxesTotals_;
    mutable std::vector<Scalar> oilFluxes_;
    mutable std::vector<Scalar> gasFluxes_;

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
    
    std::string outputFileName_;

    mutable std::mt19937 rand_gen;
    mutable std::normal_distribution<double> norm_dist{0., 1e-4 / 3.};
};
} // namespace Opm

#endif
