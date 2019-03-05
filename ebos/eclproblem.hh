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
 * \copydoc Ewoms::EclProblem
 */
#ifndef EWOMS_ECL_PROBLEM_HH
#define EWOMS_ECL_PROBLEM_HH

//#define DISABLE_ALUGRID_SFC_ORDERING 1
//#define EBOS_USE_ALUGRID 1

// make sure that the EBOS_USE_ALUGRID macro. using the preprocessor for this is slightly
// hacky...
#if EBOS_USE_ALUGRID
//#define DISABLE_ALUGRID_SFC_ORDERING 1
#if !HAVE_DUNE_ALUGRID
#warning "ALUGrid was indicated to be used for the ECL black oil simulator, but this "
#warning "requires the presence of dune-alugrid >= 2.4. Falling back to Dune::CpGrid"
#undef EBOS_USE_ALUGRID
#define EBOS_USE_ALUGRID 0
#endif
#else
#define EBOS_USE_ALUGRID 0
#endif

#if EBOS_USE_ALUGRID
#include "eclalugridvanguard.hh"
#else
//#include "eclpolyhedralgridvanguard.hh"
#include "eclcpgridvanguard.hh"
#endif
#include "eclwellmanager.hh"
#include "eclequilinitializer.hh"
#include "eclwriter.hh"
#include "ecloutputblackoilmodule.hh"
#include "ecltransmissibility.hh"
#include "eclthresholdpressure.hh"
#include "ecldummygradientcalculator.hh"
#include "eclfluxmodule.hh"
#include "eclbaseaquifermodel.hh"
#include "eclnewtonmethod.hh"
#include "ecltracermodel.hh"
#include "vtkecltracermodule.hh"

#include <ewoms/common/pffgridvector.hh>
#include <ewoms/models/blackoil/blackoilmodel.hh>
#include <ewoms/disc/ecfv/ecfvdiscretization.hh>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/thermal/EclThermalLawManager.hpp>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WetGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DeadOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityWaterPvt.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Eqldims.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/ConditionalStorage.hpp>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <opm/output/eclipse/EclipseIO.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <boost/date_time.hpp>

#include <set>
#include <vector>
#include <string>
#include <algorithm>

namespace Ewoms {
template <class TypeTag>
class EclProblem;
}

BEGIN_PROPERTIES

#if EBOS_USE_ALUGRID
NEW_TYPE_TAG(EclBaseProblem, INHERITS_FROM(EclAluGridVanguard, EclOutputBlackOil));
#else
NEW_TYPE_TAG(EclBaseProblem, INHERITS_FROM(EclCpGridVanguard, EclOutputBlackOil, VtkEclTracer));
//NEW_TYPE_TAG(EclBaseProblem, INHERITS_FROM(EclPolyhedralGridVanguard, EclOutputBlackOil));
#endif

// The class which deals with ECL wells
NEW_PROP_TAG(EclWellModel);

// Write all solutions for visualization, not just the ones for the
// report steps...
NEW_PROP_TAG(EnableWriteAllSolutions);

// The number of time steps skipped between writing two consequtive restart files
NEW_PROP_TAG(RestartWritingInterval);

// Enable the additional checks even if compiled in debug mode (i.e., with the NDEBUG
// macro undefined). Next to a slightly better performance, this also eliminates some
// print statements in debug mode.
NEW_PROP_TAG(EnableDebuggingChecks);

// if thermal flux boundaries are enabled an effort is made to preserve the initial
// thermal gradient specified via the TEMPVD keyword
NEW_PROP_TAG(EnableThermalFluxBoundaries);

// The class which deals with ECL aquifers
NEW_PROP_TAG(EclAquiferModel);

// time stepping parameters
NEW_PROP_TAG(EclMaxTimeStepSizeAfterWellEvent);
NEW_PROP_TAG(EclRestartShrinkFactor);
NEW_PROP_TAG(EclMaxFails);
NEW_PROP_TAG(EclEnableTuning);

// Set the problem property
SET_TYPE_PROP(EclBaseProblem, Problem, Ewoms::EclProblem<TypeTag>);

// Select the element centered finite volume method as spatial discretization
SET_TAG_PROP(EclBaseProblem, SpatialDiscretizationSplice, EcfvDiscretization);

//! for ebos, use automatic differentiation to linearize the system of PDEs
SET_TAG_PROP(EclBaseProblem, LocalLinearizerSplice, AutoDiffLocalLinearizer);

// Set the material law for fluid fluxes
SET_PROP(EclBaseProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Opm::ThreePhaseMaterialTraits<Scalar,
                                          /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                          /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                          /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

public:
    typedef Opm::EclMaterialLawManager<Traits> EclMaterialLawManager;

    typedef typename EclMaterialLawManager::MaterialLaw type;
};

// Set the material law for energy storage in rock
SET_PROP(EclBaseProblem, SolidEnergyLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    typedef Opm::EclThermalLawManager<Scalar, FluidSystem> EclThermalLawManager;

    typedef typename EclThermalLawManager::SolidEnergyLaw type;
};

// Set the material law for thermal conduction
SET_PROP(EclBaseProblem, ThermalConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    typedef Opm::EclThermalLawManager<Scalar, FluidSystem> EclThermalLawManager;

    typedef typename EclThermalLawManager::ThermalConductionLaw type;
};

// ebos can use a slightly faster stencil class because it does not need the normals and
// the integration points of intersections
SET_PROP(EclBaseProblem, Stencil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    typedef Ewoms::EcfvStencil<Scalar,
                               GridView,
                               /*needIntegrationPos=*/false,
                               /*needNormal=*/false> type;
};

// by default use the dummy aquifer "model"
SET_TYPE_PROP(EclBaseProblem, EclAquiferModel, Ewoms::EclBaseAquiferModel<TypeTag>);

// use the built-in proof of concept well model by default
SET_TYPE_PROP(EclBaseProblem, EclWellModel, EclWellManager<TypeTag>);

// Enable gravity
SET_BOOL_PROP(EclBaseProblem, EnableGravity, true);

// only write the solutions for the report steps to disk
SET_BOOL_PROP(EclBaseProblem, EnableWriteAllSolutions, false);

// The default for the end time of the simulation [s]
//
// By default, stop it after the universe will probably have stopped
// to exist. (the ECL problem will finish the simulation explicitly
// after it simulated the last episode specified in the deck.)
SET_SCALAR_PROP(EclBaseProblem, EndTime, 1e100);

// The default for the initial time step size of the simulation [s].
//
// The chosen value means that the size of the first time step is the
// one of the initial episode (if the length of the initial episode is
// not millions of trillions of years, that is...)
SET_SCALAR_PROP(EclBaseProblem, InitialTimeStepSize, 3600*24);

// the default for the allowed volumetric error for oil per second
SET_SCALAR_PROP(EclBaseProblem, NewtonRawTolerance, 1e-2);

// the tolerated amount of "incorrect" amount of oil per time step for the complete
// reservoir. this is scaled by the pore volume of the reservoir, i.e., larger reservoirs
// will tolerate larger residuals.
SET_SCALAR_PROP(EclBaseProblem, EclNewtonSumTolerance, 1e-4);

// set the exponent for the volume scaling of the sum tolerance: larger reservoirs can
// tolerate a higher amount of mass lost per time step than smaller ones! since this is
// not linear, we use the cube root of the overall pore volume by default, i.e., the
// value specified by the NewtonSumTolerance parameter is the "incorrect" mass per
// timestep for an reservoir that exhibits 1 m^3 of pore volume. A reservoir with a total
// pore volume of 10^3 m^3 will tolerate 10 times as much.
SET_SCALAR_PROP(EclBaseProblem, EclNewtonSumToleranceExponent, 1.0/3.0);

// set number of Newton iterations where the volumetric residual is considered for
// convergence
SET_INT_PROP(EclBaseProblem, EclNewtonStrictIterations, 8);

// set fraction of the pore volume where the volumetric residual may be violated during
// strict Newton iterations
SET_SCALAR_PROP(EclBaseProblem, EclNewtonRelaxedVolumeFraction, 0.03);

// the maximum volumetric error of a cell in the relaxed region
SET_SCALAR_PROP(EclBaseProblem, EclNewtonRelaxedTolerance, 1e9);

// Ignore the maximum error mass for early termination of the newton method.
SET_SCALAR_PROP(EclBaseProblem, NewtonMaxError, 10e9);

// set the maximum number of Newton iterations to 14 because the likelyhood that a time
// step succeeds at more than 14 Newton iteration is rather small
SET_INT_PROP(EclBaseProblem, NewtonMaxIterations, 14);

// also, reduce the target for the "optimum" number of Newton iterations to 6. Note that
// this is only relevant if the time step is reduced from the report step size for some
// reason. (because ebos first tries to do a report step using a single time step.)
SET_INT_PROP(EclBaseProblem, NewtonTargetIterations, 6);

// Disable the VTK output by default for this problem ...
SET_BOOL_PROP(EclBaseProblem, EnableVtkOutput, false);

// ... but enable the ECL output by default
SET_BOOL_PROP(EclBaseProblem, EnableEclOutput, true);

// If available, write the ECL output in a non-blocking manner
SET_BOOL_PROP(EclBaseProblem, EnableAsyncEclOutput, true);

// By default, use single precision for the ECL formated results
SET_BOOL_PROP(EclBaseProblem, EclOutputDoublePrecision, false);

// The default location for the ECL output files
SET_STRING_PROP(EclBaseProblem, OutputDir, ".");

// the cache for intensive quantities can be used for ECL problems and also yields a
// decent speedup...
SET_BOOL_PROP(EclBaseProblem, EnableIntensiveQuantityCache, true);

// the cache for the storage term can also be used and also yields a decent speedup
SET_BOOL_PROP(EclBaseProblem, EnableStorageCache, true);

// Use the "velocity module" which uses the Eclipse "NEWTRAN" transmissibilities
SET_TYPE_PROP(EclBaseProblem, FluxModule, Ewoms::EclTransFluxModule<TypeTag>);

// Use the dummy gradient calculator in order not to do unnecessary work.
SET_TYPE_PROP(EclBaseProblem, GradientCalculator, Ewoms::EclDummyGradientCalculator<TypeTag>);

// Use a custom Newton-Raphson method class for ebos in order to attain more
// sophisticated update and error computation mechanisms
SET_TYPE_PROP(EclBaseProblem, NewtonMethod, Ewoms::EclNewtonMethod<TypeTag>);

// The frequency of writing restart (*.ers) files. This is the number of time steps
// between writing restart files
SET_INT_PROP(EclBaseProblem, RestartWritingInterval, 0xffffff); // disable

// By default, we enable the debugging checks if we're compiled in debug mode
SET_BOOL_PROP(EclBaseProblem, EnableDebuggingChecks, true);

// store temperature (but do not conserve energy, as long as EnableEnergy is false)
SET_BOOL_PROP(EclBaseProblem, EnableTemperature, true);

// disable all extensions supported by black oil model. this should not really be
// necessary but it makes things a bit more explicit
SET_BOOL_PROP(EclBaseProblem, EnablePolymer, false);
SET_BOOL_PROP(EclBaseProblem, EnableSolvent, true);
SET_BOOL_PROP(EclBaseProblem, EnableEnergy, false);

// disable thermal flux boundaries by default
SET_BOOL_PROP(EclBaseProblem, EnableThermalFluxBoundaries, false);

SET_BOOL_PROP(EclBaseProblem, EnableTracerModel, false);

// By default, simulators derived from the EclBaseProblem are production simulators,
// i.e., experimental features must be explicitly enabled at compile time
SET_BOOL_PROP(EclBaseProblem, EnableExperiments, false);

// set defaults for the time stepping parameters
SET_SCALAR_PROP(EclBaseProblem, EclMaxTimeStepSizeAfterWellEvent, 3600*24*365.25);
SET_SCALAR_PROP(EclBaseProblem, EclRestartShrinkFactor, 3);
SET_INT_PROP(EclBaseProblem, EclMaxFails, 10);
SET_BOOL_PROP(EclBaseProblem, EclEnableTuning, false);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This problem simulates an input file given in the data format used by the
 *        commercial ECLiPSE simulator.
 */
template <class TypeTag>
class EclProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Stencil) Stencil;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // Grid and world dimension
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };
    enum { enablePolymer = GET_PROP_VALUE(TypeTag, EnablePolymer) };
    enum { enablePolymerMolarWeight = GET_PROP_VALUE(TypeTag, EnablePolymerMW) };
    enum { enableTemperature = GET_PROP_VALUE(TypeTag, EnableTemperature) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { enableThermalFluxBoundaries = GET_PROP_VALUE(TypeTag, EnableThermalFluxBoundaries) };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP(TypeTag, MaterialLaw)::EclMaterialLawManager EclMaterialLawManager;
    typedef typename GET_PROP(TypeTag, SolidEnergyLaw)::EclThermalLawManager EclThermalLawManager;
    typedef typename EclMaterialLawManager::MaterialLawParams MaterialLawParams;
    typedef typename EclThermalLawManager::SolidEnergyLawParams SolidEnergyLawParams;
    typedef typename EclThermalLawManager::ThermalConductionLawParams ThermalConductionLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, EclWellModel) EclWellModel;
    typedef typename GET_PROP_TYPE(TypeTag, EclAquiferModel) EclAquiferModel;

    typedef BlackOilSolventModule<TypeTag> SolventModule;
    typedef BlackOilPolymerModule<TypeTag> PolymerModule;

    typedef typename EclEquilInitializer<TypeTag>::ScalarFluidState InitialFluidState;

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    typedef EclWriter<TypeTag> EclWriterType;

    typedef EclTracerModel<TypeTag> TracerModel;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    struct RockParams {
        Scalar referencePressure;
        Scalar compressibility;
    };

public:
    /*!
     * \copydoc FvBaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
        EclWriterType::registerParameters();
        VtkEclTracerModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableWriteAllSolutions,
                             "Write all solutions to disk instead of only the ones for the "
                             "report steps");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableEclOutput,
                             "Write binary output which is compatible with the commercial "
                             "Eclipse simulator");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclOutputDoublePrecision,
                             "Tell the output writer to use double precision. Useful for 'perfect' restarts");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, RestartWritingInterval,
                             "The frequencies of which time steps are serialized to disk");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableTracerModel,
                             "Transport tracers found in the deck.");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclMaxTimeStepSizeAfterWellEvent,
                             "Maximum time step size after an well event");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclRestartShrinkFactor,
                             "Factor by which the time step is reduced after convergence failure");
        EWOMS_REGISTER_PARAM(TypeTag, int, EclMaxFails,
                             "Maximum consecutive convergence failures before termination");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclEnableTuning,
                             "Honor some aspects of the TUNING keyword from the ECL deck.");
    }

    /*!
     * \copydoc FvBaseProblem::prepareOutputDir
     */
    std::string prepareOutputDir() const
    { return this->simulator().vanguard().eclState().getIOConfig().getOutputDir(); }

    /*!
     * \copydoc FvBaseProblem::handlePositionalParameter
     */
    static int handlePositionalParameter(std::set<std::string>& seenParams,
                                         std::string& errorMsg,
                                         int argc OPM_UNUSED,
                                         const char** argv,
                                         int paramIdx,
                                         int posParamIdx OPM_UNUSED)
    {
        typedef typename GET_PROP(TypeTag, ParameterMetaData) ParamsMeta;
        Dune::ParameterTree& tree = ParamsMeta::tree();

        std::string param  = argv[paramIdx];
        size_t i = param.find('=');
        if (i != std::string::npos) {
            std::string oldParamName = param.substr(0, i);
            std::string oldParamValue = param.substr(i+1);
            std::string newParamName = "--" + oldParamName;
            for (size_t j = 0; j < newParamName.size(); ++j)
                if (newParamName[j] == '_')
                    newParamName[j] = '-';
            errorMsg =
                "The old syntax to specify parameters on the command line is no longer supported: "
                "Try replacing '"+oldParamName+"="+oldParamValue+"' with "+
                "'"+newParamName+"="+oldParamValue+"'!";
            return 0;
        }

        if (seenParams.count("EclDeckFileName") > 0) {
            errorMsg =
                "Parameter 'EclDeckFileName' specified multiple times"
                " as a command line parameter";
            return 0;
        }

        tree["EclDeckFileName"] = argv[paramIdx];
        seenParams.insert("EclDeckFileName");
        return 1;
    }

    /*!
     * \copydoc FvBaseProblem::helpPreamble
     */
    static std::string helpPreamble(int argc OPM_UNUSED,
                                    const char **argv)
    {
        std::string desc = Implementation::briefDescription();
        if (!desc.empty())
            desc = desc + "\n";

        return
            "Usage: "+std::string(argv[0]) + " [OPTIONS] [ECL_DECK_FILENAME]\n"
            + desc;
    }

    /*!
     * \copydoc FvBaseProblem::briefDescription
     */
    static std::string briefDescription()
    {
        if (briefDescription_.empty())
            return
                "The Ecl-deck Black-Oil reservoir Simulator (ebos); a hydrocarbon "
                "reservoir simulation program that processes ECL-formatted input "
                "files and is provided by the Open Porous Media project "
                "(https://opm-project.org).\n"
                "\n"
                "THE `ebos` SIMULATOR IS FOR RESEARCH PURPOSES ONLY! For industrial "
                "applications, use `flow`.";
        else
            return briefDescription_;
    }

    /*!
     * \brief Specifies the string returned by briefDescription()
     *
     * This string appears in the usage message.
     */
    static void setBriefDescription(const std::string& msg)
    { briefDescription_ = msg; }

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    EclProblem(Simulator& simulator)
        : ParentType(simulator)
        , transmissibilities_(simulator.vanguard())
        , thresholdPressures_(simulator)
        , wellModel_(simulator)
        , aquiferModel_(simulator)
        , pffDofData_(simulator.gridView(), this->elementMapper())
        , tracerModel_(simulator)
    {
        this->model().addOutputModule(new VtkEclTracerModule<TypeTag>(simulator));
        // Tell the black-oil extensions to initialize their internal data structures
        const auto& vanguard = simulator.vanguard();
        SolventModule::initFromDeck(vanguard.deck(), vanguard.eclState());
        PolymerModule::initFromDeck(vanguard.deck(), vanguard.eclState());
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableEclOutput))
            // create the ECL writer
            eclWriter_.reset(new EclWriterType(simulator));

        enableTuning_ = EWOMS_GET_PARAM(TypeTag, bool, EclEnableTuning);
        initialTimeStepSize_ = EWOMS_GET_PARAM(TypeTag, Scalar, InitialTimeStepSize);
        minTimeStepSize_ = EWOMS_GET_PARAM(TypeTag, Scalar, MinTimeStepSize);
        maxTimeStepSize_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxTimeStepSize);
        maxTimeStepAfterWellEvent_ = EWOMS_GET_PARAM(TypeTag, Scalar, EclMaxTimeStepSizeAfterWellEvent);
        restartShrinkFactor_ = EWOMS_GET_PARAM(TypeTag, Scalar, EclRestartShrinkFactor);
        maxFails_ = EWOMS_GET_PARAM(TypeTag, int, EclMaxFails);
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        auto& simulator = this->simulator();

        // the "NOGRAV" keyword from Frontsim or setting the EnableGravity to false
        // disables gravity, else the standard value of the gravity constant at sea level
        // on earth is used
        this->gravity_ = 0.0;
        const auto& deck = simulator.vanguard().deck();
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity))
            this->gravity_[dim - 1] = 9.80665;
        if (deck.hasKeyword("NOGRAV"))
            this->gravity_[dim - 1] = 0.0;

        if (enableTuning_) {
            // if support for the TUNING keyword is enabled, we get the initial time
            // steping parameters from it instead of from command line parameters
            const auto& schedule = simulator.vanguard().schedule();
            const auto& tuning = schedule.getTuning();
            initialTimeStepSize_ = tuning.getTSINIT(0);
            maxTimeStepAfterWellEvent_ = tuning.getTMAXWC(0);
            maxTimeStepSize_ = tuning.getTSMAXZ(0);
            restartShrinkFactor_ = tuning.getTSFCNV(0);
            minTimeStepSize_ = tuning.getTSMINZ(0);
        }

        // deal with DRSDT
        const auto& eclState = simulator.vanguard().eclState();
        unsigned ntpvt = eclState.runspec().tabdims().getNumPVTTables();
        maxDRs_.resize(ntpvt, 1e30);
        dRsDtOnlyFreeGas_.resize(ntpvt, false);
        size_t numDof = this->model().numGridDof();
        lastRs_.resize(numDof, 0.0);
        maxDRv_.resize(ntpvt, 1e30);
        lastRv_.resize(numDof, 0.0);
        maxOilSaturation_.resize(numDof, 0.0);

        initFluidSystem_();
        updateElementDepths_();
        readRockParameters_();
        readMaterialParameters_();
        readThermalParameters_();
        transmissibilities_.finishInit();

        const auto& timeMap = simulator.vanguard().schedule().getTimeMap();

        readInitialCondition_();
        // Set the start time of the simulation
        simulator.setStartTime(timeMap.getStartTime(/*timeStepIdx=*/0));

        // We want the episode index to be the same as the report step index to make
        // things simpler, so we have to set the episode index to -1 because it is
        // incremented inside beginEpisode(). The size of the initial time step and
        // length of the initial episode is set to zero for the same reason.
        simulator.setEpisodeIndex(-1);
        simulator.setEpisodeLength(0.0);
        simulator.setTimeStepSize(0.0);

        updatePffDofData_();

        if (GET_PROP_VALUE(TypeTag, EnablePolymer)) {
            const auto& vanguard = this->simulator().vanguard();
            const auto& gridView = vanguard.gridView();
            int numElements = gridView.size(/*codim=*/0);
            maxPolymerAdsorption_.resize(numElements, 0.0);
        }

        if (eclWriter_) {
            eclWriter_->writeInit();
            this->simulator().vanguard().releaseGlobalTransmissibilities();
        }

        tracerModel_.init();

        readBoundaryConditions_();
    }

    void prefetch(const Element& elem) const
    { pffDofData_.prefetch(elem); }

    /*!
     * \brief This method restores the complete state of the problem and its sub-objects
     *        from disk.
     *
     * The serialization format used by this method is ad-hoc. It is the inverse of the
     * serialize() method.
     *
     * \tparam Restarter The deserializer type
     *
     * \param res The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter& res)
    {
        // reload the current episode/report step from the deck
        beginEpisode(/*isOnRestart=*/true);

        // deserialize the wells
        wellModel_.deserialize(res);

        // deserialize the aquifer
        aquiferModel_.deserialize(res);
    }

    /*!
     * \brief This method writes the complete state of the problem and its subobjects to
     *        disk.
     *
     * The file format used here is ad-hoc.
     */
    template <class Restarter>
    void serialize(Restarter& res)
    {
        wellModel_.serialize(res);
        aquiferModel_.serialize(res);
    }

    /*!
     * \brief Called by the simulator before an episode begins.
     */
    void beginEpisode(bool isOnRestart = false)
    {
        // Proceed to the next report step
        Simulator& simulator = this->simulator();
        auto& eclState = this->simulator().vanguard().eclState();
        const auto& schedule = this->simulator().vanguard().schedule();
        const auto& events = schedule.getEvents();
        const auto& timeMap = schedule.getTimeMap();

        // The first thing to do in the morning of an episode is update update the
        // eclState and the deck if they need to be changed.
        int nextEpisodeIdx = simulator.episodeIndex();
        if (nextEpisodeIdx > 0 &&
            events.hasEvent(Opm::ScheduleEvents::GEO_MODIFIER, nextEpisodeIdx))
        {
            // bring the contents of the keywords to the current state of the SCHEDULE
            // section
            //
            // TODO (?): make grid topology changes possible (depending on what exactly
            // has changed, the grid may need be re-created which has some serious
            // implications on e.g., the solution of the simulation.)
            const auto& miniDeck = schedule.getModifierDeck(nextEpisodeIdx);
            eclState.applyModifierDeck(miniDeck);

            // re-compute all quantities which may possibly be affected.
            transmissibilities_.update();
            updatePorosity_();
            updatePffDofData_();
        }

        // react to TUNING changes
        if (nextEpisodeIdx > 0
            && enableTuning_
            && events.hasEvent(Opm::ScheduleEvents::TUNING_CHANGE, nextEpisodeIdx))
        {
            const auto& tuning = schedule.getTuning();
            initialTimeStepSize_ = tuning.getTSINIT(nextEpisodeIdx);
            maxTimeStepAfterWellEvent_ = tuning.getTMAXWC(nextEpisodeIdx);
            maxTimeStepSize_ = tuning.getTSMAXZ(nextEpisodeIdx);
            restartShrinkFactor_ = tuning.getTSFCNV(nextEpisodeIdx);
            minTimeStepSize_ = tuning.getTSMINZ(nextEpisodeIdx);
        }

        // Opm::TimeMap deals with points in time, so the number of time intervals (i.e.,
        // report steps) is one less!
        int numReportSteps = timeMap.size() - 1;

        // start the next episode if there are additional report steps, else finish the
        // simulation
        while (nextEpisodeIdx < numReportSteps &&
               simulator.time() >= timeMap.getTimePassedUntil(nextEpisodeIdx + 1)*(1 - 1e-10))
        {
            ++ nextEpisodeIdx;
        }

        Scalar episodeLength = timeMap.getTimeStepLength(nextEpisodeIdx);
        simulator.startNextEpisode(episodeLength);

        Scalar dt = limitNextTimeStepSize_(episodeLength);
        if (nextEpisodeIdx == 0)
            // allow the size of the initial time step to be set via an external parameter
            dt = std::min(dt, initialTimeStepSize_);

        simulator.setTimeStepSize(dt);

        const bool invalidateFromHyst = updateHysteresis_();
        const bool invalidateFromMaxOilSat = updateMaxOilSaturation_();
        const bool doInvalidate = invalidateFromHyst || invalidateFromMaxOilSat;

        if (GET_PROP_VALUE(TypeTag, EnablePolymer))
            updateMaxPolymerAdsorption_();

        // set up the wells for the next episode.
        //
        // TODO: the first two arguments seem to be unnecessary
        wellModel_.beginEpisode(this->simulator().vanguard().eclState(),
                                this->simulator().vanguard().schedule(),
                                isOnRestart);

        aquiferModel_.beginEpisode();

        if (doInvalidate)
            this->model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);
    }

    /*!
     * \brief Called by the simulator before each time integration.
     */
    void beginTimeStep()
    {
        int epsiodeIdx = this->simulator().episodeIndex();
        const auto& oilVaporizationControl = this->simulator().vanguard().schedule().getOilVaporizationProperties(epsiodeIdx);
        if (drsdtActive_())
            // DRSDT is enabled
            for (size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRs_.size(); ++pvtRegionIdx)
                maxDRs_[pvtRegionIdx] = oilVaporizationControl.getMaxDRSDT(pvtRegionIdx)*this->simulator().timeStepSize();

        if (drvdtActive_())
            // DRVDT is enabled
            for (size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRv_.size(); ++pvtRegionIdx)
                maxDRv_[pvtRegionIdx] = oilVaporizationControl.getMaxDRVDT(pvtRegionIdx)*this->simulator().timeStepSize();

        wellModel_.beginTimeStep();
        aquiferModel_.beginTimeStep();
        tracerModel_.beginTimeStep();

    }

    /*!
     * \brief Return if the storage term of the first iteration is identical to the storage
     *        term for the solution of the previous time step.
     *
     * For quite technical reasons, the storage term cannot be recycled if either DRSDT
     * or DRVDT are active in ebos.
     */
    bool recycleFirstIterationStorage() const
    { return !drsdtActive_() && !drvdtActive_(); }

    /*!
     * \brief Called by the simulator before each Newton-Raphson iteration.
     */
    void beginIteration()
    {
        wellModel_.beginIteration();
        aquiferModel_.beginIteration();
    }

    /*!
     * \brief Called by the simulator after each Newton-Raphson iteration.
     */
    void endIteration()
    {
        wellModel_.endIteration();
        aquiferModel_.endIteration();
    }

    /*!
     * \brief Called by the simulator after each time integration.
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        if (GET_PROP_VALUE(TypeTag, EnableDebuggingChecks)) {
            // in debug mode, we don't care about performance, so we check if the model does
            // the right thing (i.e., the mass change inside the whole reservoir must be
            // equivalent to the fluxes over the grid's boundaries plus the source rates
            // specified by the problem)
            this->model().checkConservativeness(/*tolerance=*/-1, /*verbose=*/true);
        }
#endif // NDEBUG

        wellModel_.endTimeStep();
        aquiferModel_.endTimeStep();
        tracerModel_.endTimeStep();

        updateCompositionChangeLimits_();
    }

    /*!
     * \brief Called by the simulator after the end of an episode.
     */
    void endEpisode()
    {
        auto& simulator = this->simulator();
        const auto& schedule = simulator.vanguard().schedule();

        int episodeIdx = simulator.episodeIndex();

        const auto& timeMap = schedule.getTimeMap();
        int numReportSteps = timeMap.size() - 1;
        if (episodeIdx + 1 >= numReportSteps) {
            simulator.setFinished(true);
            return;
        }
    }

    /*!
     * \brief Returns true if the current solution should be written
     *        to disk for visualization.
     *
     * For the ECL simulator we only write at the end of
     * episodes/report steps...
     */
    bool shouldWriteOutput() const
    {
        if (this->simulator().timeStepIndex() < 0)
            // always write the initial solution
            return true;

        if (EWOMS_GET_PARAM(TypeTag, bool, EnableWriteAllSolutions))
            return true;

        return this->simulator().episodeWillBeOver();
    }

    /*!
     * \brief Returns true if an eWoms restart file should be written to disk.
     */
    bool shouldWriteRestartFile() const
    {
        unsigned n = EWOMS_GET_PARAM(TypeTag, unsigned, RestartWritingInterval);
        unsigned i = this->simulator().timeStepIndex();
        if (i > 0 && (i%n) == 0)
            return true; // we don't write a restart file for the initial condition
        return false;
    }

    /*!
     * \brief Write the requested quantities of the current solution into the output
     *        files.
     */
    void writeOutput(bool isSubStep, bool verbose = true)
    {
        // use the generic code to prepare the output fields and to
        // write the desired VTK files.
        ParentType::writeOutput(isSubStep, verbose);

        if (!eclWriter_)
            return;

        eclWriter_->writeOutput(isSubStep);
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context,
                                           unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return transmissibilities_.permeability(globalSpaceIdx);
    }

    /*!
     * \brief This method returns the intrinsic permeability tensor
     *        given a global element index.
     *
     * Its main (only?) usage is the ECL transmissibility calculation code...
     */
    const DimMatrix& intrinsicPermeability(unsigned globalElemIdx) const
    { return transmissibilities_.permeability(globalElemIdx); }

    /*!
     * \copydoc EclTransmissiblity::transmissibility
     */
    template <class Context>
    Scalar transmissibility(const Context& context,
                            unsigned OPM_OPTIM_UNUSED fromDofLocalIdx,
                            unsigned toDofLocalIdx) const
    {
        assert(fromDofLocalIdx == 0);
        return pffDofData_.get(context.element(), toDofLocalIdx).transmissibility;
    }

    /*!
     * \copydoc EclTransmissiblity::transmissibilityBoundary
     */
    template <class Context>
    Scalar transmissibilityBoundary(const Context& elemCtx,
                                    unsigned boundaryFaceIdx) const
    {
        unsigned elemIdx = elemCtx.globalSpaceIndex(/*dofIdx=*/0, /*timeIdx=*/0);
        return transmissibilities_.transmissibilityBoundary(elemIdx, boundaryFaceIdx);
    }

    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    template <class Context>
    Scalar thermalHalfTransmissibility(const Context& context,
                                       unsigned faceIdx,
                                       unsigned timeIdx) const
    {
        const auto& face = context.stencil(timeIdx).interiorFace(faceIdx);
        unsigned toDofLocalIdx = face.exteriorIndex();
        return *pffDofData_.get(context.element(), toDofLocalIdx).thermalHalfTrans;
    }

    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    template <class Context>
    Scalar thermalHalfTransmissibilityBoundary(const Context& elemCtx,
                                               unsigned boundaryFaceIdx) const
    {
        unsigned elemIdx = elemCtx.globalSpaceIndex(/*dofIdx=*/0, /*timeIdx=*/0);
        return transmissibilities_.thermalHalfTransBoundary(elemIdx, boundaryFaceIdx);
    }

    /*!
     * \brief Return a reference to the object that handles the "raw" transmissibilities.
     */
    const EclTransmissibility<TypeTag>& eclTransmissibilities() const
    { return transmissibilities_; }

    /*!
     * \copydoc BlackOilBaseProblem::thresholdPressure
     */
    Scalar thresholdPressure(unsigned elem1Idx, unsigned elem2Idx) const
    { return thresholdPressures_.thresholdPressure(elem1Idx, elem2Idx); }

    const EclThresholdPressure<TypeTag>& thresholdPressure() const
    { return thresholdPressures_; }

    EclThresholdPressure<TypeTag>& thresholdPressure()
    { return thresholdPressures_; }

    const EclTracerModel<TypeTag>& tracerModel() const
    { return tracerModel_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return porosity_[globalSpaceIdx];
    }

    /*!
     * \brief Returns the porosity of an element
     *
     * Note that this method is *not* part of the generic eWoms problem API because it
     * would bake the assumption that only the elements are the degrees of freedom into
     * the interface.
     */
    Scalar porosity(unsigned elementIdx) const
    { return porosity_[elementIdx]; }

    /*!
     * \brief Returns the depth of an degree of freedom [m]
     *
     * For ECL problems this is defined as the average of the depth of an element and is
     * thus slightly different from the depth of an element's centroid.
     */
    template <class Context>
    Scalar dofCenterDepth(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return elementCenterDepth_[globalSpaceIdx];
    }

    /*!
     * \copydoc BlackoilProblem::rockCompressibility
     */
    template <class Context>
    Scalar rockCompressibility(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        if (rockParams_.empty())
            return 0.0;

        unsigned tableIdx = 0;
        if (!rockTableIdx_.empty()) {
            unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
            tableIdx = rockTableIdx_[globalSpaceIdx];
        }

        return rockParams_[tableIdx].compressibility;
    }

    /*!
     * \copydoc BlackoilProblem::rockReferencePressure
     */
    template <class Context>
    Scalar rockReferencePressure(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        if (rockParams_.empty())
            return 1e5;

        unsigned tableIdx = 0;
        if (!rockTableIdx_.empty()) {
            unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
            tableIdx = rockTableIdx_[globalSpaceIdx];
        }

        return rockParams_[tableIdx].referencePressure;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return materialLawParams(globalSpaceIdx);
    }

    const MaterialLawParams& materialLawParams(unsigned globalDofIdx) const
    { return materialLawManager_->materialLawParams(globalDofIdx); }

    /*!
     * \brief Return the parameters for the energy storage law of the rock
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyLawParams(const Context& context OPM_UNUSED,
                         unsigned spaceIdx OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return thermalLawManager_->solidEnergyLawParams(globalSpaceIdx);
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams &
    thermalConductionLawParams(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return thermalLawManager_->thermalConductionLawParams(globalSpaceIdx);
    }

    /*!
     * \brief Returns the ECL material law manager
     *
     * Note that this method is *not* part of the generic eWoms problem API because it
     * would force all problens use the ECL material laws.
     */
    std::shared_ptr<const EclMaterialLawManager> materialLawManager() const
    { return materialLawManager_; }

    /*!
     * \copydoc materialLawManager()
     */
    std::shared_ptr<EclMaterialLawManager> materialLawManager()
    { return materialLawManager_; }

    /*!
     * \brief Returns the initial solvent saturation for a given a cell index
     */
    Scalar solventSaturation(unsigned elemIdx) const
    {
        if (solventSaturation_.empty())
            return 0;

        return solventSaturation_[elemIdx];
    }

    /*!
     * \brief Returns the initial polymer concentration for a given a cell index
     */
    Scalar  polymerConcentration(unsigned elemIdx) const
    {
        if (polymerConcentration_.empty())
            return 0;

        return polymerConcentration_[elemIdx];
    }

    /*!
    * \brief Returns the polymer molecule weight for a given cell index
    */
    // TODO: remove this function if not called
    Scalar polymerMolecularWeight(const unsigned elemIdx) const
    {
        if (polymerMoleWeight_.empty())
            return 0.0;

        return polymerMoleWeight_[elemIdx];
    }

    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned pvtRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return pvtRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \brief Returns the index the relevant PVT region given a cell index
     */
    unsigned pvtRegionIndex(unsigned elemIdx) const
    {
        if (pvtnum_.empty())
            return 0;

        return pvtnum_[elemIdx];
    }

    const std::vector<int>& pvtRegionArray() const
    { return pvtnum_; }

    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned satnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return satnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \brief Returns the index the relevant saturation function region given a cell index
     */
    unsigned satnumRegionIndex(unsigned elemIdx) const
    {
        if (satnum_.empty())
            return 0;

        return satnum_[elemIdx];
    }

    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned miscnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return miscnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \brief Returns the index the relevant MISC region given a cell index
     */
    unsigned miscnumRegionIndex(unsigned elemIdx) const
    {
        if (miscnum_.empty())
            return 0;

        return miscnum_[elemIdx];
    }

    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned plmixnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return plmixnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \brief Returns the index the relevant PLMIXNUM (for polymer module) region given a cell index
     */
    unsigned plmixnumRegionIndex(unsigned elemIdx) const
    {
        if (plmixnum_.empty())
            return 0;

        return plmixnum_[elemIdx];
    }

    /*!
     * \brief Returns the max polymer adsorption value
     */
    template <class Context>
    Scalar maxPolymerAdsorption(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return maxPolymerAdsorption(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \brief Returns the max polymer adsorption value
     */
    Scalar maxPolymerAdsorption(unsigned elemIdx) const
    {
        if (maxPolymerAdsorption_.empty())
            return 0;

        return maxPolymerAdsorption_[elemIdx];
    }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return this->simulator().vanguard().caseName(); }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        // use the temporally constant temperature, i.e. use the initial temperature of
        // the DOF
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        if (!initialFluidStates_.empty())
            return initialFluidStates_[globalDofIdx].temperature(/*phaseIdx=*/0);
        return initialTemperature_[globalDofIdx];
    }

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * ECLiPSE uses no-flow conditions for all boundaries. \todo really?
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    {
        if (!enableEnergy || !enableThermalFluxBoundaries)
            values.setNoFlow();
        else {
            // in the energy case we need to specify a non-trivial boundary condition
            // because the geothermal gradient needs to be maintained. for this, we
            // simply assume the initial temperature at the boundary and specify the
            // thermal flow accordingly. in this context, "thermal flow" means energy
            // flow due to a temerature gradient while assuming no-flow for mass
            unsigned interiorDofIdx = context.interiorScvIndex(spaceIdx, timeIdx);
            unsigned globalDofIdx = context.globalSpaceIndex(interiorDofIdx, timeIdx);
            values.setThermalFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
        }

        if (hasFreeBoundaryConditions()) {
            unsigned indexInInside  = context.intersection(spaceIdx).indexInInside();
            unsigned interiorDofIdx = context.interiorScvIndex(spaceIdx, timeIdx);
            unsigned globalDofIdx = context.globalSpaceIndex(interiorDofIdx, timeIdx);
            unsigned pvtRegionIdx = pvtRegionIndex(context, spaceIdx, timeIdx);
            switch (indexInInside) {
            case 0:
                if (freebcXMinus_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcXMinus_[globalDofIdx], pvtRegionIdx);
                break;
            case 1:
                if (freebcX_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcX_[globalDofIdx], pvtRegionIdx);
                break;
            case 2:
                if (freebcYMinus_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcYMinus_[globalDofIdx], pvtRegionIdx);
                break;
            case 3:
                if (freebcY_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcY_[globalDofIdx], pvtRegionIdx);
                break;
            case 4:
                if (freebcZMinus_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcZMinus_[globalDofIdx], pvtRegionIdx);
                break;
            case 5:
                if (freebcZ_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcZ_[globalDofIdx], pvtRegionIdx);
                break;
            default:
                throw std::logic_error("invalid face index for boundary condition");

            }
        }
    }

    /*!
     * \copydoc FvBaseProblem::initial
     *
     * The reservoir problem uses a constant boundary condition for
     * the whole domain.
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

        values.setPvtRegionIndex(pvtRegionIndex(context, spaceIdx, timeIdx));
        values.assignNaive(initialFluidStates_[globalDofIdx]);

        if (enableSolvent)
            values[Indices::solventSaturationIdx] = solventSaturation_[globalDofIdx];

        if (enablePolymer)
            values[Indices::polymerConcentrationIdx] = polymerConcentration_[globalDofIdx];

        if (enablePolymerMolarWeight)
            values[Indices::polymerMoleWeightIdx]= polymerMoleWeight_[globalDofIdx];

        values.checkDefined();
    }

    /*!
     * \copydoc FvBaseProblem::initialSolutionApplied()
     */
    void initialSolutionApplied()
    {
        const auto& eclState = this->simulator().vanguard().eclState();

        // initialize the wells. Note that this needs to be done after initializing the
        // intrinsic permeabilities and the after applying the initial solution because
        // the well model uses these...
        wellModel_.init(eclState, this->simulator().vanguard().schedule());

        // let the object for threshold pressures initialize itself. this is done only at
        // this point, because determining the threshold pressures may require to access
        // the initial solution.
        thresholdPressures_.finishInit();

        // release the memory of the EQUIL grid since it's no longer needed after this point
        this->simulator().vanguard().releaseEquilGrid();

        updateCompositionChangeLimits_();

        aquiferModel_.initialSolutionApplied();

        const auto& initconfig = eclState.getInitConfig();
        if (initconfig.restartRequested())
            readEclRestartSolution_();
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0 everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context,
                unsigned spaceIdx,
                unsigned timeIdx) const
    {
        rate = 0.0;

        wellModel_.computeTotalRatesForDof(rate, context, spaceIdx, timeIdx);

        // convert the source term from the total mass rate of the
        // cell to the one per unit of volume as used by the model.
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
            rate[eqIdx] /= this->model().dofTotalVolume(globalDofIdx);

            Opm::Valgrind::CheckDefined(rate[eqIdx]);
            assert(Opm::isfinite(rate[eqIdx]));
        }

        aquiferModel_.addToSource(rate, context, spaceIdx, timeIdx);
    }

    /*!
     * \brief Returns the maximum value of the gas dissolution factor at the current time
     *        for a given degree of freedom.
     */
    Scalar maxGasDissolutionFactor(unsigned timeIdx, unsigned globalDofIdx) const
    {
        int pvtRegionIdx = pvtRegionIndex(globalDofIdx);
        if (!drsdtActive_() || maxDRs_[pvtRegionIdx] < 0.0)
            return std::numeric_limits<Scalar>::max()/2.0;

        // this is a bit hacky because it assumes that a time discretization with only
        // two time indices is used.
        if (timeIdx == 0)
            return lastRs_[globalDofIdx] + maxDRs_[pvtRegionIdx];
        else
            return lastRs_[globalDofIdx];
    }

    /*!
     * \brief Returns the maximum value of the oil vaporization factor at the current
     *        time for a given degree of freedom.
     */
    Scalar maxOilVaporizationFactor(unsigned timeIdx, unsigned globalDofIdx) const
    {
        int pvtRegionIdx = pvtRegionIndex(globalDofIdx);
        if (!drvdtActive_() || maxDRv_[pvtRegionIdx] < 0.0)
            return std::numeric_limits<Scalar>::max()/2.0;

        // this is a bit hacky because it assumes that a time discretization with only
        // two time indices is used.
        if (timeIdx == 0)
            return lastRv_[globalDofIdx] + maxDRv_[pvtRegionIdx];
        else
            return lastRv_[globalDofIdx];
    }

    /*!
     * \brief Returns an element's maximum oil phase saturation observed during the
     *        simulation.
     *
     * This is a bit of a hack from the conceptional point of view, but it is required to
     * match the results of the 'flow' and ECLIPSE 100 simulators.
     */
    Scalar maxOilSaturation(unsigned globalDofIdx) const
    {
        if (!vapparsActive())
            return 0.0;

        return maxOilSaturation_[globalDofIdx];
    }

    /*!
     * \brief Sets an element's maximum oil phase saturation observed during the
     *        simulation.
     *
     * This a hack on top of the maxOilSaturation() hack but it is currently required to
     * do restart externally. i.e. from the flow code.
     */
    void setMaxOilSaturation(unsigned globalDofIdx, Scalar value)
    {
        if (!vapparsActive())
            return;

        maxOilSaturation_[globalDofIdx] = value;
    }

    /*!
     * \brief Returns a reference to the ECL well manager used by the problem.
     *
     * This can be used for inspecting wells outside of the problem.
     */
    const EclWellModel& wellModel() const
    { return wellModel_; }

    EclWellModel& wellModel()
    { return wellModel_; }

    // temporary solution to facilitate output of initial state from flow
    const InitialFluidState& initialFluidState(unsigned globalDofIdx) const
    { return initialFluidStates_[globalDofIdx]; }

    const Opm::EclipseIO& eclIO() const
    { return eclWriter_->eclIO(); }

    bool vapparsActive() const
    {
        int epsiodeIdx = std::max(this->simulator().episodeIndex(), 0);
        const auto& oilVaporizationControl = this->simulator().vanguard().schedule().getOilVaporizationProperties(epsiodeIdx);
        return (oilVaporizationControl.getType() == Opm::OilVaporizationEnum::VAPPARS);
    }

    bool hasFreeBoundaryConditions() const
    { return hasFreeBoundaryConditions_; }

    /*!
     * \brief Propose the size of the next time step to the simulator.
     *
     * This method is only called if the Newton solver does converge, the simulator
     * automatically cuts the time step in half without consultating this method again.
     */
    Scalar nextTimeStepSize() const
    {
        // allow external code to do the timestepping
        if (this->nextTimeStepSize_ > 0.0)
            return this->nextTimeStepSize_;

        int episodeIdx = this->simulator().episodeIndex();

        // for the initial episode, we use a fixed time step size
        if (episodeIdx < 0)
            return initialTimeStepSize_;

        // ask the newton method for a suggestion. This suggestion will be based on how
        // well the previous time step converged. After that, apply the runtime time
        // stepping constraints.
        const auto& newtonMethod = this->model().newtonMethod();
        const auto& simulator = this->simulator();
        return limitNextTimeStepSize_(newtonMethod.suggestTimeStepSize(simulator.timeStepSize()));
    }

    /*!
     * \brief Called by Ewoms::Simulator in order to do a time integration on the model.
     */
    void timeIntegration()
    {
        Simulator& simulator = this->simulator();
        // if the time step size of the simulator is smaller than
        // the specified minimum size and we're not going to finish
        // the simulation or an episode, try with the minimum size.
        if (simulator.timeStepSize() < minTimeStepSize_ &&
            !simulator.episodeWillBeOver() &&
            !simulator.willBeFinished())
        {
            simulator.setTimeStepSize(minTimeStepSize_);
        }

        for (unsigned i = 0; i < maxFails_; ++i) {
            bool converged = this->model().update();
            if (converged)
                return;

            Scalar dt = simulator.timeStepSize();
            Scalar nextDt = dt / restartShrinkFactor_;
            if (nextDt < minTimeStepSize_)
                break; // give up: we can't make the time step smaller anymore!
            simulator.setTimeStepSize(nextDt);

            // update failed
            if (this->gridView().comm().rank() == 0)
                std::cout << "Newton solver did not converge with "
                          << "dt=" << dt << " seconds. Retrying with time step of "
                          << nextDt << " seconds\n" << std::flush;
        }

        throw std::runtime_error("Newton solver didn't converge after "
                                 +std::to_string(maxFails_)+" time-step divisions. dt="
                                 +std::to_string(double(simulator.timeStepSize())));
    }


private:
    bool drsdtActive_() const
    {
        int epsiodeIdx = std::max(this->simulator().episodeIndex(), 0);
        const auto& oilVaporizationControl = this->simulator().vanguard().schedule().getOilVaporizationProperties(epsiodeIdx);
        return (oilVaporizationControl.drsdtActive());
    }

    bool drvdtActive_() const
    {
        int epsiodeIdx = std::max(this->simulator().episodeIndex(), 0);
        const auto& oilVaporizationControl = this->simulator().vanguard().schedule().getOilVaporizationProperties(epsiodeIdx);
        return (oilVaporizationControl.drvdtActive());

    }

    Scalar cellCenterDepth(const Element& element) const
    {
        typedef typename Element::Geometry Geometry;
        static constexpr int zCoord = Element::dimension - 1;
        Scalar zz = 0.0;

        const Geometry geometry = element.geometry();
        const int corners = geometry.corners();
        for (int i=0; i < corners; ++i)
            zz += geometry.corner(i)[zCoord];

        return zz/Scalar(corners);
    }

    void updateElementDepths_()
    {
        const auto& vanguard = this->simulator().vanguard();
        const auto& gridView = vanguard.gridView();
        const auto& elemMapper = this->elementMapper();;

        int numElements = gridView.size(/*codim=*/0);
        elementCenterDepth_.resize(numElements);

        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& element = *elemIt;
            const unsigned int elemIdx = elemMapper.index(element);

            elementCenterDepth_[elemIdx] = cellCenterDepth(element);
        }
    }

    // update the parameters needed for DRSDT and DRVDT
    void updateCompositionChangeLimits_()
    {
        // update the "last Rs" values for all elements, including the ones in the ghost
        // and overlap regions
        int epsiodeIdx = std::max(this->simulator().episodeIndex(), 0);
        const auto& oilVaporizationControl = this->simulator().vanguard().schedule().getOilVaporizationProperties(epsiodeIdx);

        if (oilVaporizationControl.drsdtActive()) {
            ElementContext elemCtx(this->simulator());
            const auto& vanguard = this->simulator().vanguard();
            auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
            const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const Element& elem = *elemIt;

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& iq = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = iq.fluidState();

                typedef typename std::decay<decltype(fs)>::type FluidState;

                int pvtRegionIdx = pvtRegionIndex(compressedDofIdx);
                if (oilVaporizationControl.getOption(pvtRegionIdx) || fs.saturation(gasPhaseIdx) > freeGasMinSaturation_)
                    lastRs_[compressedDofIdx] =
                        Opm::BlackOil::template getRs_<FluidSystem,
                                                       FluidState,
                                                       Scalar>(fs, iq.pvtRegionIndex());
                else
                    lastRs_[compressedDofIdx] = std::numeric_limits<Scalar>::infinity();
            }
        }

        // update the "last Rv" values for all elements, including the ones in the ghost
        // and overlap regions
        if (drvdtActive_()) {
            ElementContext elemCtx(this->simulator());
            const auto& vanguard = this->simulator().vanguard();
            auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
            const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const Element& elem = *elemIt;

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& iq = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = iq.fluidState();

                typedef typename std::decay<decltype(fs)>::type FluidState;

                lastRv_[compressedDofIdx] =
                    Opm::BlackOil::template getRv_<FluidSystem,
                                                   FluidState,
                                                   Scalar>(fs, iq.pvtRegionIndex());
            }
        }
    }

    bool updateMaxOilSaturation_()
    {
        // we use VAPPARS
        if (vapparsActive()) {
            ElementContext elemCtx(this->simulator());
            const auto& vanguard = this->simulator().vanguard();
            auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
            const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const Element& elem = *elemIt;

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& iq = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = iq.fluidState();

                Scalar So = Opm::decay<Scalar>(fs.saturation(oilPhaseIdx));

                maxOilSaturation_[compressedDofIdx] = std::max(maxOilSaturation_[compressedDofIdx], So);
            }

            // we need to invalidate the intensive quantities cache here because the
            // derivatives of Rs and Rv will most likely have changed
            return true;
        }

        return false;
    }

    void readRockParameters_()
    {
        const auto& deck = this->simulator().vanguard().deck();
        const auto& eclState = this->simulator().vanguard().eclState();
        const auto& vanguard = this->simulator().vanguard();

        // the ROCK keyword has not been specified, so we don't need
        // to read rock parameters
        if (!deck.hasKeyword("ROCK"))
            return;

        const auto& rockKeyword = deck.getKeyword("ROCK");
        rockParams_.resize(rockKeyword.size());
        for (size_t rockRecordIdx = 0; rockRecordIdx < rockKeyword.size(); ++ rockRecordIdx) {
            const auto& rockRecord = rockKeyword.getRecord(rockRecordIdx);
            rockParams_[rockRecordIdx].referencePressure =
                rockRecord.getItem("PREF").getSIDouble(0);
            rockParams_[rockRecordIdx].compressibility =
                rockRecord.getItem("COMPRESSIBILITY").getSIDouble(0);
        }

        // check the kind of region which is supposed to be used by checking the ROCKOPTS
        // keyword. note that for some funny reason, the ROCK keyword uses PVTNUM by
        // default, *not* ROCKNUM!
        std::string propName = "PVTNUM";
        if (deck.hasKeyword("ROCKOPTS")) {
            const auto& rockoptsKeyword = deck.getKeyword("ROCKOPTS");
            std::string rockTableType =
                rockoptsKeyword.getRecord(0).getItem("TABLE_TYPE").getTrimmedString(0);
            if (rockTableType == "PVTNUM")
                propName = "PVTNUM";
            else if (rockTableType == "SATNUM")
                propName = "SATNUM";
            else if (rockTableType == "ROCKNUM")
                propName = "ROCKNUM";
            else {
                throw std::runtime_error("Unknown table type '"+rockTableType
                                         +" for the ROCKOPTS keyword given");
            }
        }

        // the deck does not specify the selected keyword, so everything uses the first
        // record of ROCK.
        if (!eclState.get3DProperties().hasDeckIntGridProperty(propName))
            return;

        const std::vector<int>& tablenumData =
            eclState.get3DProperties().getIntGridProperty(propName).getData();
        unsigned numElem = vanguard.gridView().size(0);
        rockTableIdx_.resize(numElem);
        for (size_t elemIdx = 0; elemIdx < numElem; ++ elemIdx) {
            unsigned cartElemIdx = vanguard.cartesianIndex(elemIdx);

            // reminder: Eclipse uses FORTRAN-style indices
            rockTableIdx_[elemIdx] = tablenumData[cartElemIdx] - 1;
        }
    }

    void readMaterialParameters_()
    {
        const auto& vanguard = this->simulator().vanguard();
        const auto& deck = vanguard.deck();
        const auto& eclState = vanguard.eclState();

        // the PVT and saturation region numbers
        updatePvtnum_();
        updateSatnum_();

        // the MISC region numbers (solvent model)
        updateMiscnum_();
        // the PLMIX region numbers (polymer model)
        updatePlmixnum_();

        ////////////////////////////////
        // porosity
        updatePorosity_();
        ////////////////////////////////

        ////////////////////////////////
        // fluid-matrix interactions (saturation functions; relperm/capillary pressure)
        size_t numDof = this->model().numGridDof();
        std::vector<int> compressedToCartesianElemIdx(numDof);
        for (unsigned elemIdx = 0; elemIdx < numDof; ++elemIdx)
            compressedToCartesianElemIdx[elemIdx] = vanguard.cartesianIndex(elemIdx);

        materialLawManager_ = std::make_shared<EclMaterialLawManager>();
        materialLawManager_->initFromDeck(deck, eclState, compressedToCartesianElemIdx);
        ////////////////////////////////
    }

    void readThermalParameters_()
    {
        if (!enableEnergy)
            return;

        const auto& vanguard = this->simulator().vanguard();
        const auto& deck = vanguard.deck();
        const auto& eclState = vanguard.eclState();

        // fluid-matrix interactions (saturation functions; relperm/capillary pressure)
        size_t numDof = this->model().numGridDof();
        std::vector<int> compressedToCartesianElemIdx(numDof);
        for (unsigned elemIdx = 0; elemIdx < numDof; ++elemIdx)
            compressedToCartesianElemIdx[elemIdx] = vanguard.cartesianIndex(elemIdx);

        thermalLawManager_ = std::make_shared<EclThermalLawManager>();
        thermalLawManager_->initFromDeck(deck, eclState, compressedToCartesianElemIdx);
    }

    void updatePorosity_()
    {
        const auto& vanguard = this->simulator().vanguard();
        const auto& eclState = vanguard.eclState();
        const auto& eclGrid = eclState.getInputGrid();
        const auto& props = eclState.get3DProperties();

        size_t numDof = this->model().numGridDof();

        porosity_.resize(numDof);

        const std::vector<double>& porvData =
            props.getDoubleGridProperty("PORV").getData();
        const std::vector<int>& actnumData =
            props.getIntGridProperty("ACTNUM").getData();

        int nx = eclGrid.getNX();
        int ny = eclGrid.getNY();
        for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx) {
            unsigned cartElemIdx = vanguard.cartesianIndex(dofIdx);
            Scalar poreVolume = porvData[cartElemIdx];

            // sum up the pore volume of the active cell and all inactive ones above it
            // which were disabled due to their pore volume being too small. If energy is
            // conserved, cells are not disabled due to a too small pore volume because
            // such cells still store and conduct energy.
            if (!enableEnergy && eclGrid.getMinpvMode() == Opm::MinpvMode::ModeEnum::OpmFIL) {
                const std::vector<Scalar>& minPvVector = eclGrid.getMinpvVector();
                for (int aboveElemCartIdx = static_cast<int>(cartElemIdx) - nx*ny;
                     aboveElemCartIdx >= 0;
                     aboveElemCartIdx -= nx*ny)
                {
                    if (porvData[aboveElemCartIdx] >= minPvVector[aboveElemCartIdx])
                        // the cartesian element above exhibits a pore volume which larger or
                        // equal to the minimum one
                        break;

                    Scalar aboveElemVolume = eclGrid.getCellVolume(aboveElemCartIdx);
                    if (actnumData[aboveElemCartIdx] == 0 && aboveElemVolume > 1e-3)
                        // stop at explicitly disabled elements, but only if their volume is
                        // greater than 10^-3 m^3
                        break;

                    poreVolume += porvData[aboveElemCartIdx];
                }
            }

            // we define the porosity as the accumulated pore volume divided by the
            // geometric volume of the element. Note that -- in pathetic cases -- it can
            // be larger than 1.0!
            Scalar dofVolume = this->simulator().model().dofTotalVolume(dofIdx);
            assert(dofVolume > 0.0);
            porosity_[dofIdx] = poreVolume/dofVolume;
        }
    }

    void initFluidSystem_()
    {
        const auto& deck = this->simulator().vanguard().deck();
        const auto& eclState = this->simulator().vanguard().eclState();

        FluidSystem::initFromDeck(deck, eclState);
   }

    void readInitialCondition_()
    {
        const auto& vanguard = this->simulator().vanguard();

        const auto& deck = vanguard.deck();
        if (!deck.hasKeyword("EQUIL"))
            readExplicitInitialCondition_();
        else
            readEquilInitialCondition_();

        readBlackoilExtentionsInitialConditions_();

    }

    void readEquilInitialCondition_()
    {
        // initial condition corresponds to hydrostatic conditions.
        typedef Ewoms::EclEquilInitializer<TypeTag> EquilInitializer;
        EquilInitializer equilInitializer(this->simulator(), *materialLawManager_);

        size_t numElems = this->model().numGridDof();
        initialFluidStates_.resize(numElems);
        for (size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemFluidState = initialFluidStates_[elemIdx];
            elemFluidState.assign(equilInitializer.initialFluidState(elemIdx));
        }
    }

    void readEclRestartSolution_()
    {
        // Set the start time of the simulation
        const auto& schedule = this->simulator().vanguard().schedule();
        const auto& eclState = this->simulator().vanguard().eclState();
        const auto& timeMap = schedule.getTimeMap();
        const auto& initconfig = eclState.getInitConfig();
        int episodeIdx = initconfig.getRestartStep() - 1;

        this->simulator().setStartTime(timeMap.getStartTime(/*timeStepIdx=*/0));
        this->simulator().setTime(timeMap.getTimePassedUntil(episodeIdx));
        this->simulator().setEpisodeIndex(episodeIdx);
        this->simulator().setEpisodeLength(timeMap.getTimeStepLength(episodeIdx));
        this->simulator().setTimeStepSize(eclWriter_->restartTimeStepSize());

        eclWriter_->beginRestart();

        size_t numElems = this->model().numGridDof();
        initialFluidStates_.resize(numElems);
        if (enableSolvent)
            solventSaturation_.resize(numElems, 0.0);

        if (enablePolymer)
            polymerConcentration_.resize(numElems, 0.0);

        if (enablePolymerMolarWeight) {
            const std::string msg {"Support of the RESTART for polymer molecular weight "
                                   "is not implemented yet. The polymer weight value will be "
                                   "zero when RESTART begins"};
            Opm::OpmLog::warning("NO_POLYMW_RESTART", msg);
            polymerMoleWeight_.resize(numElems, 0.0);
        }

        // this is a hack to preserve the initial fluid states
        auto tmpInitialFs = initialFluidStates_;

        for (size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemFluidState = initialFluidStates_[elemIdx];
            elemFluidState.setPvtRegionIndex(pvtRegionIndex(elemIdx));
            eclWriter_->eclOutputModule().initHysteresisParams(this->simulator(), elemIdx);
            eclWriter_->eclOutputModule().assignToFluidState(elemFluidState, elemIdx);

            processRestartSaturations_(elemFluidState);

            lastRs_[elemIdx] = elemFluidState.Rs();
            lastRv_[elemIdx] = elemFluidState.Rv();
            if (enableSolvent)
                 solventSaturation_[elemIdx] = eclWriter_->eclOutputModule().getSolventSaturation(elemIdx);
            if (enablePolymer)
                 polymerConcentration_[elemIdx] = eclWriter_->eclOutputModule().getPolymerConcentration(elemIdx);
            // if we need to restart for polymer molecular weight simulation, we need to add related here
        }

        const int epsiodeIdx = this->simulator().episodeIndex();
        const auto& oilVaporizationControl = this->simulator().vanguard().schedule().getOilVaporizationProperties(epsiodeIdx);
        if (drsdtActive_())
            // DRSDT is enabled
            for (size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRs_.size(); ++pvtRegionIdx)
                maxDRs_[pvtRegionIdx] = oilVaporizationControl.getMaxDRSDT(pvtRegionIdx)*this->simulator().timeStepSize();

        if (drvdtActive_())
            // DRVDT is enabled
            for (size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRv_.size(); ++pvtRegionIdx)
                maxDRv_[pvtRegionIdx] = oilVaporizationControl.getMaxDRVDT(pvtRegionIdx)*this->simulator().timeStepSize();

        if (tracerModel().numTracers() > 0)
            std::cout << "Warning: Restart is not implemented for the tracer model, it will initialize with initial tracer concentration" << std::endl;

        // assign the restart solution to the current solution. note that we still need
        // to compute real initial solution after this because the initial fluid states
        // need to be correct for stuff like boundary conditions.
        auto& sol = this->model().solution(/*timeIdx=*/0);
        const auto& gridView = this->gridView();
        ElementContext elemCtx(this->simulator());
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updatePrimaryStencil(elem);
            int elemIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            initial(sol[elemIdx], elemCtx, /*spaceIdx=*/0, /*timeIdx=*/0);
        }

        // make sure that the ghost and overlap entities exhibit the correct
        // solution. alternatively, this could be done in the loop above by also
        // considering non-interior elements. Since the initial() method might not work
        // 100% correctly for such elements, let's play safe and explicitly synchronize
        // using message passing.
        this->model().syncOverlap();

        // this is a hack to preserve the initial fluid states
        initialFluidStates_ = tmpInitialFs;
        // make sure that the stuff which needs to be done at the beginning of an episode
        // is run.
        this->beginEpisode(/*isOnRestart=*/true);

        eclWriter_->endRestart();
    }

    void processRestartSaturations_(InitialFluidState& elemFluidState)
    {
        // each phase needs to be above certain value to be claimed to be existing
        // this is used to recover some RESTART running with the defaulted single-precision format
        const Scalar smallSaturationTolerance = 1.e-6;
        Scalar sumSaturation = 0.0;
        for (size_t phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                if (elemFluidState.saturation(phaseIdx) < smallSaturationTolerance)
                    elemFluidState.setSaturation(phaseIdx, 0.0);

                sumSaturation += elemFluidState.saturation(phaseIdx);
            }
        }

        assert(sumSaturation > 0.0);

        for (size_t phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                const Scalar saturation = elemFluidState.saturation(phaseIdx) / sumSaturation;
                elemFluidState.setSaturation(phaseIdx, saturation);
            }
        }
    }

    void readExplicitInitialCondition_()
    {
        const auto& vanguard = this->simulator().vanguard();
        const auto& eclState = vanguard.eclState();
        const auto& eclProps = eclState.get3DProperties();

        // make sure all required quantities are enables
        if (FluidSystem::phaseIsActive(waterPhaseIdx) && !eclProps.hasDeckDoubleGridProperty("SWAT"))
            throw std::runtime_error("The ECL input file requires the presence of the SWAT keyword if "
                                     "the water phase is active");
        if (FluidSystem::phaseIsActive(gasPhaseIdx) && !eclProps.hasDeckDoubleGridProperty("SGAS"))
            throw std::runtime_error("The ECL input file requires the presence of the SGAS keyword if "
                                     "the gas phase is active");

        if (!eclProps.hasDeckDoubleGridProperty("PRESSURE"))
             throw std::runtime_error("The ECL input file requires the presence of the PRESSURE "
                                      "keyword if the model is initialized explicitly");
        if (FluidSystem::enableDissolvedGas() && !eclProps.hasDeckDoubleGridProperty("RS"))
            throw std::runtime_error("The ECL input file requires the RS keyword to be present if"
                                     " dissolved gas is enabled");
        if (FluidSystem::enableVaporizedOil() && !eclProps.hasDeckDoubleGridProperty("RV"))
            throw std::runtime_error("The ECL input file requires the RV keyword to be present if"
                                     " vaporized oil is enabled");

        size_t numDof = this->model().numGridDof();

        initialFluidStates_.resize(numDof);

        const auto& cartSize = this->simulator().vanguard().cartesianDimensions();
        size_t numCartesianCells = cartSize[0] * cartSize[1] * cartSize[2];

        std::vector<double> waterSaturationData;
        if (FluidSystem::phaseIsActive(waterPhaseIdx))
            waterSaturationData = eclProps.getDoubleGridProperty("SWAT").getData();
        else
            waterSaturationData.resize(numCartesianCells, 0.0);

        std::vector<double> gasSaturationData;
        if (FluidSystem::phaseIsActive(gasPhaseIdx))
            gasSaturationData = eclProps.getDoubleGridProperty("SGAS").getData();
        else
            gasSaturationData.resize(numCartesianCells, 0.0);

        const std::vector<double>& pressureData =
            eclProps.getDoubleGridProperty("PRESSURE").getData();
        std::vector<double> rsData;
        if (FluidSystem::enableDissolvedGas())
            rsData = eclProps.getDoubleGridProperty("RS").getData();
        std::vector<double> rvData;
        if (FluidSystem::enableVaporizedOil())
            rvData = eclProps.getDoubleGridProperty("RV").getData();
        // initial reservoir temperature
        const std::vector<double>& tempiData =
            eclState.get3DProperties().getDoubleGridProperty("TEMPI").getData();

        // make sure that the size of the data arrays is correct
#ifndef NDEBUG
        assert(waterSaturationData.size() == numCartesianCells);
        assert(gasSaturationData.size() == numCartesianCells);
        assert(pressureData.size() == numCartesianCells);
        if (FluidSystem::enableDissolvedGas())
            assert(rsData.size() == numCartesianCells);
        if (FluidSystem::enableVaporizedOil())
            assert(rvData.size() == numCartesianCells);
#endif

        // calculate the initial fluid states
        for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofFluidState = initialFluidStates_[dofIdx];

            dofFluidState.setPvtRegionIndex(pvtRegionIndex(dofIdx));
            size_t cartesianDofIdx = vanguard.cartesianIndex(dofIdx);
            assert(0 <= cartesianDofIdx);
            assert(cartesianDofIdx <= numCartesianCells);

            //////
            // set temperature
            //////
            Scalar temperature = tempiData[cartesianDofIdx];
            if (!std::isfinite(temperature) || temperature <= 0)
                temperature = FluidSystem::surfaceTemperature;
            dofFluidState.setTemperature(temperature);

            //////
            // set saturations
            //////
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))
                dofFluidState.setSaturation(FluidSystem::waterPhaseIdx,
                                            waterSaturationData[cartesianDofIdx]);
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
                dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                            gasSaturationData[cartesianDofIdx]);
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx))
                dofFluidState.setSaturation(FluidSystem::oilPhaseIdx,
                                            1.0
                                            - waterSaturationData[cartesianDofIdx]
                                            - gasSaturationData[cartesianDofIdx]);

            //////
            // set phase pressures
            //////
            Scalar oilPressure = pressureData[cartesianDofIdx];

            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            Dune::FieldVector<Scalar, numPhases> pc(0.0);
            const auto& matParams = materialLawParams(dofIdx);
            MaterialLaw::capillaryPressures(pc, matParams, dofFluidState);
            Opm::Valgrind::CheckDefined(oilPressure);
            Opm::Valgrind::CheckDefined(pc);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                dofFluidState.setPressure(phaseIdx, oilPressure + (pc[phaseIdx] - pc[oilPhaseIdx]));
            }

            if (FluidSystem::enableDissolvedGas())
                dofFluidState.setRs(rsData[cartesianDofIdx]);
            else if (Indices::gasEnabled && Indices::oilEnabled)
                dofFluidState.setRs(0.0);

            if (FluidSystem::enableVaporizedOil())
                dofFluidState.setRv(rvData[cartesianDofIdx]);
            else if (Indices::gasEnabled && Indices::oilEnabled)
                dofFluidState.setRv(0.0);

            //////
            // set invB_
            //////
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                const auto& b = FluidSystem::inverseFormationVolumeFactor(dofFluidState, phaseIdx, pvtRegionIndex(dofIdx));
                dofFluidState.setInvB(phaseIdx, b);

                const auto& rho = FluidSystem::density(dofFluidState, phaseIdx, pvtRegionIndex(dofIdx));
                dofFluidState.setDensity(phaseIdx, rho);

            }
        }
    }

    void readBlackoilExtentionsInitialConditions_()
    {
        const auto& vanguard = this->simulator().vanguard();
        const auto& eclState = vanguard.eclState();
        size_t numDof = this->model().numGridDof();

        if (enableSolvent) {
            const std::vector<double>& solventSaturationData = eclState.get3DProperties().getDoubleGridProperty("SSOL").getData();
            solventSaturation_.resize(numDof, 0.0);
            for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
                size_t cartesianDofIdx = vanguard.cartesianIndex(dofIdx);
                assert(0 <= cartesianDofIdx);
                assert(cartesianDofIdx <= solventSaturationData.size());
                solventSaturation_[dofIdx] = solventSaturationData[cartesianDofIdx];
            }
        }

        if (enablePolymer) {
            const std::vector<double>& polyConcentrationData = eclState.get3DProperties().getDoubleGridProperty("SPOLY").getData();
            polymerConcentration_.resize(numDof, 0.0);
            for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
                size_t cartesianDofIdx = vanguard.cartesianIndex(dofIdx);
                assert(0 <= cartesianDofIdx);
                assert(cartesianDofIdx <= polyConcentrationData.size());
                polymerConcentration_[dofIdx] = polyConcentrationData[cartesianDofIdx];
            }
        }

        if (enablePolymerMolarWeight) {
            const std::vector<double>& polyMoleWeightData = eclState.get3DProperties().getDoubleGridProperty("SPOLYMW").getData();
            polymerMoleWeight_.resize(numDof, 0.0);
            for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
                const size_t cartesianDofIdx = vanguard.cartesianIndex(dofIdx);
                assert(0 <= cartesianDofIdx);
                assert(cartesianDofIdx <= polyMoleWeightData.size());
                polymerMoleWeight_[dofIdx] = polyMoleWeightData[cartesianDofIdx];
            }
        }
    }

    // update the hysteresis parameters of the material laws for the whole grid
    bool updateHysteresis_()
    {
        if (!materialLawManager_->enableHysteresis())
            return false;

        // we need to update the hysteresis data for _all_ elements (i.e., not just the
        // interior ones) to avoid desynchronization of the processes in the parallel case!
        ElementContext elemCtx(this->simulator());
        const auto& vanguard = this->simulator().vanguard();
        auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
        const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            materialLawManager_->updateHysteresis(intQuants.fluidState(), compressedDofIdx);
        }
        return true;
    }

    void updateMaxPolymerAdsorption_()
    {
        // we need to update the max polymer adsoption data for all elements
        ElementContext elemCtx(this->simulator());
        const auto& vanguard = this->simulator().vanguard();
        auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
        const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);

            maxPolymerAdsorption_[compressedDofIdx] = std::max(maxPolymerAdsorption_[compressedDofIdx] , Opm::scalarValue(intQuants.polymerAdsorption()));
        }
    }

    void updatePvtnum_()
    {
        const auto& eclState = this->simulator().vanguard().eclState();
        const auto& eclProps = eclState.get3DProperties();

        if (!eclProps.hasDeckIntGridProperty("PVTNUM"))
            return;

        const auto& pvtnumData = eclProps.getIntGridProperty("PVTNUM").getData();
        const auto& vanguard = this->simulator().vanguard();

        unsigned numElems = vanguard.gridView().size(/*codim=*/0);
        pvtnum_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            unsigned cartElemIdx = vanguard.cartesianIndex(elemIdx);
            pvtnum_[elemIdx] = pvtnumData[cartElemIdx] - 1;
        }
    }

    void updateSatnum_()
    {
        const auto& eclState = this->simulator().vanguard().eclState();
        const auto& eclProps = eclState.get3DProperties();

        if (!eclProps.hasDeckIntGridProperty("SATNUM"))
            return;

        const auto& satnumData = eclProps.getIntGridProperty("SATNUM").getData();
        const auto& vanguard = this->simulator().vanguard();

        unsigned numElems = vanguard.gridView().size(/*codim=*/0);
        satnum_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            unsigned cartElemIdx = vanguard.cartesianIndex(elemIdx);
            satnum_[elemIdx] = satnumData[cartElemIdx] - 1;
        }
    }

    void updateMiscnum_()
    {
        const auto& eclState = this->simulator().vanguard().eclState();
        const auto& eclProps = eclState.get3DProperties();

        if (!eclProps.hasDeckIntGridProperty("MISCNUM"))
            return;

        const auto& miscnumData = eclProps.getIntGridProperty("MISCNUM").getData();
        const auto& vanguard = this->simulator().vanguard();

        unsigned numElems = vanguard.gridView().size(/*codim=*/0);
        miscnum_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            unsigned cartElemIdx = vanguard.cartesianIndex(elemIdx);
            miscnum_[elemIdx] = miscnumData[cartElemIdx] - 1;
        }
    }

    void updatePlmixnum_()
    {
        const auto& eclState = this->simulator().vanguard().eclState();
        const auto& eclProps = eclState.get3DProperties();

        if (!eclProps.hasDeckIntGridProperty("PLMIXNUM"))
            return;

        const auto& plmixnumData = eclProps.getIntGridProperty("PLMIXNUM").getData();
        const auto& vanguard = this->simulator().vanguard();

        unsigned numElems = vanguard.gridView().size(/*codim=*/0);
        plmixnum_.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            unsigned cartElemIdx = vanguard.cartesianIndex(elemIdx);
            plmixnum_[elemIdx] = plmixnumData[cartElemIdx] - 1;
        }
    }

    struct PffDofData_
    {
        Opm::ConditionalStorage<enableEnergy, Scalar> thermalHalfTrans;
        Scalar transmissibility;
    };

    // update the prefetch friendly data object
    void updatePffDofData_()
    {
        const auto& distFn =
            [this](PffDofData_& dofData,
                   const Stencil& stencil,
                   unsigned localDofIdx)
            -> void
        {
            const auto& elementMapper = this->model().elementMapper();

            unsigned globalElemIdx = elementMapper.index(stencil.entity(localDofIdx));
            if (localDofIdx != 0) {
                unsigned globalCenterElemIdx = elementMapper.index(stencil.entity(/*dofIdx=*/0));
                dofData.transmissibility = transmissibilities_.transmissibility(globalCenterElemIdx, globalElemIdx);

                if (enableEnergy)
                    *dofData.thermalHalfTrans = transmissibilities_.thermalHalfTrans(globalCenterElemIdx, globalElemIdx);
            }
        };

        pffDofData_.update(distFn);
    }

    void readBoundaryConditions_()
    {
        hasFreeBoundaryConditions_ = false;
        readBoundaryConditionKeyword_("FREEBCX", freebcX_);
        readBoundaryConditionKeyword_("FREEBCX-", freebcXMinus_);
        readBoundaryConditionKeyword_("FREEBCY", freebcY_);
        readBoundaryConditionKeyword_("FREEBCY-", freebcYMinus_);
        readBoundaryConditionKeyword_("FREEBCZ", freebcZ_);
        readBoundaryConditionKeyword_("FREEBCZ-", freebcZMinus_);

        const auto& vanguard = this->simulator().vanguard();

        if (vanguard.deck().hasKeyword("BCRATE")) {
            hasFreeBoundaryConditions_ = true;
            size_t numCartDof = vanguard.cartesianSize();
            unsigned numElems = vanguard.gridView().size(/*codim=*/0);
            std::vector<int> cartesianToCompressedElemIdx(numCartDof);

            for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx)
                cartesianToCompressedElemIdx[vanguard.cartesianIndex(elemIdx)] = elemIdx;

            massratebcXMinus_.resize(numElems, 0.0);
            massratebcX_.resize(numElems, 0.0);
            massratebcYMinus_.resize(numElems, 0.0);
            massratebcY_.resize(numElems, 0.0);
            massratebcZMinus_.resize(numElems, 0.0);
            massratebcZ_.resize(numElems, 0.0);

            const auto& ratebcs = vanguard.deck().getKeywordList("BCRATE");
            for (size_t listIdx = 0; listIdx < ratebcs.size(); ++listIdx) {
                const auto& ratebc = *ratebcs[listIdx];

                for (size_t record = 0; record < ratebc.size(); ++record) {

                    std::string compName = ratebc.getRecord(record).getItem("COMPONENT").getTrimmedString(0);
                    int compIdx = -999;

                    if (compName == "OIL")
                        compIdx = oilCompIdx;
                    else if (compName == "GAS")
                        compIdx = gasCompIdx;
                    else if (compName == "WATER")
                        compIdx = waterCompIdx;
                    else if (compName == "SOLVENT")
                    {
                        if (!enableSolvent)
                            throw std::logic_error("solvent is disabled and you're trying to add solvent to BCRATE");

                        compIdx = Indices::solventSaturationIdx;
                    }
                    else if (compName == "POLYMER")
                    {
                        if (!enablePolymer)
                            throw std::logic_error("polymer is disabled and you're trying to add polymer to BCRATE");

                        compIdx = Indices::polymerConcentrationIdx;
                    }
                    else
                        throw std::logic_error("invalid component name for BCRATE");
                    assert(compIdx >= 0);

                    std::string direction = ratebc.getRecord(record).getItem("DIRECTION").getTrimmedString(0);
                    std::vector<RateVector>* data = 0;
                    if (direction == "X-")
                        data = &massratebcXMinus_;
                    else if (direction == "X")
                        data = &massratebcX_;
                    else if (direction == "Y-")
                        data = &massratebcYMinus_;
                    else if (direction == "Y")
                        data = &massratebcY_;
                    else if (direction == "Z-")
                        data = &massratebcZMinus_;
                    else if (direction == "Z")
                        data = &massratebcZ_;
                    else
                        throw std::logic_error("invalid direction for BCRATE");

                    int i1 = ratebc.getRecord(record).getItem("I1").template get< int >(0) - 1;
                    int i2 = ratebc.getRecord(record).getItem("I2").template get< int >(0) - 1;
                    int j1 = ratebc.getRecord(record).getItem("J1").template get< int >(0) - 1;
                    int j2 = ratebc.getRecord(record).getItem("J2").template get< int >(0) - 1;
                    int k1 = ratebc.getRecord(record).getItem("K1").template get< int >(0) - 1;
                    int k2 = ratebc.getRecord(record).getItem("K2").template get< int >(0) - 1;
                    const Evaluation rate = ratebc.getRecord(record).getItem("RATE").getSIDouble(0);
                    for (int i = i1; i <= i2; ++i) {
                        for (int j = j1; j <= j2; ++j) {
                            for (int k = k1; k <= k2; ++k) {
                                std::array<int, 3> tmp = {i,j,k};
                                size_t elemIdx = cartesianToCompressedElemIdx[vanguard.cartesianIndex(tmp)];
                                (*data)[elemIdx][compIdx] = rate;
                            }
                        }
                    }
                }
            }
        }
    }

    void readBoundaryConditionKeyword_(const std::string& name, std::vector<bool>& compressedData)
    {
        const auto& eclProps = this->simulator().vanguard().eclState().get3DProperties();
        const auto& vanguard = this->simulator().vanguard();

        unsigned numElems = vanguard.gridView().size(/*codim=*/0);
        compressedData.resize(numElems, false);

        if (eclProps.hasDeckDoubleGridProperty(name)) {
            const std::vector<double>& data = eclProps.getDoubleGridProperty(name).getData();
            for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
                unsigned cartElemIdx = vanguard.cartesianIndex(elemIdx);
                compressedData[elemIdx] = (data[cartElemIdx] > 0);
            }
            hasFreeBoundaryConditions_ = true;
        }
    }

    // this method applies the runtime constraints specified via the deck and/or command
    // line parameters for the size of the next time step size.
    Scalar limitNextTimeStepSize_(Scalar dtNext) const
    {
        const auto& events = this->simulator().vanguard().schedule().getEvents();
        int episodeIdx = this->simulator().episodeIndex();

        // first thing in the morning, limit the time step size to the maximum size
        dtNext = std::min(dtNext, maxTimeStepSize_);

        // if TUNING is enabled, limit the time step size after a tuning event to TSINIT
        if (enableTuning_ && events.hasEvent(Opm::ScheduleEvents::TUNING_CHANGE, episodeIdx)) {
            const auto& tuning = this->simulator().vanguard().schedule().getTuning();
            return std::min(dtNext, tuning.getTSINIT(episodeIdx));
        }

        // use at least slightly more than half of the maximum time step size by default
        if (dtNext < maxTimeStepSize_ && maxTimeStepSize_ < dtNext*2)
            dtNext = 1.01 * maxTimeStepSize_/2.0;

        Scalar remainingEpisodeTime =
            this->simulator().episodeStartTime()
            + this->simulator().episodeLength()
            - this->simulator().time()
            - this->simulator().timeStepSize();

        // limit the time step size to the remaining time of the episode
        dtNext = std::min(dtNext, remainingEpisodeTime);

        // if we would have a small amount of time left over in the current episode, make
        // two equal time steps instead of a big and a small one
        if (dtNext < remainingEpisodeTime*(1.0 - 1e-5) && dtNext*1.25 > remainingEpisodeTime)
            // note: limiting to the maximum time step size here is probably not strictly
            // necessary, but it should not hurt and is more fool-proof
            dtNext = std::min(maxTimeStepSize_, remainingEpisodeTime/2.0);

        // if a well event occured, respect the limit for the maximum time step after
        // that, too
        int reportStepIdx = std::max(episodeIdx, 0);
        bool wellEventOccured =
            events.hasEvent(Opm::ScheduleEvents::NEW_WELL, reportStepIdx)
            || events.hasEvent(Opm::ScheduleEvents::PRODUCTION_UPDATE, reportStepIdx)
            || events.hasEvent(Opm::ScheduleEvents::INJECTION_UPDATE, reportStepIdx)
            || events.hasEvent(Opm::ScheduleEvents::WELL_STATUS_CHANGE, reportStepIdx);
        if (episodeIdx >= 0 && wellEventOccured)
            dtNext = std::min(dtNext, maxTimeStepAfterWellEvent_);

        return dtNext;
    }

    static std::string briefDescription_;

    std::vector<Scalar> porosity_;
    std::vector<Scalar> elementCenterDepth_;
    EclTransmissibility<TypeTag> transmissibilities_;

    std::shared_ptr<EclMaterialLawManager> materialLawManager_;
    std::shared_ptr<EclThermalLawManager> thermalLawManager_;

    EclThresholdPressure<TypeTag> thresholdPressures_;

    std::vector<int> pvtnum_;
    std::vector<unsigned short> satnum_;
    std::vector<unsigned short> miscnum_;
    std::vector<unsigned short> plmixnum_;

    std::vector<unsigned short> rockTableIdx_;
    std::vector<RockParams> rockParams_;

    std::vector<Scalar> maxPolymerAdsorption_;

    std::vector<InitialFluidState> initialFluidStates_;
    std::vector<Scalar> initialTemperature_;

    std::vector<Scalar> polymerConcentration_;
    // polymer molecular weight
    std::vector<Scalar> polymerMoleWeight_;
    std::vector<Scalar> solventSaturation_;

    std::vector<bool> dRsDtOnlyFreeGas_; // apply the DRSDT rate limit only to cells that exhibit free gas
    std::vector<Scalar> lastRs_;
    std::vector<Scalar> maxDRs_;
    std::vector<Scalar> lastRv_;
    std::vector<Scalar> maxDRv_;
    constexpr static Scalar freeGasMinSaturation_ = 1e-7;
    std::vector<Scalar> maxOilSaturation_;

    EclWellModel wellModel_;
    EclAquiferModel aquiferModel_;
    std::unique_ptr<EclWriterType> eclWriter_;

    PffGridVector<GridView, Stencil, PffDofData_, DofMapper> pffDofData_;
    TracerModel tracerModel_;

    bool hasFreeBoundaryConditions_;
    std::vector<bool> freebcX_;
    std::vector<bool> freebcXMinus_;
    std::vector<bool> freebcY_;
    std::vector<bool> freebcYMinus_;
    std::vector<bool> freebcZ_;
    std::vector<bool> freebcZMinus_;

    std::vector<RateVector> massratebcX_;
    std::vector<RateVector> massratebcXMinus_;
    std::vector<RateVector> massratebcY_;
    std::vector<RateVector> massratebcYMinus_;
    std::vector<RateVector> massratebcZ_;
    std::vector<RateVector> massratebcZMinus_;

    // time stepping parameters
    bool enableTuning_;
    Scalar initialTimeStepSize_;
    Scalar maxTimeStepAfterWellEvent_;
    Scalar maxTimeStepSize_;
    Scalar restartShrinkFactor_;
    unsigned maxFails_;
    Scalar minTimeStepSize_;
};

template <class TypeTag>
std::string EclProblem<TypeTag>::briefDescription_;

} // namespace Ewoms

#endif
