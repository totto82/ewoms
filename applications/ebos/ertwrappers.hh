/*
  Copyright (C) 2013-2014 by Andreas Lauser
  Copyright (c) 2013 by SINTEF ICT, Applied Mathematics.
  Copyright (c) 2013 by Uni Research AS

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
*/
/*!
 * \file
 *
 * \brief This file implements several wrapper classes around the
 *        opaque ERT data types.
 *
 * These classes are shamelessly ripped-off from opm-core and are
 * required to make writing ECL files exception safe...
 */
#ifndef EWOMS_ERT_WRAPPERS_HH
#define EWOMS_ERT_WRAPPERS_HH

#if HAVE_ERT

#include <ert/ecl/fortio.h>
#include <ert/ecl/ecl_endian_flip.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_sum.h>
#include <ert/ecl/ecl_util.h>
#include <ert/ecl/ecl_init_file.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_rst_file.h>

#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/material/Valgrind.hpp>

#include "eclwellmanager.hh"

namespace Ewoms {

// forward definition of the EclGridManager class. We need this for
// specialization...
template <class TypeTag>
class EclGridManager;

/// \cond 0

// required to make the compiler happy if the grid manager is not EclGridManager...

template <class GridManager>
std::string getErtCaseName__(const GridManager &gridManager)
{ OPM_THROW(std::logic_error, "You need to chose the EclGridManager to write ECL files"); }

template <class TypeTag>
std::string getErtCaseName__(const EclGridManager<TypeTag> &gridManager)
{ return gridManager.caseName(); }

template <class GridManager>
const Opm::EclipseGridConstPtr getEclGrid__(const GridManager &gridManager)
{ OPM_THROW(std::logic_error, "You need to chose the EclGridManager to write ECL files"); }

template <class TypeTag>
const Opm::EclipseGridConstPtr getEclGrid__(const EclGridManager<TypeTag> &gridManager)
{ return gridManager.eclGrid(); }

/// \endcond

class ErtBaseKeyword
{
public:
    virtual ~ErtBaseKeyword() {}
};

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This is a smart pointer class for ERT's ecl_kw_type
 *        structure.
 */
template <typename T>
class ErtKeyword : public ErtBaseKeyword
{
public:
#if HAVE_ERT
    typedef ecl_kw_type ErtHandleType;
#else
    typedef int ErtHandleType;
#endif

    // don't allow copies for objects of this class
    ErtKeyword(const ErtKeyword &) = delete;

    // Default constructor
    ErtKeyword()
        : ertHandle_(0)
    {}

    /// Initialization from single-precision array.
    ErtKeyword(const std::string& name,
               const std::vector<float>& data)
        : ertHandle_(0)
    { set(name, data); }

    /// Initialization from double-precision array.
    ErtKeyword(const std::string& name,
               const std::vector<double>& data)
        : ertHandle_(0)
    { set(name, data); }

    /// Initialization from double-precision array.
    ErtKeyword(const std::string& name,
               const std::vector<int>& data)
        : ertHandle_(0)
    { set(name, data); }

    ~ErtKeyword()
    {
#if HAVE_ERT
        if (ertHandle_)
            ecl_kw_free(ertHandle_);
#endif
    }

    template <class DataElementType>
    void set(const std::string name, const std::vector<DataElementType>& data)
    {
#if HAVE_ERT
        if(ertHandle_) {
            ecl_kw_free(ertHandle_);
        }

        name_ = name;
        ertHandle_ = ecl_kw_alloc(name.c_str(),
                                  data.size(),
                                  ertType_());

        // number of elements to take
        const int numEntries = data.size();

        // fill it with values
        T* target = static_cast<T*>(ecl_kw_get_ptr(ertHandle()));
        for (int i = 0; i < numEntries; ++i) {
            target[i] = static_cast<T>(data[i]);
        }

        Valgrind::CheckDefined(target, numEntries);
#endif
    }

    const std::string &name() const
    { return name_; }

    ErtHandleType *ertHandle() const
    { return ertHandle_; }

private:
#if HAVE_ERT
    static ecl_type_enum ertType_()
    {
        if (std::is_same<T, float>::value)
        { return ECL_FLOAT_TYPE; }
        if (std::is_same<T, double>::value)
        { return ECL_DOUBLE_TYPE; }
        if (std::is_same<T, int>::value)
        { return ECL_INT_TYPE; }

        OPM_THROW(std::logic_error,
                  "Unhandled type for data elements in ErtKeyword");
    }
#endif

    std::string name_;
    ErtHandleType *ertHandle_;
};

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This is a smart pointer class for ERT's ecl_grid_type
 *        structure.
 */
class ErtGrid
{

public:
#if HAVE_ERT
    typedef ecl_grid_type ErtHandleType;
#else
    typedef int ErtHandleType;
#endif

    ErtGrid(const ErtGrid& ) = delete;

    /*!
     * \brief Create an ERT grid based an Opm::EclipseGrid.
     */
    template <class DeckUnits>
    ErtGrid(Opm::EclipseGridConstPtr eclGrid, const DeckUnits& deckUnits)
    {
#if HAVE_ERT
        std::vector<double> mapaxesData;
        std::vector<double> coordData;
        std::vector<double> zcornData;
        std::vector<int> actnumData;

        eclGrid->exportMAPAXES(mapaxesData);
        eclGrid->exportCOORD(coordData);
        eclGrid->exportZCORN(zcornData);
        eclGrid->exportACTNUM(actnumData);

        // conversion to deck units
        deckUnits.siToDeck(mapaxesData, DeckUnits::length);
        deckUnits.siToDeck(coordData, DeckUnits::length);
        deckUnits.siToDeck(zcornData, DeckUnits::length);

        ErtKeyword<float> mapaxesKeyword("MAPAXES", mapaxesData);
        ErtKeyword<float> coordKeyword("COORD", coordData);
        ErtKeyword<float> zcornKeyword("ZCORN", zcornData);
        ErtKeyword<int> actnumKeyword("ACTNUM", actnumData);

        ertHandle_ = ecl_grid_alloc_GRDECL_kw(eclGrid->getNX(),
                                              eclGrid->getNY(),
                                              eclGrid->getNZ(),
                                              zcornKeyword.ertHandle(),
                                              coordKeyword.ertHandle(),
                                              actnumKeyword.ertHandle(),
                                              mapaxesData.size()?mapaxesKeyword.ertHandle():NULL);
#endif // HAVE_ERT && HAVE_DUNE_CORNERPOINT
    }

    ~ErtGrid()
    {
#if HAVE_ERT
        ecl_grid_free(ertHandle_);
#endif
    }


    /*!
     * \brief Save the grid to an .EGRID file.
     */
    void write(const std::string& fileName, int reportStepIdx)
    {
#if HAVE_ERT
        ecl_grid_fwrite_EGRID(ertHandle(), fileName.c_str());
#endif
    }

    ErtHandleType *ertHandle() const
    { return ertHandle_; }

private:
    ErtHandleType *ertHandle_;
};

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This is a smart pointer class for ERT's ecl_rst_file_type
 *        structure.
 */
class ErtRestartFile
{
public:
    ErtRestartFile(const ErtRestartFile &) = delete;

    template <class Simulator>
    ErtRestartFile(const Simulator &simulator, int reportStepIdx)
    {
        std::string caseName = getErtCaseName__(simulator.gridManager());

        restartFileName_ = ecl_util_alloc_filename("./",
                                                   caseName.c_str(),
                                                   /*type=*/ECL_UNIFIED_RESTART_FILE,
                                                   /*writeFormatedOutput=*/false,
                                                   reportStepIdx);

        if (reportStepIdx == 0) {
            restartFileHandle_ = ecl_rst_file_open_write(restartFileName_);
        }
        else {
            restartFileHandle_ = ecl_rst_file_open_append(restartFileName_);
        }
    }

    ~ErtRestartFile()
    {
        ecl_rst_file_close(restartFileHandle_);
        free(restartFileName_);
    }

    template <class Simulator>
    void writeHeader(const Simulator &simulator, int reportStepIdx)
    {
        const auto eclGrid = getEclGrid__(simulator.gridManager());
        const auto eclState = simulator.gridManager().eclState();
        const auto eclSchedule = eclState->getSchedule();

        double secondsElapsed = simulator.time() + simulator.timeStepSize();
        double daysElapsed = secondsElapsed/(24*60*60);

        ecl_rsthead_type rstHeader = { 0 };
        rstHeader.sim_time = simulator.startTime() + secondsElapsed;
        rstHeader.nactive = eclGrid->getNumActive();
        rstHeader.nx = eclGrid->getNX();
        rstHeader.ny = eclGrid->getNY();
        rstHeader.nz = eclGrid->getNZ();
        rstHeader.nwells = 0; // eclSchedule->numWells(reportStepIdx);
        rstHeader.niwelz = 0;
        rstHeader.nzwelz = 0;
        rstHeader.niconz = 0;
        rstHeader.ncwmax = 0;
        rstHeader.phase_sum = ECL_OIL_PHASE | ECL_WATER_PHASE | ECL_GAS_PHASE;
        rstHeader.ncwmax = 0; // eclSchedule->getMaxNumCompletionsForWells(reportStepIdx);

        static const int niwelz = 11; // Number of data elements per well in IWEL array in restart file
        static const int nzwelz = 3;  // Number of 8-character words per well in ZWEL array restart file
        static const int niconz = 14; // Number of data elements per completion in ICON array restart file
        rstHeader.niwelz = niwelz;
        rstHeader.nzwelz = nzwelz;
        rstHeader.niconz = niconz;
        rstHeader.sim_days = daysElapsed;

        ecl_rst_file_fwrite_header(restartFileHandle_,
                                   reportStepIdx,
                                   &rstHeader);
    }

    ecl_rst_file_type *ertHandle() const
    { return restartFileHandle_; }

private:
    char *restartFileName_;
    ecl_rst_file_type *restartFileHandle_;
};

/**
 * \ingroup EclBlackOilSimulator
 *
 * \brief The ErtSolution class wraps the actions that must be done to the
 *        restart file while writing solution variables; it is not a handle on
 *        its own.
 */
class ErtSolution
{
public:
    ErtSolution(const ErtSolution&) = delete;

    ErtSolution(ErtRestartFile &restartHandle)
        : restartHandle_(&restartHandle)
    {  ecl_rst_file_start_solution(restartHandle_->ertHandle()); }

    ~ErtSolution()
    { ecl_rst_file_end_solution(restartHandle_->ertHandle()); }

    template <typename T>
    void add(std::shared_ptr<const ErtKeyword<T>> ertKeyword)
    {
        attachedKeywords_.push_back(ertKeyword);
        ecl_rst_file_add_kw(restartHandle_->ertHandle(), ertKeyword->ertHandle());
    }

    ecl_rst_file_type *ertHandle() const
    { return restartHandle_->ertHandle(); }

private:
    ErtRestartFile *restartHandle_;
    std::list<std::shared_ptr<const ErtBaseKeyword>> attachedKeywords_;
};

/**
 * \ingroup EclBlackOilSimulator
 *
 * \brief The ErtSummary class wraps the actions that must be done to write ECL
 *        summary file.
 *
 * These files log the well performance, etc...
 */
template <class TypeTag>
class ErtSummary
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

    typedef EclWellManager<TypeTag> WellManager;

public:
    ErtSummary(const ErtSummary&) = delete;

    ErtSummary(const Simulator& simulator)
    {
        std::string caseName = getErtCaseName__(simulator.gridManager());

        const auto& eclGrid = simulator.gridManager().eclGrid();

        ertHandle_ = ecl_sum_alloc_writer(caseName.c_str(),
                                          /*formatted=*/false,
                                          /*unified=*/true,
                                          /*joinString=*/":",
                                          simulator.time(),
                                          eclGrid->getNX(),
                                          eclGrid->getNY(),
                                          eclGrid->getNZ());
    }

    ~ErtSummary()
    { ecl_sum_free(ertHandle_); }

    // add all wells in the well manager to the summary output and
    // write the result.
    void writeTimeStep(const WellManager& wellManager)
    { }

    ecl_sum_type *ertHandle() const
    { return ertHandle_; }

private:
    ecl_sum_type *ertHandle_;
};

/**
 * \ingroup EclBlackOilSimulator
 *
 * \brief The ErtSummary class wraps the ERT handles which are required to write a single
 *        time step to the ECL summary file.
 */
template <class TypeTag>
class ErtSummaryTimeStep
{
public:
    ErtSummaryTimeStep(ErtSummary<TypeTag>& summaryHandle,
                       double timeInSeconds,
                       int reportStepIdx)
    {
        double timeInDays = timeInSeconds / (24*60*60);
        ertHandle_ = ecl_sum_add_tstep(summaryHandle.ertHandle(), reportStepIdx, timeInDays);
    }

    // no destructor in this class as ERT takes care of freeing the
    // handle as part of freeing the summary file handle!

    ecl_sum_tstep_type *ertHandle() const
    { return ertHandle_; };

private:
    ecl_sum_tstep_type *ertHandle_;
};

} // namespace Ewoms

#endif // HAVE_ERT

#endif // EWOMS_ERT_WRAPPERS_HH
