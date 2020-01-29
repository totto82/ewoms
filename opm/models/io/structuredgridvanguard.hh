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
 * \copydoc Opm::StructuredGridVanguard
 */
#ifndef EWOMS_STRUCTURED_GRID_VANGUARD_HH
#define EWOMS_STRUCTURED_GRID_VANGUARD_HH

#include <opm/models/io/basevanguard.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include "opm/grid/polyhedralgrid.hh"
#include "opm/grid/cart_grid.h"

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#endif

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>
#include <memory>

namespace Opm
{

template <class TypeTag>
class StructuredGridVanguard;

} // namespace Opm

BEGIN_PROPERTIES

NEW_TYPE_TAG(StructuredGridVanguard);

// declare the properties required by the for the structured grid simulator vanguard
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);

NEW_PROP_TAG(DomainSizeX);
NEW_PROP_TAG(DomainSizeY);
NEW_PROP_TAG(DomainSizeZ);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);

NEW_PROP_TAG(GridGlobalRefinements);

// GRIDDIM is only set by the finger problem
#ifndef GRIDDIM
static const int dim = 3;
#else
static const int dim = 3; //GRIDDIM;
#endif

// set the Grid and Vanguard properties
#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(StructuredGridVanguard, Grid, Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming>);
#else
SET_TYPE_PROP(StructuredGridVanguard, Grid, Dune::PolyhedralGrid<dim, dim>);
#endif

SET_TYPE_PROP(StructuredGridVanguard, Vanguard, Opm::StructuredGridVanguard<TypeTag>);

END_PROPERTIES

namespace Opm
{

/*!
 * \ingroup TestProblems
 *
 * \brief Helper class for grid instantiation of the lens problem.
 */
template <class TypeTag>
class StructuredGridVanguard : public BaseVanguard<TypeTag>
{
    typedef BaseVanguard<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

    typedef std::unique_ptr<Grid> GridPointer;

    static const int dim = Grid::dimension;

public:
    /*!
     * \brief Register all run-time parameters for the structured grid simulator vanguard.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, GridGlobalRefinements,
                             "The number of global refinements of the grid "
                             "executed after it was loaded");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeX,
                             "The size of the domain in x direction");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, CellsX,
                             "The number of intervalls in x direction");
        if (dim > 1)
        {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeY,
                                 "The size of the domain in y direction");
            EWOMS_REGISTER_PARAM(TypeTag, unsigned, CellsY,
                                 "The number of intervalls in y direction");
        }
        if (dim > 2)
        {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeZ,
                                 "The size of the domain in z direction");
            EWOMS_REGISTER_PARAM(TypeTag, unsigned, CellsZ,
                                 "The number of intervalls in z direction");
        }
    }

    /*!
     * \brief Create the grid for the lens problem
     */
    StructuredGridVanguard(Simulator &simulator)
        : ParentType(simulator)
    {
        Dune::FieldVector<int, dim> cellRes;

        typedef double GridScalar;
        Dune::FieldVector<GridScalar, dim> upperRight;
        Dune::FieldVector<GridScalar, dim> lowerLeft(0);

        upperRight[0] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeX);
        upperRight[1] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeY);

        cellRes[0] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsX);
        cellRes[1] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsY);
        if (dim == 3)
        {
            upperRight[2] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeZ);
            cellRes[2] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsZ);
        }
#if 0
        if (dim == 3)
        {
            std::cout << cellRes[0] << cellRes[1] << ", " << cellRes[2] << ", " << upperRight[0] << ", " << upperRight[1] << ", " << upperRight[2] << std::endl;


            UnstructuredGrid *grid = create_grid_hexa3d(cellRes[0], cellRes[1], cellRes[2],
                                                        upperRight[0] / cellRes[0], upperRight[1] / cellRes[1], upperRight[2] / cellRes[2]);

            Grid polyGrid(*grid);
            GridPointer polygrid(new Grid(*grid));
            gridPtr_ = std::move(polygrid);
        }
        else if (dim == 2)
        {
            UnstructuredGrid *grid = create_grid_cart2d(cellRes[0], cellRes[1], upperRight[0] / cellRes[0], upperRight[1] / cellRes[1]);
            Grid polyGrid(*grid);
            GridPointer polygrid(new Grid(*grid));
            gridPtr_ = std::move(polygrid);
        }
                else
        {
            throw "Not implemented for current dimension";
        }
#else
//        std::string file_name{"/home/runar/work/porepy/src/porepy/opm_interface/ioWriter.py"};
//        std::ostringstream command;
//        command << "python " << file_name << " " << cellRes[0] << " " << cellRes[1];
//        if (dim == 3)
//            command << " " << cellRes[2];
//        command << " " << upperRight[0] / cellRes[0] << " " << upperRight[1] / cellRes[1];
//        if (dim == 3)
//            command << " " << upperRight[2] / cellRes[2];
//        command << " ./porepy_grid.txt";

//        std::cout << "Calling python: " << std::endl
//                  << std::flush;
//        system(command.str().c_str());
//        std::cout << "finished " << std::endl
//                  << std::flush;

        std::cout << " ---- dim -------- " << dim << std::endl;
        const char *c_str = dim == 2 ? "porepy_grid2d.txt" : "porepy_grid3d.txt";
        std::cout << "reading grid " << std::endl
                  << std::flush;
        UnstructuredGrid *grid = read_grid(c_str);
        std::cout << grid->dimensions << std::endl;
        std::cout << "finished " << std::endl
                  << std::flush;

        std::cout <<"printing grid"<< std::endl;
        //print_grid(grid);
 

        Grid polyGrid(*grid);
        GridPointer polygrid(new Grid(*grid));
        gridPtr_ = std::move(polygrid);
        //UnstructuredGrid *grid = create_grid_hexa3d(2, 2, 2, 1, 1, 1);
#endif

        /*       std::stringstream dgffile;
        dgffile << "DGF" << std::endl;
        dgffile << "INTERVAL" << std::endl;
        dgffile << lowerLeft  << std::endl;
        dgffile << upperRight << std::endl;
        dgffile << cellRes    << std::endl;
        dgffile << "#" << std::endl;
        dgffile << "GridParameter" << std::endl;
        dgffile << "overlap 1" << std::endl;
        dgffile << "#" << std::endl;
        dgffile << "Simplex" << std::endl;
        dgffile << "#" << std::endl;

        // use DGF parser to create a grid from interval block
        gridPtr_.reset( Dune::GridPtr< Grid >( dgffile ).release() );

        unsigned numRefinements = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);
        gridPtr_->globalRefine(static_cast<int>(numRefinements));
 */
        this->finalizeInit_();
    }

    /*!
     * \brief Return a reference to the grid object.
     */
    Grid &grid()
    {
        return *gridPtr_;
    }

    /*!
     * \brief Return a constant reference to the grid object.
     */
    const Grid &grid() const
    {
        return *gridPtr_;
    }

private:
    GridPointer gridPtr_;
};

} // namespace Opm

#endif
