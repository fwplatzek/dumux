// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Tests the grid creator class for models using facet coupling.
 */
#include <config.h>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dune/grid/uggrid.hh>
#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dumux/common/parameters.hh>
#include <dumux/mixeddimension/facet/gridcreator.hh>

#ifndef BULKGRIDTYPE // default to ug grid if not provided by CMake
#define BULKGRIDTYPE Dune::UGGrid<3>
#endif

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // parse command line argument parameters
    Dumux::Parameters::init(argc, argv);

    using BulkGrid = BULKGRIDTYPE;
    using FacetGrid = Dune::FoamGrid<2, 3>;
    using EdgeGrid = Dune::FoamGrid<1, 3>;

    using GridCreator = Dumux::FacetCouplingGridCreator<BulkGrid, FacetGrid, EdgeGrid>;
    GridCreator gridCreator;
    gridCreator.makeGrids("grid.msh");

    // check grid sizes
    const auto& bulkGridView = gridCreator.grid<0>().leafGridView();
    const auto& facetGridView = gridCreator.grid<1>().leafGridView();
    const auto& edgeGridView = gridCreator.grid<2>().leafGridView();

    if (bulkGridView.size(0) != 507)
        DUNE_THROW(Dune::InvalidStateException, "Bulk grid has " << bulkGridView.size(0) << " instead of 507 elements");
    if (bulkGridView.size(3) != 153)
        DUNE_THROW(Dune::InvalidStateException, "Bulk grid has " << bulkGridView.size(3) << " instead of 153 vertices");

    if (facetGridView.size(0) != 32)
        DUNE_THROW(Dune::InvalidStateException, "Facet grid has " << facetGridView.size(0) << " instead of 32 elements");
    if (facetGridView.size(2) != 23)
        DUNE_THROW(Dune::InvalidStateException, "Facet grid has " << facetGridView.size(2) << " instead of 23 vertices");

    if (edgeGridView.size(0) != 2)
        DUNE_THROW(Dune::InvalidStateException, "Edge grid has " << edgeGridView.size(0) << " instead of 2 elements");
    if (edgeGridView.size(1) != 3)
        DUNE_THROW(Dune::InvalidStateException, "Edge grid has " << edgeGridView.size(1) << " instead of 3 vertices");

    // all entities should have the domain index 10
    for (const auto& e : elements(bulkGridView))
    {
        const auto insertionIdx = gridCreator.gridFactory<0>().insertionIndex(e);
        const auto domainMarker = gridCreator.elementMarkerMap(0)[insertionIdx];
        if (domainMarker != 10)
            DUNE_THROW(Dune::InvalidStateException, "Bulk element has domain marker " << domainMarker << " instead of 10");
    }

    for (const auto& e : elements(facetGridView))
    {
        const auto insertionIdx = gridCreator.gridFactory<1>().insertionIndex(e);
        const auto domainMarker = gridCreator.elementMarkerMap(1)[insertionIdx];
        if (domainMarker != 10)
            DUNE_THROW(Dune::InvalidStateException, "Facet element has domain marker " << domainMarker << " instead of 10");
    }

    for (const auto& e : elements(edgeGridView))
    {
        const auto insertionIdx = gridCreator.gridFactory<2>().insertionIndex(e);
        const auto domainMarker = gridCreator.elementMarkerMap(2)[insertionIdx];
        if (domainMarker != 10)
            DUNE_THROW(Dune::InvalidStateException, "Edge element has domain marker " << domainMarker << " instead of 10");
    }

    // we have defined boundary segments for the bulk grid
    if (gridCreator.boundaryMarkerMap(0).size() != 200)
        DUNE_THROW(Dune::InvalidStateException, "Bulk grid has " << gridCreator.boundaryMarkerMap(0).size() << " instead of 200 boundary segments");
    if (gridCreator.boundaryMarkerMap(1).size() != 0)
        DUNE_THROW(Dune::InvalidStateException, "Facet grid has " << gridCreator.boundaryMarkerMap(1).size() << " instead of 0 boundary segments");
    if (gridCreator.boundaryMarkerMap(2).size() != 0)
        DUNE_THROW(Dune::InvalidStateException, "Edge grid has " << gridCreator.boundaryMarkerMap(3).size() << " instead of 0 boundary segments");

    // Check boundary boundary segments. There should be mesh
    // file-defined segments on all sides except the positive x-direction
    std::size_t negXCount, posYCount, negYCount, posZCount, negZCount;
    negXCount = posYCount = negYCount = posZCount = negZCount = 0;
    for (const auto& e : elements(bulkGridView))
    {
        for (const auto& is : intersections(bulkGridView, e))
        {
            if (is.boundary())
            {
                const auto n = is.centerUnitOuterNormal();

                if (gridCreator.gridFactory<0>().wasInserted(is))
                {
                    const auto insIdx = gridCreator.gridFactory<0>().insertionIndex(is);
                    const auto marker = gridCreator.boundaryMarkerMap(0)[insIdx];
                    if (Dune::FloatCmp::eq(n[0], -1.0, 1e-6)) // neg x-dir
                    {
                        if (marker != 4)
                            DUNE_THROW(Dune::InvalidStateException, "negative x-face should have marker 4, but has " << marker);
                        negXCount++;
                    }
                    else if (Dune::FloatCmp::eq(n[1], 1.0, 1e-6)) // pos y-dir
                    {
                        if (marker != 5)
                            DUNE_THROW(Dune::InvalidStateException, "positive y-face should have marker 5, but has " << marker);
                        posYCount++;
                    }
                    else if (Dune::FloatCmp::eq(n[1], -1.0, 1e-6)) // neg y-dir
                    {
                        if (marker != 6)
                            DUNE_THROW(Dune::InvalidStateException, "negative y-face should have marker 6, but has " << marker);
                        negYCount++;
                    }
                    else if (Dune::FloatCmp::eq(n[2], 1.0, 1e-6)) // pos z-dir
                    {
                        if (marker != 1)
                            DUNE_THROW(Dune::InvalidStateException, "positive z-face should have marker 1, but has " << marker);
                        posZCount++;
                    }
                    else if (Dune::FloatCmp::eq(n[2], -1.0, 1e-6)) // neg z-dir
                    {
                        if (marker != 2)
                            DUNE_THROW(Dune::InvalidStateException, "negative z-face should have marker 2, but has " << marker);
                        negZCount++;
                    }
                }
                else
                {
                    if (!Dune::FloatCmp::eq(n[0], 1.0, 1e-6)) // pos x-dir
                        DUNE_THROW(Dune::InvalidStateException, "Found a boundary segment not specified in .msh file that is not on pos-x face");
                }
            }
        }
    }

    //! make sure we found the right number of boundary segments
    if (negXCount != 40) DUNE_THROW(Dune::InvalidStateException, "Found " << negXCount << " instead of 0 boundary segments in neg x-direction");
    if (posYCount != 40) DUNE_THROW(Dune::InvalidStateException, "Found " << posYCount << " instead of 0 boundary segments in pos y-direction");
    if (negYCount != 40) DUNE_THROW(Dune::InvalidStateException, "Found " << negYCount << " instead of 0 boundary segments in neg y-direction");
    if (posZCount != 40) DUNE_THROW(Dune::InvalidStateException, "Found " << posZCount << " instead of 0 boundary segments in pos z-direction");
    if (negZCount != 40) DUNE_THROW(Dune::InvalidStateException, "Found " << negZCount << " instead of 0 boundary segments in neg z-direction");

    //! check if we found the right number of embeddings
    const auto& edgeEmbeddings = gridCreator.embeddedEntityMap(2);
    const auto& edgeEmbedments = gridCreator.embedmentMap(2);
    if (edgeEmbeddings.size() != 0) DUNE_THROW(Dune::InvalidStateException, "The grid with lowest dimension can't have embedded entities");
    if (edgeEmbedments.size() != 2) DUNE_THROW(Dune::InvalidStateException, "Found " << edgeEmbedments.size() << " instead of 2 edge element embedments");
    for (const auto& embedments : edgeEmbedments)
        if (embedments.second.size() != 4)
            DUNE_THROW(Dune::InvalidStateException, "edge element is embedded in " << embedments.second.size() << " facet elements instead of 4");

    const auto& facetEmbeddings = gridCreator.embeddedEntityMap(1);
    const auto& facetEmbedments = gridCreator.embedmentMap(1);
    if (facetEmbeddings.size() != 8) DUNE_THROW(Dune::InvalidStateException, "Found " << facetEmbeddings.size() << " instead of 8 embeddings in facet grid");
    if (facetEmbedments.size() != 32) DUNE_THROW(Dune::InvalidStateException, "Found " << facetEmbedments.size() << " instead of 32 facet element embedments");
    for (const auto& embedments : facetEmbedments)
        if (embedments.second.size() != 2)
            DUNE_THROW(Dune::InvalidStateException, "facet element is embedded in " << embedments.second.size() << " bulk elements instead of 2");
    for (const auto& embeddings : facetEmbeddings)
        if (embeddings.second.size() != 1)
            DUNE_THROW(Dune::InvalidStateException, "facet element has " << embeddings.second.size() << " embedded entities instead of 1");

    const auto& bulkEmbeddings = gridCreator.embeddedEntityMap(0);
    const auto& bulkEmbedments = gridCreator.embedmentMap(0);
    if (bulkEmbeddings.size() != 56) DUNE_THROW(Dune::InvalidStateException, "Found " << bulkEmbeddings.size() << " instead of 56 embeddings in bulk grid");
    if (bulkEmbedments.size() != 0) DUNE_THROW(Dune::InvalidStateException, "The grid with highest dimension can't have embedments");

    std::size_t singleEmbeddings = 0;
    std::size_t doubleEmbeddings = 0;
    for (const auto& embeddings : bulkEmbeddings)
    {
        if (embeddings.second.size() != 1 && embeddings.second.size() != 2)
            DUNE_THROW(Dune::InvalidStateException, "bulk element has " << embeddings.second.size() << " embedded entities instead of 1 or 2");
        if (embeddings.second.size() == 1) singleEmbeddings++;
        if (embeddings.second.size() == 2) doubleEmbeddings++;
    }

    if (singleEmbeddings != 48) DUNE_THROW(Dune::InvalidStateException, "Found " << singleEmbeddings << " instead of 48 bulk elements with 1 embedding");
    if (doubleEmbeddings != 8) DUNE_THROW(Dune::InvalidStateException, "Found " << doubleEmbeddings << " instead of 8 bulk elements with 2 embeddings");

    // write .vtk file for each grid
    using BulkWriter = Dune::VTKWriter<typename BulkGrid::LeafGridView>;
    BulkWriter bulkWriter(bulkGridView); bulkWriter.write("bulkgrid");

    using FacetWriter = Dune::VTKWriter<typename FacetGrid::LeafGridView>;
    FacetWriter facetWriter(facetGridView); facetWriter.write("facetgrid");

    using EdgeWriter = Dune::VTKWriter<typename EdgeGrid::LeafGridView>;
    EdgeWriter edgeWriter(edgeGridView); edgeWriter.write("edgegrid");
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception e) {
    std::cout << e << std::endl;
    return 1;
}