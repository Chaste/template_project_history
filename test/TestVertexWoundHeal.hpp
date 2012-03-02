/*

Copyright (c) 2005-2012, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTVERTEXWOUNDHEAL_HPP_
#define TESTVERTEXWOUNDHEAL_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "VolumeTrackedOffLatticeSimulation.hpp"
#include "NagaiHondaForce.hpp"
#include "WelikyOsterForce.hpp"
#include "ModifiedWelikyOsterForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

/**
 * Test suite defining a vertex dynamics simulation of a wound-healing
 * cell migration assay, similar to that developed by Nagai and Honda
 * (2006).
 */
class TestVertexWoundHeal: public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

    /**
     * Helper method. Smooth out edges of a vertex mesh.
     * 
     * @param rCellPopulation a cell population
     */
    void SmoothVertexMeshEdges(VertexBasedCellPopulation<2>& rCellPopulation)
    {
        MutableVertexMesh<DIM, DIM>& r_mesh = rCellPopulation.rGetMesh();

        for (VertexMesh<2,2>::NodeIterator node_iter = r_mesh.GetNodeIteratorBegin();
             node_iter != r_mesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            if (containing_element_indices.size() == 1)
            {
                // Get this element
                unsigned elem_index = (*containing_element_indices.begin());
                VertexElement<2,2>* p_element = r_mesh.GetElement(elem_index);

                // Remove node from this element and delete the node
                p_element->DeleteNode(p_element->GetNodeLocalIndex(node_index));
                r_mesh.DeleteNodePriorToReMesh(node_index);
            }
        }
        r_mesh.ReMesh();
    }

    /**
     * Helper method. Iterate over all cells and define the 'hole' by
     * killing those cells whose centres are located in a given region.
     * 
     * @param rCellPopulation a cell population
     * @param holeWidth the width of the hole
     * @param xMin the left boundary of the hole
     * @param xMax the right boundary of the hole
     * @param yMin the bottom boundary of the hole
     * @param yMax the top boundary of the hole
     */
    void CreateHoleInCellPopulation(AbstractCellPopulation<2>& rCellPopulation,
                                    double holeWidth,
                                    double xMin,
                                    double xMax,
                                    double yMin,
                                    double yMax)
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            // Get the coordinates of this cell centre
            c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double x = centre_of_cell[0];
            double y = centre_of_cell[1];

            if ((fabs(y-x)<holeWidth) && (x>xMin) && (x<xMax) && (y>yMin) && (y<yMax))
            {
                cell_iter->Kill();
            }
        }

        rCellPopulation.Update();
    }

public:

    void TestNagaiHondaMonolayerWithHole() throw (Exception)
    {
        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(12, 13);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->Scale(0.8, 0.8);

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        CellsGenerator<ContactInhibitionCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,p_mesh->GetNumElements(), DIFFERENTIATED);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up a cell-based simulation, output directory, time step and end time
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNagaiHondaMonolayerWithHole");
        simulator.SetDt(0.01);
        simulator.SetEndTime(10.0-0.01);
        simulator.SetSamplingTimestepMultiple(200);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);          // lambda
        p_force->SetMatureCellTargetArea(1.0);                           // A_0
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);       // beta
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(5.0);      // gamma_cell
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0); // gamma_boundary
        simulator.AddForce(p_force);

        // Create boundary condition y > 0
        c_vector<double, 2> point1 = zero_vector<double>(2);
        c_vector<double, 2> normal1 = zero_vector<double>(2);
        normal1(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition1, (&cell_population, point1, normal1));

        // Create boundary condition x > 0
        c_vector<double, 2> point2 = zero_vector<double>(2);
        c_vector<double, 2> normal2 = zero_vector<double>(2);
        normal2(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition2, (&cell_population, point2, normal2));

        // Create boundary condition y < 10
        c_vector<double, 2> point3 = zero_vector<double>(2);
        point3(1) = 10.0;
        normal1(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition3, (&cell_population, point3, normal1));

        // Create boundary condition x < 10
        c_vector<double, 2> point4 = zero_vector<double>(2);
        point4(0) = 10.0;
        normal2(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition4, (&cell_population, point4, normal2));

        // Pass the boundary conditions to the simulation
        simulator.AddCellPopulationBoundaryCondition(p_condition1);
        simulator.AddCellPopulationBoundaryCondition(p_condition2);
        simulator.AddCellPopulationBoundaryCondition(p_condition3);
        simulator.AddCellPopulationBoundaryCondition(p_condition4);

        // Smooth out edges to get nice box domain
        SmoothVertexMeshEdges(cell_population);

        // Run simulation
        simulator.Solve();

        // Now remove cells in a given region using a helper method
        CreateHoleInCellPopulation(simulator.rGetCellPopulation(), 1.75, 2.0, 8.0, 2.0, 8.0);
        SmoothVertexMeshEdges(simulator.rGetCellPopulation());

        // Reset timestep, sampling timestep and end time for simulation and run for a further duration
        simulator.SetDt(0.005);
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetEndTime(20.0);
        simulator.Solve();
    }

    void TestWelikyOsterMonolayerWithHole() throw (Exception)
    {
        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(12, 13);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->Scale(0.8, 0.8);

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        CellsGenerator<ContactInhibitionCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,p_mesh->GetNumElements(), DIFFERENTIATED);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up a cell-based simulation, output directory, time step and end time
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestWelikyOsterMonolayerWithHole");
        simulator.SetDt(0.01);
        simulator.SetEndTime(10.0-0.01);
        simulator.SetSamplingTimestepMultiple(200);

        // Create a force law and pass it to the simulation
        MAKE_PTR(WelikyOsterForce<2>, p_force);
        p_force->SetWelikyOsterAreaParameter(2.5);                // beta // Chosen so that resting size of cell with 6 vertices is about 1
        p_force->SetWelikyOsterPerimeterParameter(1.0/sqrt(2.5)); // kappa // Chosen so that resting size of cell with 6 vertices is about 1
        simulator.AddForce(p_force);

        // Create boundary condition y > 0
        c_vector<double, 2> point1 = zero_vector<double>(2);
        c_vector<double, 2> normal1 = zero_vector<double>(2);
        normal1(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition1, (&cell_population, point1, normal1));

        // Create boundary condition x > 0
        c_vector<double, 2> point2 = zero_vector<double>(2);
        c_vector<double, 2> normal2 = zero_vector<double>(2);
        normal2(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition2, (&cell_population, point2, normal2));

        // Create boundary condition y < 10
        c_vector<double, 2> point3 = zero_vector<double>(2);
        point3(1) = 10.0;
        normal1(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition3, (&cell_population, point3, normal1));

        // Create boundary condition x < 10
        c_vector<double, 2> point4 = zero_vector<double>(2);
        point4(0) = 10.0;
        normal2(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_condition4, (&cell_population, point4, normal2));

        // Pass the boundary conditions to the simulation
        simulator.AddCellPopulationBoundaryCondition(p_condition1);
        simulator.AddCellPopulationBoundaryCondition(p_condition2);
        simulator.AddCellPopulationBoundaryCondition(p_condition3);
        simulator.AddCellPopulationBoundaryCondition(p_condition4);

        // Smooth out edges to get nice box domain
        SmoothVertexMeshEdges(cell_population);

        // Run simulation
        simulator.Solve();

        // Now remove cells in a given region using a helper method
        CreateHoleInCellPopulation(simulator.rGetCellPopulation(), 1.75, 2.0, 8.0, 2.0, 8.0);
        SmoothVertexMeshEdges(simulator.rGetCellPopulation());

        // Reset timestep, sampling timestep and end time for simulation and run for a further duration
        simulator.SetDt(0.005);
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetEndTime(20.0);
        simulator.Solve();
    }

};

#endif /*TESTVERTEXMONOLAYERDEMO_HPP_*/
