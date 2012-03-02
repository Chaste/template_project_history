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

#ifndef TESTVERTEXACTIVEMIGRATION_HPP_
#define TESTVERTEXACTIVEMIGRATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"


#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "NagaiHondaForce.hpp"
#include "ModifiedWelikyOsterForce.hpp"
#include "VertexAngleForce.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "MotileMutationState.hpp"
#include "MotileCellForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

class TestVertexActiveMigration: public AbstractCellBasedTestSuite
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

public:

    void noTestNHVertexMonolayerWithActiveMigration() throw (Exception)
    {
        double box_height = 10.0;
        double box_width = 20.0;

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(24, 14);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->Scale(0.8,0.8);

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,p_mesh->GetNumElements(), DIFFERENTIATED);

        // Make this pointer first as if we move it after creating the cell population the mutation numbers arent tracked
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<MotileMutationState>());

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellVariables(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAges(true);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> centre_of_cell = cell_population.GetLocationOfCellCentre(*cell_iter);
            double x = centre_of_cell[0];
            double y = centre_of_cell[1];
            if  (fabs((x-1.5)*(x-1.5) + (y-5)*(y-5))<0.2)
            {
                (*cell_iter)->SetMutationState(p_state);
            }
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNHVertexCellMoving");
        simulator.SetDt(0.01);
        simulator.SetEndTime(1.0); //30.0
        simulator.SetSamplingTimestepMultiple(100);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(5.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        simulator.AddForce(p_force);

        // Create boundary condition y > 0
        c_vector<double, 2> point = zero_vector<double>(2);
        c_vector<double, 2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_boundary_condition_1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_1);

        // Create boundary condition y < box_height
        point(1) = box_height;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_boundary_condition_2, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_2);

        // Create boundary condition x > 0
        point(1) = 0.0;
        normal(0) = -1.0;
        normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_boundary_condition_3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_3);

        // Create boundary condition x < box_width
        point(0) = box_width;
        normal(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_boundary_condition_4, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_4);
        
        // Smooth out edges of tissue
        VertexBasedCellPopulation<2>* p_vertex_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation()));

        for (VertexMesh<2,2>::NodeIterator node_iter = p_vertex_cell_population->rGetMesh().GetNodeIteratorBegin();
             node_iter != p_vertex_cell_population->rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();

            std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
            if (containing_element_indices.size() == 1)
            {
                unsigned element_index = (*containing_element_indices.begin());

                // Remove node from element and delete node
                p_vertex_cell_population->rGetMesh().GetElement(element_index)->DeleteNode(p_vertex_cell_population->rGetMesh().GetElement(element_index)->GetNodeLocalIndex(node_index));
                p_vertex_cell_population->rGetMesh().DeleteNodePriorToReMesh(node_index);
            }
        }
        p_vertex_cell_population->rGetMesh().ReMesh();

        // Run simulation
        simulator.Solve();

        // Now flag a cell as motile
        MAKE_PTR(MotileCellForce<2>, p_motile_cell_force);
        simulator.AddForce(p_motile_cell_force);
        
        simulator.SetEndTime(2.0); //80

        simulator.SetSamplingTimestepMultiple(10);
        simulator.Solve();
    }

    void TestWOVertexMonolayerWithActiveMigration() throw (Exception)
   {
       double box_height = 5.0;
       double box_width = 10.0;

       // Create a simple 2D MutableVertexMesh
       HoneycombVertexMeshGenerator generator(12, 7);
       MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
       p_mesh->Scale(0.4,0.4);
       p_mesh->SetCheckForInternalIntersections(true);

       // Set up cells, one for each VertexElement
       std::vector<CellPtr> cells;
       CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
       cells_generator.GenerateBasicRandom(cells,p_mesh->GetNumElements(), DIFFERENTIATED);

       // Make this pointer first as if we move it after creating the cell population the mutation numbers arent tracked
       boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<MotileMutationState>());

       // Create cell population
       VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
       cell_population.SetOutputCellIdData(true);
       cell_population.SetOutputCellMutationStates(true);
       cell_population.SetOutputCellAncestors(true);
       cell_population.SetOutputCellProliferativeTypes(true);
       cell_population.SetOutputCellVariables(true);
       cell_population.SetOutputCellCyclePhases(true);
       cell_population.SetOutputCellAges(true);

       for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
       {
           c_vector<double, 2> centre_of_cell = cell_population.GetLocationOfCellCentre(*cell_iter);
           double x = centre_of_cell[0];
           double y = centre_of_cell[1];
           if  (fabs((x-1.2)*(x-1.0) + (y-1.25)*(y-1.25))<0.05)
           {
               (*cell_iter)->SetMutationState(p_state);
           }
       }

       // Set up cell-based simulation
       OffLatticeSimulation<2> simulator(cell_population);
       simulator.SetOutputDirectory("TestWOVertexCellMoving");
       simulator.SetDt(0.001);
       simulator.SetEndTime(1.0); //30.0
       simulator.SetSamplingTimestepMultiple(100);

       // Create a force law and pass it to the simulation
       MAKE_PTR(ModifiedWelikyOsterForce<2>, p_weliky_force);
       p_weliky_force->SetWelikyOsterAreaParameter(0.1);      // beta
       p_weliky_force->SetWelikyOsterPerimeterParameter(0.1); // kappa
       simulator.AddForce(p_weliky_force);

       MAKE_PTR(VertexAngleForce<2>, p_angle_force);
       p_angle_force->SetAngleRestrainingParameter(1.0);
       simulator.AddForce(p_angle_force);

       // Create boundary condition y > 0
       c_vector<double, 2> point = zero_vector<double>(2);
       c_vector<double, 2> normal = zero_vector<double>(2);
       normal(1) = -1.0;
       MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_boundary_condition_1, (&cell_population, point, normal));
       simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_1);

       // Create boundary condition y < box_height
       point(1) = box_height/2.0;
       normal(1) = 1.0;
       MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_boundary_condition_2, (&cell_population, point, normal));
       simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_2);


       // Create boundary condition x > 0
       point(1) = 0.0;
       normal(0) = -1.0;
       normal(1) = 0.0;
       MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_boundary_condition_3, (&cell_population, point, normal));
       simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_3);

       // Create boundary condition x < box_width
       point(0) = box_width/2.0;
       normal(0) = 1.0;
       MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_boundary_condition_4, (&cell_population, point, normal));
       simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_4);

//       // Smooth out edges of tissue
//       VertexBasedCellPopulation<2>* p_vertex_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation()));
//
//       for (VertexMesh<2,2>::NodeIterator node_iter = p_vertex_cell_population->rGetMesh().GetNodeIteratorBegin();
//            node_iter != p_vertex_cell_population->rGetMesh().GetNodeIteratorEnd();
//            ++node_iter)
//       {
//           unsigned node_index = node_iter->GetIndex();
//
//           std::set<unsigned> containing_element_indices = node_iter->rGetContainingElementIndices();
//           if (containing_element_indices.size() == 1)
//           {
//               unsigned element_index = (*containing_element_indices.begin());
//
//               // Remove node from element and delete node
//               p_vertex_cell_population->rGetMesh().GetElement(element_index)->DeleteNode(p_vertex_cell_population->rGetMesh().GetElement(element_index)->GetNodeLocalIndex(node_index));
//               p_vertex_cell_population->rGetMesh().DeleteNodePriorToReMesh(node_index);
//           }
//       }
//       p_vertex_cell_population->rGetMesh().ReMesh();

       // Run simulation
       simulator.Solve();

       // Now flag a cell as motile
       MAKE_PTR(MotileCellForce<2>, p_motile_cell_force);
       simulator.AddForce(p_motile_cell_force);

       simulator.SetEndTime(4.0); //80
       simulator.SetSamplingTimestepMultiple(10);
       simulator.Solve();
   }

};

#endif /*TESTVERTEXACTIVEMIGRATION_HPP_*/
