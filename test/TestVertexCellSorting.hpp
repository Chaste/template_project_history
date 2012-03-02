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

#ifndef TESTVERTEXCELLSORTING_HPP_
#define TESTVERTEXCELLSORTING_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VolumeTrackedOffLatticeSimulation.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "VertexDiffusionForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "MotileMutationState.hpp"
#include "SmartPointers.hpp"

/**
 * Test suite defining a vertex dynamics simulation of cell sorting.
 */
class TestVertexCellSorting: public AbstractCellBasedTestSuite
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

    /**
     * Simulate a population of cells exhibiting cell sorting using the 
     * vertex dynamics model proposed by T. Nagai and H. Honda ("A dynamic 
     * cell model for the formation of epithelial tissues", Philosophical 
     * Magazine Part B 81:699-719).
     * 
     * Each of the vertex dynamics model parameter member variables are
     * rescaled such that mDampingConstantNormal takes the default value 1,
     * whereas Nagai and Honda (who denote the parameter by nu) take the
     * value 0.01.
     */
    void TestVertexMonolayerCellSorting() throw (Exception)
    {
        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(10, 10);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.05);

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,p_mesh->GetNumElements(), DIFFERENTIATED);

        // Make this pointer first as if we move it after creating the cell population the mutation numbers aren't tracked
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<MotileMutationState>());

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        // Set population to output all data to results files
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellVariables(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAges(true);
        cell_population.SetOutputCellVolumes(true);

        // Ran
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (RandomNumberGenerator::Instance()->ranf() < 0.5)
            {
                (*cell_iter)->SetMutationState(p_state);
            }
        }

        // Set up cell-based simulation and output directory
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexCellSorting");

        // Set time step and end time for simulation
        simulator.SetDt(0.005);
        simulator.SetEndTime(50.0);

        // Only record results every 100 time steps
        simulator.SetSamplingTimestepMultiple(200);

        // Set up force law and pass it to the simulation
        MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(550.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(5.0); //5
        p_force->SetNagaiHondaLabeledCellCellAdhesionEnergyParameter(30.0); //30
        p_force->SetNagaiHondaLabeledCellLabeledCellAdhesionEnergyParameter(15.0); //15
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(60.0); //60
        p_force->SetNagaiHondaLabeledCellBoundaryAdhesionEnergyParameter(40.0); //40
        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(VertexDiffusionForce<2>, p_random_force, (0.5));
        simulator.AddForce(p_random_force);
        
        // Run simulation
        simulator.Solve();
   }
};

#endif /*TESTVERTEXCELLSORTING_HPP_*/
