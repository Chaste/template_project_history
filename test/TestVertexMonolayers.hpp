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

#ifndef TESTVERTEXMONOLAYERS_HPP_
#define TESTVERTEXMONOLAYERS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VolumeTrackedOffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "NagaiHondaForce.hpp"
#include "WelikyOsterForce.hpp"
#include "ModifiedWelikyOsterForce.hpp"
#include "VertexAngleForce.hpp"

#include "AbstractCellBasedTestSuite.hpp"

/**
 * Test suite defining a vertex dynamics simulation of a growing monolayer.
 */
class TestVertexMonolayers: public AbstractCellBasedTestSuite
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

        // Call Destroy() on CellwiseData singleton class to avoid memory leaks
        CellwiseData<2>::Destroy();

        AbstractCellBasedTestSuite::tearDown();
    }

    /**
     * Helper method. Populates a vector of cells of specified size.
     * Each cell has a simple stochastic area-dependent cell-cycle
     * model: it enters quiescence if its area drops below a specified
     * threshold. 
     */
    void SetUpCells(std::vector<CellPtr>& rCells, unsigned numCells, double quiescentVolumeFraction)
    {
        rCells.clear();
        rCells.reserve(numCells);

        for (unsigned i=0; i<numCells; i++)
        {
            // Create a cell-cycle model and set contact inhibition parameter
        	ContactInhibitionCellCycleModel* p_cell_cycle_model = new ContactInhibitionCellCycleModel;
            p_cell_cycle_model->SetDimension(2);
            p_cell_cycle_model->SetCellProliferativeType(TRANSIT);
            p_cell_cycle_model->SetQuiescentVolumeFraction(quiescentVolumeFraction);
            p_cell_cycle_model->SetSDuration(1e-10);
            p_cell_cycle_model->SetG2Duration(1e-10);

            // Create a cell using the cell-cycle model and a 'wild-type' cell mutation state
            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            // Give the cell a random birth time
            double birth_time = -p_cell_cycle_model->GetAverageTransitCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
            p_cell->SetBirthTime(birth_time);

            rCells.push_back(p_cell);
        }
    }

public:

    /**
     * Simulate a monolayer using the vertex dynamics model proposed by
     * T. Nagai and H. Honda ("A dynamic cell model for the formation of
     * epithelial tissues", Philosophical Magazine Part B 81:699-719).
     * 
     * Each of the vertex dynamics model parameter member variables are
     * rescaled such that mDampingConstantNormal takes the default value 1,
     * whereas Nagai and Honda (who denote the parameter by nu) take the
     * value 0.01.
     */
    void TestNagaiHondaMonolayer() throw (Exception)
    {
        // Create a simple 2D vertex-based mesh with a single element
        HoneycombVertexMeshGenerator mesh_generator(1, 1);
        MutableVertexMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Set up a vector of cells
        std::vector<CellPtr> cells;

        // Create 1 cell, with stochastic area-dependent cell-cycle model with phi=0.1
        SetUpCells(cells, 1, 0.1);

        // Create a vertex-based cell population
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

        // Set up a cell-based simulation, output directory, time step and end time
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNagaiHondaMonolayer");
        simulator.SetDt(0.001);
        simulator.SetEndTime(200);

        // Only record results every 100 time steps
        simulator.SetSamplingTimestepMultiple(100);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);          // lambda
        p_force->SetMatureCellTargetArea(1.0);                           // A_0
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);       // beta
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(5.0);      // gamma_cell
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0); // gamma_boundary
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();
    }

    /**
     * Simulate a monolayer using a modified version of the vertex dynamics
     * model proposed by M. Weliky and G. Oster ("The mechanical basis of
     * cell rearrangement. I. Epithelial morphogenesis during Fundulus epiboly",
     * Development 109:373-386).
     *
     * The default values for the two model parameter member variables are
     * our own best estimates, since they are not given in the Weliky and
     * Oster paper.
     */
    void TestModifiedWelikyOsterMonolayer() throw (Exception)
    {
        // Create a simple 2D vertex-based mesh with a single element
        HoneycombVertexMeshGenerator mesh_generator(1, 1);
        MutableVertexMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Set up a vector of cells
        std::vector<CellPtr> cells;

        // Create 1 cell, with stochastic area-dependent cell-cycle model with phi=0.1
        SetUpCells(cells, 1, 0.1);

        // Create a vertex-based cell population
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

        // Set up a cell-based simulation, output directory, time step and end time
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestModifiedWelikyOsterMonolayer");
        simulator.SetDt(0.005);
        simulator.SetEndTime(300);

        // Only record results every 100 time steps
        simulator.SetSamplingTimestepMultiple(100);

        // Create a force law and pass it to the simulation
        MAKE_PTR(ModifiedWelikyOsterForce<2>, p_force);
        p_force->SetWelikyOsterAreaParameter(0.02);      // beta
        p_force->SetWelikyOsterPerimeterParameter(0.02); // kappa
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();
    }

    /**
     * Simulate a monolayer using a modified version of the vertex dynamics
     * model proposed by M. Weliky and G. Oster ("The mechanical basis of
     * cell rearrangement. I. Epithelial morphogenesis during Fundulus epiboly",
     * Development 109:373-386), including an additional membrane stiffness
     * force.
     *
     * The default values for the two model parameter member variables are
     * our own best estimates, since they are not given in the Weliky and
     * Oster paper.
     */
    void TestModifiedWelikyOsterMonolayerWithVertexAngleForce() throw (Exception)
    {
        // Create a simple 2D vertex-based mesh with a single element
        HoneycombVertexMeshGenerator mesh_generator(1, 1);
        MutableVertexMesh<2,2>* p_mesh = mesh_generator.GetMesh();

        // Set up a vector of cells
        std::vector<CellPtr> cells;

        // Create 1 cell, with stochastic area-dependent cell-cycle model with phi=0.1
        SetUpCells(cells, 1, 0.1);

        // Create a vertex-based cell population
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

        // Set up a cell-based simulation, output directory, time step and end time
        VolumeTrackedOffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestModifiedWelikyOsterMonolayerWithVertexAngleForce");
        simulator.SetDt(0.005);
        simulator.SetEndTime(200);

        // Only record results every 100 time steps
        simulator.SetSamplingTimestepMultiple(100);

        // Create a force law and pass it to the simulation
        MAKE_PTR(ModifiedWelikyOsterForce<2>, p_force);
        p_force->SetWelikyOsterAreaParameter(0.1);      // beta
        p_force->SetWelikyOsterPerimeterParameter(0.1); // kappa
        simulator.AddForce(p_force);

        // Create another force law and pass it to the simulation
        MAKE_PTR(VertexAngleForce<2>, p_force2);
        p_force2->SetAngleRestrainingParameter(1.0); // gamma
        simulator.AddForce(p_force2);

        // Run simulation
        simulator.Solve();
    }
};

#endif /*TESTVERTEXMONOLAYERS_HPP_*/
