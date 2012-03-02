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

#ifndef TESTVERTEXRELAXATION_HPP_
#define TESTVERTEXRELAXATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "OffLatticeSimulationWithEquilibriumStoppingEvent.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "WelikyOsterForce.hpp"
#include "ModifiedWelikyOsterForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "Warnings.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"

class TestVertexRelaxation : public AbstractCellBasedTestSuite
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

    void TestVertexRelaxationWelikyOster() throw (Exception)
    {
        // Create output file
        OutputFileHandler output_file_handler("TestVertexRelaxation", false);
        out_stream p_output_file = output_file_handler.OpenOutputFile("WelikyAreaData.dat");

        // Repeat simulation for increasing number of nodes
        for (unsigned num_nodes=3; num_nodes<21; num_nodes++)
        {
            // Reset the simulation time
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Construct a 2D vertex mesh comprising a single regular polygonal element
            std::vector<Node<2>*> nodes;

            double radius = pow(0.25/sqrt((double)num_nodes*tan(M_PI/(double)num_nodes))/sin(M_PI/(double)num_nodes),1.0/3.0)/sqrt(M_PI);
            for (unsigned i=0; i<num_nodes; i++)
            {
                double theta = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
                nodes.push_back(new Node<2>(i, true, radius*cos(theta), radius*sin(theta)));
            }

            std::vector<VertexElement<2,2>*> elements;
            elements.push_back(new VertexElement<2,2>(0, nodes));

            double cell_swap_threshold = 1e-5; // so vertices never merge
            MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, 1e-4);

            // Create a non-proliferating cell corresponding to this element
            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, 1, std::vector<unsigned>(), DIFFERENTIATED);
            cells[0]->SetBirthTime(-1.0);

            // Create a cell population
            VertexBasedCellPopulation<2> cell_population(mesh, cells);

            // Set up a cell-based simulation
            OffLatticeSimulationWithEquilibriumStoppingEvent<2> simulator(cell_population);

            std::stringstream stream;
            stream << num_nodes;
            std::string output_directory = "TestVertexRelaxationWelikyOster_" + stream.str();
            simulator.SetOutputDirectory(output_directory);

            simulator.SetSamplingTimestepMultiple(500);
            simulator.SetEndTime(1000);

            // Create a force law and pass it to the simulation
            MAKE_PTR(WelikyOsterForce<2>, p_force);
            simulator.AddForce(p_force);

            // Run simulation
            simulator.Solve();

            TS_ASSERT_EQUALS(num_nodes, cell_population.rGetMesh().GetElement(0)->GetNumNodes());

            // Write results to output file
            double end_time = SimulationTime::Instance()->GetTime() - simulator.GetDt()*simulator.GetNumberOfStaticTimesThreshold();
            *p_output_file << num_nodes << "\t" << cell_population.rGetMesh().GetVolumeOfElement(0) << "\t" << end_time << "\n";
        }

        // Close output file and tidy up
        p_output_file->close();
        Warnings::QuietDestroy();
    }

    void TestVertexRelaxationNagaiHonda() throw (Exception)
    {
        // Create output file
        OutputFileHandler output_file_handler("TestVertexRelaxation", false);
        out_stream p_output_file = output_file_handler.OpenOutputFile("NagaiHondaAreaData.dat");

        // Repeat simulation for increasing number of nodes
        for (unsigned num_nodes=3; num_nodes<21; num_nodes++)
        {
            // Reset the simulation time
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Construct a 2D vertex mesh comprising a single regular polygonal element
            std::vector<Node<2>*> nodes;

            double radius = 1.0/sqrt(M_PI);
            for (unsigned i=0; i<num_nodes; i++)
            {
                double theta = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
                nodes.push_back(new Node<2>(i, true, radius*cos(theta), radius*sin(theta)));
            }

            std::vector<VertexElement<2,2>*> elements;
            elements.push_back(new VertexElement<2,2>(0, nodes));

            double cell_swap_threshold = 1e-5; // so vertices never merge
            MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, 0.2);

            // Create a non-proliferating cell corresponding to this element
            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, 1, std::vector<unsigned>(), DIFFERENTIATED);
            cells[0]->SetBirthTime(-1.0);

            // Create a cell population
            VertexBasedCellPopulation<2> cell_population(mesh, cells);

            // Set up a cell-based simulation
            OffLatticeSimulationWithEquilibriumStoppingEvent<2> simulator(cell_population);

            std::stringstream stream;
            stream << num_nodes;
            std::string output_directory = "TestVertexRelaxationNagaiHonda_" + stream.str();
            simulator.SetOutputDirectory(output_directory);

            simulator.SetSamplingTimestepMultiple(500);
            simulator.SetEndTime(1000);

            // Create a force law and pass it to the simulation
            MAKE_PTR(NagaiHondaForce<2>, p_force);
            simulator.AddForce(p_force);

            // Run simulation
            simulator.Solve();

            TS_ASSERT_EQUALS(num_nodes, cell_population.rGetMesh().GetElement(0)->GetNumNodes());

            // Write results to output file
            double end_time = SimulationTime::Instance()->GetTime() - simulator.GetDt()*simulator.GetNumberOfStaticTimesThreshold();
            *p_output_file << num_nodes << "\t" << cell_population.rGetMesh().GetVolumeOfElement(0) << "\t" << end_time << "\n";
        }

        // Close output file and tidy up
        p_output_file->close();
        Warnings::QuietDestroy();
    }

    void TestVertexRelaxationModifiedWelikyOster() throw (Exception)
    {
        // Create output file
        OutputFileHandler output_file_handler("TestVertexRelaxation", false);
        out_stream p_output_file = output_file_handler.OpenOutputFile("ModifiedWelikyAreaData.dat");

        // Repeat simulation for increasing number of nodes
        for (unsigned num_nodes=3; num_nodes<21; num_nodes++)
        {
            // Reset the simulation time
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Construct a 2D vertex mesh comprising a single regular polygonal element
            std::vector<Node<2>*> nodes;

            double radius = 1.0/sqrt(M_PI);
            for (unsigned i=0; i<num_nodes; i++)
            {
                double theta = M_PI+2.0*M_PI*(double)(i)/(double)(num_nodes);
                nodes.push_back(new Node<2>(i, true, radius*cos(theta), radius*sin(theta)));
            }

            std::vector<VertexElement<2,2>*> elements;
            elements.push_back(new VertexElement<2,2>(0, nodes));

            double cell_swap_threshold = 1e-5; // so vertices never merge
            MutableVertexMesh<2,2> mesh(nodes, elements, cell_swap_threshold, 0.2);

            // Create a non-proliferating cell corresponding to this element
            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, 1, std::vector<unsigned>(), DIFFERENTIATED);
            cells[0]->SetBirthTime(-1.0);

            // Create a cell population
            VertexBasedCellPopulation<2> cell_population(mesh, cells);

            // Set up a cell-based simulation
            OffLatticeSimulationWithEquilibriumStoppingEvent<2> simulator(cell_population);

            std::stringstream stream;
            stream << num_nodes;
            std::string output_directory = "TestVertexRelaxationModifiedWelikyOster_" + stream.str();
            simulator.SetOutputDirectory(output_directory);

            simulator.SetSamplingTimestepMultiple(500);
            simulator.SetEndTime(1000);

            // Create a force law and pass it to the simulation
            MAKE_PTR(ModifiedWelikyOsterForce<2>, p_force);
            simulator.AddForce(p_force);

            // Run simulation
            simulator.Solve();

            TS_ASSERT_EQUALS(num_nodes, cell_population.rGetMesh().GetElement(0)->GetNumNodes());

            // Write results to output file
            double end_time = SimulationTime::Instance()->GetTime() - simulator.GetDt()*simulator.GetNumberOfStaticTimesThreshold();
            *p_output_file << num_nodes << "\t" << cell_population.rGetMesh().GetVolumeOfElement(0) << "\t" << end_time << "\n";
        }

        // Close output file and tidy up
        p_output_file->close();
        Warnings::QuietDestroy();
    }
};

#endif /*TESTVERTEXRELAXATION_HPP_*/
