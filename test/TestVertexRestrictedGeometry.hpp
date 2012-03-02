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

#ifndef TESTVERTEXRESTRICTEDGEOMETRY_HPP_
#define TESTVERTEXRESTRICTEDGEOMETRY_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "NagaiHondaForce.hpp"
#include "ModifiedWelikyOsterForce.hpp"
#include "VertexAngleForce.hpp"
#include "PlaneBasedCellKiller.hpp"

//#include "PipeBoundaryCondition.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "ObstructionBoundaryCondition.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

class TestVertexRestrictedGeometry: public AbstractCellBasedTestSuite
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

//    void noTestVertexInPipe() throw (Exception)
//    {
//        // Create a simple 2D MutableVertexMesh
//        HoneycombVertexMeshGenerator generator(3, 7);
//        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
//        p_mesh->Translate(0.5, 1.0);
//        p_mesh->SetCellRearrangementThreshold(0.05);
//
//        // Set up cells, one for each VertexElement
//        std::vector<CellPtr> cells;
//        CellsGenerator<SimpleWntCellCycleModel, 2> cells_generator;
//        cells_generator.GenerateBasicRandom(cells,p_mesh->GetNumElements(), TRANSIT);
//
//        // Create cell population
//        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
//
//        // Set up cell-based simulation
//        CellBasedSimulation<2> simulator(cell_population);
//        simulator.SetOutputDirectory("TestVertexPipe");
//        simulator.SetDt(0.005);
//        simulator.SetEndTime(1.0);
//        simulator.SetSamplingTimestepMultiple(100);
//
//        // Create a force law and pass it to the simulation
//        NagaiHondaForce<2> nagai_honda_force;
//        nagai_honda_force.SetNagaiHondaDeformationEnergyParameter(55.0);
//        nagai_honda_force.SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
//        nagai_honda_force.SetNagaiHondaCellCellAdhesionEnergyParameter(5.0);
//        nagai_honda_force.SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
//        simulator.AddForce(&nagai_honda_force);
//
//        // Create cell killer for x = 0 and pass to the simulation
//        c_vector<double, 2> point = zero_vector<double>(2);
//        point(1) = 12.0;
//        c_vector<double, 2> normal = zero_vector<double>(2);
//        normal(1) = -0.1;
//        PlaneBasedCellKiller<2> cell_killer(&cell_population, point, normal);
//        simulator.AddCellKiller(&cell_killer);
//
//        // Create some boundary conditions and pass them to the simulation
//        PipeBoundaryCondition<2> boundary_condition(&cell_population);
//        simulator.AddCellPopulationBoundaryCondition(&boundary_condition);
//
//        // Create an instance of a Wnt concentration
//        WntConcentration<2>::Instance()->SetType(LINEAR);
//        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
//        WntConcentration<2>::Instance()->SetCryptLength(12);
//
//        // Run simulation
//        simulator.Solve();
//
//        // Tidy up
//        WntConcentration<2>::Instance()->Destroy();
//    }

    void TestVertexWithObstuction() throw (Exception)
    {
        // Create a simple 2D MutableVertexMesh
        CylindricalHoneycombVertexMeshGenerator generator(6, 1, true, 0.6);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();
        p_mesh->SetCellRearrangementThreshold(0.01);

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        CellsGenerator<SimpleWntCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,p_mesh->GetNumElements(), TRANSIT);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexObstruction");
        simulator.SetDt(0.001);
        simulator.SetEndTime(50.0); // 200
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetOutputNodeVelocities(true);

        // Create a force law and pass it to the simulation
//        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
//        p_nagai_honda_force->SetNagaiHondaDeformationEnergyParameter(55.0);
//        p_nagai_honda_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
//        p_nagai_honda_force->SetNagaiHondaCellCellAdhesionEnergyParameter(5.0);
//        p_nagai_honda_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
//        simulator.AddForce(p_nagai_honda_force);

        // Create a force law and pass it to the simulation
        MAKE_PTR(ModifiedWelikyOsterForce<2>, p_weliky_force);
        p_weliky_force->SetWelikyOsterAreaParameter(0.1);      // beta
        p_weliky_force->SetWelikyOsterPerimeterParameter(0.1); // kappa
        simulator.AddForce(p_weliky_force);

        MAKE_PTR(VertexAngleForce<2>, p_angle_force);
        p_angle_force->SetAngleRestrainingParameter(1.0);
        simulator.AddForce(p_angle_force);

        // Create cell killer for y = 10 and pass to the simulation
        c_vector<double, 2> point1 = zero_vector<double>(2);
        point1(1) = 6.0; //12
        c_vector<double, 2> normal1 = zero_vector<double>(2);
        normal1(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, point1, normal1));
        simulator.AddCellKiller(p_killer);



        // Create boundary condition y > 0
        c_vector<double, 2> point2 = zero_vector<double>(2);
        c_vector<double, 2> normal2 = zero_vector<double>(2);
        normal2(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_boundary_condition_1, (&cell_population, point2, normal2));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_1);

        // Create obstruction boundary condition
        MAKE_PTR_ARGS(ObstructionBoundaryCondition<2>, p_boundary_condition_2, (&cell_population));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_2);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(50); // So only cells at very bottom divide

        // Run simulation
        simulator.Solve();

        // Tidy up
        WntConcentration<2>::Instance()->Destroy();
    }
};

#endif /*TESTVERTEXRESTRICTEDGEOMETRY_HPP_*/
