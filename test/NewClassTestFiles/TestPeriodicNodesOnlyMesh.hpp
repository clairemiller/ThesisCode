/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_
#define TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "Periodic3dNodesOnlyMesh.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "LogFile.hpp"
#include "WildTypeCellMutationState.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "CellVolumesWriter.hpp"
#include "FileComparison.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellPopulationAreaWriter.hpp"
#include "NodeVelocityWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestOffLatticeSimulationWithNodeBasedCellPopulation : public AbstractCellBasedWithTimingsTestSuite
{
public:
    void TestSimplePeriodicMonolayer() throw (Exception)
    {
        // Create a simple periodic mesh
        unsigned num_cells_width = 7;
        double periodic_width = 6.0;

        // Generate random initial locations
        std::vector<Node<2>*> nodes(num_cells_width);
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        unsigned i = 0;
        for ( std::vector<Node<2>*>::iterator it = nodes.begin(); it!=nodes.end(); it++ )
        {
            double x = (p_gen->ranf())*periodic_width;
            *it = new Node<2>( i++, false, x, 0.0 );
        }

        // Reseed the random number generator
        p_gen->Reseed(0);

        // Convert this to a Cylindrical2dNodesOnlyMesh 
        Periodic3dNodesOnlyMesh<2> mesh(periodic_width);
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithPeriodicNodeBasedCellPopulation");

        // Run for long enough to see the periodic bounday influencing the cells
        simulator.SetEndTime(10.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't outside the domain
        // for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        // {
        //     TS_ASSERT_LESS_THAN_EQUALS(0,simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0]);
        //     TS_ASSERT_LESS_THAN_EQUALS(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0],periodic_width);
        // }

        // Now run the simulation again with the periodic boundary in a different place and check its the same

        // First reset the singletons
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        p_gen->Reseed(0);

        // Convert this to a Cylindrical2dNodesOnlyMesh
        Periodic3dNodesOnlyMesh<2> mesh_2(periodic_width);
        mesh_2.ConstructNodesWithoutMesh(nodes, 1.5);

        // Add an offset
        double x_offset = periodic_width/2.0;
        mesh_2.Translate(-x_offset,0.0);

        // Create cells
        std::vector<CellPtr> cells_2;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator_2;
        cells_generator_2.GenerateBasicRandom(cells_2, mesh_2.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population_2(mesh_2, cells_2);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator_2(node_based_cell_population_2);
        simulator_2.SetOutputDirectory("TestOffLatticeSimulationWith2ndPeriodicNodeBasedCellPopulation");

        // Run for long enough to see the periodic boundary influencing the cells
        simulator_2.SetEndTime(10.0);

        // Pass the same force law to the simulation
        simulator_2.AddForce(p_linear_force);

        simulator_2.Solve();

        // Check with a different interaction distance
        // First reset the singletons
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        p_gen->Reseed(0);

        // Convert this to a Cylindrical2dNodesOnlyMesh
        Periodic3dNodesOnlyMesh<2> mesh_3(periodic_width);
        mesh_3.ConstructNodesWithoutMesh(nodes, 2.0);

        // Create cells
        std::vector<CellPtr> cells_3;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator_3;
        cells_generator_3.GenerateBasicRandom(cells_3, mesh_3.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population_3(mesh_3, cells_3);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator_3(node_based_cell_population_3);
        simulator_3.SetOutputDirectory("TestOffLatticeSimulationWith3rdPeriodicNodeBasedCellPopulation");

        // Run for long enough to see the periodic boundary influencing the cells
        simulator_3.SetEndTime(10.0);

        // Pass the same force law to the simulation
        simulator_3.AddForce(p_linear_force);

        simulator_3.Solve();


        // Check that nothing's gone badly wrong by testing that nodes aren't outside the domain
        // for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        // {
        //     double x_1 = simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[0];
        //     double x_2 = simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[0];
        //     double x_3 = simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[0];

        //     if (x_1 < x_offset)
        //     {
        //         TS_ASSERT_DELTA(x_1+x_offset, x_2, 1e-6)
        //     }
        //     else
        //     {
        //         TS_ASSERT_DELTA(x_1-x_offset, x_2, 1e-6)
        //     }

        //     TS_ASSERT_DELTA(x_1,x_3,1e-6);

        //     TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[1],simulator_2.rGetCellPopulation().GetNode(i)->rGetLocation()[1],1e-6);
        //     TS_ASSERT_DELTA(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()[1],simulator_3.rGetCellPopulation().GetNode(i)->rGetLocation()[1],1e-6);
        // }
    }

};

#endif /*TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_*/
