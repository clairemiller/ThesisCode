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

#ifndef TESTGRADIENTCELLKILLER_HPP_
#define TESTGRADIENTCELLKILLER_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>
#include "AbstractCellBasedTestSuite.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "StemCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"

#include "GradientCellKiller.hpp"
#include "NoNeighbourCellKiller.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * This class contains tests for the cell killer for cells with no neighbours
 */
class TestNoNeighbourCellKiller : public AbstractCellBasedTestSuite
{
    public:
    void TestCellKiller() throw(Exception)
    {
        // Set the node position
        std::vector<Node<2>*> nodes(2);
        nodes[0] = new Node<2>(0, false, 0.0, 0.0);
        nodes[1] = new Node<2>(1, false, 0.0,10.0);
        c_vector<double,2> initial_location_node2 = nodes[1]->rGetLocation();

        // Create the mesh and cell population
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        // Get a reference to the cells held in cell population
        std::list<CellPtr>& r_cells = cell_population.rGetCells();

        // Assign neighbours
        r_cells.front()->GetCellData()->SetItem("n_neighbours",0);
        r_cells.back()->GetCellData()->SetItem("n_neighbours",1);

        // Create the cell killer
        NoNeighbourCellKiller<2> cell_killer(&cell_population);

        // Check that apoptosis happens for the first cell only
        cell_killer.CheckAndLabelCellsForApoptosisOrDeath();
        TS_ASSERT(r_cells.front()->HasApoptosisBegun());
        TS_ASSERT( !(r_cells.back()->HasApoptosisBegun()) );

        // Set up singleton classes
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(32.0, 32);
        double death_time = p_simulation_time->GetTime() + r_cells.front()->GetApoptosisTime();
        
        // Increment time to a time after cell death
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_DELTA(p_simulation_time->GetTime(), 1.0, 1e-3);
        TS_ASSERT_DELTA(death_time, 0.25, 1e-3);
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(death_time+1.0, 1);
        p_simulation_time->IncrementTimeOneStep();

        // Remove dead cells
        cell_population.RemoveDeadCells();
        TS_ASSERT( !(r_cells.front()->IsDead()) );

        // Check only one cell at correct location remains
        TS_ASSERT_EQUALS((int) r_cells.size(),1);
        Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(r_cells.front());
        c_vector<double, 2> location = p_node->rGetLocation();
        TS_ASSERT_EQUALS(initial_location_node2[1], location[1]);
    }
};




/**
 * This class contains tests for methods on the gradient cell killer class
 */
class TestGradientCellKiller : public AbstractCellBasedTestSuite
{
private:
    static double constScalingKillRateFunction(c_vector<double,2> loc)
    {
        return (1.0);
    }
    static double boolScalingKillRateFunction(c_vector<double,2> loc)
    {
        return ( (double) (loc[1] > 1.5) );
    }
public:
    void TestGradientCellKillerClassFunctions() throw(Exception)
    {
        // Set up singleton classes
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(32.0, 32);

        // Set the node position
        std::vector<Node<2>*> nodes(3);
        for ( int i = 0; i < 3; i++ )
        {
            nodes[i] = new Node<2>(i, false, 0.0, (double) i);
        }
        c_vector<double,2> initial_location_node2 = nodes[1]->rGetLocation();

        // Create the mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
 
        double death_time = p_simulation_time->GetTime() + cells[0]->GetApoptosisTime();


        // Create the population
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        // Get a reference to the cells held in cell population
        std::list<CellPtr>& r_cells = cell_population.rGetCells();

        // Create the cell killer
        GradientCellKiller<2> grad_cell_killer(&cell_population,this->constScalingKillRateFunction);

        // Get a reference to the cells
        std::list<CellPtr>::iterator cell_it = r_cells.begin();

        // Check that apoptosis happens using the const cell killer for the first cell
        grad_cell_killer.CheckAndLabelSingleCellForApoptosis(*r_cells.begin());
        TS_ASSERT((*cell_it)->HasApoptosisBegun());

        // Check apoptosis happens for the boolean cell killer (cells 2 should be in apoptosis, cell 3 should not)
        GradientCellKiller<2> grad_cell_killer_2(&cell_population,this->boolScalingKillRateFunction);
        grad_cell_killer_2.CheckAndLabelCellsForApoptosisOrDeath();
        ++cell_it;
        TS_ASSERT(!(*cell_it)->HasApoptosisBegun());
        ++cell_it;
        TS_ASSERT((*cell_it)->HasApoptosisBegun());
        ++cell_it;
        TS_ASSERT(cell_it == r_cells.end());

        // Increment time to a time after cell death
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_DELTA(p_simulation_time->GetTime(), 1.0, 1e-3);
        TS_ASSERT_DELTA(death_time, 0.25, 1e-3);
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(death_time+1.0, 1);
        p_simulation_time->IncrementTimeOneStep();

        // Remove dead cells
        cell_population.RemoveDeadCells();
        cell_it = r_cells.begin();
        TS_ASSERT(!(*cell_it)->IsDead());

        // Check only one cell at correct location remains
        Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_it);
        c_vector<double, 2> location = p_node->rGetLocation();
        TS_ASSERT_EQUALS(initial_location_node2[1], location[1]);
        cell_it++;
        TS_ASSERT(cell_it == r_cells.end());
    }

    // void TestArchivingOfGradientCellKiller() throw (Exception)
    // {
       
    // }
};

#endif /*TESTGRADIENTCELLKILLER_HPP_*/
