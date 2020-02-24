
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>
#include "AbstractCellBasedTestSuite.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"

#include "StochasticExponentialDecayPalssonAdhesionForce.hpp"
#include "ExponentialDecayPalssonAdhesionForce.hpp"
#include "Debug.hpp"

#include "FakePetscSetup.hpp"

class TestStochasticDecayPalssonAdhesionForce : public AbstractCellBasedTestSuite
{
public:
    void TestStochasticExponentialDecayPalssonAdhesionForce() throw(Exception)
    {
        // Set the node positions
        std::vector<Node<2>*> nodes(2);
        nodes[0] = new Node<2>( 0, false, 0.0, 0.0 );
        nodes[1] = new Node<2>( 1, false, 0.0, 1.15 );

        // Construct the mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create the cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_diff_type); 

        // Create the population
        NodeBasedCellPopulation<2> cell_population(mesh,cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        // Add the adhesion force
        // Set up simulator
        OffLatticeSimulation<2> simulator(cell_population);

        // Add the cell-cell force
        MAKE_PTR_ARGS( StochasticExponentialDecayPalssonAdhesionForce<2>, p_cell_force, (1.0,1.0,1.0,2.0) );
        simulator.AddForce(p_cell_force);

        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Need to run CalculateForceBetweenNodes to initialise
        c_vector<double, 2> force = p_cell_force->CalculateForceBetweenNodes(0,1,cell_population);

        // // Get the decay rates
        NodeBasedCellPopulation<2>::Iterator cell_it = cell_population.Begin();
        CellPtr pCell0 = *(cell_it);
        ++cell_it;
        CellPtr pCell1 = *(cell_it);
        double decayrate0 = pCell0->GetCellData()->GetItem("AdhesionDecayRate");
        double decayrate1 = pCell1->GetCellData()->GetItem("AdhesionDecayRate");

        // Check they are within the range and not equal
        TS_ASSERT_LESS_THAN_EQUALS(decayrate0,2.0);
        TS_ASSERT_LESS_THAN_EQUALS(decayrate1,2.0);
        TS_ASSERT_LESS_THAN_EQUALS(1.0,decayrate0);
        TS_ASSERT_LESS_THAN_EQUALS(1.0,decayrate1);
        TS_ASSERT_DIFFERS(decayrate0,decayrate1);

        // Test the member scaling function
        pCell0->SetBirthTime(-2.0);
        double scale_factor = p_cell_force->ScalingFunction(pCell0);
        double scale_factor_comp = std::exp(-decayrate0);
        TS_ASSERT_EQUALS(scale_factor, scale_factor_comp);

        pCell0->SetBirthTime(-5.0);
        scale_factor = p_cell_force->ScalingFunction(pCell0);
        TS_ASSERT_EQUALS(scale_factor,std::exp(-decayrate0*4.0));

        // Test the force calculation
        pCell1->SetBirthTime(-2.0);
        double scale_factor_cell1 = std::exp(-decayrate1);
        force = p_cell_force->CalculateForceBetweenNodes(0,1,cell_population);
        TS_ASSERT_DELTA(0.5*(scale_factor+scale_factor_cell1),force[1],1.0e-3);
        TS_ASSERT_DELTA(0,force[0],1.0e-6);
    }

    void TestComparisonToExponential() throw(Exception)
    {
        // Set the node positions
        double y_pos = 1.1;
        std::vector<Node<2>*> nodes(2);
        nodes[0] = new Node<2>( 0, false, 1.0, 0.0 );
        nodes[1] = new Node<2>( 1, false, 1.03, y_pos );

        std::vector<Node<2>*> nodes_comp(2);
        nodes_comp[0] = new Node<2>( 0, false, 1.0, 0.0 );
        nodes_comp[1] = new Node<2>( 1, false, 1.03, y_pos );
        
        // Construct the mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        NodesOnlyMesh<2> mesh_comp;
        mesh_comp.ConstructNodesWithoutMesh(nodes_comp,1.5);

        // Create the cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_diff_type); 

        std::vector<CellPtr> cells_comp;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type_comp);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator_comp;
        cells_generator_comp.GenerateBasicRandom(cells_comp,mesh_comp.GetNumNodes(), p_diff_type_comp); 

        // Create the population
        NodeBasedCellPopulation<2> cell_population(mesh,cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        NodeBasedCellPopulation<2> cell_population_comp(mesh_comp,cells_comp);
        cell_population_comp.SetAbsoluteMovementThreshold(1.5);

        // Add the adhesion force
        // Set up simulator
        OffLatticeSimulation<2> simulator(cell_population);

        OffLatticeSimulation<2> simulator_comp(cell_population_comp);

        // Add the cell-cell force
        MAKE_PTR_ARGS( StochasticExponentialDecayPalssonAdhesionForce<2>, p_cell_force, (1.0,1.0,1.5,1.5) );
        simulator.AddForce(p_cell_force);

        MAKE_PTR_ARGS( ExponentialDecayPalssonAdhesionForce<2>, p_cell_force_comp, (1.0,1.0,1.5) );
        simulator_comp.AddForce(p_cell_force_comp);

        // Set up the simulation time
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Need to run CalculateForceBetweenNodes to initialise
        c_vector<double, 2> force = p_cell_force->CalculateForceBetweenNodes(0,1,cell_population);
        c_vector<double, 2> force_comp = p_cell_force_comp->CalculateForceBetweenNodes(0,1,cell_population_comp);

        // Get the decay rates for the stochastic cells
        NodeBasedCellPopulation<2>::Iterator cell_it = cell_population.Begin();
        CellPtr pCell0 = *(cell_it);
        ++cell_it;
        CellPtr pCell1 = *(cell_it);
        double decayrate0 = pCell0->GetCellData()->GetItem("AdhesionDecayRate");
        double decayrate1 = pCell1->GetCellData()->GetItem("AdhesionDecayRate");

        // Check they are within the range and not equal
        TS_ASSERT_EQUALS(decayrate0,1.5);
        TS_ASSERT_EQUALS(decayrate1,1.5);

        // Test the member scaling function
        NodeBasedCellPopulation<2>::Iterator cell_it_comp = cell_population_comp.Begin();
        CellPtr pCell0_comp = *(cell_it_comp);
        ++cell_it_comp;
        CellPtr pCell1_comp = *(cell_it_comp);
        pCell0->SetBirthTime(-4.0);
        pCell0_comp->SetBirthTime(-4.0);
        double scale_factor = p_cell_force->ScalingFunction(pCell0);
        double scale_factor_comp = p_cell_force_comp->ScalingFunction(pCell0_comp);
        TS_ASSERT_EQUALS(scale_factor, scale_factor_comp);

        pCell1->SetBirthTime(-3.0);
        pCell1_comp->SetBirthTime(-3.0);
        scale_factor = p_cell_force->ScalingFunction(pCell1);
        scale_factor_comp = p_cell_force_comp->ScalingFunction(pCell1_comp);
        TS_ASSERT_EQUALS(scale_factor,scale_factor_comp);

        // Test the force calculation
        pCell1->SetBirthTime(-2.0);
        pCell1_comp->SetBirthTime(-2.0);
        force = p_cell_force->CalculateForceBetweenNodes(0,1,cell_population);
        force_comp = p_cell_force_comp->CalculateForceBetweenNodes(0,1,cell_population_comp);
        TS_ASSERT_EQUALS(force[0],force_comp[0]);
        TS_ASSERT_EQUALS(force[1],force_comp[1]);
    }

};