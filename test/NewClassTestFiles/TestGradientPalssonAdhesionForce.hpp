
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

#include "ExponentialDecayPalssonAdhesionForce.hpp"
#include "Debug.hpp"

#include "FakePetscSetup.hpp"

class TestGradientPalssonAdhesionForce : public AbstractCellBasedTestSuite
{
public:
    void TestExponentialDecayPalssonAdhesionForce() throw(Exception)
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
        MAKE_PTR(DifferentiatedCellProliferativeType, p_stem_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_stem_type); 

        // Create the population
        NodeBasedCellPopulation<2> cell_population(mesh,cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        // Add the adhesion force
        // Set up simulator
        OffLatticeSimulation<2> simulator(cell_population);

        // Add the cell-cell force
        MAKE_PTR_ARGS( ExponentialDecayPalssonAdhesionForce<2>, p_cell_force, (1.0,10.0,0.5) );
        simulator.AddForce(p_cell_force);

        // Test the member scaling function
        CellPtr pCell0 = *(cell_population.Begin());
        pCell0->SetBirthTime(-10.0);
        double scale_factor = p_cell_force->ScalingFunction(pCell0);
        TS_ASSERT_EQUALS(scale_factor, 1.0);

        pCell0->SetBirthTime(-11.0);
        scale_factor = p_cell_force->ScalingFunction(pCell0);
        TS_ASSERT_EQUALS(scale_factor,std::exp(-0.5));

        // Test the force calulation
        CellPtr pCell1 = *(++cell_population.Begin());
        pCell1->SetBirthTime(-1.0);
        c_vector<double,2> force = p_cell_force->CalculateForceBetweenNodes(0,1,cell_population);
        TS_ASSERT_DELTA(force[0],0.0,1.0e-6);
        TS_ASSERT_DELTA(force[1],0.5*(1.0+std::exp(-0.5)),1.0e-3);
    }
};