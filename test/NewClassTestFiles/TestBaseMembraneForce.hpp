


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>

// To run in parallel
//#include "PetscSetupAndFinalize.hpp"
// When run in serial
#include "FakePetscSetup.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "BaseMembraneForce.hpp"


class Test2DSimulation : public AbstractCellBasedTestSuite
{
public:
    void DNRTestForceCalculation() throw(Exception)
    {
        // Set the node positions
        unsigned n = 10;
        std::vector<Node<2>*> nodes(n);
        for ( unsigned i = 0; i < n; i++)
        {
            nodes[i] = new Node<2>( i, false, i*3, -0.5 + (0.1*i) );
        }

        // Construct the mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create the cells
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_stem_type); 

        // Force birth time to be 0
        for ( std::vector<CellPtr>::iterator p_cell = cells.begin(); p_cell != cells.end(); p_cell++ )
        {
            (*p_cell)->SetBirthTime(0.0);
            // Set the duration high to avoid division during simulation
            dynamic_cast<FixedG1GenerationalCellCycleModel*>( (*p_cell)->GetCellCycleModel())->SetSDuration(100);
        }

        // Create the population
        NodeBasedCellPopulation<2> cell_population(mesh,cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        // Set up the simulator
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestStablePoint");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(1.0/120.0);

        // Add the membrane adhesion/repulsion
        MAKE_PTR_ARGS(BaseMembraneForce<2>,p_base_force,(500.0,500.0,0.0,15.0,1.5));
        simulator.AddForce(p_base_force);

        // Run simulation
        simulator.Solve();

        // Check the force that is calculated
        printf("node,y0,force,\n");
        for ( unsigned i=0;i < n; i++)
        {
            double force = simulator.rGetCellPopulation().GetNode(i)->rGetAppliedForce()[1];
            printf("%i,%f,%f,\n",i,-0.5 + (0.1*i),force);
        }
    }

	void DNRTestStablePoint() throw(Exception)
	{
		// We want a setup with a stem, transit and differentiated cell all far enough away they feel no effects from each other.
		int n=3;

		// Set the node positions
		std::vector<Node<2>*> nodes(n);
		for ( int i = 0; i < n; i++)
		{
			nodes[i] = new Node<2>( i, false, i*3, -0.5 );
		}

		// Construct the mesh
		NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes,1.5);

		// Create the cells
    	std::vector<CellPtr> cells;
    	MAKE_PTR(StemCellProliferativeType, p_stem_type);
    	CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
    	cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_stem_type); 

    	// Assign the other 2 proliferative types
    	MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        cells[1]->SetCellProliferativeType(p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        cells[2]->SetCellProliferativeType(p_diff_type);

        // Force birth time to be 0
        for ( std::vector<CellPtr>::iterator p_cell = cells.begin(); p_cell != cells.end(); p_cell++ )
        {
        	(*p_cell)->SetBirthTime(0.0);
            // Set the duration high to avoid division during simulation
            dynamic_cast<FixedG1GenerationalCellCycleModel*>( (*p_cell)->GetCellCycleModel())->SetSDuration(100);
        }

        // Create the population
        NodeBasedCellPopulation<2> cell_population(mesh,cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        // Set up the simulator
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestStablePoint");
        simulator.SetSamplingTimestepMultiple(24);
        simulator.SetEndTime(12);

        // Add the membrane adhesion/repulsion
        MAKE_PTR_ARGS(BaseMembraneForce<2>,p_base_force,(500.0,500.0,0.0,15.0,1.5));
        //MAKE_PTR_ARGS(BaseMembraneForce<2>,p_base_force,(0.0,0.0,0.0,15.0,1.5));
        simulator.AddForce(p_base_force);

        // Run simulation
        simulator.Solve();

        // Get the cell locations and check they are as expected
        double y = simulator.rGetCellPopulation().GetNode(0)->rGetLocation()[1];
        TS_ASSERT_DELTA( y, 0.0, 5e-3 );
        y = simulator.rGetCellPopulation().GetNode(1)->rGetLocation()[1];
        TS_ASSERT_DELTA( y, 0.0, 5e-3 );
        y = simulator.rGetCellPopulation().GetNode(2)->rGetLocation()[1];
        TS_ASSERT_DELTA( y, 0.3, 5e-3 );
    }
    
    void TestSaveLoadBaseMembraneForce() throw(Exception)
    {
                // Set the node positions
                unsigned n = 10;
                std::vector<Node<2>*> nodes(n);
                for ( unsigned i = 0; i < n; i++)
                {
                    nodes[i] = new Node<2>( i, false, i*3, -0.5 + (0.1*i) );
                }
        
                // Construct the mesh
                NodesOnlyMesh<2> mesh;
                mesh.ConstructNodesWithoutMesh(nodes,1.5);
        
                // Create the cells
                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, p_stem_type);
                CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_stem_type); 
        
                // Force birth time to be 0
                for ( std::vector<CellPtr>::iterator p_cell = cells.begin(); p_cell != cells.end(); p_cell++ )
                {
                    (*p_cell)->SetBirthTime(0.0);
                    // Set the duration high to avoid division during simulation
                    dynamic_cast<FixedG1GenerationalCellCycleModel*>( (*p_cell)->GetCellCycleModel())->SetSDuration(100);
                }
        
                // Create the population
                NodeBasedCellPopulation<2> cell_population(mesh,cells);
                cell_population.SetAbsoluteMovementThreshold(1.5);
        
                // Set up the simulator
                OffLatticeSimulation<2> simulator(cell_population);
                simulator.SetOutputDirectory("TestSaveLoadBaseMembraneForce");
                simulator.SetDt(1.0/120.0);
                simulator.SetSamplingTimestepMultiple(1);
                simulator.SetEndTime(1.0/120.0);
        
                // Add the membrane adhesion/repulsion
                MAKE_PTR_ARGS(BaseMembraneForce<2>,p_base_force,(500.0,500.0,0.0,15.0,1.5));
                simulator.AddForce(p_base_force);
        
                // Run simulation
                simulator.Solve();

                // Save simulation
                CellBasedSimulationArchiver<2,OffLatticeSimulation<2>, 2>::Save(&simulator);

                // Reload
                OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2>,2 >::Load("TestSaveLoadBaseMembraneForce", 1.0/120.0);
                delete p_simulator;
    }
};