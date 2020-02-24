

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>
#include "AbstractCellBasedTestSuite.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"

#include "PalssonAdhesionForce.hpp"
#include "LowNeighbourForce.hpp"
#include "RotationalDivisionForce.hpp"
#include "Debug.hpp"

#include "FakePetscSetup.hpp"


class TestPalssonAdhesionForce : public AbstractCellBasedTestSuite
{
public:
    void DNRTestMaxForceCalculation() throw(Exception)
    {
        PalssonAdhesionForce<2> force = PalssonAdhesionForce<2>(1.0);

        TS_ASSERT_EQUALS((int)(force.CalculatePalssonAdhesionForce(0.153)*10000),267);

    }
    void DNRTestAllClassFunctions() throw(Exception)
    {
        // DNRTest constructor
        PalssonAdhesionForce<2> force = PalssonAdhesionForce<2>(5.5);
        TS_ASSERT_EQUALS(force.GetPeakForce(),5.5);

        // Set the node positions
        std::vector<Node<2>*> nodes(3);
        nodes[0] = new Node<2>( 0, false, 0.0, 1.0 );
        nodes[1] = new Node<2>( 1, false, 0.0, 2.1 );
        nodes[2] = new Node<2>( 2, false, 1.0, 2.1 );

        // Construct the mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5); // Want all nodes on process 0

        // Create the cells
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_stem_type); 

        // Create the population
        NodeBasedCellPopulation<2> cell_population(mesh,cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);
        
        // Check the force calculation
        c_vector<double, 2> calculated_force_01 = force.CalculateForceBetweenNodes(0,1,cell_population);
        double c1 = std::sqrt(1.0/(2.0*7.0));
        double c2 = c1*std::exp(-7.0*c1*c1);
        double b01 = 0.1/0.5;
        double F01 = -5.5*( (b01 + c1)*std::exp(-7.0*std::pow(b01+c1,2.0)) - c2*std::exp(-7.0*b01*b01) );
        TS_ASSERT_DELTA(calculated_force_01[1],F01,1.0e-5);
        TS_ASSERT_DELTA(calculated_force_01[0],0.0,1.0e-5);
        c_vector<double, 2> calculated_force_12 = force.CalculateForceBetweenNodes(1,2,cell_population);
        TS_ASSERT_DELTA(calculated_force_12[0],0.0,1e-10);
        TS_ASSERT_DELTA(calculated_force_12[1],0.0,1e-10);
        c_vector<double, 2> calculated_force_02 = force.CalculateForceBetweenNodes(0,2,cell_population);
        double b02 = (std::sqrt(2.21) - 1.0)/0.5;
        double F02_magnitude = -5.5*( (b02 + c1)*std::exp(-7.0*std::pow(b02+c1,2.0)) - c2*std::exp(-7.0*b02*b02) );
        TS_ASSERT_DELTA(calculated_force_02[0],F02_magnitude*1.0/sqrt(2.21),1e-10);
        TS_ASSERT_DELTA(calculated_force_02[1],F02_magnitude*1.1/sqrt(2.21),1e-10);

        // DNRTest add force contribution method
        // First need to run update to calculate the node pairs
        cell_population.Update(false);
        force.AddForceContribution(cell_population);
        TS_ASSERT_EQUALS(cell_population.GetNode(0)->rGetAppliedForce()[0],calculated_force_01[0] + calculated_force_02[0]);
        TS_ASSERT_EQUALS(cell_population.GetNode(0)->rGetAppliedForce()[1],calculated_force_01[1] + calculated_force_02[1]);
        TS_ASSERT_EQUALS(cell_population.GetNode(1)->rGetAppliedForce()[0],-1.0*calculated_force_01[0]);
        TS_ASSERT_EQUALS(cell_population.GetNode(1)->rGetAppliedForce()[1],-1.0*calculated_force_01[1]);
        TS_ASSERT_EQUALS(cell_population.GetNode(2)->rGetAppliedForce()[0],-1.0*calculated_force_02[0]);
        TS_ASSERT_EQUALS(cell_population.GetNode(2)->rGetAppliedForce()[1],-1.0*calculated_force_02[1]);
    }
};



class TestLowNeighbourForce : public AbstractCellBasedTestSuite
{
    public:
    void DNRTestAllClassFunctions() throw(Exception)
    {
        // DNRTest constructor
        c_vector<double,2> vec_2d_force = scalar_vector<double>(2,0.0);
        vec_2d_force[1] = 1.0;
        LowNeighbourForce<2> force_2d(vec_2d_force,-0.1);
        c_vector<double,3> vec_3d_force = scalar_vector<double>(3,0.0);
        LowNeighbourForce<3> force_3d(vec_3d_force,-0.1);

        TS_ASSERT(force_2d.GetNeighbourThreshold() == 4);
        TS_ASSERT(force_3d.GetNeighbourThreshold() == 10);

        // DNRTest the force
        //----------------
        // Set the node positions
        int nNodes = 3;
        std::vector<Node<2>*> nodes(nNodes);
        for ( int i = 0; i < nNodes; i++ )
        {
            nodes[i] = new Node<2>( i, false, i, i );
        }

        // Construct the mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create the cell population
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_stem_type);
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        // Assign numbers of neighbours to each cell
        int i = 0;
        for (NodeBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin(); cell_iter != cell_population.End(); ++cell_iter, i++ )
        {
            cell_iter->GetCellData()->SetItem("n_neighbours",i*4);
        }

        // Run the force contribution calc
        force_2d.AddForceContribution(cell_population);

        // Check that only the first 2 cells experience the force
        TS_ASSERT_EQUALS( cell_population.GetNode(0)->rGetAppliedForce()[1], 1.0 );
        TS_ASSERT_EQUALS( cell_population.GetNode(1)->rGetAppliedForce()[1], 1.0 );
        TS_ASSERT_EQUALS( cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0 );
    }
};
