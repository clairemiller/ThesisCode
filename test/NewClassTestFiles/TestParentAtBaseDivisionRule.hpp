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

#include "ParentAtBaseDivisionRule.hpp"



//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"
class TestParentAtBaseDivisionRule : public AbstractCellBasedTestSuite
{
public:
    void TestCalculateDivisionVector() throw(Exception)
    {
        // Set the node positions
        std::vector<Node<2>*> nodes(2);
        nodes[0] = new Node<2>( 0, false, 0.0, 1.0 );
        nodes[1] = new Node<2>( 1, false, 10.0, 5.0);

        // Construct the mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create the cells
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_stem_type); 

        // Create the population
        NodeBasedCellPopulation<2> cell_population(mesh,cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        // Add the division rule
        MAKE_PTR_ARGS(ParentAtBaseDivisionRule<2>, p_div_rule, (1));
        cell_population.SetCentreBasedDivisionRule(p_div_rule);

        // Check the get function
        TS_ASSERT_EQUALS(p_div_rule->GetVerticalDirection(),(unsigned) 1);

        // Check the calculate division vector function
        std::pair<c_vector<double, 2>, c_vector<double, 2> > div_vec;
        for (int i = 0; i < 10; i++)
        {
            div_vec = p_div_rule->CalculateCellDivisionVector(cell_population.GetCellUsingLocationIndex(0),cell_population);
            TS_ASSERT( div_vec.first[1] < div_vec.second[1] );

            div_vec = p_div_rule->CalculateCellDivisionVector(cell_population.GetCellUsingLocationIndex(1),cell_population);
            TS_ASSERT( div_vec.first[1] < div_vec.second[1] );
        }
        
    }
};