#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>
#include "AbstractCellBasedTestSuite.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "NeighbourTrackingModifier.hpp"
#include "TopOfTissueTrackingModifier.hpp"

#include "StemCellProliferativeType.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

#include "FakePetscSetup.hpp"
#include "Debug.hpp"


class TestNeighbourTrackingModifier : public AbstractCellBasedTestSuite
{
public:
    void DNRTestAllClassFunctions() throw(Exception)
    {
        // Test the constructor and get methods
        NeighbourTrackingModifier<2> modifier_1 = NeighbourTrackingModifier<2>();
        NeighbourTrackingModifier<2> modifier_2 = NeighbourTrackingModifier<2>(1.0);
        TS_ASSERT_EQUALS(modifier_1.GetNeighbourhoodRadius(),1.5);
        TS_ASSERT_EQUALS(modifier_2.GetNeighbourhoodRadius(),1.0);

        // Construct the cell population
        std::vector<Node<2>*> nodes(3);
        nodes[0] = new Node<2>(0, false, 0.0, 0.0);
        nodes[1] = new Node<2>(1, false, 0.0, 0.9);
        nodes[2] = new Node<2>(2, false, 0.0, 2.3);
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_stem_type); 
        NodeBasedCellPopulation<2> cell_population(mesh,cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        // Check the neighbour allocation
        modifier_1.UpdateCellData(cell_population);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(0)->GetCellData()->GetItem("n_neighbours"),1);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(1)->GetCellData()->GetItem("n_neighbours"),2);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(2)->GetCellData()->GetItem("n_neighbours"),1);
        modifier_2.UpdateCellData(cell_population);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(0)->GetCellData()->GetItem("n_neighbours"),1);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(1)->GetCellData()->GetItem("n_neighbours"),1);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(2)->GetCellData()->GetItem("n_neighbours"),0);
    }
};


class TestTopOfTissueTrackingModifier : public AbstractCellBasedTestSuite
{
public:
    void TestSimple2dSystem() throw(Exception)
    {
        // Set up the cells, with some detached
        unsigned width = 10;
        std::vector<Node<2>*> nodes(width*width);
        for ( unsigned j=0; j < width; j++ )
        {
            for (unsigned i = 0; i < width; i++)
            {
                unsigned id = j*width + i;
                double x = i;
                double y = j;
                nodes[id] = new Node<2>( id,false, x,y);
            }
        }
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        // Run through the modifier
        TopOfTissueTrackingModifier<2> topoftissue_modifier = TopOfTissueTrackingModifier<2>();
        topoftissue_modifier.UpdateCellData(cell_population);

        // Now check the correct cells are labelled
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned cell_id = cell_population.GetLocationIndexUsingCell(*cell_iter);
            if (cell_id < (width-1)*width)
            {
                assert((*cell_iter)->GetMutationState()->IsType<WildTypeCellMutationState>());
            }
            else
            {
                assert((*cell_iter)->GetMutationState()->IsType<TopOfTissueCellMutationState>());
            }
        }
    }

    void DNRTestSimple3dSystem() throw(Exception)
    {
        // Set up the cells, with some detached
        unsigned width = 10;
        std::vector<Node<3>*> nodes(width*width*width);
        for ( unsigned k=0; k < width; k++ )
        {
            for ( unsigned j=0; j < width; j++ )
            {
                for (unsigned i = 0; i < width; i++)
                {
                    unsigned id = k*width*width + j*width + i;
                    double x = i;
                    double y = j;
                    double z = k;
                    nodes[id] = new Node<3>( id,false, x,y, z);
                }
            }
        }
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        // Run through the modifier
        TopOfTissueTrackingModifier<3> topoftissue_modifier = TopOfTissueTrackingModifier<3>();
        topoftissue_modifier.UpdateCellData(cell_population);

        // Now check the correct cells are labelled
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned cell_id = cell_population.GetLocationIndexUsingCell(*cell_iter);
            if (cell_id < (width-1)*width*width)
            {
                assert((*cell_iter)->GetMutationState()->IsType<WildTypeCellMutationState>());
            }
            else
            {
                assert((*cell_iter)->GetMutationState()->IsType<TopOfTissueCellMutationState>());
            }
        }
    }
};