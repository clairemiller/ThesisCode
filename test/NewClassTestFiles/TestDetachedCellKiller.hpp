
#ifndef TESTDETACHEDCELLKILLER_HPP_
#define TESTDETACHEDCELLKILLER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "Debug.hpp"

#include "CellsGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"

#include "DetachedCellKiller.hpp"
#include "AttachedCellMutationState.hpp"

#include <chrono>

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestDetachedCellKiller : public AbstractCellBasedTestSuite
{
public:

    void TestSimple2dSystem() throw(Exception)
    {
        // Set up the cells, with some detached
        unsigned width = 10;
        std::vector<Node<2>*> nodes(width*width);
        for ( unsigned i=0; i < width; i++ )
        {
            for (unsigned j = 0; j < width; j++)
            {
                unsigned id = j*width + i;
                double x = i;
                double y = j;
                if (j==(width-1))
                {
                    y += 2.0;
                }
                nodes[id] = new Node<2>( id,false, x,y);
            }
        }
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
        // Add the attached mutation state to the lower width cells
        for (unsigned id = 0; id < width; id++)
        {
            boost::shared_ptr<AbstractCellMutationState> p_state(new AttachedCellMutationState);
            cells[id]->SetMutationState(p_state);
        }
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);
        cell_population.Update(); // To set up node neighbours

        // Create cell killer
        DetachedCellKiller<2> killer = DetachedCellKiller<2>(&cell_population,1.0);

        // Now run the label for apoptosis method
        killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Now check the correct cells are labelled
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned cell_id = cell_population.GetLocationIndexUsingCell(*cell_iter);
            if (cell_id < (width-1)*width)
            {
                assert(!(*cell_iter)->HasApoptosisBegun());
            }
            else
            {
                assert((*cell_iter)->IsDead());
            }
        }
    }

    void TestSimple3dSystem() throw(Exception)
    {
        // Set up the cells, with some detached
        unsigned width = 10;
        std::vector<Node<3>*> nodes(width*width*width);
        for ( unsigned i=0; i < width; i++ )
        {
            for (unsigned j = 0; j < width; j++)
            {
                for (unsigned k = 0; k < width; k++)
                {
                    unsigned id = k*width*width + j*width + i;
                    double x = i;
                    double y = j;
                    double z = k;
                    if (k==(width-1))
                    {
                        z += 2.0;
                    }
                    nodes[id] = new Node<3>( id,false, x,y, z);
                }   
            }
        }

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
        // Add the attached mutation state to the lower width cells
        for (unsigned id = 0; id < width*width; id++)
        {
            boost::shared_ptr<AbstractCellMutationState> p_state(new AttachedCellMutationState);
            cells[id]->SetMutationState(p_state);
        }
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);
        cell_population.Update(); // To set up node neighbours

        // Create cell killer
        DetachedCellKiller<3> killer = DetachedCellKiller<3>(&cell_population,1.0);

        // Now run the label for apoptosis method
        killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Now check the correct cells are labelled
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned cell_id = cell_population.GetLocationIndexUsingCell(*cell_iter);
            if (cell_id < (width-1)*width*width)
            {
                assert(!(*cell_iter)->HasApoptosisBegun());
            }
            else
            {
                assert((*cell_iter)->HasApoptosisBegun());
            }
        }
    }

    void TestSpeedSetDifference() throw(Exception)
    {
        // Set up the cells, with some attached
        unsigned width = 15;
        std::vector<Node<3>*> nodes(width*width*width);
        for ( unsigned i=0; i < width; i++ )
        {
            for (unsigned j = 0; j < width; j++)
            {
                for (unsigned k = 0; k < width; k++)
                {
                    unsigned id = k*width*width + j*width + i;
                    double x = i;
                    double y = j;
                    double z = k;
                    if (k==(width-1))
                    {
                        z += 2.0;
                    }
                    nodes[id] = new Node<3>( id,false, x,y, z);
                }   
            }
        }
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
        // Add the attached mutation state to the lower width cells
        for (unsigned id = 0; id < width*width; id++)
        {
            boost::shared_ptr<AbstractCellMutationState> p_state(new AttachedCellMutationState);
            cells[id]->SetMutationState(p_state);
        }
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);
        cell_population.Update(); // To set up node neighbours
           


        // Test using current implementation (sorted vectors)
        //------------------------------------------------------------
        std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now(); 
        std::chrono::high_resolution_clock::time_point t1 = t0; 
        std::vector<unsigned> tissue_set;
        for (typename AbstractCellPopulation<3,3>::Iterator cell_iter = cell_population.Begin();
        cell_iter != cell_population.End(); ++cell_iter)
        {
            if ((*cell_iter)->GetMutationState()->template IsType<AttachedCellMutationState>())
            {
                const unsigned index = cell_population.GetLocationIndexUsingCell(*cell_iter);
                tissue_set.push_back(index);
            }
        }
        std::sort(tissue_set.begin(), tissue_set.end());
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        TRACE("Speed: current find attached cells:")
        PRINT_VARIABLE(time_span.count());

        // Iteratively find the attached sets neighbours
        t1 = std::chrono::high_resolution_clock::now(); 
        std::vector<unsigned> set_neighbours = tissue_set; // The neighbours added during the previous loop
        unsigned iterator_check = 0;
        while (!set_neighbours.empty())
        {
            std::set<unsigned> new_neighbours;
            for ( std::vector<unsigned>::iterator it = set_neighbours.begin(); it != set_neighbours.end(); it++ )
            {
                std::set<unsigned> neighbours = cell_population.GetNodesWithinNeighbourhoodRadius(*it,1.5);
                std::vector<unsigned> added_indices;
                // The set additions
                std::set_difference(neighbours.begin(), neighbours.end(), // Returns neighbours in this set...
                                tissue_set.begin(), tissue_set.end(), //... that are not in this set
                                std::back_inserter(added_indices));
                // Merge with the other new editions
                new_neighbours.insert(added_indices.begin(), added_indices.end());
            }
            // Merge the new neighbours into the full set
            tissue_set.insert(tissue_set.end(),new_neighbours.begin(), new_neighbours.end());
            std::sort(tissue_set.begin(), tissue_set.end());
            set_neighbours.resize(new_neighbours.size());
            std::copy(new_neighbours.begin(), new_neighbours.end(), set_neighbours.begin());
            iterator_check++;
            assert(iterator_check < cell_population.GetNumNodes());
        }
        t2 = std::chrono::high_resolution_clock::now(); 
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        TRACE("Speed: current find the tissue set");
        PRINT_VARIABLE(time_span.count());

        t1 = std::chrono::high_resolution_clock::now(); 
        // Finally, determine those cells outside of the tissue
        std::vector<unsigned> all_indices = (cell_population.rGetMesh()).GetAllNodeIndices();
        std::sort(all_indices.begin(), all_indices.end());
        std::vector<unsigned> detached_cells;
        std::set_difference(all_indices.begin(), all_indices.end(), 
            tissue_set.begin(), tissue_set.end(), 
            std::back_inserter(detached_cells));
        t2 = std::chrono::high_resolution_clock::now(); 
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        TRACE("Speed: current find the detached cells")
        PRINT_VARIABLE(time_span.count());

        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t0);
        TRACE("Speed: total using current implementation");
        PRINT_VARIABLE(time_span.count());
        TRACE("--------------------------------------------");


        // New test case
        //------------------------------------------------------------
        t0 = std::chrono::high_resolution_clock::now(); 
        t1 = t0; 
        std::vector<unsigned> tissue_set2;
        for (typename AbstractCellPopulation<3,3>::Iterator cell_iter = cell_population.Begin();
        cell_iter != cell_population.End(); ++cell_iter)
        {
            if ((*cell_iter)->GetMutationState()->template IsType<AttachedCellMutationState>())
            {
                const unsigned index = cell_population.GetLocationIndexUsingCell(*cell_iter);
                tissue_set2.push_back(index);
            }
        }
        std::sort(tissue_set2.begin(), tissue_set2.end());
        t2 = std::chrono::high_resolution_clock::now(); 
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        TRACE("Speed: new find attached cells:")
        PRINT_VARIABLE(time_span.count());

        // Iteratively find the attached sets neighbours
        t1 = std::chrono::high_resolution_clock::now(); 
        std::vector<unsigned> set_neighbours2 = tissue_set2; // The neighbours added during the previous loop
        iterator_check = 0;
        while (!set_neighbours2.empty())
        {
            std::vector<unsigned> new_neighbours;
            for ( std::vector<unsigned>::iterator it = set_neighbours2.begin(); it != set_neighbours2.end(); it++ )
            {
                std::set<unsigned> neighbours = cell_population.GetNodesWithinNeighbourhoodRadius(*it,1.5);
                std::vector<unsigned> added_indices;
                // The set additions
                std::set_difference(neighbours.begin(), neighbours.end(), // Returns neighbours in this set...
                                tissue_set2.begin(), tissue_set2.end(), //... that are not in this set
                                std::back_inserter(added_indices));
                // Merge with the other new editions
                std::vector<unsigned> tmp;
                std::set_union(new_neighbours.begin(), new_neighbours.end(),
                                added_indices.begin(), added_indices.end(),
                                std::back_inserter(tmp));
                new_neighbours = tmp;
            }
            // Merge the new neighbours into the full set
            tissue_set2.insert(tissue_set2.end(),new_neighbours.begin(), new_neighbours.end());
            std::sort(tissue_set2.begin(), tissue_set2.end());
            set_neighbours2.resize(new_neighbours.size());
            std::copy(new_neighbours.begin(), new_neighbours.end(), set_neighbours2.begin());
            iterator_check++;
            assert(iterator_check < cell_population.GetNumNodes());
        }
        t2 = std::chrono::high_resolution_clock::now(); 
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        TRACE("Speed: new find the tissue set");
        PRINT_VARIABLE(time_span.count());

        t1 = std::chrono::high_resolution_clock::now(); 
        // Finally, determine those cells outside of the tissue
        std::vector<unsigned> all_indices2 = (cell_population.rGetMesh()).GetAllNodeIndices();
        std::sort(all_indices2.begin(), all_indices2.end());
        std::vector<unsigned> detached_cells2;
        std::set_difference(all_indices2.begin(), all_indices2.end(), 
            tissue_set2.begin(), tissue_set2.end(), 
            std::back_inserter(detached_cells2));
        t2 = std::chrono::high_resolution_clock::now(); 
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        TRACE("Speed: new find the detached cells")
        PRINT_VARIABLE(time_span.count());

        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t0);
        TRACE("Speed: total using new implementation");
        PRINT_VARIABLE(time_span.count());
        
    }

    void TestIterationSkipping() throw(Exception)
    {
        // Set up the cells, with some detached
        unsigned width = 10;
        std::vector<Node<2>*> nodes(width*width);
        for ( unsigned i=0; i < width; i++ )
        {
            for (unsigned j = 0; j < width; j++)
            {
                unsigned id = j*width + i;
                double x = i;
                double y = j;
                if (j==(width-1))
                {
                    y += 2.0;
                }
                nodes[id] = new Node<2>( id,false, x,y);
            }
        }
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
        // Add the attached mutation state to the lower width cells
        for (unsigned id = 0; id < width; id++)
        {
            boost::shared_ptr<AbstractCellMutationState> p_state(new AttachedCellMutationState);
            cells[id]->SetMutationState(p_state);
        }
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);
        cell_population.Update(); // To set up node neighbours

        // Create cell killer
        DetachedCellKiller<2> killer = DetachedCellKiller<2>(&cell_population,1.0);

        // Now run the label for apoptosis method, should run the first time (checked in previous tests)
        killer.CheckAndLabelCellsForApoptosisOrDeath();

        // Now move the next layer of cells up
        for (unsigned id = (width-2)*width; id < (width-1)*width; id++)
        {
            c_vector<double,2>& cell_loc = mesh.GetNode(id)->rGetModifiableLocation();
            cell_loc[1] += 2.0;
        }

        // The following iterations should not label anything
        unsigned n_skipped = killer.GetIterationsToSkip();
        for (unsigned i = 0; i < n_skipped; i++)
        {
            assert(killer.mIterationsSkipped == i);
            killer.CheckAndLabelCellsForApoptosisOrDeath();
            // Check that nothing has been labelled
            bool cellsStillInPopulation = false;
            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
            {
                unsigned cell_id = cell_population.GetLocationIndexUsingCell(*cell_iter);
                if (cell_id < (width-1)*width)
                {
                    assert(!(*cell_iter)->IsDead());
                    if (cell_id > (width-2)*width)
                    {
                        cellsStillInPopulation = true; // Cells could be removed from iterator if killed
                    }
                }
            }
            assert(cellsStillInPopulation);
        }
        
        // Now one more should label the next row
        killer.CheckAndLabelCellsForApoptosisOrDeath();
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            unsigned cell_id = cell_population.GetLocationIndexUsingCell(*cell_iter);
            if (cell_id < (width-2)*width)
            {
                assert(!(*cell_iter)->IsDead());
            }
            else
            {
                assert((*cell_iter)->IsDead());
            }
        }
    }
};

#endif /*TESTDETACHEDCELLKILLERHPP_*/
