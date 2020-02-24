
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

#include "DetachedCellKillerWithWriter.hpp"
#include "CellAgeAtDeathWriter.hpp"
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

        // Create the writer
        boost::shared_ptr<CellAgeAtDeathWriter<2,2> > p_cell_writer(new CellAgeAtDeathWriter<2,2>());
        cell_population.AddPopulationWriter(p_cell_writer);

        // Create cell killer
        DetachedCellKillerWithWriter<2> killer = DetachedCellKillerWithWriter<2>(&cell_population,1.0,p_cell_writer);

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

        // Check that these are stored in the cell writer
        PRINT_VECTOR(p_cell_writer->GetDeadCellIDs());
        PRINT_VECTOR(p_cell_writer->GetDeadCellAges());
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

        boost::shared_ptr<CellAgeAtDeathWriter<3,3> > p_cell_writer(new CellAgeAtDeathWriter<3,3>());
        cell_population.AddPopulationWriter(p_cell_writer);

        // Create cell killer
        DetachedCellKillerWithWriter<3> killer = DetachedCellKillerWithWriter<3>(&cell_population,1.0,p_cell_writer);

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

        PRINT_VECTOR(p_cell_writer->GetDeadCellIDs());
    }

      void TestArchiving()
    {
        // Set up singleton classes
        OutputFileHandler handler("archive", false);    // don't erase contents of folder
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "detached_withwriter_killer.arch";

        {
            std::vector<Node<2>*> nodes(0);
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes,1.5);
            std::vector<CellPtr> cells(0);
            NodeBasedCellPopulation<2> cell_population(mesh, cells);

            // Create an output archive
            boost::shared_ptr<CellAgeAtDeathWriter<2,2> > p_cell_writer(new CellAgeAtDeathWriter<2,2>());

            // Create cell killer
            DetachedCellKillerWithWriter<2> cell_killer = DetachedCellKillerWithWriter<2>(&cell_population,1.0,p_cell_writer);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            DetachedCellKillerWithWriter<2>* const p_cell_killer = &cell_killer;
            output_arch << p_cell_killer;

            //TS_ASSERT_DELTA(p_cell_killer->GetDeathProbabilityInAnHour(), 0.134, 1e-9);
       }

       {
            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            DetachedCellKillerWithWriter<2>* p_cell_killer;

            // Restore from the archive
            input_arch >> p_cell_killer;

            // Test we have restored the probability correctly
            //TS_ASSERT_DELTA(p_cell_killer->GetDeathProbabilityInAnHour(), 0.134, 1e-9);

            delete p_cell_killer;
        }
    }
};

#endif /*TESTDETACHEDCELLKILLERHPP_*/
