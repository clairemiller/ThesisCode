// Test run
//#include "ScratchAssayCellKiller.hpp"

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "HoneycombMeshGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "SmartPointers.hpp"
#include "Cylindrical2dNodesOnlyMesh.hpp"
#include "UniformlyDistributedCellCycleModel.hpp"
#include "PalssonAdhesionForce.hpp"
#include "RepulsionForce.hpp"
#include "PlaneBoundaryCondition.hpp"

#include "NNeighbourInhibitedStochasticCellCycleModel.hpp"
#include "NNeighbourTrackingModifier.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellDivisionPopulationWriter.hpp"

#include "CellProliferativePhasesWriter.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"
class TestCreatingAndUsingANewCellKillerTutorial : public AbstractCellBasedTestSuite
{
public:

    // void DNRTestMyCellKiller() throw(Exception)
    // {

    // 	double killer_width = 8.0;
    // 	double killer_centre = 10.0;

    //     HoneycombMeshGenerator generator(20, 20, 0);
    //     MutableMesh<2,2>* p_mesh = generator.GetMesh();

    //     std::vector<CellPtr> cells;
    //     CellsGenerator<FixedDurationGenerationBasedFixedG1GenerationalCellCycleModel, 2> cells_generator;
    //     cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

    //     MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

    //     ColumnCellKiller<2> my_cell_killer(&cell_population, killer_width, killer_centre);

    //     my_cell_killer.CheckAndLabelCellsForApoptosisOrDeath();


    //     double xMin = killer_centre-(killer_width/2.0);
    //     double xMax = killer_centre+(killer_width/2.0);
    //     for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
    //          cell_iter != cell_population.End();
    //          ++cell_iter)
    //     {
    //         double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];

    //         if ( x >= xMin && x <= xMax )
    //         {
    //             TS_ASSERT_EQUALS(cell_iter->IsDead(), true);
    //         }
    //         else
    //         {
    //             TS_ASSERT_EQUALS(cell_iter->IsDead(), false);
    //         }
    //     }

    //     cell_population.RemoveDeadCells();

    //     for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
    //          cell_iter != cell_population.End();
    //          ++cell_iter)
    //     {
    //         double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
    //         double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

    //         TS_ASSERT_LESS_THAN_EQUALS(pow(x/20, 2) + pow(y/10, 2) > 1.0, 1.0);
    //     }

    //     OutputFileHandler handler("archive", false);
    //     std::string archive_filename = handler.GetOutputDirectoryFullPath() + "my_cell_killer.arch";

    //     {
    //         AbstractCellKiller<2>* const p_cell_killer = new ColumnCellKiller<2>(NULL,0.0,0.0);

    //         std::ofstream ofs(archive_filename.c_str());
    //         boost::archive::text_oarchive output_arch(ofs);

    //         output_arch << p_cell_killer;
    //         delete p_cell_killer;
    //     }

    //     {
    //         std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
    //         boost::archive::text_iarchive input_arch(ifs);

    //         AbstractCellKiller<2>* p_cell_killer;

    //         input_arch >> p_cell_killer;
    //         delete p_cell_killer;
    //     }
    // }

    // void DNRTestOffLatticeSimulationWithMyCellKiller() throw(Exception)
    // {
    //     HoneycombMeshGenerator generator(20, 20, 0);
    //     MutableMesh<2,2>* p_mesh = generator.GetMesh();

    //     std::vector<CellPtr> cells;
    //     CellsGenerator<FixedDurationGenerationBasedFixedG1GenerationalCellCycleModel, 2> cells_generator;
    //     cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

    //     MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

    //     MAKE_PTR_ARGS(ColumnCellKiller<2>, p_killer, (&cell_population,8.0,10.0));

    //     OffLatticeSimulation<2> simulator(cell_population);
    //     simulator.SetOutputDirectory("TestOffLatticeSimulationWithMyCellKiller");
    //     simulator.SetEndTime(1.0);

    //     MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
    //     p_linear_force->SetCutOffLength(3);
    //     simulator.AddForce(p_linear_force);

    //     simulator.AddCellKiller(p_killer);

    //     simulator.Solve();
    // }

    void TestSimple() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        unsigned nInitialCells = 4;
        // Set up the node positions
        std::vector<Node<2>*> nodes(nInitialCells);
        // Initialise the node locations using a random number generator
        nodes[0] = new Node<2>(0, false, 3.0,0.0);
        nodes[1] = new Node<2>(1,false,4.0,0.0);
        nodes[2] = new Node<2>(2,false,3.5,1.0);
        nodes[3] = new Node<2>(3,false,4.5,1.0);

        // Create the mesh
        NodesOnlyMesh<2> mesh; 
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create the cells
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<NNeighbourInhibitedStochasticCellCycleModel, 2> cells_generator;
        std::vector<CellPtr> cells;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);

        // Modify the cell cycle to set the division probabilities
        for ( std::vector<CellPtr>::iterator iter = cells.begin(); iter != cells.end(); ++iter )
        {
            dynamic_cast<NNeighbourInhibitedStochasticCellCycleModel*>((*iter)->GetCellCycleModel())->SetDivisionProbability(0.6);
        }

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        // Add the output of the cell ancestors

        // Set up simulator
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Test");
        simulator.SetSamplingTimestepMultiple(24);
        simulator.SetEndTime(100);

        // Output the division information
        simulator.SetOutputDivisionLocations(true);

        // Add the extra adhesion to base membrane force for stem and transit
        MAKE_PTR_ARGS( PalssonAdhesionForce<2>, p_base_force, (500.0, 500.0, 0.2, 2.0) );
        //p_base_force->ReduceForceWithGeneration(true);
        simulator.AddForce(p_base_force);

        MAKE_PTR(RepulsionForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR(NNeighbourTrackingModifier<2>,p_modifier);
        simulator.AddSimulationModifier(p_modifier);


        // Add the fixed boundary conditions
        // Note that the Adhesion to base membrane force is equivalent to a y=0 plane boundary condition
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        // At y=0
        normal(1)=-1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        p_bc1->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        // Run solve
        simulator.Solve();
    }

};