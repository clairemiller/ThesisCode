


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>

// To run in parallel
// #include "PetscSetupAndFinalize.hpp"
// When run in serial
#include "FakePetscSetup.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellAncestorWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellDataItemWriter.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"

// Includes to use command line options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

// Claire's new classes
#include "PeriodicNdNodesOnlyMesh.hpp"
#include "NeighbourTrackingModifier.hpp"
#include "TopOfTissueForce.hpp"
#include "FixedDirectionCentreBasedDivisionRule.hpp"
#include "ParentAtBaseDivisionRule.hpp"
#include "TestFillTissueProject2.hpp"
#include "DetachedCellKiller.hpp"
#include "TissueHeightWriter.hpp"

#include "Debug.hpp"

class TestCupScenario : public AbstractCellBasedTestSuite
{
public:
    void Test3dProtectedSkinScenario() throw(Exception)
    {  
        // Get input info
        std::pair<std::string, double> fill_info = getFillMainDirNameStreamAndSimTime(10.0);

        // Timing parameters
        double force_off = fill_info.second + 60.0*24.0;
        double force_on = force_off + 3.0*24.0;
        double end_sim = force_on + 40.0*24.0;

        // Other simulation parameters
        unsigned base_width = 10;
        unsigned output_frequency = 120;

        // Set up the reporting
        CellBasedEventHandler::Reset();
        
        // Input simulation results
        std::stringstream input_file_stream;
        input_file_stream << fill_info.first + "/Seed";
        unsigned seed = getSeedAndAppendToFolder(input_file_stream,2);
        RandomNumberGenerator::Instance()->Reseed(seed);
        
        // Load the simulation results
        OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3>,3 >::Load(input_file_stream.str(), fill_info.second);
        
        // Create the output folder
        std::stringstream output_folder_stream;
        output_folder_stream << "Paper2/TestProtectedSkinScenario" << "/Seed";
        zeroFill(output_folder_stream,seed,2);

        // Add the force at the top
        c_vector<double,3> force = zero_vector<double>(3);
        force[2] = gRemovalForce;
        MAKE_PTR_ARGS(TopOfTissueForce<3>, p_neighbour_force, (force,60.0,1.5));
        p_simulator->AddForce(p_neighbour_force);

        // Add the detached cell killer
        MAKE_PTR_ARGS(DetachedCellKiller<3>, p_detached_killer, (&(p_simulator->rGetCellPopulation()),1.5));
        p_simulator->AddCellKiller(p_detached_killer);

        // Add the sloughing at the top (to speed up computation)
        c_vector<double,3> pt = zero_vector<double>(3);
        c_vector<double,3> nml = zero_vector<double>(3);
        pt[2] = 20.0;
        nml[2] = 1.0;
        AbstractCellPopulation<3>& r_population = (p_simulator->rGetCellPopulation());
        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>,p_killer,(&(p_simulator->rGetCellPopulation()),pt,nml));
        p_simulator->AddCellKiller(p_killer);

        // Changes to the simulator
        p_simulator->SetOutputDirectory(output_folder_stream.str());
        p_simulator->SetSamplingTimestepMultiple(output_frequency);

        // Add the cell height writer
        r_population.template AddCellPopulationCountWriter<TissueHeightWriter>();

        // Add output of division information
        p_simulator->SetOutputDivisionLocations(true);

        // Initial force on
        p_simulator->SetEndTime(force_off);
        p_simulator->Solve();

        // Turn force off for a bit
        c_vector<double,3> zero_force = zero_vector<double>(3);
        p_neighbour_force->SetAppliedForce(zero_force);
        p_simulator->SetEndTime(force_on);
        p_simulator->Solve();

        // Then turn it back on again
        p_neighbour_force->SetAppliedForce(force);
        p_simulator->SetEndTime(end_sim);
        p_simulator->Solve();

        // Reporting
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        // Reset the singletons
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};
