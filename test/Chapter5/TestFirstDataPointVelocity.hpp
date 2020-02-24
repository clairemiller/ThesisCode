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

#include "CellBasedEventHandler.hpp"
#include "RotationalDivisionForce.hpp"
#include "VerticallyFixedStemCellBoundaryCondition.hpp"

#include "../CommonFunctions.hpp"
#include "TestFillTissuePaper1.hpp" // Include for accessing required global variables
#include <boost/assign/list_of.hpp>

#include "Debug.hpp"

class TestFirstDataPointVelocity : public AbstractCellBasedTestSuite
{
public:
    void TestPinned() throw(Exception)
    {
	double continued_sim_length = 0.1;

	// Get appropriate folder paths
	std::string input_file = GetFillTissueOutputFolderAndReseed(3.0);
	std::stringstream output_file_stream;
        output_file_stream << "Paper1/SecondSubmission/DataPointTApprox0/PinnedCase/Seed";
	getSeedAndAppendToFolder(output_file_stream,2);

        // Load the simulation results
        OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3>,3 >::Load(input_file, gTissueFillSimLength);
	
	// Add the stem cell restriction
	MAKE_PTR_ARGS(VerticallyFixedStemCellBoundaryCondition<3>, p_sc_bc, (&(p_simulator->rGetCellPopulation())));
        p_simulator->AddCellPopulationBoundaryCondition(p_sc_bc);	

	// Set up the reporting
	CellBasedEventHandler::Reset();
	
	// Change the simulator setup
	p_simulator->SetOutputDirectory(output_file_stream.str());
        p_simulator->SetEndTime(gTissueFillSimLength + continued_sim_length);
        p_simulator->SetSamplingTimestepMultiple(1);

	// Run solver
	p_simulator->Solve();

	// Reportings
	TRACE("Pinned case");
	CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

	// Delete pointer
	delete p_simulator;
    }


    void TestBaseCase() throw(Exception)
    {
        // Only want to run minimal length
        double continued_sim_length = 0.1;

        // Get appropriate folder paths
        std::string input_file = GetFillTissueOutputFolderAndReseed(3.0);
        std::stringstream output_file_stream;
        output_file_stream << "Paper1/SecondSubmission/DataPointTApprox0/BaseCase/Seed";
        getSeedAndAppendToFolder(output_file_stream,2);

        // Load the simulation results
        OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3>,3 >::Load(input_file, gTissueFillSimLength);

        // Set up the reporting
        CellBasedEventHandler::Reset();

        // Change the simulator setup
        p_simulator->SetOutputDirectory(output_file_stream.str());
        p_simulator->SetEndTime(gTissueFillSimLength + continued_sim_length);
	p_simulator->SetSamplingTimestepMultiple(1);

        // Run solver
        p_simulator->Solve();

        // Reportings
        TRACE("Base case");
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        // Delete pointer
        delete p_simulator;
    }

    void TestRotationalCase() throw(Exception)
    {
        std::vector<double> v_rot_constant = boost::assign::list_of(1.0)(5.0)(10.0)(20.0);
        for (unsigned i = 0; i < v_rot_constant.size(); i++)
        {
            RunRotationalCase(v_rot_constant[i]);
        }
    }

    void RunRotationalCase(double spring_constant) throw(Exception)
    {
        // Only want to run minimal length
        double continued_sim_length = 0.1;

        // Get appropriate folder paths
        std::string input_file = GetFillTissueOutputFolderAndReseed(3.0);
        std::stringstream output_file_stream;
        output_file_stream << "Paper1/SecondSubmission/DataPointTApprox0/Rotational/SpringConstant";
        output_file_stream << (unsigned)(spring_constant) << "/Seed";
        getSeedAndAppendToFolder(output_file_stream,2);

        // Load the simulation results
        OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3>,3 >::Load(input_file, gTissueFillSimLength);

        // Set up the reporting
        CellBasedEventHandler::Reset();

        // Change the simulator setup
        p_simulator->SetOutputDirectory(output_file_stream.str());
        p_simulator->SetEndTime(gTissueFillSimLength + continued_sim_length);
	p_simulator->SetSamplingTimestepMultiple(1);

        // Add the rotational division force
        MAKE_PTR_ARGS(RotationalDivisionForce<3>, p_rot_force, (spring_constant));
        p_simulator->AddForce(p_rot_force);

        // Run solver
        p_simulator->Solve();

        // Reportings
        TRACE("Rotational case");
	PRINT_VARIABLE(spring_constant);
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        // Delete pointer
        delete p_simulator;
    }
};
