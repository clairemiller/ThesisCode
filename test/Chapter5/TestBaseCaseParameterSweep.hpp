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
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "PalssonAdhesionForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "PeriodicNdNodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "FlatBaseMembraneBoundaryCondition.hpp"
#include "UndulatingBaseMembraneAdhesionForce.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "CellIdWriter.hpp"
#include "HeightDependentDivisionModifier.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "ParentAtBaseDivisionRule.hpp"
#include "FixedDirectionCentreBasedDivisionRule.hpp"
#include "CellHeightAtMPhaseCellModifier.hpp"
#include "CellDataItemWriter.hpp"
#include "PlaneBasedCellKiller.hpp"

#include "../CommonFunctions.hpp"
#include "TestFillTissuePaper1.hpp" // Include for accessing required global variables
#include <boost/assign/list_of.hpp>

#include "Debug.hpp"

const double gContinuedSimLength = 500.0*24.0;

class TestFromFullTissue : public AbstractCellBasedTestSuite
{
private: 
    std::string GetOutputFolder(double spring_power, double adhesion_factor);
    void RunSimulation(double spring_power, double adhesion_factor) throw(Exception);
public:
    void TestParameterSweepBaseCase() throw(Exception)
    {
	    double attachment_factor = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-att_index");
        double spring_length_i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-spr_index");

        assert(attachment_factor < 4);
	    assert(spring_length_i < 3);

	    double spring_power = 1.0 + spring_length_i;
	    RunSimulation(spring_power,attachment_factor);
    }

    void TestContinuedSimulations() throw(Exception)
    {
        // Desired setup
        double attachment_factor = 2.0;
        double spring_length_i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-spr_index");
	    assert(spring_length_i < 2);
        double spring_power = 2.0 + spring_length_i;

        // Required simulation info
        std::string sim_dir = GetOutputFolder(spring_power, attachment_factor);
        double original_sim_length = gTissueFillSimLength + gContinuedSimLength;
        double new_sim_length = original_sim_length + 2.0*gContinuedSimLength;

        // Load the simulation results
        OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3>,3 >::Load(sim_dir, original_sim_length);

        // Set up the reporting
        CellBasedEventHandler::Reset();

        // Change the simulator setup
        p_simulator->SetEndTime(new_sim_length);
        
	// Run solver
        p_simulator->Solve();

        // Reportings
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        // Save the results
        CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(p_simulator);

        // Delete pointer
        delete p_simulator;
    }
};

std::string TestFromFullTissue::GetOutputFolder(double spring_power, double adhesion_factor)
{
    std::stringstream output_file_stream;
    output_file_stream << "Paper1/SecondSubmission/BaseSetupParameterSweep/";
    output_file_stream << "LogSpringLength";
    zeroFill(output_file_stream, (unsigned)(spring_power), 2);
    output_file_stream << "_AdhesionMultiplier";
    zeroFill(output_file_stream, (unsigned)(adhesion_factor), 2);
    output_file_stream << "/Seed";
    getSeedAndAppendToFolder(output_file_stream,2);
    return(output_file_stream.str());
}

void TestFromFullTissue::RunSimulation(double spring_power, double adhesion_factor) throw(Exception)
{
    // Get appropriate folder paths
    std::string input_dir = GetFillTissueOutputFolderAndReseed(spring_power);
    std::string output_dir = GetOutputFolder(spring_power, adhesion_factor);

    // Load the simulation results
    OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3>,3 >::Load(input_dir, gTissueFillSimLength);

    // Set up the reporting
    CellBasedEventHandler::Reset();

    // Change the simulator setup
    p_simulator->SetOutputDirectory(output_dir);
    p_simulator->SetEndTime(gTissueFillSimLength + gContinuedSimLength);

    // Set the adhesion strength multiplier
    const std::vector<boost::shared_ptr<AbstractForce<3,3> > >& rForces = p_simulator->rGetForceCollection();
    for (unsigned i = 0; i < rForces.size(); i++)
    {
        boost::shared_ptr<UndulatingBaseMembraneAdhesionForce<3> > pBaseForce = boost::dynamic_pointer_cast< UndulatingBaseMembraneAdhesionForce<3>  >(rForces[i]);
        if ( pBaseForce )
        {
            double base_value = pBaseForce->GetAlphaProlifCell();
            double diff_value = pBaseForce->GetAlphaDiffCell();
            pBaseForce->SetAdhesionParameters(base_value*adhesion_factor,diff_value);
        }
    }

    // Add the cell modifier
    MAKE_PTR(CellHeightAtMPhaseCellModifier<3>, p_modifier);
    p_simulator->AddSimulationModifier(p_modifier);

    // Add the output of the cell modifier cell data
    NodeBasedCellPopulation<3>& r_population = dynamic_cast<NodeBasedCellPopulation<3>& >(p_simulator->rGetCellPopulation());
    boost::shared_ptr<CellDataItemWriter<3,3> > p_cell_data_item_writer(new CellDataItemWriter<3,3>("HeightAtDivision"));
    r_population.AddCellWriter(p_cell_data_item_writer);

    // Run solver
    p_simulator->Solve();

    // Reportings
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();

    // Save the results
    CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(p_simulator);

    // Delete pointer
    delete p_simulator;
    
}
