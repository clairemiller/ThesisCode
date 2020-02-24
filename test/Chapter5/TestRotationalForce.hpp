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
#include "PlaneBasedCellKiller.hpp"
#include "CellIdWriter.hpp"
#include "HeightDependentDivisionModifier.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "ParentAtBaseDivisionRule.hpp"
#include "RotationalDivisionForce.hpp"
#include "FixedDirectionCentreBasedDivisionRule.hpp"
#include "CellHeightAtMPhaseCellModifier.hpp"
#include "CellDataItemWriter.hpp"
#include "UndulatingBaseMembraneAdhesionForce.hpp"

#include "../CommonFunctions.hpp"
#include "TestFillTissuePaper1.hpp" // Include for accessing required global variables
#include <boost/assign/list_of.hpp>

#include "Debug.hpp"

class TestRotationalForce : public AbstractCellBasedTestSuite
{
private:
    template<unsigned DIM>    
    void RunSimulation(double spring_constant, double spring_power, double adhesion_factor) throw(Exception);
public:
    void TestRotationalForce3d() throw(Exception)
    {
        std::vector<double> v_attachment_force = boost::assign::list_of(0.0)(0.1)(0.2)(0.4)(1.0);
        std::vector<double> v_rot_constant = boost::assign::list_of(1.0)(5.0)(10.0);

	unsigned attachment_i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-att_index");
	unsigned rotational_i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-rot_index");
	double spring_length_i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-spr_index");

	assert(attachment_i < 5);
	assert(rotational_i < 3);
	assert(spring_length_i < 3);
	
	double spring_power = 1.0+spring_length_i;
	RunSimulation<3>(v_rot_constant[rotational_i], spring_power, v_attachment_force[attachment_i]);

    }
};

template<unsigned DIM>
void TestRotationalForce::RunSimulation(double spring_constant, double spring_power, double adhesion_factor) throw(Exception)
{
    // New setup
    double continued_sim_length = 500.0*24.0;

    // Get appropriate folder paths
    std::string input_file = GetFillTissueOutputFolderAndReseed(spring_power);
    std::stringstream output_file_stream;
    output_file_stream << "Paper1/SecondSubmission/RotationalForce/";
    output_file_stream << "LogSpringLength";
    zeroFill(output_file_stream, (unsigned)(spring_power*10), 2);
    output_file_stream << "e-1_AdhesionMultiplier";
    zeroFill(output_file_stream, (unsigned)(adhesion_factor*10), 2);
    output_file_stream << "e-1_RotationalSpringConstant";
    zeroFill(output_file_stream, (unsigned)(spring_constant*100), 4);
    output_file_stream << "e-2/Seed";
    getSeedAndAppendToFolder(output_file_stream,2);

    // Load the simulation results
    OffLatticeSimulation<DIM>* p_simulator = CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM>,DIM >::Load(input_file, gTissueFillSimLength);

    // Set up the reporting
    CellBasedEventHandler::Reset();

    // Change the simulator setup
    p_simulator->SetOutputDirectory(output_file_stream.str());
    p_simulator->SetEndTime(gTissueFillSimLength + continued_sim_length);

    // Change the spring length and set the adhesion strength multiplier
    const std::vector<boost::shared_ptr<AbstractForce<DIM,DIM> > >& rForces = p_simulator->rGetForceCollection();
    for (unsigned i = 0; i < rForces.size(); i++)
    {
        boost::shared_ptr<UndulatingBaseMembraneAdhesionForce<DIM> > pBaseForce = boost::dynamic_pointer_cast< UndulatingBaseMembraneAdhesionForce<DIM> >(rForces[i]);
        if ( pBaseForce )
        {
            double base_value = pBaseForce->GetAlphaProlifCell();
            double diff_value = pBaseForce->GetAlphaDiffCell();
            pBaseForce->SetAdhesionParameters(base_value*adhesion_factor,diff_value);
        }
    }

    // Add the rotational division force
    MAKE_PTR_ARGS(RotationalDivisionForce<DIM>, p_rot_force, (spring_constant));
    p_simulator->AddForce(p_rot_force);

    // Add the cell modifier
    MAKE_PTR(CellHeightAtMPhaseCellModifier<DIM>, p_modifier);
    p_simulator->AddSimulationModifier(p_modifier);

    // Add the output of the cell modifier cell data
    NodeBasedCellPopulation<DIM>& p_population = dynamic_cast<NodeBasedCellPopulation<DIM>& >(p_simulator->rGetCellPopulation());
    boost::shared_ptr<CellDataItemWriter<DIM,DIM> > p_cell_data_item_writer(new CellDataItemWriter<DIM,DIM>("HeightAtDivision"));
    p_population.AddCellWriter(p_cell_data_item_writer);

    // Run solver
    p_simulator->Solve();

    // Reportings
    TRACE("\nFinished simulation:");
    PRINT_VARIABLE(spring_constant);
    PRINT_VARIABLE(spring_power);
    PRINT_VARIABLE(adhesion_factor);
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();

    // Save the results
    CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(p_simulator);

    // Delete pointer
    delete p_simulator;
    
}
