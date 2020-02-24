

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
#include "CellMutationStatesWriter.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellDataItemWriter.hpp"
#include "PlaneBasedCellKiller.hpp"

// Includes to use command line options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

// Claire's new classes
#include "PeriodicNdNodesOnlyMesh.hpp"
#include "NeighbourTrackingModifier.hpp"
#include "ExponentialDecayPalssonAdhesionForce.hpp"
#include "TopOfTissueForce.hpp"
#include "FixedDirectionCentreBasedDivisionRule.hpp"
#include "ParentAtBaseDivisionRule.hpp"
#include "TestFillTissueProject2.hpp"
#include "DetachedCellKiller.hpp"
#include "TissueHeightWriter.hpp"
#include "AttachmentPointCellWriter.hpp"

#include "Debug.hpp"

class TestVaryingDecayRateAndCellCycle : public AbstractCellBasedTestSuite
{
    void RunSimulation(double lambda_pow, unsigned cycle_i) throw(Exception);
public:
    // This test case takes the parameter combinations as input to run only one setup
    void DNRTestBaseParameterCombination() throw(Exception)
    {
	    RunSimulation(3,2);
    }

    void TestIndividualParameterCombinations3D() throw(Exception)
    {
        // Get indices for the setup
        unsigned lambda_i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-lambda_i");
        unsigned cycle_i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-cycle_i");
        // Check they are within the indices desired
        assert(lambda_i < 7);
        assert(cycle_i < 5);

        // Run simulation
	    RunSimulation(lambda_i,cycle_i);
    }
};

void TestVaryingDecayRateAndCellCycle::RunSimulation(double lambda_pow, unsigned cycle_i) throw(Exception)
{
    // Calculate parameters
    double lambda = 5.0*std::pow(2.0,lambda_pow)*0.001;
    unsigned cell_cycle = 2*(cycle_i+3);

    // Determine desired input height
    double height = 10.0;
    if ( lambda < 0.025 )
    {
        if ( cell_cycle == 6 )
        {
            height += 10.0;
            if ( lambda == 0.01 )
            {
                height += 10.0;
            }
        }
        if ( lambda < 0.015 & cell_cycle < 12 )
        {
            height += 10.0;
        }
        if ( lambda == 0.005 )
        {
            height += 10.0;
            if ( cell_cycle < 10)
            {
                height += 10.0;
            }
        }
        if (cell_cycle < 10)
        {
            height += 10;
        }
    }
    else if (cell_cycle == 6 & lambda < 0.3)
    {
        height += 10.0;
    }

    PRINT_VARIABLE(lambda);
    PRINT_VARIABLE(cell_cycle);
    PRINT_VARIABLE(height);

    // Determine input folder
    std::pair<std::string, double> fill_info = getFillMainDirNameStreamAndSimTime(height);
    std::stringstream input_file_stream;
    input_file_stream << fill_info.first << "/Seed";
    unsigned seed = getSeedAndAppendToFolder(input_file_stream,2);
    RandomNumberGenerator::Instance()->Reseed(seed);

    // Load the simulation results
    OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3>,3 >::Load(input_file_stream.str(), fill_info.second);

    // Set up output folder
    std::stringstream output_folder_stream;
    output_folder_stream << "Paper2/VaryingDecayRateAndCellCycle/DecayRate";
    output_folder_stream << (unsigned)(lambda*1000.0) << "e-3";
    output_folder_stream << "_CellCycle" << cell_cycle;   
    output_folder_stream << "/Seed";
    zeroFill(output_folder_stream,seed,2);
    p_simulator->SetOutputDirectory(output_folder_stream.str());

    // Set up the reporting
    CellBasedEventHandler::Reset();

    // The parameters for the sim   
    unsigned base_width = 10;
    double sim_length = 60.0*24.0;
    unsigned output_frequency = 120*24;

    // Setup timing parameters
    p_simulator->SetSamplingTimestepMultiple(output_frequency);
    p_simulator->SetEndTime(fill_info.second + sim_length);

    // Add the force at the top
    c_vector<double,3> force = zero_vector<double>(3);
    force[2] = gRemovalForce;
    MAKE_PTR_ARGS(TopOfTissueForce<3>, p_top_force, (force,60.0,1.5));
    p_simulator->AddForce(p_top_force);

    // Get a reference to the cell population
    NodeBasedCellPopulation<3>& r_population = dynamic_cast<NodeBasedCellPopulation<3>& >(p_simulator->rGetCellPopulation());

    // Add the detached cell killer
    MAKE_PTR_ARGS(DetachedCellKiller<3>, p_detached_killer, (&r_population,0.7));
    p_simulator->AddCellKiller(p_detached_killer);

    // Change the adhesion reduction rate
    const std::vector<boost::shared_ptr<AbstractForce<3, 3> > >& r_forces = p_simulator->rGetForceCollection();
    for ( unsigned i = 0; i < r_forces.size(); i++ )
    {
        boost::shared_ptr<ExponentialDecayPalssonAdhesionForce<3> > p_gradient_force = boost::dynamic_pointer_cast< ExponentialDecayPalssonAdhesionForce<3> >(r_forces[i]);
        if (p_gradient_force)
        {
            p_gradient_force->SetDecayRate(lambda);
        }
    }

    // Change the cell cycle lengths (not this only needs to be done for the stem cells)
    double cycle_duration = (double)(cell_cycle);
    for (typename NodeBasedCellPopulation<3>::Iterator cell_iter = (r_population.Begin()); 
        cell_iter != (r_population.End()); ++cell_iter)
    {
        if ( (*cell_iter)->GetCellProliferativeType()->template IsType<StemCellProliferativeType>() )
        {
            // Both cell cycle models inherit from AbstractPhaseBased, where the cell cycle lengths are set 
            AbstractPhaseBasedCellCycleModel* p_cell_cycle_model = dynamic_cast<AbstractPhaseBasedCellCycleModel*>((*cell_iter)->GetCellCycleModel());
            if (cycle_duration >= 10.0)
            {
                p_cell_cycle_model->SetStemCellG1Duration(cycle_duration-8.0);
            }
            else
            {
                assert(cycle_duration>2.0);
                double ratio = (cycle_duration-1)/13.0; // Fraction of the old durations we want to change to, excluding mitotic phase
                p_cell_cycle_model->SetStemCellG1Duration(ratio*4.0); 
                p_cell_cycle_model->SetSDuration(ratio*5.0);
                p_cell_cycle_model->SetG2Duration(ratio*4.0);
            }
            PRINT_VARIABLE(p_cell_cycle_model->GetAverageStemCellCycleTime());
            // Reset the birth times (otherwise they will group for longer simulations and do weird things at the start for shorter)
            double birth_time = (p_cell_cycle_model->GetAverageStemCellCycleTime())*(RandomNumberGenerator::Instance()->ranf());
            (*cell_iter)->SetBirthTime(SimulationTime::Instance()->GetTime() - birth_time);
        }
    }

    // Add the cell height and attachment point writers
    r_population.AddCellPopulationCountWriter<TissueHeightWriter>();
    r_population.AddCellWriter<AttachmentPointCellWriter>();

    // Add the forces writer and initialise cell population
    boost::shared_ptr<CellDataItemWriter<3,3> > p_cell_data_item_writer(new CellDataItemWriter<3,3>("membrane_force_magnitude"));
    r_population.AddCellWriter(p_cell_data_item_writer);

    for ( typename AbstractCellPopulation<3>::Iterator cell_iter = r_population.Begin(); cell_iter!= r_population.End(); ++cell_iter )
    {
        (*cell_iter)->GetCellData()->SetItem("membrane_force_magnitude",-1.0);
    }

    // Run solver
    p_simulator->Solve();

    // Reporting
    TRACE("\nFinished simulation:");
    PRINT_2_VARIABLES(lambda, cell_cycle);
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();

    // Save the simulation at this time point and delete pointer
    CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(p_simulator);
    delete p_simulator;
}
