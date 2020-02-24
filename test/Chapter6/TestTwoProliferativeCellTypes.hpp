

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
#include "TopOfTissueForce.hpp"
#include "FixedDirectionCentreBasedDivisionRule.hpp"
#include "ParentAtBaseDivisionRule.hpp"
#include "TestFillTissueProject2.hpp"
#include "DetachedCellKiller.hpp"
#include "TissueHeightWriter.hpp"

#include "Debug.hpp"

class TestTwoProliferativeCellTypes : public AbstractCellBasedTestSuite
{
    void RunSimulation(double short_cycle_length, double long_cycle_length, std::string output_folder) throw(Exception);

public:

    void RunUsingArithmeticMean(unsigned i) throw(Exception)
    {
        std::stringstream output_folder_stream;
        output_folder_stream << "Paper2/TwoProliferativeCellTypes/ArithmeticMean/CycleDiff";
        
        double cycle_diff = ((double) i) + 1.0;
        double short_cycle_length = 12.0-cycle_diff;
        double long_cycle_length = 12.0 + cycle_diff;
	    output_folder_stream << (unsigned) cycle_diff << "hr_ShortCycle";
        output_folder_stream << (unsigned) short_cycle_length << "hr_LongCycle";
        output_folder_stream << (unsigned) long_cycle_length << "hr/";

        RunSimulation(short_cycle_length, long_cycle_length, output_folder_stream.str());
    }

    void RunUsingHarmonicMean(double deltaT) throw(Exception)
    {
        double short_cycle_length = 0.5*(12.0-deltaT) + 0.5*std::sqrt(std::pow(deltaT-12.0,2.0) + 24.0*deltaT);
        double long_cycle_length = short_cycle_length + deltaT;

        std::stringstream output_folder_stream;
        output_folder_stream << "Paper2/TwoProliferativeCellTypes/HarmonicMean/";
        output_folder_stream << "CycleDiff" << (unsigned) deltaT << "hr_";
        output_folder_stream << "ShortCycle" << (unsigned) (short_cycle_length*100) << "e-2hr_";
        output_folder_stream << "LongCycle" << (unsigned) (long_cycle_length*100) << "e-2hr/";

        RunSimulation(short_cycle_length,long_cycle_length,output_folder_stream.str());
    }

    void TestRunArithmeticMeanFromInput() throw(Exception)
    {
        // Get the index for the simulation to run
        unsigned i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-index");
        assert(i<5);
        RunUsingArithmeticMean(i);
    }

    void TestRunHarmonicMeanFromInput() throw(Exception)
    {
        // Get the index for the simulation to run
        unsigned i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-index");
        assert(i<4);
        double deltaT = (double) std::pow(2,i);
        RunUsingArithmeticMean(deltaT);
    }

    void DNRTestRunAll() throw(Exception)
    {
        for (unsigned i = 0; i < 5; i++)
        {
            RunUsingArithmeticMean(i);
        }
    }
};

void TestTwoProliferativeCellTypes::RunSimulation(double short_cycle_length, double long_cycle_length, std::string output_folder) throw(Exception)
{   
    double decay_rate = 0.04;

    // Fix up the cell cycle lengths to the minimum of the stochastics
    short_cycle_length = short_cycle_length - 2.0;
    long_cycle_length = long_cycle_length - 2.0;

    // Get input info
    std::pair<std::string, double> fill_info = getFillMainDirNameStreamAndSimTime(10.0);
    std::stringstream input_file_stream;
    input_file_stream << fill_info.first + "/Seed";
    unsigned seed = getSeedAndAppendToFolder(input_file_stream,2);
    RandomNumberGenerator::Instance()->Reseed(seed);

    // Set up the reporting
    CellBasedEventHandler::Reset();

    // The parameters for the sim   
    unsigned base_width = 10;
    double sim_length = 60.0*24.0;
    unsigned output_frequency = 120;

    // Load the simulation results
    OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3>,3 >::Load(input_file_stream.str(), fill_info.second);

    // Setup output
    std::stringstream output_folder_stream;
    output_folder_stream << output_folder << "/Seed";
    zeroFill(output_folder_stream,seed,2);
    p_simulator->SetOutputDirectory(output_folder_stream.str());
    p_simulator->SetSamplingTimestepMultiple(output_frequency);
    p_simulator->SetEndTime(fill_info.second + sim_length);

    // Add output of division information
    p_simulator->SetOutputDivisionLocations(true);

    // Add the force at the top
    c_vector<double,3> force = zero_vector<double>(3);
    force[2] = gRemovalForce;
    MAKE_PTR_ARGS(TopOfTissueForce<3>, p_top_force, (force,60.0,1.5));
    p_simulator->AddForce(p_top_force);

    // Add the top of tissue tracking modifier
    MAKE_PTR(TopOfTissueTrackingModifier<3>,p_topoftissue);
    p_simulator->AddSimulationModifier(p_topoftissue);

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
            p_gradient_force->SetDecayRate(decay_rate);
        }
    }

    // Change the cell cycle lengths (not this only needs to be done for the stem cells)
    // Ref: Sada et al. (2016). Defining the cellular lineage hierarchy in the interfollicular epidermis of adult skin
    // The study shows 65% fast cycling cells dividing every 2 days, and 35% slow cycling dividing every 3 days
    unsigned n_stem_cells = 100;
    unsigned n_fast_cycling = 50;
    unsigned n_slow_cycling = 50;
    unsigned count = 0;
    for (typename NodeBasedCellPopulation<3>::Iterator cell_iter = (r_population.Begin()); 
        cell_iter != (r_population.End()); ++cell_iter)
    {
        if ( (*cell_iter)->GetCellProliferativeType()->template IsType<StemCellProliferativeType>() )
        {
            // Both cell cycle models inherit from AbstractPhaseBased, where the cell cycle lengths are set 
            AbstractPhaseBasedCellCycleModel* p_cell_cycle_model = dynamic_cast<AbstractPhaseBasedCellCycleModel*>((*cell_iter)->GetCellCycleModel());
            double cycle_duration = p_cell_cycle_model->GetAverageStemCellCycleTime();
	    if (count < n_fast_cycling)
            {
	 	// Short cycle in [5,10]
                p_cell_cycle_model->SetStemCellG1Duration(short_cycle_length-4.0); 
                p_cell_cycle_model->SetSDuration(2.0);
                p_cell_cycle_model->SetG2Duration(1.0);
                count++;
		TS_ASSERT_DELTA(p_cell_cycle_model->GetAverageStemCellCycleTime(),short_cycle_length,1e-6);
            }
            else
            {
                assert(count<n_stem_cells); // Check, should never hit
                p_cell_cycle_model->SetStemCellG1Duration(long_cycle_length-8.0); // Always greater than 10
                count++;
		TS_ASSERT_DELTA(p_cell_cycle_model->GetAverageStemCellCycleTime(),long_cycle_length,1e-6);
            }
        }
    }
    assert(count==n_stem_cells);

    // Add the cell height writer
    r_population.template AddCellPopulationCountWriter<TissueHeightWriter>();

    // Run solver
    p_simulator->Solve();

    // Reporting
    TRACE("\nSimulation completed:");
    PRINT_2_VARIABLES(short_cycle_length, long_cycle_length);
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();

    // Delete pointer and reset the singletons
    delete p_simulator;
    SimulationTime::Destroy();
    SimulationTime::Instance()->SetStartTime(0.0);
    RandomNumberGenerator::Destroy();
}
