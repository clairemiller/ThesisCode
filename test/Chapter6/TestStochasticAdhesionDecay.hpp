

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
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "CellDataItemWriter.hpp"
#include "CellIdWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellAgesWriter.hpp"


// Includes to use command line options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

// Claire's new classes
#include "PeriodicNdNodesOnlyMesh.hpp"
#include "NeighbourTrackingModifier.hpp"
#include "StochasticExponentialDecayPalssonAdhesionForce.hpp"
#include "TopOfTissueForce.hpp"
#include "FixedDirectionCentreBasedDivisionRule.hpp"
#include "ParentAtBaseDivisionRule.hpp"
#include "DetachedCellKiller.hpp"
#include "CellDeathWriter.hpp"
#include "VerticallyFixedStemCellBoundaryCondition.hpp"
#include "CellCountThresholderSimulationModifier.hpp"
#include "HeightDependentDifferentiationCellModifier.hpp"
#include "TopOfTissueTrackingModifier.hpp"
#include "TissueHeightWriter.hpp"
#include "FlatBaseMembraneBoundaryCondition.hpp"
#include "UndulatingBaseMembraneAdhesionForce.hpp"
#include "RotationalDivisionForce.hpp"
#include "CellDataItemWriter.hpp"

#include "../CommonFunctions.hpp"
#include "Debug.hpp"

const double gFillLength = 21.0*24.0;
const double gRemovalForce = 0.5;

class TestStochasticAdhesionDecay : public AbstractCellBasedTestSuite
{
    void RunFill(double min_decay_rate, double max_decay_rate, std::string output_folder) throw(Exception);
    void RunSimulation(std::string input_folder, std::string output_folder) throw(Exception);
public:

    void RunFromVariationPercent(unsigned lambda_percent) throw(Exception)
    {
	    double lambda = 0.04;
        double lambda_var = lambda*(((double)lambda_percent)/100.0);
        
        std::stringstream parameters_stream;
        parameters_stream << "Variation_" <<  lambda_percent << "Percent/";
        std::string fill_output_folder = "Paper2/FilledTissue/StochasticDecay/" + (parameters_stream.str());
        RunFill(lambda-lambda_var, lambda+lambda_var, fill_output_folder);

        std::string simulation_output_folder = "Paper2/StochasticDecay/" + (parameters_stream.str());
        RunSimulation(fill_output_folder, simulation_output_folder);
    }

    void DNRTestDifferentLambdaSequentially() throw(Exception)
    {
        double lambda = 0.04;

        for (unsigned i = 0; i < 5; i++)
        {
            unsigned lambda_percent = i*20;
            RunFromVariationPercent(lambda_percent);
        }
    }

    void TestDifferentLambdaIndividually() throw(Exception)
    {
        unsigned i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-index");
        assert(i<9);
        
        unsigned lambda_percent = i*10;
        RunFromVariationPercent(lambda_percent);
    }
};


void TestStochasticAdhesionDecay::RunSimulation(std::string input_folder, std::string output_folder) throw(Exception)
{
    // Parameters
    double sim_length = 60.0*24.0;
    double output_freq = 120.0*24.0;

    // Get the input folder
    std::stringstream input_file_stream;
    input_file_stream << input_folder + "/Seed";
    unsigned seed = getSeedAndAppendToFolder(input_file_stream,2);
    RandomNumberGenerator::Instance()->Reseed(seed);

     // Set up the reporting
     CellBasedEventHandler::Reset();

     // Load the simulation results
     OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3>,3 >::Load(input_file_stream.str(), gFillLength);

    // Setup the output folder
    std::stringstream output_folder_stream;
    output_folder_stream << output_folder << "/Seed";
    zeroFill(output_folder_stream,seed,2);

    // Change relevant properties of the simulator
    p_simulator->SetOutputDirectory(output_folder_stream.str());
    p_simulator->SetSamplingTimestepMultiple(output_freq);
    p_simulator->SetEndTime(gFillLength + sim_length);

    // Add the force at the top
    c_vector<double,3> force = zero_vector<double>(3);
    force[2] = gRemovalForce;
    MAKE_PTR_ARGS(TopOfTissueForce<3>, p_top_force, (force,60.0,1.5));
    p_simulator->AddForce(p_top_force);

    // Add the detached cell killer
    NodeBasedCellPopulation<3>& r_population = dynamic_cast<NodeBasedCellPopulation<3>& >(p_simulator->rGetCellPopulation());
    MAKE_PTR_ARGS(DetachedCellKiller<3>, p_detached_killer, (&r_population,0.7));
    p_simulator->AddCellKiller(p_detached_killer);

    // Add the tissue height output
    r_population.AddCellPopulationCountWriter<TissueHeightWriter>();

    // Run solver
    p_simulator->Solve();

    // Save the results
    CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(p_simulator);

    // Reporting
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();

    // Delete pointer and reset the singletons
    delete p_simulator;
    SimulationTime::Destroy();
    SimulationTime::Instance()->SetStartTime(0.0);
    RandomNumberGenerator::Destroy();
}

void TestStochasticAdhesionDecay::RunFill(double min_decay_rate, double max_decay_rate, std::string output_folder) throw(Exception)
{   
    // Get the seed and assign folder
    std::stringstream output_folder_stream;
    output_folder_stream << output_folder;
    output_folder_stream << "/Seed";
    unsigned seed = getSeedAndAppendToFolder(output_folder_stream,2);
    RandomNumberGenerator::Instance()->Reseed(seed);

    // Set up the reporting
    CellBasedEventHandler::Reset();

    // The parameters for the sim 
    unsigned base_width = 10;
    double height = (double) base_width;
    unsigned n_stem = base_width*base_width;
    unsigned output_freq = 24*120;
    double spring_length = std::pow(10.0,-3.0);
    double peak_adhesion = gRemovalForce*2.0; 
    double minimum_age = 80.0;
    double simulation_length = gFillLength;

    // Set up the nodes in standard setup
    std::vector<Node<3>*> nodes(n_stem);
    for ( unsigned i=0; i < base_width; i++ )
    {
        for (unsigned j = 0; j < base_width; j++)
        {
            unsigned id = j*base_width + i;
            double x = i;
            double y = j;
            double z = 0.0;
            nodes[id] = new Node<3>( id,false, x,y,z);
        }
    }

    // Construct the mesh
    std::vector<double> periodic_widths(2, (double) base_width);
    PeriodicNdNodesOnlyMesh<3> mesh(periodic_widths,true,true,false);
    mesh.ConstructNodesWithoutMesh(nodes,2.0);

    // Create the cells
    std::vector<CellPtr> cells;
    MAKE_PTR(StemCellProliferativeType, p_stem_type);
    for (unsigned i=0; i<n_stem; i++)
    {
	    UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel;
        p_cell_cycle_model->SetDimension(3);
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell_cycle_model->SetMaxTransitGenerations(0);
        // Note: durations are S=5 and M=1. So this gives a total of 10 hours
        p_cell_cycle_model->SetG2Duration(2);
        p_cell_cycle_model->SetStemCellG1Duration(2); 

        double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
        p_cell->SetBirthTime(birth_time);

        cells.push_back(p_cell);
    }
    
    // Create cell population
    NodeBasedCellPopulation<3> cell_population(mesh, cells);
    cell_population.SetAbsoluteMovementThreshold(1.5);

    // Set up simulator
    OffLatticeSimulation<3> simulator(cell_population);
    simulator.SetOutputDirectory(output_folder_stream.str());
    simulator.SetSamplingTimestepMultiple(output_freq);
    simulator.SetEndTime(simulation_length);

    // Add output of division information
    simulator.SetOutputDivisionLocations(true);

    // Add any extra output
    cell_population.template AddCellWriter<CellIdWriter>();
    cell_population.template AddCellWriter<CellMutationStatesWriter>();
    cell_population.template AddPopulationWriter<NodeVelocityWriter>();
    cell_population.template AddCellPopulationCountWriter<CellDeathWriter>();
    cell_population.template AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
    cell_population.template AddCellWriter<CellAgesWriter>();
    boost::shared_ptr<CellDataItemWriter<3,3> > p_lambda_writer(new CellDataItemWriter<3,3>("AdhesionDecayRate"));
    cell_population.AddCellWriter(p_lambda_writer);

    // Add the adhesive cell-cell force and the repulsion force
    MAKE_PTR_ARGS( StochasticExponentialDecayPalssonAdhesionForce<3>, p_adh_force, (peak_adhesion,minimum_age,min_decay_rate,max_decay_rate));
    p_adh_force->SetMeinekeSpringGrowthDuration(1.0);
    p_adh_force->SetMeinekeDivisionRestingSpringLength(spring_length);
    simulator.AddForce(p_adh_force);
    MAKE_PTR( RepulsionForce<3>, p_rep_force);
    p_rep_force->SetMeinekeDivisionRestingSpringLength(spring_length);
    simulator.AddForce(p_rep_force);

    // Add the sloughing at the top
    c_vector<double,3> pt = zero_vector<double>(3);
    c_vector<double,3> nml = zero_vector<double>(3);
    pt[2] = height;
    nml[2] = 1.0;
    MAKE_PTR_ARGS(PlaneBasedCellKiller<3>,p_killer,(&cell_population,pt,nml));
    simulator.AddCellKiller(p_killer);

    // Add the bottom boundaries
    MAKE_PTR_ARGS(FlatBaseMembraneBoundaryCondition<3>, p_base_bc, (&cell_population));
    p_base_bc->SetUseJiggledBottomCells(true);
    simulator.AddCellPopulationBoundaryCondition(p_base_bc);
    MAKE_PTR(UndulatingBaseMembraneAdhesionForce<3>, p_base_force);
    // Alpha* = 50 when alpha = 0.02
    // Therefore alpha* = (50/0.02)*alpha = 50*50*alpha
    // Here, alpha = peak_adhesion / 0.0267131
    double peak_base_force = 50.0/ 0.0267131;
    p_base_force->SetAdhesionParameters(peak_base_force,0.0);
    simulator.AddForce(p_base_force);

    // Add the restriction on the stem cells for the fill
    MAKE_PTR_ARGS(VerticallyFixedStemCellBoundaryCondition<3>, p_sc_bc, (&cell_population));
    simulator.AddCellPopulationBoundaryCondition(p_sc_bc);

    // Add the rotational force
    MAKE_PTR_ARGS(RotationalDivisionForce<3>,p_rot_force, (10.0));
    simulator.AddForce(p_rot_force);

    // Add the modifier to change stem cells to differentiated at a certain height
    MAKE_PTR_ARGS(HeightDependentDifferentiationCellModifier<3>, p_diff_modifier, (2.0));
    simulator.AddSimulationModifier(p_diff_modifier);

    // Add in the division direction
    c_vector<double,3> div_vec = zero_vector<double>(3);
    div_vec[2] = std::pow(10.0,-3.0);
    MAKE_PTR_ARGS(FixedDirectionCentreBasedDivisionRule<3>,p_div_rule,(div_vec));
    cell_population.SetCentreBasedDivisionRule(p_div_rule);

    // Add the top of tissue tracking
    MAKE_PTR(TopOfTissueTrackingModifier<3>, p_top_modifier);
    simulator.AddSimulationModifier(p_top_modifier);

    // Run solver
    simulator.Solve();

    // Remove the sloughing
    simulator.RemoveAllCellKillers();

    // Remove the extra boundary condition (need to re-add the base membrane)
    simulator.RemoveAllCellPopulationBoundaryConditions();
    p_base_bc->SetUseJiggledBottomCells(true);
    p_base_bc->SetIsInitialised(true);
    simulator.AddCellPopulationBoundaryCondition(p_base_bc);

    // Save the results
    CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(&simulator);

    // Reporting
    TRACE("\nFill simulation complete for seed");
    PRINT_VARIABLE(seed);
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();

    // Destroy the singletons
    SimulationTime::Destroy();
    SimulationTime::Instance()->SetStartTime(0.0);
    RandomNumberGenerator::Destroy();
}
