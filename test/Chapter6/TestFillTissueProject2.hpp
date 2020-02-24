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
#include "ExponentialDecayPalssonAdhesionForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "PeriodicNdNodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "FlatBaseMembraneBoundaryCondition.hpp"
#include "UndulatingBaseMembraneAdhesionForce.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "CellIdWriter.hpp"
#include "HeightDependentDivisionModifier.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellDeathWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "LiRepulsionForce.hpp"
#include "CellCountThresholderSimulationModifier.hpp"
#include "VerticallyFixedStemCellBoundaryCondition.hpp"
#include "RotationalDivisionForce.hpp"
#include "FixedDirectionCentreBasedDivisionRule.hpp"
#include "TopOfTissueTrackingModifier.hpp"
#include "HeightDependentDifferentiationCellModifier.hpp"

#include "../CommonFunctions.hpp"
#include "Debug.hpp"

// Global variables that are needed for subsequent simulations
const double gRemovalForce = 0.5;

std::pair<std::string,double> getFillMainDirNameStreamAndSimTime(double height)
{
    // Main directory name
    std::stringstream fill_folder_stream;
    fill_folder_stream << "Paper2/FilledTissue/TissueHeight" << (int) height;
    
    // Simulation time
    double sim_time = 45.0*24.0;

    return(std::make_pair(fill_folder_stream.str(), sim_time));
}

class TestInitialTransient : public AbstractCellBasedTestSuite
{
    void RunSimulation(double lambda, double height) throw(Exception);
    void RunWithReducedHeight(double new_height, double old_height) throw(Exception);

public:
    void TestFillTissueBaseSetupHeight() throw(Exception)
    {
        RunSimulation(0.04,10.0);
        RunWithReducedHeight(10.0,10.0);
    }

    void Test3dFillTissue() throw(Exception)
    {
        // Want to run heights: 10, 20, 30, 40, 50
        // We run to 50 first (for 40 days) and then change to the other ones
        RunSimulation(0.04, 50.0);

        // Now run for another 5 days with a reduced height
        for (unsigned height = 10; height <= 50; height += 10)
        {
            PRINT_VARIABLE(height);
            RunWithReducedHeight((double) height, 50.0);
        }
    }
};

void TestInitialTransient::RunSimulation(double lambda, double height) throw(Exception)
{

    // Get the simulation info based on height
    std::pair<std::string, double> sim_info = getFillMainDirNameStreamAndSimTime(height);
    // Get the seed and assign folder
    std::stringstream output_folder_stream;
    output_folder_stream << sim_info.first << "/Seed";
    unsigned seed = getSeedAndAppendToFolder(output_folder_stream,2);
    RandomNumberGenerator::Instance()->Reseed(seed);

    // Set up the reporting
    CellBasedEventHandler::Reset();

    // The parameters for the sim 
    unsigned base_width = 10;
    unsigned n_stem = base_width*base_width;
    unsigned output_freq = 5*24*120;
    double spring_length = std::pow(10.0,-3.0);
    double peak_adhesion = gRemovalForce*2.0; // alpha = peak_adhesion / 0.0267131
    double minimum_age = 80.0;
    double simulation_length = sim_info.second-(5.0*24.0);

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

    // Add any extra output
    cell_population.template AddCellWriter<CellIdWriter>();
    cell_population.template AddCellWriter<CellMutationStatesWriter>();
    cell_population.template AddPopulationWriter<NodeVelocityWriter>();
    cell_population.template AddCellPopulationCountWriter<CellDeathWriter>();
    cell_population.template AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
    cell_population.template AddCellWriter<CellAgesWriter>();

    // Add the adhesive cell-cell force and the repulsion force
    MAKE_PTR_ARGS( ExponentialDecayPalssonAdhesionForce<3>, p_adh_force, (peak_adhesion,minimum_age,lambda));
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

    // Save the results
    CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(&simulator);

    // Reporting
    TRACE("\nFill simulation complete for:");
    PRINT_VARIABLE(height);
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();
}

void TestInitialTransient::RunWithReducedHeight(double new_height, double old_height) throw(Exception)
{
    // Get the old simulation
    std::pair<std::string, double> fill_info = getFillMainDirNameStreamAndSimTime(old_height);
    std::stringstream input_file_stream;
    input_file_stream << fill_info.first << "/Seed";
    unsigned seed = getSeedAndAppendToFolder(input_file_stream,2);
    double initial_fill_length = fill_info.second - (5.0*24.0);
    OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3>,3 >::Load(input_file_stream.str(), initial_fill_length);

    // Change the folder and end time
    std::stringstream output_file_stream;
    output_file_stream << (getFillMainDirNameStreamAndSimTime(new_height)).first << "/Seed";
    zeroFill(output_file_stream,seed,2);

    p_simulator->SetEndTime(fill_info.second);
    p_simulator->SetOutputDirectory(output_file_stream.str());

    // Change the sloughing height
    p_simulator->RemoveAllCellKillers();
    c_vector<double,3> pt = zero_vector<double>(3);
    c_vector<double,3> nml = zero_vector<double>(3);
    pt[2] = new_height;
    nml[2] = 1.0;
    NodeBasedCellPopulation<3>& r_population = dynamic_cast<NodeBasedCellPopulation<3>& >(p_simulator->rGetCellPopulation());
    MAKE_PTR_ARGS(PlaneBasedCellKiller<3>,p_killer,(&r_population,pt,nml));
    p_simulator->AddCellKiller(p_killer);

    // Run solver
    p_simulator->Solve();

    // Remove the population boundary conditions and re-add the base membrane
    p_simulator->RemoveAllCellPopulationBoundaryConditions();
    MAKE_PTR_ARGS(FlatBaseMembraneBoundaryCondition<3>, p_base_bc, (&r_population));
    p_base_bc->SetUseJiggledBottomCells(true);
    p_base_bc->SetIsInitialised(true);
    p_simulator->AddCellPopulationBoundaryCondition(p_base_bc);
    
    // Remove the cell killer
    p_simulator->RemoveAllCellKillers();
    
    // Change the output frequency back to daily
    p_simulator->SetSamplingTimestepMultiple(24*120);

    // Save the results and delete pointer
    CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(p_simulator);
    delete p_simulator;

}

