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
#include "Debug.hpp"

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
#include "NodeVelocityWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellDeathWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "LiRepulsionForce.hpp"
#include "CellCountThresholderSimulationModifier.hpp"
#include "VerticallyFixedStemCellBoundaryCondition.hpp"
#include "CellDataItemWriter.hpp"
#include "FixedDirectionCentreBasedDivisionRule.hpp"

#include "../CommonFunctions.hpp"

// Global variables that are needed for subsequent simulations
const double gTissueFillSimLength = 21.0*24.0;

std::string GetFillTissueOutputFolderAndReseed(double log_spring_length)
{
    assert(std::abs(log_spring_length) >= 1.0);
    std::stringstream output_file_stream;
    output_file_stream << "Paper1/FilledTissue/";
    output_file_stream << "AbsLogSpringLength0" << (unsigned)std::abs(log_spring_length) << "/Seed";
    unsigned seed = getSeedAndAppendToFolder(output_file_stream,2);
    RandomNumberGenerator::Instance()->Reseed(seed);
    return(output_file_stream.str());
}

class TestInitialTransient : public AbstractCellBasedTestSuite
{
    void RunFill(double log_spring_length) throw(Exception);
    
public:
    void TestFillTissue() throw(Exception)
    {
	    // Run for each spring length
        double spring_length_index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-spr_index");
        assert(spring_length_index < 3.0);
        double log_sl = 1.0 + spring_length_index;
        RunFill(-1.0*log_sl);
    }
};

void TestInitialTransient::RunFill(double log_spring_length) throw(Exception)
{
    // Get the seed and assign folder
    std::string output_folder = GetFillTissueOutputFolderAndReseed(log_spring_length);

    // Set up the reporting
    CellBasedEventHandler::Reset();

    // The parameters for the sim   
    unsigned base_width = 10;
    double height = base_width;
    unsigned n_stem = base_width*base_width;
    unsigned output_freq = 5*24*120;
    double spring_length = std::pow(10.0,log_spring_length);

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
    CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
    MAKE_PTR(StemCellProliferativeType, p_stem_type);
    cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
    for ( std::vector<CellPtr>::iterator cell_it = cells.begin(); cell_it != cells.end(); cell_it++ )
    {
        UniformG1GenerationalCellCycleModel* p_model = dynamic_cast<UniformG1GenerationalCellCycleModel*>((*cell_it)->GetCellCycleModel());
        p_model->SetMaxTransitGenerations(0);
        p_model->SetStemCellG1Duration(3);
        // Add membrane_force_magnitude to cell data for output
        (*cell_it)->GetCellData()->SetItem("membrane_force_magnitude",-1.0);
    }

    // Create cell population
    NodeBasedCellPopulation<3> cell_population(mesh, cells);
    cell_population.SetAbsoluteMovementThreshold(1.5);

    // Set up simulator
    OffLatticeSimulation<3> simulator(cell_population);
    simulator.SetOutputDirectory(output_folder);
    simulator.SetSamplingTimestepMultiple(output_freq);
    simulator.SetEndTime(gTissueFillSimLength);

    // Add any extra output
    cell_population.AddCellWriter<CellIdWriter>();
    cell_population.AddCellWriter<CellMutationStatesWriter>();
    cell_population.AddPopulationWriter<NodeVelocityWriter>();
    cell_population.AddCellPopulationCountWriter<CellDeathWriter>();
    cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
    boost::shared_ptr<CellDataItemWriter<3,3> > p_cell_data_item_writer(new CellDataItemWriter<3,3>("membrane_force_magnitude"));
    cell_population.AddCellWriter(p_cell_data_item_writer);

    // Add the adhesive cell-cell force and the repulsion force
    MAKE_PTR( PalssonAdhesionForce<3>, p_adhesion_force );
    MAKE_PTR(RepulsionForce<3>, p_rep_force);
    // Change the spring lengths
    p_adhesion_force->SetMeinekeDivisionRestingSpringLength(spring_length);
    p_rep_force->SetMeinekeDivisionRestingSpringLength(spring_length);
    simulator.AddForce(p_adhesion_force);
    simulator.AddForce(p_rep_force);

    // Add the sloughing at the top
    c_vector<double,3> pt = zero_vector<double>(3);
    c_vector<double,3> nml = zero_vector<double>(3);
    pt[2] = height;
    nml[2] = 1.0;
    MAKE_PTR_ARGS(PlaneBasedCellKiller<3>,p_killer,(&cell_population,pt,nml));
    simulator.AddCellKiller(p_killer);

    // Add the bottom boundary
    MAKE_PTR_ARGS(FlatBaseMembraneBoundaryCondition<3>,p_membrane_bc,(&cell_population));
    p_membrane_bc->SetUseJiggledBottomCells(true);
    simulator.AddCellPopulationBoundaryCondition(p_membrane_bc);
    MAKE_PTR(UndulatingBaseMembraneAdhesionForce<3>, p_membrane_force);
    simulator.AddForce(p_membrane_force);
    
    // Add the restriction on the stem cells
    MAKE_PTR_ARGS(VerticallyFixedStemCellBoundaryCondition<3>, p_sc_bc, (&cell_population));
    simulator.AddCellPopulationBoundaryCondition(p_sc_bc);

    // Add the simulation modifier to check that the output of the simulation is not too large
    MAKE_PTR_ARGS(CellCountThresholderSimulationModifier<3>, p_count_modifier, (base_width*base_width*10*2));
    simulator.AddSimulationModifier(p_count_modifier);

    // Add in the division direction
    c_vector<double,3> div_vec = zero_vector<double>(3);
    div_vec[2] = spring_length;
    MAKE_PTR_ARGS(FixedDirectionCentreBasedDivisionRule<3>,p_div_rule,(div_vec));
    cell_population.SetCentreBasedDivisionRule(p_div_rule);

    // Run solver
    simulator.Solve();

    // Remove the population boundary condition
    simulator.RemoveAllCellPopulationBoundaryConditions();

    // Re-add the flat boundary
    simulator.AddCellPopulationBoundaryCondition(p_membrane_bc);

    // Save the results
    CellBasedSimulationArchiver<3,OffLatticeSimulation<3>, 3>::Save(&simulator);

    // Reporting
    TRACE("\nCompleted fill for spring length: ");
    PRINT_VARIABLE(spring_length);
    CellBasedEventHandler::Headings();
    CellBasedEventHandler::Report();

    // Reset the singletons
    SimulationTime::Destroy();
    SimulationTime::Instance()->SetStartTime(0.0);
    RandomNumberGenerator::Destroy();
}
