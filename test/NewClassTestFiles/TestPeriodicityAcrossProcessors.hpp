

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>

// To run in parallel
#include "PetscSetupAndFinalize.hpp"
// When run in serial
//#include "FakePetscSetup.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellsGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellAncestorWriter.hpp"

// Includes to use command line options
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

// Claire's new classes
#include "PeriodicNdNodesOnlyMesh.hpp"
#include "FixedDirectionCentreBasedDivisionRule.hpp"


class Test3DSimulation : public AbstractCellBasedTestSuite
{
private:
    std::vector<Node<3>*> GenerateMesh( unsigned nx, unsigned ny, unsigned nz )
    {
        std::vector<Node<3>*> nodes(nx*ny*nz);
        for ( unsigned k = 0; k < nz; k++ )
        {
            for ( unsigned j = 0; j < ny; j++ )
            {
                for ( unsigned i = 0; i < nx; i++ )
                {
                    double x = (double)i + 0.5*(double)(j%2) + 0.5*(double)(k%2) - (double)( ((j+2)*k)%2);
                    double y = sqrt(3.0)/2.0 * (double)j;
                    double z = sqrt(3.0)/2.0 * (double)k;
                    nodes[ k*nx*ny + j*nx + i ] = new Node<3>(i+j*nx+k*nx*ny, false, x, y, z );
                }
            }
        }
        return nodes;
    }

public:
    void TestPeriodicXZ() throw(Exception)
    {
        // Set up the node positions
	    std::vector<Node<3>*> nodes = GenerateMesh(6,1,16);

	    // Create the mesh
        std::vector<double> periodic_width(2);
        periodic_width[0] = 6.0;
        periodic_width[1] = 16.0;
	    PeriodicNdNodesOnlyMesh<3> mesh(periodic_width, true, false, true); 
	    mesh.ConstructNodesWithoutMesh(nodes,1.0);

        // Create cells
        std::vector<CellPtr> cells;
        PeriodicNdNodesOnlyMesh<3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
        //PeriodicNdNodesOnlyMesh<3>::NodeIterator node_iter_end = mesh.GetNodeIteratorEnd();
        for (unsigned i=0; i<mesh.GetNumNodes(); i++, ++node_iter)
        {
            FixedG1GenerationalCellCycleModel* p_cell_cycle_model = new FixedG1GenerationalCellCycleModel;
            p_cell_cycle_model->SetDimension(3);

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

            c_vector<double,3> node_location = (*node_iter).rGetLocation();
            double birth_time = - ( 0.2*(1+node_location[1])*(1+node_location[2]) );

            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create a node-based cell population
        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);
        c_vector<double,3> division_direction = zero_vector<double>(3);
        division_direction(1) = 1.0;
        MAKE_PTR_ARGS(FixedDirectionCentreBasedDivisionRule<3>, p_div_rule, (division_direction));
        node_based_cell_population.SetCentreBasedDivisionRule(p_div_rule);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(node_based_cell_population);
        char foldername[8]; 
        sprintf(foldername, "%02iProcs", PetscTools::GetNumProcs());
        simulator.SetOutputDirectory("TestPeriodicAcrossProcs"+std::string(foldername));

        // Run for long enough to see the periodic bounday influencing the cells
        simulator.SetSamplingTimestepMultiple(3);
        simulator.SetEndTime(30.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();
    }
};