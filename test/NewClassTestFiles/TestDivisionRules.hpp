#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>

// To run in parallel
//#include "PetscSetupAndFinalize.hpp"
// When run in serial
#include "FakePetscSetup.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "FlatBaseMembraneForce.hpp"
#include "SinusoidalBaseMembraneForce.hpp"
#include "ApproxGaussianBaseMembraneForce.hpp"
#include "RotationalDivisionForce.hpp"
#include "RotationalDivisionForceOnUndulatingMembrane.hpp"
#include "VerticalToMembraneDivisionRule.hpp"

#include <boost/assign/list_of.hpp>


class TestRotationalDivisionForce : public AbstractCellBasedTestSuite
{
public:
    void DNRTestVerticalToMembraneDivisionRuleFlat() throw(Exception)
    {
        // Sinusoidal
        MAKE_PTR(FlatBaseMembraneForce<2>,p_force_flat);
        MAKE_PTR_ARGS(RotationalDivisionForceOnUndulatingMembrane<2>, p_rot_force, (p_force_flat,1.0));

        // Set up the node positions
        int n_nodes = 10;
        std::vector<Node<2>*> nodes(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            nodes[i] = new Node<2>(i,false,(double)i,0.0);
        }
        // Create the cell population
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_transit_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);
        for (int i = 0; i < n_nodes; i++)
        {
            cells[i]->SetBirthTime(-23.9);
            dynamic_cast<FixedG1GenerationalCellCycleModel* >(cells[i]->GetCellCycleModel())->SetMaxTransitGenerations(0);
        }
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        // Add the division rule
        MAKE_PTR_ARGS(VerticalToMembraneDivisionRule<2>,p_div_rule,(p_force_flat,0.5));
        cell_population.SetCentreBasedDivisionRule(p_div_rule);

        // Now we solve
        OffLatticeSimulation<2> simulator(cell_population);
        // Add the forces
        simulator.AddForce(p_force_flat);
        simulator.AddForce(p_rot_force);
        simulator.SetOutputDirectory("TestVerticalDivisionFlat");
        simulator.SetEndTime(0.2);
        simulator.Solve();
    }


    void DNRTest2DVerticalToMembraneDivisionRuleSine() throw(Exception)
    {
        // Sinusoidal
        MAKE_PTR(SinusoidalBaseMembraneForce<2>,p_force_sine);
        MAKE_PTR_ARGS(RotationalDivisionForceOnUndulatingMembrane<2>, p_rot_force, (p_force_sine,1.0));
        double period = p_force_sine->GetPeriod();
        double amplitude = p_force_sine->GetAmplitude();

        // Set up the node positions
        int n_nodes = period + 1;
        std::vector<Node<2>*> nodes(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            nodes[i] = new Node<2>(i,false,i,amplitude*(1.0+sin(2.0*(double)(i)*M_PI/period)));
        }
        // Create the cell population
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_transit_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);
        for (int i = 0; i < n_nodes; i++)
        {
            cells[i]->SetBirthTime(-23.9);
            dynamic_cast<FixedG1GenerationalCellCycleModel* >(cells[i]->GetCellCycleModel())->SetMaxTransitGenerations(0);
        }
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        // Add the division rule
        MAKE_PTR_ARGS(VerticalToMembraneDivisionRule<2>,p_div_rule,(p_force_sine,0.9));
        cell_population.SetCentreBasedDivisionRule(p_div_rule);

        // Now we solve
        OffLatticeSimulation<2> simulator(cell_population);
        // Add the forces
        simulator.AddForce(p_force_sine);
        simulator.AddForce(p_rot_force);
        simulator.SetOutputDirectory("Test2DVerticalDivisionSine");
        simulator.SetEndTime(0.2);
        simulator.Solve();
    }

    void Test3DVerticalToMembraneDivisionRuleSine() throw(Exception)
    {
        // Sinusoidal
        MAKE_PTR(SinusoidalBaseMembraneForce<3>,p_force_sine);
        MAKE_PTR_ARGS(RotationalDivisionForceOnUndulatingMembrane<3>, p_rot_force, (p_force_sine,1.0));
        double period = p_force_sine->GetPeriod();
        double amplitude = p_force_sine->GetAmplitude();

        // Set up the node positions
        int n_nodes = 3*(int)period;
        std::vector<Node<3>*> nodes(n_nodes);
        std::vector<double> y_vals =  boost::assign::list_of(0)(period/2.0)(period*5.0/4.0);
        for (int i = 0; i < (int)period; i++)
        {
            double x = (double) i;
            for (unsigned j = 0; j < 3; j++)
            {
                unsigned index = i*3 + j;
                double z = amplitude/2.0*(2.0+sin(2.0*x*M_PI/period) + sin(2.0*y_vals[j]*M_PI/period));
                nodes[index] = new Node<3>(index,false,x,y_vals[j],z);

            }
        }
        // Create the cell population
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_transit_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);
        for (int i = 0; i < n_nodes; i++)
        {
            cells[i]->SetBirthTime(-23.9);
            dynamic_cast<FixedG1GenerationalCellCycleModel* >(cells[i]->GetCellCycleModel())->SetMaxTransitGenerations(0);
        }
        NodeBasedCellPopulation<3> cell_population(mesh,cells);

        // Add the division rule
        MAKE_PTR_ARGS(VerticalToMembraneDivisionRule<3>,p_div_rule,(p_force_sine,0.1));
        cell_population.SetCentreBasedDivisionRule(p_div_rule);

        // Now we solve
        OffLatticeSimulation<3> simulator(cell_population);
        // Add the forces
        simulator.AddForce(p_force_sine);
        //simulator.AddForce(p_rot_force);
        simulator.SetOutputDirectory("Test3DVerticalDivisionSine");
        simulator.SetEndTime(0.2);
        simulator.Solve();
    }

    void DNRTestVerticalToMembraneDivisionRuleGaussian() throw(Exception)
    {
        // Sinusoidal
        MAKE_PTR(ApproxGaussianBaseMembraneForce<2>,p_force_gaussian);
        MAKE_PTR_ARGS(RotationalDivisionForceOnUndulatingMembrane<2>, p_rot_force, (p_force_gaussian,1.0));
        double period = p_force_gaussian->GetPeriod();
        double amplitude = p_force_gaussian->GetAmplitude();

        // Set up the node positions
        int n_nodes = period + 1;
        std::vector<Node<2>*> nodes(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            nodes[i] = new Node<2>(i,false,i,amplitude*std::pow( 0.5*(1.0+sin((double)i*2.0*_PI/(period))), 4.0));
        }
        // Create the cell population
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_transit_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);
        for (int i = 0; i < n_nodes; i++)
        {
            cells[i]->SetBirthTime(-23.9);
            dynamic_cast<FixedG1GenerationalCellCycleModel* >(cells[i]->GetCellCycleModel())->SetMaxTransitGenerations(0);
        }
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        // Add the division rule
        MAKE_PTR_ARGS(VerticalToMembraneDivisionRule<2>,p_div_rule,(p_force_gaussian,0.3));
        cell_population.SetCentreBasedDivisionRule(p_div_rule);

        // Now we solve
        OffLatticeSimulation<2> simulator(cell_population);
        // Add the forces
        simulator.AddForce(p_force_gaussian);
        simulator.AddForce(p_rot_force);
        simulator.SetOutputDirectory("TestVerticalDivisionGaussian");
        simulator.SetEndTime(0.2);
        simulator.Solve();
    }


};