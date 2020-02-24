
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
#include "ParentAtBaseDivisionRule.hpp"


class TestRotationalDivisionForce : public AbstractCellBasedTestSuite
{
public:
    void TestConstructor() throw(Exception)
    {
        RotationalDivisionForce<2> rot_force(1.0);
        TS_ASSERT_EQUALS(rot_force.GetTorsionCoefficient(),1.0);
        c_vector<double,2> test = zero_vector<double>(2);
        test[1] = 1.0;
        
        TS_ASSERT_EQUALS(rot_force.GetNormalVector(test)[0],0.0);
        TS_ASSERT_EQUALS(rot_force.GetNormalVector(test)[1],1.0);
        TS_ASSERT_EQUALS(rot_force.GetGrowthDuration(),1.0);
    }

    std::pair<c_vector<double,2>,c_vector<double,2> > RunForceCalc(NodeBasedCellPopulation<2>& cell_population, 
        RotationalDivisionForce<2>& rot_force) throw(Exception)
    {
        cell_population.rGetMesh().GetNode(0)->ClearAppliedForce();
        cell_population.rGetMesh().GetNode(1)->ClearAppliedForce();        
        rot_force.AddForceContribution(cell_population);

        // Get the forces
        std::pair<c_vector<double,2>,c_vector<double,2> > applied_forces;
        applied_forces.first = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
        applied_forces.second = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

        return(applied_forces);
    }

    void TestCellPairForces() throw(Exception)
    {
        // Set up force with vertical division
        RotationalDivisionForce<2> rot_force(1.0);

        // Create the cells and mesh
        std::vector<Node<2>*> nodes(2);
        for ( int i = 0; i < 2; i++ )
        {
            nodes[i] = new Node<2>( i, false, (double)i, 0.0 );
        }
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells(2);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells,mesh.GetNumNodes());
        cells[1]->SetBirthTime(0);

        // Create the cell population and mark the springs
        NodeBasedCellPopulation<2> cell_population(mesh,cells);
        std::list<CellPtr>& r_pop_cells = cell_population.rGetCells();
        TS_ASSERT_EQUALS((int)r_pop_cells.size(),2);
        std::pair<CellPtr,CellPtr> cell_pair = std::make_pair(r_pop_cells.front(),r_pop_cells.back());
        cell_population.MarkSpring(cell_pair);

        // Add the node pair
        std::vector<std::pair<Node<2>*,Node<2>*> > node_pairs;
        node_pairs.push_back( std::make_pair(cell_population.GetNode(0),cell_population.GetNode(1)) );
        cell_population.rGetNodePairs() = node_pairs;

        // Get references to the nodes
        c_vector<double,2>& node_a_location = cell_population.rGetMesh().GetNode(0)->rGetModifiableLocation();
        c_vector<double,2>& node_b_location = cell_population.rGetMesh().GetNode(1)->rGetModifiableLocation();

        // Output forces container
        std::pair<c_vector<double,2>,c_vector<double,2> > applied_forces;

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetEndTime(1.0);
        
        /* Scenario 1: 
         *      cell A is stem, cell B is differentiated
         *      cell B is clockwise 90 degrees from desired division angle (vertical)
        */
        r_pop_cells.front()->SetCellProliferativeType(p_stem_type);
        r_pop_cells.back()->SetCellProliferativeType(p_diff_type);  
        
        applied_forces = RunForceCalc(cell_population,rot_force);
        // TS_ASSERT_DELTA(applied_forces.first[0],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[1],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[0],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[1],_PI/2.0,1.0e-10);

        // /* Scenario 2: 
        //  *      cell A is stem, cell B is differentiated
        //  *      cell B is clockwise 45 degrees from desired division angle (vertical)
        // */
        // node_a_location = scalar_vector<double>(2,1.0);
        // node_b_location = scalar_vector<double>(2,2.0);
        
        // applied_forces = RunForceCalc(cell_population,rot_force);
        // TS_ASSERT_DELTA(applied_forces.first[0],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[1],0,1.0e-10);
        // double calc_force = _PI/(4.0*2.0);
        // TS_ASSERT_DELTA(applied_forces.second[0],-calc_force,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[1],calc_force,1.0e-10);

        // /* Scenario 3: 
        //  *      cell A is stem, cell B is differentiated
        //  *      cell B is anti-clockwise 30 degrees from desired division angle (vertical)
        // */
        // node_a_location = scalar_vector<double>(2,0.0);
        // node_b_location[0] = cos(4.0*_PI/6.0);
        // node_b_location[1] = sin(4.0*_PI/6.0);
        
        // applied_forces = RunForceCalc(cell_population,rot_force);
        // TS_ASSERT_DELTA(applied_forces.first[0],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[1],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[0],(_PI/6.0)*(sqrt(3.0)/2.0),1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[1],(_PI/6.0)*0.5,1.0e-10);

        // /* Scenario 4: 
        //  *      cell A is stem, cell B is differentiated
        //  *      cell B is anti-clockwise 150 degrees from desired division angle (vertical)
        // */

        // node_a_location = scalar_vector<double>(2,0.0);
        // node_b_location[0] = cos(4.0*_PI/3.0);
        // node_b_location[1] = sin(4.0*_PI/3.0);
        
        // applied_forces = RunForceCalc(cell_population,rot_force);
        // TS_ASSERT_DELTA(applied_forces.first[0],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[1],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[0],-(5.0*_PI/6.0)*(sqrt(3.0)/2.0),1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[1],(5.0*_PI/6.0)*0.5,1.0e-10);

        // /* Scenario 5: 
        //  *      cell A is stem, cell B is differentiated
        //  *      cell B is clockwise 120 degrees from desired division angle (vertical)
        // */
        // node_b_location[0] = cos(-_PI/6.0);
        // node_b_location[1] = sin(-_PI/6.0);
        
        // applied_forces = RunForceCalc(cell_population,rot_force);
        // TS_ASSERT_DELTA(applied_forces.first[0],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[1],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[0],(2.0*_PI/3.0)*0.5,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[1],(2.0*_PI/3.0)*(sqrt(3.0)/2.0),1.0e-10);

        // /* Scenario 6: 
        //  *      cell B is stem, cell A is differentiated
        //  *      cell A is clockwise 45 degrees from desired division angle (vertical)
        // */
        // r_pop_cells.front()->SetCellProliferativeType(p_diff_type);
        // r_pop_cells.back()->SetCellProliferativeType(p_stem_type);
        // node_a_location = scalar_vector<double>(2,1.0);
        // node_b_location = scalar_vector<double>(2,0.0);
        
        // applied_forces = RunForceCalc(cell_population,rot_force);
        // TS_ASSERT_DELTA(applied_forces.second[0],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[1],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[0],-_PI/(4.0*2.0),1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[1],_PI/(4.0*2.0),1.0e-10);

        // /* Scenario 7: 
        //  *      cell B is stem, cell A is differentiated
        //  *      cell A is anti-clockwise 30 degrees from desired division angle (vertical)
        // */
        // node_a_location[0] = cos(4.0*_PI/6.0);
        // node_a_location[1] = sin(4.0*_PI/6.0);
        
        // applied_forces = RunForceCalc(cell_population,rot_force);
        // TS_ASSERT_DELTA(applied_forces.second[0],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[1],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[0],(_PI/6.0)*(sqrt(3.0)/2.0),1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[1],(_PI/6.0)*0.5,1.0e-10);

        // /* Scenario 8: 
        //  *      cell B is stem, cell A is differentiated
        //  *      cell A is anti-clockwise 150 degrees from desired division angle (vertical)
        // */
        // node_a_location[0] = cos(4.0*_PI/3.0);
        // node_a_location[1] = sin(4.0*_PI/3.0);
        
        // applied_forces = RunForceCalc(cell_population,rot_force);
        // TS_ASSERT_DELTA(applied_forces.second[0],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[1],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[0],-(5.0*_PI/6.0)*(sqrt(3.0)/2.0),1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[1],(5.0*_PI/6.0)*0.5,1.0e-10);

        // /* Scenario 9: 
        //  *      cell B is stem, cell A is differentiated
        //  *      cell A is clockwise 120 degrees from desired division angle (vertical)
        // */
        // node_a_location[0] = cos(-_PI/6.0);
        // node_a_location[1] = sin(-_PI/6.0);
        // node_b_location = scalar_vector<double>(2,0.0);
        
        // applied_forces = RunForceCalc(cell_population,rot_force);
        // TS_ASSERT_DELTA(applied_forces.second[0],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.second[1],0,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[0],(2.0*_PI/3.0)*0.5,1.0e-10);
        // TS_ASSERT_DELTA(applied_forces.first[1],(2.0*_PI/3.0)*(sqrt(3.0)/2.0),1.0e-10);
    }
};

class TestRotationalDivisionForceOnUndulatingMembrane : public AbstractCellBasedTestSuite
{
public:
    // void DNRTestDivisionAngleCalculation() throw(Exception)
    // {
    //     // First easy check: the flat membrane
    //     MAKE_PTR(FlatBaseMembraneForce<2>, p_force_flat);
    //     RotationalDivisionForceOnUndulatingMembrane<2> rot_force_flat(p_force_flat,1.0);
    //     c_vector<double, 2> p = zero_vector<double>(2);
    //     p[0] = 12.3;
    //     TS_ASSERT_DELTA(rot_force_flat.GetDivisionAngle(p),M_PI/2.0,1.0e-6);


    //     // Now the sinusoidal
    //     MAKE_PTR(SinusoidalBaseMembraneForce<2>,p_force_sine);
    //     RotationalDivisionForceOnUndulatingMembrane<2> rot_force_sine(p_force_sine,1.0);
    //     p[0] = 3.0*7.0/4.0;
    //     TS_ASSERT_DELTA(rot_force_sine.GetDivisionAngle(p),M_PI/2.0,1.0e-6);
    // }

    void TestNewlyDivideCellFlat() throw(Exception)
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
        MAKE_PTR(ParentAtBaseDivisionRule<2>,p_div_rule);
        cell_population.SetCentreBasedDivisionRule(p_div_rule);

        // Now we solve
        OffLatticeSimulation<2> simulator(cell_population);
        // Add the forces
        simulator.AddForce(p_force_flat);
        simulator.AddForce(p_rot_force);
        simulator.SetOutputDirectory("TestRotationalDivisionForceFlat");
        simulator.SetEndTime(1.2);
        simulator.Solve();
    }

    void TestNewlyDivideCellSine() throw(Exception)
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
        MAKE_PTR(ParentAtBaseDivisionRule<2>,p_div_rule);
        cell_population.SetCentreBasedDivisionRule(p_div_rule);

        // Now we solve
        OffLatticeSimulation<2> simulator(cell_population);
        // Add the forces
        simulator.AddForce(p_force_sine);
        simulator.AddForce(p_rot_force);
        simulator.SetOutputDirectory("TestRotationalDivisionForceSine");
        simulator.SetEndTime(1.2);
        simulator.Solve();
    }

    void TestNewlyDivideCellGaussian() throw(Exception)
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
        MAKE_PTR(ParentAtBaseDivisionRule<2>,p_div_rule);
        cell_population.SetCentreBasedDivisionRule(p_div_rule);

        // Now we solve
        OffLatticeSimulation<2> simulator(cell_population);
        // Add the forces
        simulator.AddForce(p_force_gaussian);
        simulator.AddForce(p_rot_force);
        simulator.SetOutputDirectory("TestRotationalDivisionForceGaussian");
        simulator.SetEndTime(1.2);
        simulator.Solve();
    }

    void TestNewlyDivideCellFlat3d() throw(Exception)
    {
        // Sinusoidal
        MAKE_PTR(FlatBaseMembraneForce<3>,p_force_flat);
        MAKE_PTR_ARGS(RotationalDivisionForceOnUndulatingMembrane<3>, p_rot_force, (p_force_flat,1.0));

        // Set up the node positions
        int n_nodes = 10*10;
        std::vector<Node<3>*> nodes(n_nodes);
        for (int i = 0; i < 10; i++)
        {
            for (int j = 0; j < 10; j++)
            {
                nodes[i*10+j] = new Node<3>(i*10+j,false,(double)j,(double)i,0.0);
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
        MAKE_PTR(ParentAtBaseDivisionRule<3>,p_div_rule);
        cell_population.SetCentreBasedDivisionRule(p_div_rule);

        // Now we solve
        OffLatticeSimulation<3> simulator(cell_population);
        // Add the forces
        simulator.AddForce(p_force_flat);
        simulator.AddForce(p_rot_force);
        simulator.SetOutputDirectory("Test3dRotationalDivisionForceFlat");
        simulator.SetEndTime(1.2);
        simulator.Solve();
    }

    void TestNewlyDivideCellSine3d() throw(Exception)
    {
        // Sinusoidal
        MAKE_PTR(SinusoidalBaseMembraneForce<3>,p_force_sine);
        MAKE_PTR_ARGS(RotationalDivisionForceOnUndulatingMembrane<3>, p_rot_force, (p_force_sine,1.0));
        double period = p_force_sine->GetPeriod();
        double amplitude = p_force_sine->GetAmplitude();

        // Set up the node positions
        int n_nodes = (int)(period) * (int)(period);
        std::vector<Node<3>*> nodes(n_nodes);
        for (int i = 0; i < (int)(period); i++)
        {
            for (int j = 0; j < (int)(period); j++)
            {
                double z = amplitude*(1.0+sin((double)(i+j)*M_PI/period) * cos((double)(i-j)*M_PI/period));
                nodes[i*(int)(period)+j] = new Node<3>(i*(int)(period)+j,false,(double)i,(double)j,z);
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
        MAKE_PTR(ParentAtBaseDivisionRule<3>,p_div_rule);
        cell_population.SetCentreBasedDivisionRule(p_div_rule);

        // Now we solve
        OffLatticeSimulation<3> simulator(cell_population);
        // Add the forces
        simulator.AddForce(p_force_sine);
        simulator.AddForce(p_rot_force);
        simulator.SetOutputDirectory("Test3dRotationalDivisionForceSine");
        simulator.SetEndTime(1.2);
        simulator.Solve();
    }
};