
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>
#include "AbstractCellBasedTestSuite.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "DifferentiatedCellProliferativeType.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"

#include "KLKSrnModel.hpp"
#include "KLKOdeSystem.hpp"
#include "KLKOdeModelCellModifier.hpp"
#include "KLKDrivenPalssonAdhesionForce.hpp"
#include "AttachedCellMutationState.hpp"

#include "Debug.hpp"
#include "FakePetscSetup.hpp"

class TestKLKSystem : public AbstractCellBasedTestSuite
{
    double CalculateSteadyState(double z, double eT, double iT)
    {
        double pH, kp3, km3, a, b, c, e_ana;
        pH = 6.8482 - 0.3765*z - 5.1663*z*z + 3.1792*z*z*z;
        kp3 = (15.6*pH-58.5)*std::pow(10.0,7.0);
        km3 = 6.90e6*exp(-3*pH);
        a = kp3;
        b = (iT/eT-1)*eT*kp3+km3;
        c = -km3*eT;
        e_ana = -0.5*b/a + 0.5*sqrt(b*b-4*a*c)/a;
        return(e_ana/eT);
    }
public:
    void TestSteadyState() throw(Exception)
    {
        /* Test that we can construct a KLKSrnModel object: */
        TS_ASSERT_THROWS_NOTHING(KLKSrnModel srn_model);
        /* Now we construct and initialise a cell with a KLKSrnModel.*/
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel();
        double eT = 0.5e-9;
        // Set the s0 really low so we can determine steady state
        KLKSrnModel* p_srn_model = new KLKSrnModel(1e-18,eT,eT);
        // Initialise cell
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model, p_srn_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->InitialiseCellCycleModel();
        // Initialise required cell data
        double z = 0.3;
        p_cell->GetCellData()->SetItem("ZLocation",z);
        p_cell->InitialiseSrnModel();

        /* Now increment time and check the ODE in KLKSrnModel is solved correctly. */
        double end_time = 24.0;
        unsigned num_steps = 1e5;
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();

            double current_time = SimulationTime::Instance()->GetTime();

            /* Check that the ODE system is solved correctly */
            p_srn_model->SimulateToCurrentTime();
        }

        // Test converged to steady state
        double e_ana = CalculateSteadyState(z,eT,eT);
        double e_num = p_srn_model->GetKLKEnzyme();
        TS_ASSERT_DELTA(e_ana,e_num,1.0e-3);

        // Check that at higher z the steady state is higher
        double new_z = 0.9;
        assert(new_z > z);
        TS_ASSERT_LESS_THAN(z,new_z);
        p_cell->GetCellData()->SetItem("ZLocation",new_z);
        SimulationTime::Destroy();
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);
        for (unsigned i=0; i<num_steps; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();

            double current_time = SimulationTime::Instance()->GetTime();

            /* Check that the ODE system is solved correctly */
            p_srn_model->SimulateToCurrentTime();
        }
        double e_num_new = p_srn_model->GetKLKEnzyme();
        TS_ASSERT_LESS_THAN(e_num, e_num_new);
        e_ana = CalculateSteadyState(new_z,eT,eT);
        TS_ASSERT_DELTA(e_ana,e_num_new,1.0e-3);
    }

    void TestInitialConditions() throw(Exception)
    {
        /* Now we construct and initialise a cell with a KLKSrnModel.*/
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel();
        
        // Change the substrate concentrations
        double eT = 0.5e-9;
        double iT = 0.1e-9;
        double s0 = 1.0e-6;
        KLKSrnModel* p_srn_model = new KLKSrnModel(s0,eT,iT);

        // Initialise cell
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model, p_srn_model));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->InitialiseCellCycleModel();
        // Initialise required cell data
        double z = 0.4;
        p_cell->GetCellData()->SetItem("ZLocation",z);
        p_cell->InitialiseSrnModel();

        // Test that it has the expected initial condition
        double e_num = p_srn_model->GetKLKEnzyme();
        double e_ic = (eT-iT)/eT;
        double k2 = 10.97e3;
        double kp1 = 23.87e7;
        e_ic = k2*e_ic/(kp1*s0+k2);
        TS_ASSERT_DELTA(e_num, e_ic, 1.0e-4)

        // Now create another cell with different initial conditions
        UniformG1GenerationalCellCycleModel* p_cell_cycle_model_2 = new UniformG1GenerationalCellCycleModel();
        
        // Change the substrate concentrations
        double iT_2 = 0.4e-9;
        KLKSrnModel* p_srn_model_2 = new KLKSrnModel(s0,eT,iT_2);

        // Initialise cell
        CellPtr p_cell_2(new Cell(p_state, p_cell_cycle_model_2, p_srn_model_2));
        p_cell_2->SetCellProliferativeType(p_diff_type);
        p_cell_2->InitialiseCellCycleModel();
        // Initialise required cell data
        p_cell_2->GetCellData()->SetItem("ZLocation",z);
        p_cell_2->InitialiseSrnModel();

        /* Now increment time and check the total amount of enzyme is correct. */
        double end_time = 100.0;
        unsigned num_steps = 1e5;
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();

            double current_time = SimulationTime::Instance()->GetTime();

            /* Check that the ODE system is solved correctly */
            p_srn_model->SimulateToCurrentTime();
            p_srn_model_2->SimulateToCurrentTime();
        }

        e_num = p_srn_model->GetKLKEnzyme();
        double e_ana = CalculateSteadyState(z,eT,iT);
        TS_ASSERT_DELTA(e_ana,e_num,1.0e-4);
        e_num = p_srn_model_2->GetKLKEnzyme();
        e_ana = CalculateSteadyState(z,eT,iT_2);
        TS_ASSERT_DELTA(e_ana,e_num,1.0e-4);
    }

    void TestDifferentCellTypesAndCellDivision() throw(Exception)
    {
        /* Construct and initialise 2 cell types.*/
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        FixedG1GenerationalCellCycleModel* p_cell_cycle_model_diff = new FixedG1GenerationalCellCycleModel();
        FixedG1GenerationalCellCycleModel* p_cell_cycle_model_stem = new FixedG1GenerationalCellCycleModel();
        KLKSrnModel* p_srn_model_diff = new KLKSrnModel(1e-6,0.5e-9,0.2e-9);
        KLKSrnModel* p_srn_model_stem = new KLKSrnModel(1e-6,0.5e-9,0.2e-9);
        
        CellPtr p_cell_diff(new Cell(p_state, p_cell_cycle_model_diff, p_srn_model_diff));
        p_cell_diff->SetCellProliferativeType(p_diff_type);
        p_cell_diff->InitialiseCellCycleModel();
        p_cell_diff->SetBirthTime(0.0);
        
        CellPtr p_cell_stem(new Cell(p_state, p_cell_cycle_model_stem, p_srn_model_stem));
        p_cell_stem->SetCellProliferativeType(p_stem_type);
        p_cell_stem->InitialiseCellCycleModel();
        p_cell_stem->SetBirthTime(-24.0);
        
        // Initialise required cell data
        p_cell_diff->GetCellData()->SetItem("ZLocation",0.0);
        p_cell_diff->InitialiseSrnModel();
        p_cell_stem->GetCellData()->SetItem("ZLocation",0.0);
        p_cell_stem->InitialiseSrnModel();

        // Set up the simulation time in order to run the division and solve
        double end_time = 2.0;
        unsigned num_steps = 240;
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        // Divide cell
        SimulationTime::Instance()->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_stem->ReadyToDivide(), true);
        CellPtr p_daughter_cell = p_cell_stem->Divide();
        KLKSrnModel* p_daughter_srn_model = dynamic_cast<KLKSrnModel*>(p_daughter_cell->GetSrnModel());
        assert(p_daughter_srn_model);
        p_cell_stem->SetBirthTime(24.0);

        /* Now increment time and check the ODE in KLKSrnModel is solved correctly. */
        for (unsigned i=0; i<1; i++) // (num_steps-1); i++)
        {
            // Move the cells
            double new_loc = std::min(1.0, (double)i * 0.05); 
            p_cell_diff->GetCellData()->SetItem("ZLocation", new_loc);
            p_cell_stem->GetCellData()->SetItem("ZLocation", new_loc);
            p_daughter_cell->GetCellData()->SetItem("ZLocation", new_loc);

            // Now move the ode forward
            SimulationTime::Instance()->IncrementTimeOneStep();
            double current_time = SimulationTime::Instance()->GetTime();

            /* Check that the ODE system is solved correctly */
            p_srn_model_diff->SimulateToCurrentTime();
            p_srn_model_stem->SimulateToCurrentTime();
            p_daughter_srn_model->SimulateToCurrentTime();

            TS_ASSERT_DELTA(p_srn_model_diff->GetKLKEnzyme(), p_srn_model_stem->GetKLKEnzyme(), 1.0e-4);
            TS_ASSERT_DELTA(p_srn_model_diff->GetKLKEnzyme(), p_daughter_srn_model->GetKLKEnzyme(), 1.0e-4);
        }
    }

    void TestWithModifier() throw(Exception)
    {
        // Set up the cells
        unsigned n_nodes = 1;
        std::vector<Node<2>*> nodes(n_nodes);
        for ( unsigned i=0; i < n_nodes; i++ )
        {
            nodes[i] = new Node<2>( i,false, 0, (double)(2*i)+1);
        }
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);
        for (unsigned i = 0; i < n_nodes; i++)
        {
            KLKSrnModel* p_srn_model = new KLKSrnModel;
            cells[i]->SetSrnModel(p_srn_model);
        }
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up simulator
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestKLKSubcellularNetwork");
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(100.0);

        // Add the modifier
        MAKE_PTR_ARGS(KLKOdeModelCellModifier<2>, p_cell_modifier, (10.0,10.0));
        simulator.AddSimulationModifier(p_cell_modifier);

        // Solve
        simulator.Solve();

        // Check the ode state variables
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End(); ++cell_iter)
        {
            unsigned id = (*cell_iter)->GetCellId();
            
            KLKSrnModel* p_srn_model = dynamic_cast<KLKSrnModel*>((*cell_iter)->GetSrnModel());
            assert(p_srn_model);

            if ( id < 5)
            {
                TS_ASSERT_DELTA(p_srn_model->GetKLKEnzyme(), 0.0, 1.0e-3);
            }
            else
            {
                double z = (*cell_iter)->GetCellData()->GetItem("ZLocation");
                TS_ASSERT_EQUALS(z,(double)(2*id-9)/10.0);
                double s_ana = CalculateSteadyState(((double)(2*id-9))/10.0, 0.1e-9, 0.1e-9);
                TS_ASSERT_DELTA(p_srn_model->GetKLKEnzyme(), s_ana, 1.0e-3);
            }
            
        }
    }

    void TestForce() throw(Exception)
    {
        MAKE_PTR_ARGS(KLKDrivenPalssonAdhesionForce<2>,p_cell_force,(1.0));
        assert(p_cell_force);

        // Set the node positions
        std::vector<Node<2>*> nodes(2);
        nodes[0] = new Node<2>( 0, false, 0.0, 0.0 );
        nodes[1] = new Node<2>( 1, false, 0.0, 1.15 );

        // Construct the mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create the cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_stem_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_stem_type); 

        // Assign the adhesion level
        for (unsigned i = 0; i < 2; i++)
        {
            cells[i]->GetCellData()->SetItem("AdhesiveProteinLevel",0.5*(double)(i+1));
        }

        // Create the population
        NodeBasedCellPopulation<2> cell_population(mesh,cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        // Add the adhesion force
        // Set up simulator
        OffLatticeSimulation<2> simulator(cell_population);

        // Add the force
        simulator.AddForce(p_cell_force);

        // Test the member scaling function
        CellPtr pCell0 = *(cell_population.Begin());
        double scale_factor = p_cell_force->ScalingFunction(pCell0);
        TS_ASSERT_DELTA(scale_factor, 0.5,1e-6);

        // Test the force calculation
        c_vector<double,2> force = p_cell_force->CalculateForceBetweenNodes(0,1,cell_population);
        TS_ASSERT_DELTA(force[0],0.0,1.0e-6);
        TS_ASSERT_DELTA(force[1],0.75,1.0e-3);
    }

    void TestArchiveOdeSystem() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "KLK_ode.arch";
        {
            std::vector<double> state_variables;
            state_variables.push_back(0.1);
            state_variables.push_back(0.2);
            state_variables.push_back(0.3);
            state_variables.push_back(0.4);
            state_variables.push_back(0.5);

            KLKOdeSystem ode_system(state_variables);
            std::vector<double> initial_conditions = ode_system.GetInitialConditions();

            // These are the initial conditions hard-coded in the constructor.
            TS_ASSERT_EQUALS(initial_conditions.size(), 5u);
            TS_ASSERT_DELTA(initial_conditions[0], 1.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[1], 0.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[2], 0.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[3], 0.0, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[4], 1.0, 1e-6);
            ode_system.SetParameter("z_location", 0.8);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem* const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }
        {
            // Create an input archive and restore the ode system from the archive
            AbstractOdeSystem* p_ode_system;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_ode_system;
            TS_ASSERT_DELTA(p_ode_system->rGetConstStateVariables()[1], 0.2, 1e-6);
        }
    }

    void TestArchiveCell() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "KLK_srn_model.arch";
        double end_time = 40.0;
        unsigned num_steps = 2e6;

        /* Create an output archive. */
        {
            /* Destroy the current instance of {{{SimulationTime}}} and create another instance.
             * Set the start time, end time and number of time steps. */
            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

            /* Create a cell with associated srn and cell-cycle model. */
            UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel();
            KLKSrnModel* p_srn_model = new KLKSrnModel(1.0e-6,0.5e-9,0.8e-9);
            
            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->GetCellData()->SetItem("ZLocation",0.35);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();

            /* Move forward two time steps. */
            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();
            /* Solve the SRN. */
            p_srn_model->SimulateToCurrentTime();

            /* Now archive the cell-cycle model through its cell. */
            CellPtr const p_const_cell = p_cell;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_const_cell;
        }

        /* Now create an input archive. Begin by again destroying the current
         * instance of {{{SimulationTime}}} and creating another instance. Set
         * the start time, end time and number of time steps. note that this is
         * overwritten when you load the archive.
         */
        {
            SimulationTime::Destroy();
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), 0.0, 1e-4);

            /* Create a pointer to a cell. */
            CellPtr p_cell;

            /* Create an input archive and restore the cell from the archive. */
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_cell;

            /* Run to steady state */
            KLKSrnModel* p_srn_model = dynamic_cast<KLKSrnModel*>(p_cell->GetSrnModel());
            assert(p_srn_model);
            for (unsigned i=2; i<num_steps; i++)
            {
                SimulationTime::Instance()->IncrementTimeOneStep();
                double current_time = SimulationTime::Instance()->GetTime();
                p_srn_model->SimulateToCurrentTime();
            }
            TS_ASSERT_DELTA(p_simulation_time->GetTime(), end_time, 1e-4);
            /* Check it's moved on OK */
            double s_ana = CalculateSteadyState(0.35, 0.5e-9, 0.8e-9);
            TS_ASSERT_DELTA(p_srn_model->GetKLKEnzyme(), s_ana, 1.0e-3);
        }
    }

    void TestArchiveModifier() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "klk_modifier.arch";
        {

            KLKOdeModelCellModifier<2> cell_modifier = KLKOdeModelCellModifier<2>(10.0,10.0);
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractCellBasedSimulationModifier<2>* const p_const_cell_modifier = &cell_modifier;
            output_arch << p_const_cell_modifier;
        }
        {
            // Create an input archive and restore the ode system from the archive
            AbstractCellBasedSimulationModifier<2>* p_cell_modifier;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_cell_modifier;

            // Check that it is correctly loaded
            std::vector<Node<2>*> nodes(1);
            nodes[0] = new Node<2>(0,false,0,12.0);
            std::vector<CellPtr> cells;
            CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            cells_generator.GenerateBasicRandom(cells, 1, p_diff_type);
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes,1.5);
            NodeBasedCellPopulation<2> cell_population(mesh, cells);
            p_cell_modifier->UpdateAtEndOfTimeStep(cell_population);
            double stored_z_loc = (*cell_population.Begin())->GetCellData()->GetItem("ZLocation");
            TS_ASSERT_DELTA(stored_z_loc,0.2, 1.0e-6);
        }
    }

    void TestArchiveForce() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "klk_modifier.arch";
        {

            KLKDrivenPalssonAdhesionForce<2> force = KLKDrivenPalssonAdhesionForce<2>(1.0);
            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractForce<2>* const p_const_force = &force;
            output_arch << p_const_force;
        }
        {
            // Create an input archive and restore the ode system from the archive
            AbstractForce<2>* p_force;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_force;

            // Cast back to the KLK force
            KLKDrivenPalssonAdhesionForce<2>* p_klk_force = dynamic_cast<KLKDrivenPalssonAdhesionForce<2>*>(p_force);
            assert(p_klk_force);

            // Check the parameters have been stored
            TS_ASSERT_DELTA(p_klk_force->GetPeakForce(), 1.0, 1.0e-6);
        }
    }
};