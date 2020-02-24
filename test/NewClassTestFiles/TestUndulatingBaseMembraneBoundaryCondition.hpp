#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>
#include "AbstractCellBasedTestSuite.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

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
#include "RandomNumberGenerator.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "FlatBaseMembraneBoundaryCondition.hpp"
#include "SinusoidalBaseMembraneBoundaryCondition.hpp"
#include "ApproxGaussianBaseMembraneBoundaryCondition.hpp"
#include "UndulatingBaseMembraneAdhesionForce.hpp"
#include "AttachedCellMutationState.hpp"
#include "Debug.hpp"


class TestAbstractUndulatingBaseMembraneBoundaryCondition : public AbstractCellBasedTestSuite
{
public:
    void TestMembraneCalculations() throw(Exception)
    {

        // Create an empty cell population
        NodesOnlyMesh<2> mesh;        
        std::vector<CellPtr> cells;
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        // Construct the classes
        FlatBaseMembraneBoundaryCondition<2> flat_force(&cell_population);
        SinusoidalBaseMembraneBoundaryCondition<2> sine_force(&cell_population);
        // CLAIRETODO: Gaussian force
        // ApproxGaussianBaseMembraneBoundaryCondition<2> gauss_force; 

        // Check the normal
        c_vector<double,2> p;
        c_vector<double,2> n;
        p[0] = 7.0;
        p[1] = 0.0;
        n = flat_force.CalculateDerivativesAtPoint(p);
        TS_ASSERT_DELTA(n[0],0.0,1.0e-6);
        TS_ASSERT_DELTA(n[1],1.0,1.0e-6);
        n = sine_force.CalculateDerivativesAtPoint(p);
        TS_ASSERT_DELTA(n[0],4.0*M_PI/7.0,1e-6);
        TS_ASSERT_DELTA(n[1],1.0,1e-6);

        // CLAIRETODO: Check the membrane calculations
        c_vector<double,2> pm;
        pm = flat_force.DetermineClosestMembraneLocation(p);
        TS_ASSERT_DELTA(pm[0],7.0,1.0e-6);
        TS_ASSERT_DELTA(pm[1],0.0,1.0e-6);
    }

    void TestPeriodicSimulation() throw(Exception)
    {
        // Set up the cell population
        std::vector<Node<2>*> nodes(0);
        
        // Create the cell population
        std::vector<double> widths(1);
        widths[0] = 6.0;
        PeriodicNdNodesOnlyMesh<2> mesh(widths,true,false);        
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        // Initialise the boundary condition
        // Default period: 7
        // Should break
        TS_ASSERT_THROWS_THIS(MAKE_PTR_ARGS(SinusoidalBaseMembraneBoundaryCondition<2>, new_bc, (&cell_population)),
            "The periodic boundary width must be a multiple of membrane period in periodic dimensions.");

        // But if we change the period, this should now be fine
        MAKE_PTR_ARGS(SinusoidalBaseMembraneBoundaryCondition<2>, new_bc, (&cell_population, 1, 2.0, 6.0));
    }

    void TestAdhesion2D() throw(Exception)
    {
        // With y as vertical
        // Set up the node positions
        int n_nodes = 10;
        std::vector<Node<2>*> nodes(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            nodes[i] = new Node<2>(i,false,i,0.0);
        }
        
        // Create the cell population
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        // Initialise the boundary condition
        MAKE_PTR_ARGS(FlatBaseMembraneBoundaryCondition<2>,p_bc, (&cell_population));
        std::map<Node<2>*, c_vector<double, 2> > emptymap;
        p_bc->ImposeBoundaryCondition(emptymap);

        // Cells should all be attached
        for (NodeBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin(); 
            cell_iter != cell_population.End(); ++cell_iter)
            {
                TS_ASSERT((*cell_iter)->GetMutationState()->IsType<AttachedCellMutationState>());
            }

        // Now move them upwards and calculate the force
        for (int i = 0; i < n_nodes; i++)
        {
            cell_population.GetNode(i)->rGetModifiableLocation()[1] = ((double)i)/10.0;
        }        
        UndulatingBaseMembraneAdhesionForce<2> membrane_force;
        membrane_force.AddForceContribution(cell_population);
        
        // Check the force calculation
        double alpha_star = 50.0;
        double lambda = 7.0;
        double c = std::sqrt(0.5/lambda);
        for (int i = 0; i < n_nodes; i++)
        {
            double y = ((double)i)/10.0;
            double F = (y+c)*std::exp(-lambda*std::pow(y+c,2.0)) - c*std::exp(-lambda*(y*y+c*c));
            F *= alpha_star;
            c_vector<double,2> v_applied_force = cell_population.GetNode(i)->rGetAppliedForce();
            TS_ASSERT(v_applied_force[1] <= 0.0);
            TS_ASSERT_DELTA(v_applied_force[1], F, 1.0e-6);
            TS_ASSERT_DELTA(v_applied_force[0], 0.0, 1.0e-6);
            // Clear the force for the next check
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Change them to differentiated and make sure they're 0
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (NodeBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin(); 
            cell_iter != cell_population.End(); ++cell_iter)
            {
                (*cell_iter)->SetCellProliferativeType(p_diff_type);
            }
        membrane_force.AddForceContribution(cell_population);
        for (int i = 0; i < n_nodes; i++)
        {
            c_vector<double,2> v_applied_force = cell_population.GetNode(i)->rGetAppliedForce();
            TS_ASSERT_DELTA(v_applied_force[1], 0.0, 1.0e-6);
            TS_ASSERT_DELTA(v_applied_force[0], 0.0, 1.0e-6);
        }
        

    }

    void TestAdhesion3D() throw(Exception)
    {
        // With z as vertical
        // Set up the node positions
        int n_nodes = 10;
        std::vector<Node<3>*> nodes(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            double x = RandomNumberGenerator::Instance()->ranf();
            double y = RandomNumberGenerator::Instance()->ranf();
            nodes[i] = new Node<3>(i,false,x*5,y*10,0.0);
        }
        
        // Create the cell population
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
        NodeBasedCellPopulation<3> cell_population(mesh,cells);

        // Initialise the boundary condition
        MAKE_PTR_ARGS(FlatBaseMembraneBoundaryCondition<3>,p_bc, (&cell_population));
        std::map<Node<3>*, c_vector<double, 3> > emptymap;
        p_bc->ImposeBoundaryCondition(emptymap);

        // Cells should all be attached
        for (NodeBasedCellPopulation<3>::Iterator cell_iter = cell_population.Begin(); 
            cell_iter != cell_population.End(); ++cell_iter)
            {
                TS_ASSERT((*cell_iter)->GetMutationState()->IsType<AttachedCellMutationState>());
            }

        // Now move them upwards and calculate the force
        for (int i = 0; i < n_nodes; i++)
        {
            cell_population.GetNode(i)->rGetModifiableLocation()[2] = ((double)i)/10.0;
        }   
        double alpha_star = 50.0/0.0267131;    
        UndulatingBaseMembraneAdhesionForce<3> membrane_force;
        membrane_force.SetAdhesionParameters(alpha_star, 0.0);
        membrane_force.AddForceContribution(cell_population);
        
        // Check the force calculation
        double lambda = 7.0;
        double c = std::sqrt(0.5/lambda);
        for (int i = 0; i < n_nodes; i++)
        {
            double y = ((double)i)/10.0;
            double F = (y+c)*std::exp(-lambda*std::pow(y+c,2.0)) - c*std::exp(-lambda*(y*y+c*c));
            F *= alpha_star;
            c_vector<double,3> v_applied_force = cell_population.GetNode(i)->rGetAppliedForce();
            TS_ASSERT(v_applied_force[2] <= 0.0);
            TS_ASSERT_DELTA(v_applied_force[2], F, 1.0e-6);
            TS_ASSERT_DELTA(v_applied_force[0], 0.0, 1.0e-6);
            TS_ASSERT_DELTA(v_applied_force[1], 0.0, 1.0e-6);
            // Clear the force for the next check
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        // Change them to differentiated and make sure they're 0
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (NodeBasedCellPopulation<3>::Iterator cell_iter = cell_population.Begin(); 
            cell_iter != cell_population.End(); ++cell_iter)
            {
                (*cell_iter)->SetCellProliferativeType(p_diff_type);
            }
        membrane_force.AddForceContribution(cell_population);
        for (int i = 0; i < n_nodes; i++)
        {
            c_vector<double,3> v_applied_force = cell_population.GetNode(i)->rGetAppliedForce();
            TS_ASSERT_DELTA(v_applied_force[2], 0.0, 1.0e-6);
            TS_ASSERT_DELTA(v_applied_force[1], 0.0, 1.0e-6);
            TS_ASSERT_DELTA(v_applied_force[0], 0.0, 1.0e-6);
        }
        

    }

    void TestSinusoidalBoundaryCondition2d() throw(Exception)
    {
        // Set up the node positions
        int n_nodes = 10;
        std::vector<Node<2>*> nodes(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            nodes[i] = new Node<2>(i,false,i,2*(1+sin(2.0/7.0*M_PI*(double)i))-1.0);
        }
        
        // Create the cell population
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        // Initialise the boundary conditions
        MAKE_PTR_ARGS(SinusoidalBaseMembraneBoundaryCondition<2>, new_bc, (&cell_population));

        // Now we solve
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestUndulatingMembraneBoundaryCondition");
        simulator.SetEndTime(0.1);
        simulator.AddCellPopulationBoundaryCondition(new_bc);
        simulator.Solve();

        // Check they are all on the membrane
        // First need to get the membrane parameters
        double A = new_bc->GetAmplitude();
        double P = new_bc->GetPeriod();
        // CLAIRETODO: A couple manually calculated
        for ( int i = 0; i < n_nodes; i++ )
        {
            c_vector<double,2> cell_loc = cell_population.rGetMesh().GetNode(i)->rGetLocation();
            double yPos = A * (1.0+sin(cell_loc[0]*2*_PI/P));
            TS_ASSERT_DELTA(cell_loc[1],yPos,1e-5);
        }
    }

    void TestFlatBoundaryCondition3d() throw(Exception)
    {
        // Set up the node positions
        int n_nodes = 10;
        std::vector<Node<3>*> nodes(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            nodes[i] = new Node<3>(i,false,i,i,-((double)i)/10.0);
        }
        
        // Create the cell population
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
        NodeBasedCellPopulation<3> cell_population(mesh,cells);

        // Initialise the boundary conditions
        MAKE_PTR_ARGS(FlatBaseMembraneBoundaryCondition<3>, new_bc, (&cell_population));

        // Now we solve
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestFlatMembraneBoundaryCondition");
        simulator.SetEndTime(0.1);
        simulator.AddCellPopulationBoundaryCondition(new_bc);
        simulator.Solve();

        // Check they are all on the membrane
        for ( int i = 0; i < n_nodes; i++ )
        {
            c_vector<double,3> cell_loc = cell_population.rGetMesh().GetNode(i)->rGetLocation();
            TS_ASSERT_DELTA(cell_loc[0],(double)i,1e-5);
            TS_ASSERT_DELTA(cell_loc[1],(double)i,1e-5);
            TS_ASSERT_DELTA(cell_loc[2],0.0,1e-5);
        }
    }

    void TestInvertedSim() throw(Exception)
    {
        // With y as vertical
        // Set up the node positions
        int n_nodes = 10;
        std::vector<Node<2>*> nodes_y(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            nodes_y[i] = new Node<2>(i,false,i,0.0);
        }
        
        // Create the cell population
        NodesOnlyMesh<2> mesh_y;
        mesh_y.ConstructNodesWithoutMesh(nodes_y,1.5);
        std::vector<CellPtr> cells_y;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        cells_generator.GenerateBasicRandom(cells_y, mesh_y.GetNumNodes(), p_transit_type);
        NodeBasedCellPopulation<2> cell_population_y(mesh_y,cells_y);

        // Initialise the boundary condition
        MAKE_PTR_ARGS(SinusoidalBaseMembraneBoundaryCondition<2>,new_bc_y, (&cell_population_y));

        // Now we solve
        OffLatticeSimulation<2> simulator(cell_population_y);
        simulator.SetOutputDirectory("TestSinusoidalBaseMembraneBoundaryConditionY");
        simulator.SetEndTime(1.5);
        simulator.AddCellPopulationBoundaryCondition(new_bc_y);
        simulator.Solve();

        // Check they are all on the membrane
        for ( int i = 0; i < n_nodes; i++ )
        {
            c_vector<double,2> cell_loc = cell_population_y.rGetMesh().GetNode(i)->rGetLocation();
            double y = new_bc_y->BaseShapeFunction(cell_loc);
            TS_ASSERT_DELTA(cell_loc[1],y,1e-5);
        }

        // First reset the singletons
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(0);

        // With x as vertical
        // Set up the node positions
        std::vector<Node<2>*> nodes_x(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            nodes_x[i] = new Node<2>(i,false,0.0,i);
        }
        
        // Create the cell population
        NodesOnlyMesh<2> mesh_x;
        mesh_x.ConstructNodesWithoutMesh(nodes_x,1.5);
        std::vector<CellPtr> cells_x;
        cells_generator.GenerateBasicRandom(cells_x, mesh_x.GetNumNodes(), p_transit_type);
        NodeBasedCellPopulation<2> cell_population_x(mesh_x,cells_x);

        // Initialise the boundary condition
        MAKE_PTR_ARGS(SinusoidalBaseMembraneBoundaryCondition<2>,new_bc_x,(&cell_population_x,0));

        // Now we solve
        OffLatticeSimulation<2> simulator_x(cell_population_x);
        simulator_x.SetOutputDirectory("TestFlatBaseMembraneBoundaryConditionX");
        simulator_x.SetEndTime(1.5);
        simulator_x.AddCellPopulationBoundaryCondition(new_bc_x);
        simulator_x.Solve();

        // Check that they match
        for ( int i=0; i < n_nodes; i++ )
        {
            c_vector<double,2> cell_loc_y = cell_population_y.rGetMesh().GetNode(i)->rGetLocation();
            c_vector<double,2> cell_loc_x = cell_population_x.rGetMesh().GetNode(i)->rGetLocation();
            TS_ASSERT_DELTA(cell_loc_y[0],cell_loc_x[1],1e-10);
            TS_ASSERT_DELTA(cell_loc_y[1],cell_loc_x[0],1e-10);
        }

    }

    void TestMembraneAttachmentProperty() throw(Exception)
    {
        int n_nodes = 30;
        std::vector<Node<2>*> nodes(n_nodes);
        std::vector<c_vector<double,2> > c_loc(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            c_loc[i][0] = (i % 10); 
            c_loc[i][1] = 2.0*(1.0 + sin(c_loc[i][0]*2.0*M_PI/7.0));
            nodes[i] = new Node<2>(i,false,c_loc[i][0],c_loc[i][1]);
        }
        
        // Create the cell population
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(AttachedCellMutationState, p_property);
        // Set birth times
        for ( int i = 0; i < n_nodes; i++ )
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            // p_model->SetDimension(3);

            CellPtr p_cell(new Cell(p_property, p_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->SetBirthTime(0.0);

            cells.push_back(p_cell);
        }
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        MAKE_PTR_ARGS(SinusoidalBaseMembraneBoundaryCondition<2>,new_bc, (&cell_population));

        // Now we solve
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.AddCellPopulationBoundaryCondition(new_bc);
        simulator.SetOutputDirectory("TestAttachmentMutation");
        simulator.SetEndTime(2.0);
        simulator.Solve();

        // Check the output
        for ( int i=0; i < n_nodes; i++ )
        {
            CellPtr cell = cell_population.GetCellUsingLocationIndex(i);
            boost::shared_ptr<AbstractCellMutationState> p_mutation = cell->GetMutationState();
            TS_ASSERT(p_mutation->IsType<AttachedCellMutationState>());
            c_vector<double,2> end_loc = cell_population.rGetMesh().GetNode(i)->rGetLocation();
            TS_ASSERT_DELTA(c_loc[i][1],cell->GetCellData()->GetItem("attachment_point_y"),1e-10);
            TS_ASSERT_DELTA(c_loc[i][1],end_loc[1],5e-3);
        }
    }


    void TestSaveLoadBaseMembraneForces() throw(Exception)
    {
        // With y as vertical
        // Set up the node positions

        int n_nodes = 10;
        std::vector<Node<2>*> nodes(n_nodes);
        for (int i = 0; i < n_nodes; i++)
        {
            nodes[i] = new Node<2>(i,false,i,6.5);
        }
        
        // Create the cell population
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);
        NodeBasedCellPopulation<2> cell_population(mesh,cells);

        // Initialise the boundary condition
        MAKE_PTR_ARGS(FlatBaseMembraneBoundaryCondition<2>,p_bc, (&cell_population));

        // Now we solve
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestSaveLoadUndulatingMembraneBoundaryCondition");
        simulator.SetEndTime(1.5);
        simulator.AddCellPopulationBoundaryCondition(p_bc);        
        simulator.Solve();

        // Fix attachment points
        for (typename NodeBasedCellPopulation<2>::Iterator cell_iter = (cell_population.Begin()); 
        cell_iter != (cell_population.End()); ++cell_iter)
        {
            (*cell_iter)->GetCellData()->SetItem("attachment_point_x",0.0);
            (*cell_iter)->GetCellData()->SetItem("attachment_point_y",0.0);
        }

        // Save simulation
        CellBasedSimulationArchiver<2,OffLatticeSimulation<2>, 2>::Save(&simulator);
        
        // Reload
        OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2>,2 >::Load("TestSaveLoadUndulatingMembraneBoundaryCondition", 1.5);
        NodeBasedCellPopulation<2>& r_population = dynamic_cast<NodeBasedCellPopulation<2>& >(p_simulator->rGetCellPopulation());
        for (typename NodeBasedCellPopulation<2>::Iterator cell_iter = (r_population.Begin()); 
        cell_iter != (r_population.End()); ++cell_iter)
        {
            double ay = (*cell_iter)->GetCellData()->GetItem("attachment_point_y");
            TS_ASSERT_DELTA(ay,0.0,1e-8);
        }

        delete p_simulator;
    }


    void TestArchivingOfBaseMembraneClasses() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "MembraneBoundaryConditions.arch";
        ArchiveLocationInfo::SetMeshFilename("MembraneBoundaryConditions");

        // Create data structures to store variables to Test for equality here
        {
            // Create cell populations
            std::vector<Node<3>*> nodes(10);
            for ( unsigned i=0; i < 10; i++ )
            {
                nodes[i] = new Node<3>( i,false, 0.0,(double)(i)/10.0,i);
            }
            NodesOnlyMesh<3> mesh;
            mesh.ConstructNodesWithoutMesh(nodes,1.5);
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_stem_type);
            NodeBasedCellPopulation<3> cell_population(mesh,cells);

            // Set up boundary conditions and force classes
            FlatBaseMembraneBoundaryCondition<3> flat_bc(&cell_population,1);
            flat_bc.SetUseJiggledBottomCells(true);
            SinusoidalBaseMembraneBoundaryCondition<3> sine_bc(&cell_population);
            sine_bc.SetCurveParameters(3.0,4.0);
            UndulatingBaseMembraneAdhesionForce<3> membrane_force;  
            membrane_force.SetAdhesionParameters(30.0,1.0);              

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

                
            // Create base class pointers
            AbstractCellPopulationBoundaryCondition<3>* const p_flat_bc = &flat_bc;
            AbstractCellPopulationBoundaryCondition<3>* const p_sine_bc = &sine_bc;
            AbstractForce<3,3>* const p_membrane_force = &membrane_force;

            // Archive from base class pointers
            (*p_arch) << p_flat_bc;
            (*p_arch) << p_sine_bc;
            (*p_arch) << p_membrane_force;
        }
    
        {
            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            AbstractCellPopulationBoundaryCondition<3>* p_flat_bc;
            AbstractCellPopulationBoundaryCondition<3>* p_sine_bc;
            AbstractForce<3,3>* p_membrane_force;
            (*p_arch) >> p_flat_bc;
            (*p_arch) >> p_sine_bc;
            (*p_arch) >> p_membrane_force;

            // Cast up to appropriate class
            FlatBaseMembraneBoundaryCondition<3>* p_flat_cast = dynamic_cast<FlatBaseMembraneBoundaryCondition<3>*>(p_flat_bc);
            SinusoidalBaseMembraneBoundaryCondition<3>* p_sine_cast = dynamic_cast<SinusoidalBaseMembraneBoundaryCondition<3>*>(p_sine_bc);
            UndulatingBaseMembraneAdhesionForce<3>* p_membrane_force_cast = dynamic_cast<UndulatingBaseMembraneAdhesionForce<3>*>(p_membrane_force);
    
            // Check the casting is correct 
            TS_ASSERT(p_flat_cast);
            TS_ASSERT(p_sine_cast);
            TS_ASSERT(p_membrane_force_cast);
    
            // Check things in the data structures
            TS_ASSERT(p_flat_cast->mVert == 1);          
            TS_ASSERT(p_flat_cast->mUseJiggledBottomCells == true);

            TS_ASSERT(p_sine_cast->GetAmplitude() == 3.0);
            TS_ASSERT(p_sine_cast->GetPeriod() == 4.0);

            TS_ASSERT(p_membrane_force_cast->GetAlphaProlifCell() == 30.0);
            TS_ASSERT(p_membrane_force_cast->GetAlphaDiffCell() == 1.0);
        } 
    }
    
};