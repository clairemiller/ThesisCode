#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>
#include "AbstractCellBasedTestSuite.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "TissueHeightWriter.hpp"
#include "TopOfTissueCellMutationState.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"

#include "Debug.hpp"
#include "FakePetscSetup.hpp"


class TestMyCellWriters : public AbstractCellBasedTestSuite
{
public:
    void TestTissueHeightWriter() throw(Exception)
    {   
        // Set the node positions
        std::vector<Node<2>*> nodes(4);
        nodes[0] = new Node<2>( 0, false, 0.0, 1.0 );
        nodes[1] = new Node<2>( 1, false, 10.0, 5.0);
        nodes[2] = new Node<2>( 2, false, 1.0, 4.0);
        nodes[3] = new Node<2>( 3, false, 3.0, 8.0);

        // Construct the mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes,1.5);

        // Create the cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,mesh.GetNumNodes(), p_diff_type); 

        // Change nodes 2 and 3 to top of tissue mutation state
        MAKE_PTR(TopOfTissueCellMutationState, p_state);
        cells[2]->SetMutationState(p_state);
        cells[3]->SetMutationState(p_state);

        // Create the population
        NodeBasedCellPopulation<2> cell_population(mesh,cells);
        cell_population.SetAbsoluteMovementThreshold(1.5);

        // Create an output directory for the writer
        std::string output_directory = "TestTissueHeightWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a TissueHeightWriter and generate output
        TissueHeightWriter<2> height_writer;
        height_writer.OpenOutputFile(output_file_handler);
        height_writer.WriteTimeStamp();
        height_writer.Visit(&cell_population);
        height_writer.WriteNewline();
        height_writer.CloseFile();

        // Output should say:
        // 0	0	3	5	6	6	8
    }
};