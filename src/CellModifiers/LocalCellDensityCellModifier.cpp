/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
     this list of conditions and the following disclaimer in the documentation
     and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
     contributors may be used to endorse or promote products derived from this
     software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "LocalCellDensityCellModifier.hpp"
#include "Debug.hpp"

template<unsigned DIM>
LocalCellDensityCellModifier<DIM>::LocalCellDensityCellModifier(double gridSize)
 : AbstractCellBasedSimulationModifier<DIM,DIM>(), mGridSize(gridSize), mLocalDensityGrid()
{
        assert(gridSize > 0.0);
        std::fill(mIsDimPeriodic.begin(),mIsDimPeriodic.end(),false);
}

template<unsigned DIM>
LocalCellDensityCellModifier<DIM>::~LocalCellDensityCellModifier()
{
}

template<unsigned DIM>
unsigned& LocalCellDensityCellModifier<DIM>::AccessGridCell(c_vector<unsigned, DIM> indices)
{
    unsigned grid_index = indices[0];
    if ( DIM > 1)
    {
        grid_index += indices[1]*mDomainSize[0];
        if ( DIM == 3)
        {
            grid_index += indices[2]*mDomainSize[1]*mDomainSize[0];
        }
    }
    assert(grid_index < mLocalDensityGrid.size());
    return( mLocalDensityGrid[grid_index] );
}

template<unsigned DIM>
double LocalCellDensityCellModifier<DIM>::GetCellDensity(c_vector<double,DIM> cell_location)
{
    // Get the node location
    cell_location /= mGridSize;
    c_vector<unsigned, DIM> node_index(cell_location);
    double radius = 1.0;
    int grid_radius = 1; //(unsigned) std::round(radius/mGridSize);
    double density = 0;
    for (unsigned d = 0; d < DIM; d++)
    {
        for (int i = -grid_radius; i <= grid_radius; i++)
        {
            c_vector<int, DIM> new_index = node_index;
            new_index[d] = node_index[d] + i;
            if ( new_index[d] >= 0 && new_index[d] < (int)mDomainSize[d] )
            {
                density += (double) AccessGridCell(new_index);
            }
            else if (mIsDimPeriodic[d])
            {
                new_index[d] = (mDomainSize[d] + new_index[d]) % mDomainSize[d];
                density += (double) AccessGridCell(new_index);
            }
        }
    }
    density = density / std::pow(2.0*radius*mGridSize, (double) DIM-1);
    return(density);
}

template<unsigned DIM>
void LocalCellDensityCellModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Get the mesh
    AbstractMesh<DIM,DIM>& rMesh = rCellPopulation.rGetMesh();

    // Get the domain
    unsigned n_grid_cells = 1;
    for ( unsigned i = 0; i < DIM; i++ )
    {
        mDomainSize[i] = (unsigned) std::max( std::ceil(rCellPopulation.GetWidth(i)/mGridSize), 1.0 );
        n_grid_cells *= mDomainSize[i];
    }

    // Check that the vector is still the correct size
    if ( mLocalDensityGrid.size() != n_grid_cells )
    {
        mLocalDensityGrid.resize(n_grid_cells);
    }
    // Reset the vector to 0
    std::fill(mLocalDensityGrid.begin(), mLocalDensityGrid.end(), 0);

    // Iterate over the population and add any proliferative cells
    for ( typename AbstractCellPopulation<DIM>::Iterator cell_it = rCellPopulation.Begin(); cell_it != rCellPopulation.End(); ++cell_it )
    {
        boost::shared_ptr<AbstractCellProliferativeType> pCellType = cell_it->GetCellProliferativeType();
        if ( pCellType->IsType<StemCellProliferativeType>() || pCellType->IsType<TransitCellProliferativeType>() )
        {
            c_vector<double,DIM> cell_loc = rCellPopulation.GetLocationOfCellCentre(*cell_it);
            cell_loc /= mGridSize;
            c_vector<unsigned, DIM> cell_index(cell_loc);
            AccessGridCell(cell_index)++;
        }
    }

    // Now iterate over the population, and reset cell data for the stem cells
    for ( typename AbstractCellPopulation<DIM>::Iterator cell_it = rCellPopulation.Begin(); cell_it != rCellPopulation.End(); ++cell_it )
    {
        boost::shared_ptr<AbstractCellProliferativeType> pCellType = cell_it->GetCellProliferativeType();
        if ( pCellType->IsType<StemCellProliferativeType>() || pCellType->IsType<TransitCellProliferativeType>() )
        {
            c_vector<double,DIM> cell_loc = rCellPopulation.GetLocationOfCellCentre(*cell_it);
            cell_it->GetCellData()->SetItem("LocalDensity", GetCellDensity(cell_loc));
            // TRACE("Setting proliferative local density");
            // PRINT_VARIABLE(GetCellDensity(cell_loc));
        }
    }
}

template<unsigned DIM>
void LocalCellDensityCellModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Iterate over the population, and set the data for differentiated cells
    for ( typename AbstractCellPopulation<DIM>::Iterator cell_it = rCellPopulation.Begin(); cell_it != rCellPopulation.End(); ++cell_it )
    {
        boost::shared_ptr<AbstractCellProliferativeType> pCellType = cell_it->GetCellProliferativeType();
        if ( pCellType->IsType<DifferentiatedCellProliferativeType>() )
        {
            c_vector<double,DIM> cell_loc = rCellPopulation.GetLocationOfCellCentre(*cell_it);
            cell_it->GetCellData()->SetItem("LocalDensity", GetCellDensity(cell_loc));
            // TRACE("Setting differentiated local density");
            // PRINT_VARIABLE(GetCellDensity(cell_loc));
        }
    }
}

template<unsigned DIM>
void LocalCellDensityCellModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Check for periodicity
    // Get the mesh
    PeriodicNdNodesOnlyMesh<DIM>* pMesh = dynamic_cast<PeriodicNdNodesOnlyMesh<DIM>* >( &(rCellPopulation.rGetMesh()) );
    if (pMesh)
    {
        std::vector<unsigned> periodicDims = pMesh->GetPeriodicDimensions();
        for (unsigned i = 0; i < periodicDims.size(); i++ )
        {
            mIsDimPeriodic[i] = true;
        }
    }
    // Check that it is a nodes only mesh
    if ( !dynamic_cast<NodesOnlyMesh<DIM>* >( &(rCellPopulation.rGetMesh()) ) )
    {
        EXCEPTION("Error: Local density cell modifier only implemented for node based meshes.");
    }

    // Now run a first iteration to fill the cell data
    UpdateAtEndOfTimeStep(rCellPopulation);
    UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
void LocalCellDensityCellModifier<DIM>::SetAllCellDataToValue(AbstractCellPopulation<DIM,DIM>& rCellPopulation, double p)
{
    assert(p > 0.0);
    // Iterate over the population, and set the data for differentiated cells
    for ( typename AbstractCellPopulation<DIM>::Iterator cell_it = rCellPopulation.Begin(); cell_it != rCellPopulation.End(); ++cell_it )
    {
        c_vector<double,DIM> cell_loc = rCellPopulation.GetLocationOfCellCentre(*cell_it);
        cell_it->GetCellData()->SetItem("LocalDensity", p);
    }
}

template<unsigned DIM>
void LocalCellDensityCellModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<GridSize>" << mGridSize<< "</GridSize>\n";
	// Call parent class method
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned DIM>
double LocalCellDensityCellModifier<DIM>::GetGridSize() const
{
    return mGridSize;
}

template<unsigned DIM>
void LocalCellDensityCellModifier<DIM>::SetGridSize(double gridSize)
{
    mGridSize = gridSize;
    assert(gridSize> 0.0);
}


// Explicit instantiation
template class LocalCellDensityCellModifier<1>;
template class LocalCellDensityCellModifier<2>;
template class LocalCellDensityCellModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LocalCellDensityCellModifier)
