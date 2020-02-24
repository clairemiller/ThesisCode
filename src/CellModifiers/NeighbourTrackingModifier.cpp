/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "NeighbourTrackingModifier.hpp"
#include "MeshBasedCellPopulation.hpp"

template<unsigned DIM>
NeighbourTrackingModifier<DIM>::NeighbourTrackingModifier(double neighbourhoodRadius, bool onlyAboveNeighbours)
    : AbstractCellBasedSimulationModifier<DIM>(), mNeighbourhoodRadius(neighbourhoodRadius), mOnlyAboveNeighbours(onlyAboveNeighbours)
{
}

template<unsigned DIM>
NeighbourTrackingModifier<DIM>::~NeighbourTrackingModifier()
{
}

template<unsigned DIM>
void NeighbourTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void NeighbourTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void NeighbourTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // Check we are a node based cell populatoin
    NodeBasedCellPopulation<DIM>* pPopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    if ( ! (bool) pPopulation )
    {
        EXCEPTION("Neighbour tracking modifier only implemented for node based cell populations.\n");
    }

    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the number of neighbours
        unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        std::set<unsigned> neighbours = pPopulation->GetNodesWithinNeighbourhoodRadius(index,mNeighbourhoodRadius);
        // Store the cell's neighbour count in CellData
        cell_iter->GetCellData()->SetItem("n_neighbours", (double) neighbours.size());
        // Calculate and store the above neighbours
        if ( mOnlyAboveNeighbours )
        {
            unsigned nNeighbours = 0;
            double height = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[DIM-1];
            for ( std::set<unsigned>::iterator nbr_it = neighbours.begin(); nbr_it != neighbours.end(); nbr_it++ )
            {
                CellPtr cell_neighbour = rCellPopulation.GetCellUsingLocationIndex(*nbr_it);
                double height_neighbour = rCellPopulation.GetLocationOfCellCentre(cell_neighbour)[DIM-1];
                if ( height_neighbour > height )
                {
                    nNeighbours++;
                }
            }
            cell_iter->GetCellData()->SetItem("n_above_neighbours", (double) nNeighbours);
        }

    }
}

template<unsigned DIM>
double NeighbourTrackingModifier<DIM>::GetNeighbourhoodRadius() const
{
    return mNeighbourhoodRadius;
}

template<unsigned DIM>
void NeighbourTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NeighbourhoodRadius>" << mNeighbourhoodRadius << "</NeighbourhoodRadius>\n";
    *rParamsFile << "\t\t\t<OnlyAboveNeighbours>" << mOnlyAboveNeighbours << "</OnlyAboveNeighbours>\n";
    // Now call parent method
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class NeighbourTrackingModifier<1>;
template class NeighbourTrackingModifier<2>;
template class NeighbourTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NeighbourTrackingModifier)

