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

#include "CellCountThresholderSimulationModifier.hpp"

template<unsigned DIM>
CellCountThresholderSimulationModifier<DIM>::CellCountThresholderSimulationModifier(unsigned thresholdCellCount)
 : AbstractCellBasedSimulationModifier<DIM,DIM>(), mThresholdCellCount(thresholdCellCount)
{
  assert(thresholdCellCount > 0);
}

template<unsigned DIM>
CellCountThresholderSimulationModifier<DIM>::~CellCountThresholderSimulationModifier()
{
}

template<unsigned DIM>
void CellCountThresholderSimulationModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
  unsigned nCells = rCellPopulation.GetNumNodes();
  if ( nCells > mThresholdCellCount )
  {
    EXCEPTION("Cell count threshold simulation modifier: the total number of cells (" + std::to_string(nCells) + ") has exceeded the specified threshold (" + std::to_string(mThresholdCellCount) + ") at time " + (std::to_string(SimulationTime::Instance()->GetTime()) + " hours.") );
  }
}

template<unsigned DIM>
void CellCountThresholderSimulationModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{}

template<unsigned DIM>
void CellCountThresholderSimulationModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ThresholdCellCount>" << mThresholdCellCount << "</ThresholdCellCount>\n";
	// Call parent class method
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned DIM>
unsigned CellCountThresholderSimulationModifier<DIM>::GetThresholdCellCount() const
{
  return mThresholdCellCount;
}

template<unsigned DIM>
void CellCountThresholderSimulationModifier<DIM>::SetThresholdCellCount(unsigned thresholdCellCount)
{
  mThresholdCellCount = thresholdCellCount;
  assert(thresholdCellCount > 0);
}


// Explicit instantiation
template class CellCountThresholderSimulationModifier<1>;
template class CellCountThresholderSimulationModifier<2>;
template class CellCountThresholderSimulationModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellCountThresholderSimulationModifier)
