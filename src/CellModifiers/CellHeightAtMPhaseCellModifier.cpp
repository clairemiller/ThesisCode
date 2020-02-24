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

#include "CellHeightAtMPhaseCellModifier.hpp"

template<unsigned DIM>
CellHeightAtMPhaseCellModifier<DIM>::CellHeightAtMPhaseCellModifier(double MPhaseLength, unsigned vertical)
 : AbstractCellBasedSimulationModifier<DIM,DIM>(), mMPhaseLength(MPhaseLength), mVertical(vertical)
{
    assert(MPhaseLength > 0);
    assert(vertical<DIM);
}

template<unsigned DIM>
CellHeightAtMPhaseCellModifier<DIM>::~CellHeightAtMPhaseCellModifier()
{
}

template<unsigned DIM>
void CellHeightAtMPhaseCellModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CellHeightAtMPhaseCellModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Now iterate over the population, and add data for any cells that have completed their M phase this time step
    double dt = SimulationTime::Instance()->GetTimeStep();
    for ( typename AbstractCellPopulation<DIM>::Iterator cell_it = rCellPopulation.Begin(); cell_it != rCellPopulation.End(); ++cell_it )
    {
        double age = cell_it->GetAge();
        if ( age > mMPhaseLength && (age-dt) < mMPhaseLength )
        {
            boost::shared_ptr<AbstractCellProliferativeType> pCellType = cell_it->GetCellProliferativeType();
            if ( pCellType->IsType<DifferentiatedCellProliferativeType>() )
            {
                c_vector<double,DIM> cell_loc = rCellPopulation.GetLocationOfCellCentre(*cell_it);
                cell_it->GetCellData()->SetItem("HeightAtDivision", cell_loc[mVertical]);
            }
        }
    }
}

template<unsigned DIM>
void CellHeightAtMPhaseCellModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // First set all cells to max double
    for ( typename AbstractCellPopulation<DIM>::Iterator cell_it = rCellPopulation.Begin(); cell_it != rCellPopulation.End(); ++cell_it )
    {
        cell_it->GetCellData()->SetItem("HeightAtDivision", DBL_MAX);
    }
    // Just in case, also check for end of division times
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CellHeightAtMPhaseCellModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MPhaseLength>" << mMPhaseLength << "</MPhaseLength>\n";
	  *rParamsFile << "\t\t\t<Vertical>" << mVertical << "</Vertical>\n";

    // Call parent class method
	  AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned DIM>
double CellHeightAtMPhaseCellModifier<DIM>::GetMPhaseLength() const
{
    return mMPhaseLength;
}

template<unsigned DIM>
void CellHeightAtMPhaseCellModifier<DIM>::SetMPhaseLength(double MPhaseLength)
{
    mMPhaseLength = MPhaseLength;
    assert(MPhaseLength > 0);
}

template<unsigned DIM>
unsigned CellHeightAtMPhaseCellModifier<DIM>::GetVertical() const
{
    return mVertical;
}

template<unsigned DIM>
void CellHeightAtMPhaseCellModifier<DIM>::SetVertical(unsigned vertical)
{
    mVertical = vertical;
    assert(vertical < DIM);
}


// Explicit instantiation
template class CellHeightAtMPhaseCellModifier<1>;
template class CellHeightAtMPhaseCellModifier<2>;
template class CellHeightAtMPhaseCellModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellHeightAtMPhaseCellModifier)
