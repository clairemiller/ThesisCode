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

#include "KLKOdeModelCellModifier.hpp"

template<unsigned DIM>
KLKOdeModelCellModifier<DIM>::KLKOdeModelCellModifier(double heightStartSC, double expectedSCHeight, unsigned vertical)
 : AbstractCellBasedSimulationModifier<DIM,DIM>(), mHeightStartSC(heightStartSC), mExpectedSCHeight(expectedSCHeight), mVertical(vertical)
{
    assert(heightStartSC > 0.0);
    assert(expectedSCHeight > 0.0);
    assert(vertical<DIM);
}

template<unsigned DIM>
KLKOdeModelCellModifier<DIM>::~KLKOdeModelCellModifier()
{
}

template<unsigned DIM>
void KLKOdeModelCellModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void KLKOdeModelCellModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Now iterate over the population, and add data for any cells that have completed their M phase this time step
    for ( typename AbstractCellPopulation<DIM>::Iterator cell_it = rCellPopulation.Begin(); cell_it != rCellPopulation.End(); ++cell_it )
    {
        double cell_loc = (rCellPopulation.GetLocationOfCellCentre(*cell_it))[mVertical];
        double z_location = (cell_loc - mHeightStartSC) / mExpectedSCHeight; // -ve zs are dealt with in ode system code
        cell_it->GetCellData()->SetItem("ZLocation", z_location);
    }
}

template<unsigned DIM>
void KLKOdeModelCellModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void KLKOdeModelCellModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<HeightStartSC>" << mHeightStartSC << "</HeightStartSC>\n";
    *rParamsFile << "\t\t\t<ExpectedSCHeight>" << mExpectedSCHeight << "</ExpectedSCHeight>\n";
    *rParamsFile << "\t\t\t<Vertical>" << mVertical << "</Vertical>\n";

    // Call parent class method
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template<unsigned DIM>
double KLKOdeModelCellModifier<DIM>::GetHeightStartSC() const
{
    return mHeightStartSC;
}

template<unsigned DIM>
double KLKOdeModelCellModifier<DIM>::GetExpectedSCHeight() const
{
    return mExpectedSCHeight;
}

template<unsigned DIM>
unsigned KLKOdeModelCellModifier<DIM>::GetVertical() const
{
    return mVertical;
}


// Explicit instantiation
template class KLKOdeModelCellModifier<1>;
template class KLKOdeModelCellModifier<2>;
template class KLKOdeModelCellModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(KLKOdeModelCellModifier)
