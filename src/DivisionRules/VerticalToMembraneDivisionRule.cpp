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

#include "VerticalToMembraneDivisionRule.hpp"
#include "Debug.hpp"

template<unsigned DIM>
VerticalToMembraneDivisionRule<DIM>::VerticalToMembraneDivisionRule(boost::shared_ptr<AbstractUndulatingBaseMembraneForce<DIM> > pMembrane, double divisionSpringLength)
: AbstractCentreBasedDivisionRule<DIM,DIM>(), mDivisionSpringLength(divisionSpringLength)
{
    assert(divisionSpringLength > 0.0);
    assert(divisionSpringLength <= 1.0);
    mpMembrane = pMembrane;
    mVertical = pMembrane->GetVerticalDirection();
    assert(DIM>1);
}

template<unsigned DIM>
std::pair<c_vector<double, DIM>, c_vector<double, DIM> > VerticalToMembraneDivisionRule<DIM>::CalculateCellDivisionVector(
    CellPtr pParentCell,
    AbstractCentreBasedCellPopulation<DIM>& rCellPopulation)
{
    // Get the normal vector
    c_vector<double, DIM> parent_position = rCellPopulation.GetLocationOfCellCentre(pParentCell);
    c_vector<double, DIM> division_vector = mpMembrane->CalculateDerivativesAtPoint(parent_position);
    // Convert to unit vector and negate the non-vertical for correct direction
    double vec_length = norm_2(division_vector);
    division_vector *= -1.0*mDivisionSpringLength/vec_length;
    division_vector[mVertical] *= -1.0;
    // Determine daughter position
    c_vector<double, DIM> daughter_position = parent_position + division_vector;
    std::pair<c_vector<double, DIM>, c_vector<double, DIM> > positions(parent_position, daughter_position);

    return positions;
}

template<unsigned DIM>
const boost::shared_ptr<AbstractUndulatingBaseMembraneForce<DIM> > VerticalToMembraneDivisionRule<DIM>::GetMembrane() const
{
    return mpMembrane;
}

template<unsigned DIM>
double VerticalToMembraneDivisionRule<DIM>::GetDivisionSpringLength() const
{
    return mDivisionSpringLength;
}

// Explicit instantiation
template class VerticalToMembraneDivisionRule<1>;
template class VerticalToMembraneDivisionRule<2>;
template class VerticalToMembraneDivisionRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VerticalToMembraneDivisionRule)
