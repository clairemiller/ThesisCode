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

#include "Debug.hpp"
#include "RotationalDivisionForceOnUndulatingMembrane.hpp"

template<unsigned DIM>
RotationalDivisionForceOnUndulatingMembrane<DIM>::RotationalDivisionForceOnUndulatingMembrane(boost::shared_ptr< AbstractUndulatingBaseMembraneForce<DIM> > membrane, double torsionCoefficient, double growthDuration)
    : RotationalDivisionForce<DIM>(torsionCoefficient,growthDuration), mpMembrane(membrane)
{
    unsigned vertical = mpMembrane->GetVerticalDirection();
    if (vertical != DIM-1)
    {
        EXCEPTION("Rotational force only valid for vertical direction = DIM-1");
    }
}

template<unsigned DIM>
c_vector<double,DIM> RotationalDivisionForceOnUndulatingMembrane<DIM>::GetNormalVector(c_vector<double,DIM> p) const
{
    unsigned vertical = mpMembrane->GetVerticalDirection();
    c_vector<double,DIM> normal = -1.0*(mpMembrane->CalculateDerivativesAtPoint(p));
    normal[vertical] = -1.0*normal[vertical];
    return(normal);
}

template<unsigned DIM>
double RotationalDivisionForceOnUndulatingMembrane<DIM>::GetDistanceToBasePlane(c_vector<double,DIM> p) const
{
    c_vector<double,DIM> membrane_location = mpMembrane->DetermineClosestMembraneLocation(p);
    c_vector<double,DIM> distance_vector = p-membrane_location;
    return( norm_2(distance_vector) );
}
    

template<unsigned DIM>
void RotationalDivisionForceOnUndulatingMembrane<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Nothing to output
    RotationalDivisionForce<DIM>::OutputForceParameters(rParamsFile);
}

template<unsigned DIM>
void RotationalDivisionForceOnUndulatingMembrane<DIM>::WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile)
{
    // Nothing to output
}



// Explicit instantiation
template class RotationalDivisionForceOnUndulatingMembrane<1>;
template class RotationalDivisionForceOnUndulatingMembrane<2>;
template class RotationalDivisionForceOnUndulatingMembrane<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RotationalDivisionForceOnUndulatingMembrane)