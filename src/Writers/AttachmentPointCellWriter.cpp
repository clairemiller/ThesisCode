/*

Copyright (c) 2005-2018, University of Oxford.
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

#include "AttachmentPointCellWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AttachmentPointCellWriter<ELEMENT_DIM, SPACE_DIM>::AttachmentPointCellWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellattachmentpoints.dat")
{
    this->mVtkCellDataName = "Dist to BM attachment";
}

// Function to extract the membrane point from CellData
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,SPACE_DIM> AttachmentPointCellWriter<ELEMENT_DIM, SPACE_DIM>::GetVectorFromCellData(CellPtr p_cell)
{
	c_vector<double,SPACE_DIM> vec;
	switch(SPACE_DIM) {
		case 1: NEVER_REACHED;
		case 3:
			vec[2] = p_cell->GetCellData()->GetItem("attachment_point_z");
		case 2:
			vec[1] = p_cell->GetCellData()->GetItem("attachment_point_y");
			vec[0] = p_cell->GetCellData()->GetItem("attachment_point_x");
	}
	return vec;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AttachmentPointCellWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // As we are collecting a vector, we instead use the distance to the attachment point
    c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    c_vector<double, SPACE_DIM> att_location = GetVectorFromCellData(pCell);
    c_vector<double, SPACE_DIM> dir_vec = pCellPopulation->rGetMesh().GetVectorFromAtoB(att_location,cell_location);
    return norm_2(dir_vec);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AttachmentPointCellWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // Write location index corresponding to cell
    *this->mpOutStream << pCellPopulation->GetLocationIndexUsingCell(pCell) << " ";

    // Including the cell location seems to be an issue in terms of max line length
    // // Write cell location
    // c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    // for (unsigned i=0; i<SPACE_DIM; i++)
    // {
    //     *this->mpOutStream << cell_location[i] << " ";
    // }

    // Write cell attachment point
    c_vector<double, SPACE_DIM> att_location = GetVectorFromCellData(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << att_location[i] << " ";
    }
}

// Explicit instantiation
template class AttachmentPointCellWriter<1,1>;
template class AttachmentPointCellWriter<1,2>;
template class AttachmentPointCellWriter<2,2>;
template class AttachmentPointCellWriter<1,3>;
template class AttachmentPointCellWriter<2,3>;
template class AttachmentPointCellWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AttachmentPointCellWriter)
