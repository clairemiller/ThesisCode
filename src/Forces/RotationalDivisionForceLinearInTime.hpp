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

#ifndef ROTATIONALDIVISIONFORCELINEARINTIME_HPP_
#define ROTATIONALDIVISIONFORCELINEARINTIME_HPP_

#include "RotationalDivisionForce.hpp"
#include "AbstractUndulatingBaseMembraneForce.hpp"

/**
 * A class to apply a rotational force to two cells involved in the mitotic division phase from a single cell
 */
template<unsigned  DIM>
class RotationalDivisionForceLinearInTime : public RotationalDivisionForce<DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<RotationalDivisionForce<DIM> >(*this);
    }

public:

    /**
     * Constructor.
     */
    RotationalDivisionForceLinearInTime(double torsionCoefficient, double divisionAngle = M_PI/2.0, double growthDuration=1.0);

    /**
     * Get method for division angle
     * 
     * @return the value of mTorsionCoefficient
     */
    virtual double GetTorsionCoefficient(double age) const;

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);

    /**
     * Overridden WriteDataToVisualizerSetupFile() method.
     * Write any data necessary to a visualization setup file.
     * Used by AbstractCellBasedSimulation::WriteVisualizerSetupFile().
     * 
     * @param pVizSetupFile a visualization setup file
     */
    virtual void WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile);
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RotationalDivisionForceLinearInTime)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const RotationalDivisionForceLinearInTime<DIM> * t, const unsigned int file_version)
		{
            double torsion_coefficient, division_spring_legnth, division_angle, growth_duration;
            growth_duration = t->GetGrowthDuration();
            // In order to get the max coefficent we need to use the growth duration as the age
            torsion_coefficient = t->GetTorsionCoefficient(growth_duration); 
            division_angle = t->GetDivisionAngle(zero_vector<double>(DIM));

            ar << torsion_coefficient;
            ar << division_angle;
            ar << growth_duration;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, RotationalDivisionForceLinearInTime<DIM> * t, const unsigned int file_version)
		{
            double torsion_coefficient, division_angle, growth_duration;
            ar >> torsion_coefficient;
            ar >> division_angle;
            ar >> growth_duration;

            ::new(t)RotationalDivisionForceLinearInTime<DIM>(torsion_coefficient,division_angle,growth_duration);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*ROTATIONALDIVISIONFORCELINEARINTIME_HPP_*/
