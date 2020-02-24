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

#ifndef NEIGHBOURTRACKINGMODIFIER_HPP_
#define NEIGHBOURTRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "NodeBasedCellPopulation.hpp"

/**
 * A modifier class which at each simulation time step calculates the volume of each cell
 * and stores it in in the CellData property as "volume". To be used in conjunction with
 * contact inhibition cell cycle models.
 */
template<unsigned DIM>
class NeighbourTrackingModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
        archive & mNeighbourhoodRadius;
        archive & mOnlyAboveNeighbours;
    }

protected:
    /*
     * Member variable defining the relevant radius
     */
    double mNeighbourhoodRadius;

    /*
     * Member variable, whether to only count above neighbours
     */
    bool mOnlyAboveNeighbours;

public:

    /**
     * Default constructor.
     */
    NeighbourTrackingModifier(double neighbourhoodRadius = 1.5, bool onlyAboveNeighbours = false);

    /**
     * Destructor.
     */
    virtual ~NeighbourTrackingModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to compute the volume of each cell in the population and store these in the CellData.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Method to get the value of member variable mNeighbourhoodRadius
    */
    double GetNeighbourhoodRadius() const;

    /**
     * Method to get the value of member variable mOnlyAboveNeighbours
     */
    bool GetOnlyAboveNeighbours() const;

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NeighbourTrackingModifier)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const NeighbourTrackingModifier<DIM> * t, const unsigned int file_version)
		{
		     double neighbourhood_radius = t->GetNeighbourhoodRadius();
		     ar << neighbourhood_radius;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, NeighbourTrackingModifier<DIM> * t, const unsigned int file_version)
		{
            double neighbourhood_radius;
		    ar >> neighbourhood_radius;
		    ::new(t)NeighbourTrackingModifier<DIM>(neighbourhood_radius);
        }
	} // End namespace serialization
} // End namespace boost


#endif /*NEIGHBOURTRACKINGMODIFIER_HPP_*/
