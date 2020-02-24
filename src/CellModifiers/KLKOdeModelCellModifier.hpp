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

#ifndef KLKODEMODELCELLMODIFIER
#define KLKODEMODELCELLMODIFIER

#include "ChasteSerialization.hpp"

#include "AbstractCellPopulation.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * An abstract modifier class (to implement setup, update and finalise methods), for use in cell-based simulations.
 */
template<unsigned  DIM>
class KLKOdeModelCellModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
      archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM> >(*this);
      archive & mHeightStartSC;
      archive & mExpectedSCHeight;
      archive & mVertical;
    }
    
    // Member parameters
    double mHeightStartSC; // The height at which we enter the SC (i.e. where z=0 for the ode system)
    double mExpectedSCHeight; // The expected height of the SC
    unsigned mVertical; // The vertical dimension (normally DIM-1 unless parallelised)

public:

    /**
     * Default constructor.
     */
    KLKOdeModelCellModifier(double heightStartSC, double expectedSCHeight, unsigned vertical=DIM-1);

    /**
     * Destructor.
     */
    virtual ~KLKOdeModelCellModifier();

    // Override UpdateCellData
	void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Specify what to do in the simulation at the end of each timestep.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Specify what to do in the simulation before the start of the time loop.
     * Required override
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Output any simulation modifier parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);

    /**
     * Get method for the starting height of the stratum corneum
     * @return the value of mHeightStartSC
     */
    double GetHeightStartSC() const;

    /**
     * Get method for the expected height of the stratum corneum
     * @return the value of mExpectedSCHeight
     */
    double GetExpectedSCHeight() const;

    /**
     * Get method for the vertical direction member variable
     * @return the value of mVertical
     */
    unsigned GetVertical() const;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(KLKOdeModelCellModifier)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const KLKOdeModelCellModifier<DIM> * t, const unsigned int file_version)
		{
        double heightStartSC = t->GetHeightStartSC();
        double expectedSCHeight = t->GetExpectedSCHeight();
        unsigned vertical = t->GetVertical();

        ar << heightStartSC;
        ar << expectedSCHeight;
        ar << vertical;
    }

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, KLKOdeModelCellModifier<DIM> * t, const unsigned int file_version)
		{
      double heightStartSC, expectedSCHeight;
      unsigned vertical;

      ar >> heightStartSC;
      ar >> expectedSCHeight;
      ar >> vertical;

			::new(t)KLKOdeModelCellModifier<DIM>(heightStartSC, expectedSCHeight, vertical);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*KLKODEMODELCELLMODIFIER*/
