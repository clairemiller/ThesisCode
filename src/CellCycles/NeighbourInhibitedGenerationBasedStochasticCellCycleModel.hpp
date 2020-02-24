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

#ifndef NEIGHBOURINHIBITEDGENERATIONBASEDSTOCHASTICCELLCYCLEMODEL_HPP_
#define NEIGHBOURINHIBITEDGENERATIONBASEDSTOCHASTICCELLCYCLEMODEL_HPP_

#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"

class NeighbourInhibitedGenerationBasedStochasticCellCycleModel : public FixedG1GenerationalCellCycleModel
{
private:

    friend class boost::serialization::access;

    /**
     * Boost Serialization method for archiving/checkpointing
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<FixedG1GenerationalCellCycleModel>(*this);
        archive & mLikelihoodOfDivision;
        archive & mNumNeighboursForInhibition;
        archive & mCurrentQuiescentDuration;
        archive & mCurrentQuiescentOnsetTime;
    }

protected:
    /*
     * The likelihood that the cell will divide after the G1 phase
    */
    double mLikelihoodOfDivision;

    /**
     * The number of neighbours at which a cell becomes quiescent
     */
    unsigned mNumNeighboursForInhibition;

    /**
     * The time when the current period of quiescence began.
     */
    double mCurrentQuiescentOnsetTime;

    /**
     * How long the current period of quiescence has lasted.
     * Has units of hours.
     */
    double mCurrentQuiescentDuration;

    /**
     * Whether the cell division was symmetric or asymmetric
     * Added for use with the division vector calculation in AbstractOffLatticeCellPopulation
     */
    bool mSymmetricDivision;

    /**
     * Protected copy-constructor for use by CreateCellCycleModel.
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    NeighbourInhibitedGenerationBasedStochasticCellCycleModel(const NeighbourInhibitedGenerationBasedStochasticCellCycleModel& rModel);

    /**
     * Overridden virtual method for setting the G1 duration
    */
    void SetG1Duration();

public:

    /**
     * Constructor.
     */
    NeighbourInhibitedGenerationBasedStochasticCellCycleModel();

    /**
     * Overridden UpdateCellCyclePhase() method.
     */
    void UpdateCellCyclePhase();

    /**
     * Overridden builder method to create new instances of
     * the cell-cycle model.
     *
     * @return new cell-cycle model
     *
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * @param likelihoodOfDivision
    */
    void SetLikelihoodOfDivision(double likelihoodOfDivision);

    /**
     * @return mLikelihoodOfDivision
    */
    double GetLikelihoodOfDivision() const;

    /**
     * @param numNeighboursForInhibition
     */
    void SetNumNeighboursForInhibition(unsigned numNeighboursForInhibition);

    /**
     * @return mNumNeighboursForInhibition
     */
    unsigned GetNumNeighboursForInhibition() const;

    /**
     * @return mCurrentQuiescentDuration
     */
    double GetCurrentQuiescentDuration() const;

    /**
     * @return mCurrentQuiescentOnsetTime
     */
    double GetCurrentQuiescentOnsetTime() const;

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(NeighbourInhibitedGenerationBasedStochasticCellCycleModel)

#endif // NEIGHBOURINHIBITEDGENERATIONBASEDSTOCHASTICCELLCYCLEMODEL_HPP_
