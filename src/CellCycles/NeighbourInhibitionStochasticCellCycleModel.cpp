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

#include "NeighbourInhibitionStochasticCellCycleModel.hpp"

NeighbourInhibitionStochasticCellCycleModel::NeighbourInhibitionStochasticCellCycleModel()
    : PopulationAsymmetryCellCycleModel(),
      mLikelihoodOfDivision(DOUBLE_UNSET),
      mNumNeighboursForInhibition(UNSIGNED_UNSET),
      mCurrentQuiescentOnsetTime(SimulationTime::Instance()->GetTime()),
      mCurrentQuiescentDuration(0.0)
{
}

NeighbourInhibitionStochasticCellCycleModel::NeighbourInhibitionStochasticCellCycleModel(const NeighbourInhibitionStochasticCellCycleModel& rModel)
    : PopulationAsymmetryCellCycleModel(rModel),
      mLikelihoodOfDivision(rModel.mLikelihoodOfDivision),
      mNumNeighboursForInhibition(rModel.mNumNeighboursForInhibition),
      mCurrentQuiescentOnsetTime(rModel.mCurrentQuiescentOnsetTime),
      mCurrentQuiescentDuration(rModel.mCurrentQuiescentDuration)
{
    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables will already
     * have been correctly initialized in its constructor or parent classes.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     *
     */
}

void NeighbourInhibitionStochasticCellCycleModel::SetG1Duration()
{
    // Determine the number of G1 cycles according to a geometric distribution
    // We need to check that prob!=1.0
    unsigned n_g1_cycles = 1;
    if ( mLikelihoodOfDivision != 1.0 )
    {
        // Get a uniformly distributed random number
        double v_rand = RandomNumberGenerator::Instance()->ranf();
        // We need to ensure v_rand is not 0.0
        while( v_rand == 0.0 )
        {
            v_rand = RandomNumberGenerator::Instance()->ranf();
        }
        // Convert uniformly distributed to geometrically distributed
        n_g1_cycles = 1+std::floor( std::log(v_rand)/std::log(1.0-mLikelihoodOfDivision) );
    }

    // Now determine G1 cycle length time based on cell type
    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mG1Duration = GetStemCellG1Duration()*n_g1_cycles;     
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        mG1Duration = GetTransitCellG1Duration()*n_g1_cycles;
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}

void NeighbourInhibitionStochasticCellCycleModel::UpdateCellCyclePhase()
{
    if ( mLikelihoodOfDivision == DOUBLE_UNSET)
    {
        EXCEPTION("The likelihood of division for the cell cycle model has not yet been set.\n");
    } 
    else if (mNumNeighboursForInhibition == UNSIGNED_UNSET )
    {
        EXCEPTION("The number of neighbours for inhibition in the cell cycle model has not been set yet.");
    }

    // Get cell volume
    unsigned n_neighbours = mpCell->GetCellData()->GetItem("n_neighbours");

    // Removes the cell label
    mpCell->RemoveCellProperty<CellLabel>();

    if (mCurrentCellCyclePhase == G_ONE_PHASE)
    {
        // Update G1 duration based on cell volume
        double dt = SimulationTime::Instance()->GetTimeStep();

        // If we have enough neighbours we become quiescent
        if ( n_neighbours >= mNumNeighboursForInhibition )
        {
            // Update the duration of the current period of contact inhibition.
            mCurrentQuiescentDuration = SimulationTime::Instance()->GetTime() - mCurrentQuiescentOnsetTime;
            mG1Duration += dt;

            /*
             * This method is usually called within a CellBasedSimulation, after the CellPopulation
             * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
             * CellPropertyRegistry::Instance() here when adding the CellLabel, we would be creating
             * a new CellPropertyRegistry. In this case the CellLabel's cell count would be incorrect.
             * We must therefore access the CellLabel via the cell's CellPropertyCollection.
             */
            boost::shared_ptr<AbstractCellProperty> p_label =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<CellLabel>();
            // Set the value of the cell type dependent on proliferative type
            if ( mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>() )
            {
              p_label.reset(new CellLabel(6));
            }
            else if ( mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() )
            {
              p_label.reset(new CellLabel(7));
            }
            mpCell->AddCellProperty(p_label);
        }
        else
        {
            // Reset the cell's quiescent duration and update the time at which the onset of quiescent occurs
            mCurrentQuiescentDuration = 0.0;
            mCurrentQuiescentOnsetTime = SimulationTime::Instance()->GetTime();
        }
    }

    double time_since_birth = GetAge();
    assert(time_since_birth >= 0);

    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if ( time_since_birth < GetMDuration() )
    {
        mCurrentCellCyclePhase = M_PHASE;
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration)
    {
        mCurrentCellCyclePhase = G_ONE_PHASE;
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;
    }
    else if ( time_since_birth < GetMDuration() + mG1Duration + GetSDuration() + GetG2Duration())
    {
        mCurrentCellCyclePhase = G_TWO_PHASE;
    }
}

AbstractCellCycleModel* NeighbourInhibitionStochasticCellCycleModel::CreateCellCycleModel()
{
    return new NeighbourInhibitionStochasticCellCycleModel(*this);
}

void NeighbourInhibitionStochasticCellCycleModel::SetLikelihoodOfDivision(double likelihoodOfDivision)
{
    mLikelihoodOfDivision = likelihoodOfDivision;
}

double NeighbourInhibitionStochasticCellCycleModel::GetLikelihoodOfDivision() const
{
    return mLikelihoodOfDivision;
}

void NeighbourInhibitionStochasticCellCycleModel::SetNumNeighboursForInhibition(unsigned numNeighboursForInhibition)
{
    mNumNeighboursForInhibition = numNeighboursForInhibition;
}

unsigned NeighbourInhibitionStochasticCellCycleModel::GetNumNeighboursForInhibition() const
{
    return mNumNeighboursForInhibition;
}

double NeighbourInhibitionStochasticCellCycleModel::GetCurrentQuiescentDuration() const
{
    return mCurrentQuiescentDuration;
}

double NeighbourInhibitionStochasticCellCycleModel::GetCurrentQuiescentOnsetTime() const
{
    return mCurrentQuiescentOnsetTime;
}

void NeighbourInhibitionStochasticCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<LikelihoodOfDivision>" << mLikelihoodOfDivision << "</LikelihoodOfDivision>\n";
    *rParamsFile << "\t\t\t<NumNeighboursForInhibition>" << mNumNeighboursForInhibition << "</NumNeighboursForInhibition>\n";

    // Call method on direct parent class
    PopulationAsymmetryCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(NeighbourInhibitionStochasticCellCycleModel)
