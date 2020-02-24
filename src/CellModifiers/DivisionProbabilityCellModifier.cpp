
#include "DivisionProbabilityCellModifier.hpp"

template<unsigned DIM>
// Default constructor
DivisionProbabilityCellModifier<DIM>::DivisionProbabilityCellModifier(double probNonContactInhibited, double probInhibitedStem, double probInhibitedTransit)
: mProb_NonContactInhibited(probNonContactInhibited), mProb_InhibitedStem(probInhibitedStem), mProb_InhibitedTransit(probInhibitedTransit)
{}

// Override UpdatAtEndOfTimeStep method
template<unsigned DIM>
void DivisionProbabilityCellModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	rCellPopulation.Update();

	// Get the simulation time
	double simTime = SimulationTime :: Instance()->GetTime();

	// Get the time step
	double dt = SimulationTime::Instance()->GetTimeStep();

	// If the cell cycle has passed the check point: Determine if division or not. If not, reset the birth time to re-start G1 phase
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
	{
		double cellAge = cell_iter->GetAge();
		AbstractSimplePhaseBasedCellCycleModel* cellcycle = dynamic_cast<AbstractSimplePhaseBasedCellCycleModel*>(cell_iter->GetCellCycleModel());
		if ( cellcycle == NULL ) 
		{
			EXCEPTION("DivisionProbabilityCellModifier only accepts cell cycles which inherit from AbstractSimplePhaseBasedCellCycleModel.\n");
		}
		double timeG1Checkpoint = cellcycle->GetMDuration() + cellcycle->GetG1Duration();
		CellCyclePhase currentPhase = cellcycle->GetCurrentCellCyclePhase();
		if ( currentPhase == G_ONE_PHASE && (cellAge+dt) >= timeG1Checkpoint )
		{
			// Determine the number of neighbours and consequently the probability
			unsigned cell_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
			int nNeighbours = (rCellPopulation.GetNeighbouringNodeIndices(cell_index)).size();
			// Add in the base membrane as a neighbour
			if ( (rCellPopulation.GetLocationOfCellCentre(*cell_iter))(DIM-1) <= 0.5 )
			{
				nNeighbours += 1;
			}
			double mProbDivision=0.0;
			// Contact inhibited
			if ( nNeighbours > 3 )
			{
				// Stem Cell contact inhibited
				if ( cell_iter->GetCellProliferativeType()->template IsType<StemCellProliferativeType>() )
				{
					mProbDivision = mProb_InhibitedStem;
				}
				// Transit cell contact inhibited
				else if ( cell_iter->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>() )
				{
					mProbDivision = mProb_InhibitedTransit;
				}
			}
			// Non contact inhibited
			else if ( nNeighbours >= 0 )
			{
				mProbDivision = mProb_NonContactInhibited;
			}
			else
			{
				NEVER_REACHED;
			}

			// Now determine if it divides
			if ( (RandomNumberGenerator::Instance()->ranf()) > mProbDivision ) 
    		{
    			// Cell doesn't divide so reset the cell age
    			cell_iter->SetBirthTime( simTime - cellcycle->GetMDuration() );
    		}
		}
	}
}

// Get member variables functions
template<unsigned DIM>
double DivisionProbabilityCellModifier<DIM>::GetProbNonContactInhibited() const
{
	return mProb_NonContactInhibited;
}

template<unsigned DIM>
double DivisionProbabilityCellModifier<DIM>::GetProbInhibitedStem() const
{
	return mProb_InhibitedStem;
}

template<unsigned DIM>
double DivisionProbabilityCellModifier<DIM>::GetProbInhibitedTransit() const
{
	return mProb_InhibitedTransit;
}

// Set member functions
template<unsigned DIM>
void DivisionProbabilityCellModifier<DIM>::SetProbNonContactInhibited(double probNonContactInhibited)
{
	mProb_NonContactInhibited = probNonContactInhibited;
}
template<unsigned DIM>
void DivisionProbabilityCellModifier<DIM>::SetProbInhibitedStem(double probInhibitedStem)
{
	mProb_InhibitedStem = probInhibitedStem;
}
template<unsigned DIM>
void DivisionProbabilityCellModifier<DIM>::SetProbInhibitedTransit(double probInhibitedTransit)
{
	mProb_InhibitedTransit = probInhibitedTransit;
}

// Override output method
template<unsigned DIM>
void DivisionProbabilityCellModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<ProbabilityDivisionNonContactInhibited>" << mProb_NonContactInhibited << "</ProbabilityDivisionNonContactInhibited>\n";
	*rParamsFile << "\t\t\t<ProbabilityDivisionInhibitedStem" << mProb_InhibitedStem << "</ProbabilityDivisionInhibitedStem\n";
	*rParamsFile << "\t\t\t<ProbabilityDivisionInhibitedTransit" << mProb_InhibitedTransit << "</ProbabilityDivisionInhibitedTransit\n>";

	// Call parent class method
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template class DivisionProbabilityCellModifier<1>;
template class DivisionProbabilityCellModifier<2>;
template class DivisionProbabilityCellModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DivisionProbabilityCellModifier)