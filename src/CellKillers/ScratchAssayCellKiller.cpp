
#include "ScratchAssayCellKiller.hpp"

template<unsigned DIM>
ScratchAssayCellKiller<DIM>::ScratchAssayCellKiller(AbstractCellPopulation<DIM,DIM>* pCellPopulation, double colWidth, double centre, double height, double killTime) 
: AbstractCellKiller<DIM>(pCellPopulation), mColWidth(colWidth), mCentre(centre), mHeight(height), mKillTime(killTime)
{
	assert(DIM!=(unsigned)1);
}

template<unsigned DIM>
void ScratchAssayCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{

	double simTime = SimulationTime :: Instance()->GetTime();
	double simTime_Next = simTime + SimulationTime :: Instance() ->GetTimeStep();

	// Exit if not within the time
	if ( mKillTime > simTime_Next || mKillTime < simTime )
		return;

	// Find the cells within the kill zone
	double xMin = mCentre-(mColWidth/2.0), xMax = mCentre+(mColWidth/2.0);
	for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
            cell_iter != this->mpCellPopulation->End();
            ++cell_iter)
	{
		c_vector<double, 2> xCurrent = ( this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter) );
		if ( xCurrent[0] >= xMin && xCurrent[0] <= xMax && xCurrent[DIM-1] >= mHeight )
		{
			cell_iter -> Kill();
		}
	}
}

template<unsigned DIM>
double ScratchAssayCellKiller<DIM>::GetColWidth() const
{
	return mColWidth;
}

template<unsigned DIM>
double ScratchAssayCellKiller<DIM>::GetCentre() const
{
	return mCentre;
}

template<unsigned DIM>
double ScratchAssayCellKiller<DIM>::GetHeight() const
{
	return mHeight;
}

template<unsigned DIM>
double ScratchAssayCellKiller<DIM>::GetKillTime() const
{
	return mKillTime;
}

template<unsigned DIM>
void ScratchAssayCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<ColWidth>" << mColWidth << "</ColWidth>\n";
	*rParamsFile << "\t\t\t<Centre>" << mCentre << "</Centre>\n";
	*rParamsFile << "\t\t\t<KillTime>" << mKillTime << "</KillTime>\n";
	
	// Call parent class method
	AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}


// Explicit instantiation
template class ScratchAssayCellKiller<2>;
template class ScratchAssayCellKiller<3>;

// Serialization
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ScratchAssayCellKiller)