
#include "HeightDependentDivisionModifier.hpp"

template<unsigned DIM>
// Default constructor
HeightDependentDivisionModifier<DIM>::HeightDependentDivisionModifier(double height)
: mHeight(height), mRevertAfterG1Phase(false), mVertical(DIM-1), mUseMembrane(false), mMaxStemCells(UINT_MAX)
{}

template<unsigned DIM>
void HeightDependentDivisionModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
	/*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void HeightDependentDivisionModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);
}
	
// Override UpdatAtEndOfTimeStep method
template<unsigned DIM>
void HeightDependentDivisionModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	rCellPopulation.Update();

	// Get the simulation time
	double simTime = SimulationTime :: Instance()->GetTime();

	// If the cell is below a certain height, it becomes a stem cell
	MAKE_PTR(StemCellProliferativeType, p_stem_type);
	MAKE_PTR(DifferentiatedCellProliferativeType,p_diff_type);
	unsigned stem_counter = 0;
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
	{
		c_vector<double,DIM> cell_loc = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
		double cell_height = cell_loc[mVertical];
		if (mUseMembrane)
		{
			c_vector<double,DIM> membrane_loc = mpMembrane->DetermineClosestMembraneLocation(cell_loc);
			cell_height = norm_2(cell_loc - membrane_loc);
		}
		if ( cell_height < mHeight && stem_counter < mMaxStemCells )
		{
			stem_counter++;
			if ( (*cell_iter)->GetCellProliferativeType()->template IsType<DifferentiatedCellProliferativeType>() )
			{
				(*cell_iter)->SetCellProliferativeType(p_stem_type);
				(*cell_iter)->SetBirthTime(simTime);
				// We need to run initialise to reset the G1 duration
				AbstractSimplePhaseBasedCellCycleModel* p_cell_cycle = dynamic_cast<AbstractSimplePhaseBasedCellCycleModel*>((*cell_iter)->GetCellCycleModel());
				if (!p_cell_cycle)
				{
					EXCEPTION("Height dependent division modifier only valid for simple phase based cell cycle models.");
				}
				p_cell_cycle->Initialise();
			}
		}
		else if ( (*cell_iter)->GetCellProliferativeType()->template IsType<StemCellProliferativeType>() && ((*cell_iter)->GetAge() < 1.0))
		{
			(*cell_iter)->SetCellProliferativeType(p_diff_type);
		}
		else if ( (*cell_iter)->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>() )
		{
			EXCEPTION("Height dependent division modifier does not yet include transit cells.");
		}
	}
}

// Get member variables functions
template<unsigned DIM>
double HeightDependentDivisionModifier<DIM>::GetHeight() const
{
	return mHeight;
}

template<unsigned DIM>
bool HeightDependentDivisionModifier<DIM>::GetRevertAfterG1Phase() const
{
	return mRevertAfterG1Phase;
}

template<unsigned DIM>
unsigned HeightDependentDivisionModifier<DIM>::GetVertical() const
{
	return mVertical;
}

template<unsigned DIM>
bool HeightDependentDivisionModifier<DIM>::GetUseMembrane() const
{
	return mUseMembrane;
}

// Set member functions
template<unsigned DIM>
void HeightDependentDivisionModifier<DIM>::SetHeight(double height)
{
	mHeight = height;
}
template<unsigned DIM>
void HeightDependentDivisionModifier<DIM>::SetRevertAfterG1Phase(bool revertAfterG1Phase)
{
	mRevertAfterG1Phase = revertAfterG1Phase;
}
template<unsigned DIM>
void HeightDependentDivisionModifier<DIM>::SetVertical(unsigned vertical)
{
	mVertical = vertical;
	assert(vertical < DIM);
}
template<unsigned DIM>
void HeightDependentDivisionModifier<DIM>::SetMembrane(boost::shared_ptr<AbstractUndulatingBaseMembraneBoundaryCondition<DIM> > pMembrane)
{
	mpMembrane = pMembrane;
	mUseMembrane = true;
}
template<unsigned DIM>
void HeightDependentDivisionModifier<DIM>::SetMaxStemCells(unsigned maxCells)
{
	mMaxStemCells = maxCells;
}

// Override output method
template<unsigned DIM>
void HeightDependentDivisionModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<Height>" << mHeight << "</Height>\n";
	*rParamsFile << "\t\t\t<RevertAfterG1Phase>" << mRevertAfterG1Phase << "</RevertAfterG1Phase>\n";
	*rParamsFile << "\t\t\t<Vertical>" << mVertical << "</Vertical>\n";
	*rParamsFile << "\t\t\t<UseMembrane>" << mUseMembrane << "</UseMembrane>\n";

	// Call parent class method
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template class HeightDependentDivisionModifier<1>;
template class HeightDependentDivisionModifier<2>;
template class HeightDependentDivisionModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HeightDependentDivisionModifier)
