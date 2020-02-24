
#include "FixedDurationDifferentiationCellModifier.hpp"

// Default constructor
template<unsigned DIM>
FixedDurationDifferentiationCellModifier<DIM>::FixedDurationDifferentiationCellModifier(double transitTime) 
: AbstractCellBasedSimulationModifier<DIM,DIM>(), mTransitTime(transitTime)
{}

// Override of UpdateAtEndOfTimeStep method
// This sets the colour of the differentiated cells
template<unsigned DIM>
void FixedDurationDifferentiationCellModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	rCellPopulation.Update();

	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
	{

		if ( cell_iter->GetCellProliferativeType()->template IsType<DifferentiatedCellProliferativeType>() )
		{
			double cellAge = cell_iter->GetAge();
			if ( cellAge >= mTransitTime )
			{
				MAKE_PTR_ARGS(AbstractCellProliferativeType, p_abstract_type, (3));
				cell_iter->SetCellProliferativeType( p_abstract_type );
			}

		}
		else if ( cell_iter->GetCellProliferativeType()->GetColour() == 3 && cell_iter->GetAge() >= 2.0*mTransitTime )
		{
			MAKE_PTR_ARGS(AbstractCellProliferativeType, p_abstract_type, (4));
			cell_iter->SetCellProliferativeType( p_abstract_type );
		}
		else if ( cell_iter->GetCellProliferativeType()->GetColour() == 3 && cell_iter->GetAge() >= 3.0*mTransitTime )
		{
			cell_iter -> StartApoptosis();
		}
	}
}


// 'Get' function. Mainly useful for the archiving the parameters
template<unsigned DIM>
double FixedDurationDifferentiationCellModifier<DIM>::GetTransitTime() const
{
	return mTransitTime;
}

// 'Set' function.
template<unsigned DIM>
void FixedDurationDifferentiationCellModifier<DIM>::SetTransitTime(double transitTime)
{
	mTransitTime = transitTime;
}

// Override the 'OutputSimulationModifierParameters' method to include the transit time variable
template<unsigned DIM>
void FixedDurationDifferentiationCellModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<TransitTime>" << mTransitTime << "</TransitTime>\n";

	// Call parent class method
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class FixedDurationDifferentiationCellModifier<1>;
template class FixedDurationDifferentiationCellModifier<2>;
template class FixedDurationDifferentiationCellModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixedDurationDifferentiationCellModifier)
