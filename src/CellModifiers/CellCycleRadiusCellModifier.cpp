
#include "CellCycleRadiusCellModifier.hpp"

template<unsigned DIM>
// Default constructor
CellCycleRadiusCellModifier<DIM>::CellCycleRadiusCellModifier(double radiusAtDivision)
: mRadiusAtDivision(radiusAtDivision), mRadiusAfterDivision(-1.0)
{
	assert(radiusAtDivision>0.0);
}

// Override UpdatAtEndOfTimeStep method
template<unsigned DIM>
void CellCycleRadiusCellModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	rCellPopulation.Update();

	assert(mRadiusAfterDivision>0.0);

	// Iterate over population
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
	{
		boost::shared_ptr<AbstractCellProliferativeType> cell_type = cell_iter->GetCellProliferativeType();
		Node<DIM>* p_node = rCellPopulation.GetNode(rCellPopulation.GetLocationIndexUsingCell(*cell_iter));
		// If the cell is about to divide, reset the parent radius to one
		// The division will occur next time step before the forces are recalculated
		if ( cell_iter->ReadyToDivide() )
		{
			p_node->SetRadius(mRadiusAfterDivision);
		}
		// If it is a dividing cell, we grow the radius
		else if (cell_type->template IsType<StemCellProliferativeType>() || cell_type->template IsType<TransitCellProliferativeType>() )
		{
			AbstractPhaseBasedCellCycleModel* cell_model = dynamic_cast<AbstractPhaseBasedCellCycleModel*>(cell_iter->GetCellCycleModel());
			if ( !cell_model )
			{
				EXCEPTION("CellCycleRadiusCellModifier only works with cell cycles derived from the AbstractPhaseBasedCellCycleModel class.");
			}
			double interphase_time = cell_model->GetG1Duration() + cell_model->GetSDuration() + cell_model->GetG2Duration();
			double cell_age = cell_model->GetAge();
			double new_radius = mRadiusAfterDivision + (mRadiusAtDivision-mRadiusAfterDivision)*std::min(cell_age,interphase_time)/interphase_time;
			p_node->SetRadius(new_radius);
		}
	}

}

template<unsigned DIM>
void CellCycleRadiusCellModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory)
{
	// Check the initial radius of the cells to check they are of uniform radius
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter != rCellPopulation.End(); ++cell_iter)
	{
		Node<DIM>* p_node = rCellPopulation.GetNode(rCellPopulation.GetLocationIndexUsingCell(*cell_iter));
		if ( mRadiusAfterDivision < 0.0 )
		{
			mRadiusAfterDivision = p_node->GetRadius();
		}
		else if ( p_node->GetRadius() != mRadiusAfterDivision )
		{
			EXCEPTION("CellCycleRadiusModifier is only written for cell populations of homogeneous radius.");
		}
	}
}

// Get member variables functions
template<unsigned DIM>
double CellCycleRadiusCellModifier<DIM>::GetRadiusAtDivision() const
{
	return mRadiusAtDivision;
}

// Set member functions
template<unsigned DIM>
void CellCycleRadiusCellModifier<DIM>::SetRadiusAtDivision(double radiusAtDivision)
{
	mRadiusAtDivision = radiusAtDivision;
}


// Override output method
template<unsigned DIM>
void CellCycleRadiusCellModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<RadiusAtDivision>" << mRadiusAtDivision << "</RadiusAtDivision>\n";

	// Call parent class method
	AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template class CellCycleRadiusCellModifier<1>;
template class CellCycleRadiusCellModifier<2>;
template class CellCycleRadiusCellModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellCycleRadiusCellModifier)