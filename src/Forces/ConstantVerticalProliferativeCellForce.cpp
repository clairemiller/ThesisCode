
#include "ConstantVerticalProliferativeCellForce.hpp"
#include "MeshBasedCellPopulation.hpp"

// Default constructor
template<unsigned DIM>
ConstantVerticalProliferativeCellForce<DIM>::ConstantVerticalProliferativeCellForce(double forceMagnitude, double cutOffHeight, unsigned verticalDirection)
 		: AbstractForce<DIM,DIM>(), mForceMagnitude(forceMagnitude), mCutOffHeight(cutOffHeight), mVerticalDirection(verticalDirection)
{
	assert(mCutOffHeight > 0.0);
	assert(verticalDirection<DIM);
} 

// Overridden AddForceContribution method
template<unsigned DIM>
void ConstantVerticalProliferativeCellForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
	// Check that it is a node based cell population
	NodeBasedCellPopulation<DIM>* pNodeBasedPopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
	if ( pNodeBasedPopulation == NULL )
	{
		EXCEPTION("Currently only valid for node based population.\n");
	}

	// Loop over cell population and determine the force
	for ( typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter!= rCellPopulation.End(); ++cell_iter )
	{
		// Create a container for the force
		// We only ever use the DIM-1 element
		c_vector<double,DIM> force = zero_vector<double>(DIM);
		
		// Calculate cell height
		double height = (rCellPopulation.GetLocationOfCellCentre(*cell_iter))(mVerticalDirection);
		
		// Get the proliferative type
		boost::shared_ptr<AbstractCellProliferativeType> p_prolif_type = cell_iter->GetCellProliferativeType();
		bool isProlif = ( p_prolif_type->template IsType<StemCellProliferativeType>() ) || 
							( p_prolif_type->template IsType<TransitCellProliferativeType>() );
		if ( isProlif && height < mCutOffHeight  )
		{
			force[mVerticalDirection] = mForceMagnitude;
			// Find the index of the cell
			unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
			// Add force
			rCellPopulation.GetNode(index)->AddAppliedForceContribution(force); 
		}
	}
}


// Get methods for the member variables
template<unsigned DIM>
double ConstantVerticalProliferativeCellForce<DIM>::GetForceMagnitude() const
{
	return mForceMagnitude;
} 

template<unsigned DIM>
double ConstantVerticalProliferativeCellForce<DIM>::GetCutOffHeight() const
{
	return mCutOffHeight;
}

template<unsigned DIM>
unsigned ConstantVerticalProliferativeCellForce<DIM>::GetVerticalDirection() const
{
	return mVerticalDirection;
}


// Overridden OutputForceParameters method
template<unsigned DIM>
void ConstantVerticalProliferativeCellForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<ForceMagnitude>" << mForceMagnitude << "</ForceMagnitude>\n";
	*rParamsFile << "\t\t\t<CutOffHeight>" << mCutOffHeight << "</CutOffHeight>\n";
	*rParamsFile << "\t\t\t<VerticalDirection>" << mVerticalDirection << "</VerticalDirection>\n";
	

	AbstractForce<DIM,DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class ConstantVerticalProliferativeCellForce<1>;
template class ConstantVerticalProliferativeCellForce<2>;
template class ConstantVerticalProliferativeCellForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ConstantVerticalProliferativeCellForce)