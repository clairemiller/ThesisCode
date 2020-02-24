
#include "AbstractGradientPalssonAdhesionForce.hpp"

// Default constructor
template<unsigned DIM>
AbstractGradientPalssonAdhesionForce<DIM>::AbstractGradientPalssonAdhesionForce(double peakForce, double decayStartValue)
 : PalssonAdhesionForce<DIM>(), mPeakForce(peakForce), mDecayStartValue(decayStartValue)
{
    assert(decayStartValue >= 0.0);
    this->mAlpha = peakForce/0.02669;
}

// Overridden CalculateForceBetweenNodes method
// NOTE: Currently only implemented for node based cell populations
template<unsigned DIM>
c_vector<double, DIM> AbstractGradientPalssonAdhesionForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    c_vector<double,DIM> normal_force = PalssonAdhesionForce<DIM>::CalculateForceBetweenNodes(nodeAGlobalIndex, nodeBGlobalIndex, rCellPopulation);

    // Get the cell pointers
    CellPtr cellA = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr cellB = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    // Check if either cell is proliferative
    bool isProlif = (cellA->GetCellProliferativeType()->IsType<StemCellProliferativeType>() || cellA->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() );
    isProlif = isProlif || (cellB->GetCellProliferativeType()->IsType<StemCellProliferativeType>() || cellB->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() );
    if ( isProlif )
    {
        return (normal_force);
    }

    // Get the scaling coefficient
	double scaling_val_a = ScalingFunction(cellA);
	double scaling_val_b = ScalingFunction(cellB);

    // Take the average of the 2
    double scaling_val = 0.5*(scaling_val_a + scaling_val_b);
    assert(scaling_val >= 0.0);
    c_vector<double,DIM> force = scaling_val*normal_force;
    return (force);
}

template<unsigned DIM>
double AbstractGradientPalssonAdhesionForce<DIM>::GetPeakForce() const
{
    return(mPeakForce);
}

template<unsigned DIM>
double AbstractGradientPalssonAdhesionForce<DIM>::GetDecayStartValue() const
{
    return(mDecayStartValue);
}

template<unsigned DIM>
void AbstractGradientPalssonAdhesionForce<DIM>::SetPeakForce(double peakForce)
{
    mPeakForce = peakForce;
    this->mAlpha = peakForce/0.02669;
}

template<unsigned DIM>
void AbstractGradientPalssonAdhesionForce<DIM>::SetDecayStartValue(double decayStartValue)
{
    assert(decayStartValue >= 0.0);
    mDecayStartValue = decayStartValue;
}


// Overridden OutputForceParameters method
template<unsigned DIM>
void AbstractGradientPalssonAdhesionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<PeakForce>" << mPeakForce << "</PeakForce\n";
	*rParamsFile << "\t\t\t<DecayStartValue>" << mDecayStartValue << "</DecayStartValue\n";

	PalssonAdhesionForce<DIM>::OutputForceParameters(rParamsFile);
}

// Overridden WriteDataToVisualisedSetupFile method
template<unsigned DIM>
void AbstractGradientPalssonAdhesionForce<DIM>::WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile)
{
	PalssonAdhesionForce<DIM>::WriteDataToVisualizerSetupFile(pVizSetupFile);
}

// Explicit instantiation
template class AbstractGradientPalssonAdhesionForce<1>;
template class AbstractGradientPalssonAdhesionForce<2>;
template class AbstractGradientPalssonAdhesionForce<3>;