
#include "StochasticExponentialDecayPalssonAdhesionForce.hpp"

// Default constructor
template<unsigned DIM>
StochasticExponentialDecayPalssonAdhesionForce<DIM>::StochasticExponentialDecayPalssonAdhesionForce(double peakForce, double decayStartValue, double minDecayRate, double maxDecayRate)
 : AbstractGradientPalssonAdhesionForce<DIM>(peakForce, decayStartValue), mMinDecayRate(minDecayRate), mMaxDecayRate(maxDecayRate), mIsInitialised(false)
{
    assert(minDecayRate > 0.0);
    assert(maxDecayRate >= minDecayRate);
}

template<unsigned DIM>
void StochasticExponentialDecayPalssonAdhesionForce<DIM>::AssignCellDecayRate(CellPtr pCell)
{
    // Get a uniformly generated number between the two decay rate
    double decay_rate = mMinDecayRate + (RandomNumberGenerator::Instance()->ranf()*(mMaxDecayRate-mMinDecayRate));
    pCell->GetCellData()->SetItem("AdhesionDecayRate",decay_rate);
}

template<unsigned DIM>
c_vector<double, DIM> StochasticExponentialDecayPalssonAdhesionForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    if (!mIsInitialised)
    {
        for ( typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin(); cell_iter!= rCellPopulation.End(); ++cell_iter )
	    {
            AssignCellDecayRate(*cell_iter);
        }
        mIsInitialised = true;
    }
    // Run parent class
    c_vector<double, DIM> force = AbstractGradientPalssonAdhesionForce<DIM>::CalculateForceBetweenNodes(nodeAGlobalIndex,nodeBGlobalIndex,rCellPopulation);
    return(force);
}

template<unsigned DIM>
double StochasticExponentialDecayPalssonAdhesionForce<DIM>::ScalingFunction(CellPtr cellA)
{
    double age = cellA->GetAge();
    // if it's a new cell or the start of the simulation, calculate decay rate
    if (age <= SimulationTime::Instance()->GetTimeStep())
    {
        AssignCellDecayRate(cellA);
    }
    // Only applies if age is greater than 2.5 days, or 60 hours
    if ( age < (this->mDecayStartValue) )
    {
        return(1.0);
    }
    // Otherwise return exponential decay scaling factor
    double decay_rate = cellA->GetCellData()->GetItem("AdhesionDecayRate");
    double scaling_factor = std::exp(-decay_rate*(age-(this->mDecayStartValue)));
    return(scaling_factor);
}

template<unsigned DIM>
double StochasticExponentialDecayPalssonAdhesionForce<DIM>::GetMinDecayRate() const
{
    return(mMinDecayRate);
}

template<unsigned DIM>
void StochasticExponentialDecayPalssonAdhesionForce<DIM>::SetMinDecayRate(double minDecayRate)
{
    assert(minDecayRate > 0.0);
    mMinDecayRate = minDecayRate;
}

template<unsigned DIM>
double StochasticExponentialDecayPalssonAdhesionForce<DIM>::GetMaxDecayRate() const
{
    return(mMaxDecayRate);
}

template<unsigned DIM>
void StochasticExponentialDecayPalssonAdhesionForce<DIM>::SetMaxDecayRate(double maxDecayRate)
{
    assert(maxDecayRate > 0.0);
    mMaxDecayRate = maxDecayRate;
}

template<unsigned DIM>
bool StochasticExponentialDecayPalssonAdhesionForce<DIM>::IsInitialised() const
{
    return mIsInitialised;
}

template<unsigned DIM>
void StochasticExponentialDecayPalssonAdhesionForce<DIM>::SetInitialised(bool isInitialised)
{
    mIsInitialised = isInitialised;
}


// Explicit instantiation
template class StochasticExponentialDecayPalssonAdhesionForce<1>;
template class StochasticExponentialDecayPalssonAdhesionForce<2>;
template class StochasticExponentialDecayPalssonAdhesionForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(StochasticExponentialDecayPalssonAdhesionForce)