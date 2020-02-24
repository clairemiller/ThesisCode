
#include "ExponentialDecayPalssonAdhesionForce.hpp"

// Default constructor
template<unsigned DIM>
ExponentialDecayPalssonAdhesionForce<DIM>::ExponentialDecayPalssonAdhesionForce(double peakForce, double decayStartValue, double decayRate)
 : AbstractGradientPalssonAdhesionForce<DIM>(peakForce, decayStartValue), mDecayRate(decayRate)
{
    assert(decayRate > 0.0);
}

template<unsigned DIM>
double ExponentialDecayPalssonAdhesionForce<DIM>::ScalingFunction(CellPtr cellA)
{
    double age = cellA->GetAge();
    // Only applies if age is greater than 2.5 days, or 60 hours
    if ( age < (this->mDecayStartValue) )
    {
        return(1.0);
    }
    // Otherwise return exponential decay scaling factor
    double scaling_factor = std::exp(-mDecayRate*(age-(this->mDecayStartValue)));
    return(scaling_factor);
}

template<unsigned DIM>
double ExponentialDecayPalssonAdhesionForce<DIM>::GetDecayRate() const
{
    return(mDecayRate);
}

template<unsigned DIM>
void ExponentialDecayPalssonAdhesionForce<DIM>::SetDecayRate(double decayRate)
{
    assert(decayRate > 0.0);
    mDecayRate = decayRate;
}


// Explicit instantiation
template class ExponentialDecayPalssonAdhesionForce<1>;
template class ExponentialDecayPalssonAdhesionForce<2>;
template class ExponentialDecayPalssonAdhesionForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExponentialDecayPalssonAdhesionForce)