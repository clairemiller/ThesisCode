
#include "LinearGradientPalssonAdhesionForce.hpp"

// Default constructor
template<unsigned DIM>
LinearGradientPalssonAdhesionForce<DIM>::LinearGradientPalssonAdhesionForce(double peakForce, double decayStartValue, double decaySlope)
 : AbstractGradientPalssonAdhesionForce<DIM>(peakForce, decayStartValue), mDecayRate(decayRate)
{
    assert(decayRate > 0.0);
}

template<unsigned DIM>
double LinearGradientPalssonAdhesionForce<DIM>::ScalingFunction(CellPtr cellA)
{
    double age = cellA->GetAge();
    // Only applies if age is greater than 2.5 days, or 60 hours
    if ( age < (this->mDecayStartValue) )
    {
        return(1.0);
    }
    // Otherwise return exponential decay scaling factor
    double scaling_factor = 1.0-decaySlope*(age-(this->mDecayStartValue));
    return(scaling_factor);
}

template<unsigned DIM>
double LinearGradientPalssonAdhesionForce<DIM>::GetDecayRate() const
{
    return(mDecayRate);
}

template<unsigned DIM>
void LinearGradientPalssonAdhesionForce<DIM>::SetDecayRate(double decayRate)
{
    assert(decayRate > 0.0);
    mDecayRate = decayRate;
}


// Explicit instantiation
template class LinearGradientPalssonAdhesionForce<1>;
template class LinearGradientPalssonAdhesionForce<2>;
template class LinearGradientPalssonAdhesionForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LinearGradientPalssonAdhesionForce)