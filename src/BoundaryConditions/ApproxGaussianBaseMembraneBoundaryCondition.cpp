#include "ApproxGaussianBaseMembraneBoundaryCondition.hpp"


template<unsigned DIM>
ApproxGaussianBaseMembraneBoundaryCondition<DIM>::ApproxGaussianBaseMembraneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation, unsigned verticalDirection)
    : SinusoidalBaseMembraneBoundaryCondition<DIM>(pCellPopulation, verticalDirection) 
{
    if ( DIM == 3 )
	{
		EXCEPTION("Undulating membrane forces not yet implemented for 3 dimensions");
	}
}

template<unsigned DIM>
double ApproxGaussianBaseMembraneBoundaryCondition<DIM>::BaseShapeFunction(c_vector<double,DIM> p)
{
    double z = 0.0;
    for (unsigned i=0; i < DIM; i++)
    {
        if (DIM != this->mVert)
        {
            double x = p[i];
            z += std::pow( 0.5*(1.0+sin(x*2.0*_PI/(this->mPeriod))), 4.0);
        }
    }
    z *= (this->mAmplitude)/((double)(DIM-1)); // Should never have DIM 1 due to assert in Abstract class
    return(z);
}

template<unsigned DIM>
c_vector<double,DIM> ApproxGaussianBaseMembraneBoundaryCondition<DIM>::CalculateDerivativesAtPoint(c_vector<double,DIM> p)
{
    // Store variables locally for convenience
    double A = this->mAmplitude;
    double P = this->mPeriod;
    
    c_vector<double,DIM> dzd;
    for (unsigned int i = 0; i < DIM; i++)
    {
        if ( i==(this->mVert) )
        {
            dzd[i] = 1.0;
        }
        else
        {
            dzd[i] = ( 0.5*_PI*A/P )/((double)(DIM-1)) * std::pow(1.0+sin(2.0*_PI*p[i]/P),3.0) * cos(2.0*_PI*p[i]/P);
        }
    }
    return(dzd);
}


// Explicit instantiation
template class ApproxGaussianBaseMembraneBoundaryCondition<1>;
template class ApproxGaussianBaseMembraneBoundaryCondition<2>;
template class ApproxGaussianBaseMembraneBoundaryCondition<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ApproxGaussianBaseMembraneBoundaryCondition)