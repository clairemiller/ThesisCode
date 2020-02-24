#include "SinusoidalBaseMembraneBoundaryCondition.hpp"
#include "Debug.hpp"

template<unsigned DIM>
SinusoidalBaseMembraneBoundaryCondition<DIM>::SinusoidalBaseMembraneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation, unsigned verticalDirection, double amplitude, double period)
    : AbstractUndulatingBaseMembraneBoundaryCondition<DIM>(pCellPopulation, verticalDirection), 
        mAmplitude(amplitude), mPeriod(period)
{
    // If it's a periodic mesh, we need to check the width suits the function (x0=xN)
	PeriodicNdNodesOnlyMesh<DIM>* rMesh = dynamic_cast<PeriodicNdNodesOnlyMesh<DIM>* >( &(pCellPopulation->rGetMesh()) );
	if (rMesh)
	{
		std::vector<unsigned> periodic_dims = rMesh->GetPeriodicDimensions();
		for (unsigned i = 0; i < DIM; i++)
		{
            if ( i != this->mVert)
            {
                std::vector<unsigned>::iterator it_dim = std::find(periodic_dims.begin(),periodic_dims.end(),i);
                if ( it_dim != periodic_dims.end() )
                {
                    double width = rMesh->GetPeriodicWidths()[it_dim-periodic_dims.begin()];
                    if ( fmod(width,mPeriod) > 1.0e-10 )
                    {
                        EXCEPTION("The periodic boundary width must be a multiple of membrane period in periodic dimensions.");
                    }
                }   
            }         
		}
	}
}

template<unsigned DIM>
double SinusoidalBaseMembraneBoundaryCondition<DIM>::BaseShapeFunction(c_vector<double,DIM> p)
{
    double z = 0.0;
    for ( unsigned i = 0; i < DIM; i++ )
    {
        if (i != this->mVert)
        {
            double x = p[i];
            z += 1.0 + sin(x*2.0*M_PI/mPeriod);
        }
    }
    z *= mAmplitude/((double)(DIM-1)); // Should never have DIM 1 due to assert in Abstract class

    return(z);
}

template<unsigned DIM>
c_vector<double,DIM> SinusoidalBaseMembraneBoundaryCondition<DIM>::CalculateDerivativesAtPoint(c_vector<double,DIM> p)
{
    c_vector<double,DIM> dzd;
    for (unsigned int i = 0; i < DIM; i++)
    {
        if ( i==(this->mVert) )
        {
            dzd[i] = 1.0;
        }
        else
        {
            dzd[i] = ( 2.0*_PI*mAmplitude / (mPeriod * (double)(DIM-1) ) ) * cos(2.0*_PI*p[i]/mPeriod);   
        }
    }
    return(dzd);
}

template<unsigned DIM>
void SinusoidalBaseMembraneBoundaryCondition<DIM>::SetCurveParameters(double amplitude, double period)
{
    mAmplitude = amplitude;
    mPeriod = period;
}

template<unsigned DIM>
double SinusoidalBaseMembraneBoundaryCondition<DIM>::GetAmplitude() const
{
    return mAmplitude;
}

template<unsigned DIM>
double SinusoidalBaseMembraneBoundaryCondition<DIM>::GetPeriod() const
{
    return mPeriod;
}

// Explicit instantiation
template class SinusoidalBaseMembraneBoundaryCondition<1>;
template class SinusoidalBaseMembraneBoundaryCondition<2>;
template class SinusoidalBaseMembraneBoundaryCondition<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SinusoidalBaseMembraneBoundaryCondition)
