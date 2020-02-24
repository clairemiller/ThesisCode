#ifndef APPROXGAUSSIANBASEMEMBRANEBOUNDARYCONDITION_HPP_
#define APPROXGAUSSIANBASEMEMBRANEBOUNDARYCONDITION_HPP_

#include "FakePetscSetup.hpp"
#include "SinusoidalBaseMembraneBoundaryCondition.hpp"


template<unsigned DIM>
class ApproxGaussianBaseMembraneBoundaryCondition : public SinusoidalBaseMembraneBoundaryCondition<DIM>
{
private:
	// Add archiving functions
	friend class boost::serialization::access;
	
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<SinusoidalBaseMembraneBoundaryCondition<DIM> >(*this);
	}

public:
	/**
	 * Default constructor
	 */
	ApproxGaussianBaseMembraneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation, unsigned verticalDirection = DIM-1);

	/**
	 * Default deconstructor
	*/
	~ApproxGaussianBaseMembraneBoundaryCondition(){}

	/** 
	 * A function for the membrane shape
	 * @param p an x,y(,z) location to determine the appropriate height given a point location
	*/
	virtual double BaseShapeFunction(c_vector<double,DIM> p);

	/**
	 * Calculates the derivative to the membrane at a given point
	 * @param the point
	 * @return the vector of derivatives
	 */
	virtual c_vector<double,DIM> CalculateDerivativesAtPoint(c_vector<double,DIM> p);	
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ApproxGaussianBaseMembraneBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a ApproxGaussianBaseMembraneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const ApproxGaussianBaseMembraneBoundaryCondition<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
}

/**
 * De-serialize constructor parameters and initialize a ApproxGaussianBaseMembraneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, ApproxGaussianBaseMembraneBoundaryCondition<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)ApproxGaussianBaseMembraneBoundaryCondition<DIM>(p_cell_population);
}
}
} // namespace ...

#endif /*APPROXGAUSSIANBASEMEMBRANEBOUNDARYCONDITION_HPP_*/
