#ifndef SINUSOIDALBASEMEMBRANEBOUNDARYCONDITION_HPP_
#define SINUSOIDALBASEMEMBRANEBOUNDARYCONDITION_HPP_

#include "FakePetscSetup.hpp"
#include "AbstractUndulatingBaseMembraneBoundaryCondition.hpp"

// Define pi for the purposes of the class
#ifndef _PI
#define _PI 3.14159265359
#endif


template<unsigned DIM>
class SinusoidalBaseMembraneBoundaryCondition : public AbstractUndulatingBaseMembraneBoundaryCondition<DIM>
{
private:
	// Add archiving functions
	friend class boost::serialization::access;
	friend class TestAbstractUndulatingBaseMembraneBoundaryCondition;
	
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<AbstractUndulatingBaseMembraneBoundaryCondition<DIM> >(*this);
		archive & mAmplitude;
		archive & mPeriod;
	}

protected:
	// Member variable for the magnitude
	double mAmplitude;
	double mPeriod;
public:
	/**
	 * Default constructor
	 * We are taking the amplitude and period (scaled from micro-metre to cells) from:
     * Niels Grabe and Karsten Neuber 2005, A Multicellular Systems Biology Model Predicts Epidermal Mor- phology, Kinetics and Ca2+ Flow
	 */
	SinusoidalBaseMembraneBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation, unsigned verticalDirection = DIM-1, double amplitude = 2.0, double period = 7.0);

	/**
	 * Default deconstructor
	*/
	~SinusoidalBaseMembraneBoundaryCondition(){}

	/** 
	 * A function for the membrane shape
	 * @param p an x,y(,z) location to determine the appropriate height given a point location
	*/
	virtual double BaseShapeFunction(c_vector<double,DIM> p);

	/**
	 * Calculates the derivatives to the membrane at a given point
	 * @param the point
	 * @return the vector of derivatives
	 */
	virtual c_vector<double,DIM> CalculateDerivativesAtPoint(c_vector<double,DIM> p);	

	/**
	 * A function to define the parameters for the sinusoidal function
	 * @param amplitude the amplitude of the sine curve
	 * @param period the period of the sine curve
	*/
	void SetCurveParameters(double amplitude, double period);

	/**
	 * Get function for amplitude
	 * @return the value of mAmplitude
	 */
	double GetAmplitude() const;

	/**
	 * Get function for period
	 * @return the value of mPeriod
	 */
	double GetPeriod() const;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SinusoidalBaseMembraneBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a SinusoidalBaseMembraneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SinusoidalBaseMembraneBoundaryCondition<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
	double amplitude = t->GetAmplitude();
	double period = t->GetPeriod();

	ar << p_cell_population;
	ar << amplitude;
	ar << period;
}

/**
 * De-serialize constructor parameters and initialize a SinusoidalBaseMembraneBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, SinusoidalBaseMembraneBoundaryCondition<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;
	double amplitude, period;
	ar >> amplitude;
	ar >> period;

    // Invoke inplace constructor to initialise instance
    ::new(t)SinusoidalBaseMembraneBoundaryCondition<DIM>(p_cell_population);
	t->SetCurveParameters(amplitude,period);
}
}
} // namespace ...

#endif /*SINUSOIDALBASEMEMBRANEBOUNDARYCONDITION_HPP_*/
