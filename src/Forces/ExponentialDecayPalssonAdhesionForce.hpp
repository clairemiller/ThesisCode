
#ifndef EXPONENTIALDECAYPALSSONADHESIONFORCE_HPP_
#define EXPONENTIALDECAYPALSSONADHESIONFORCE_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractGradientPalssonAdhesionForce.hpp"


#include "FakePetscSetup.hpp"


template<unsigned DIM>
class ExponentialDecayPalssonAdhesionForce : public AbstractGradientPalssonAdhesionForce<DIM>
{
private:
	// Member parameters
	double mDecayRate;

	friend class TestExponentialDecayPalssonAdhesionForce;
	// Add archiving functions
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<AbstractGradientPalssonAdhesionForce<DIM> >(*this);
		archive & mDecayRate;
	}

public:
	/**
	 * Default constructor
	 */
	ExponentialDecayPalssonAdhesionForce(double peakForce, double decayStartAge, double decayRate);

    /**
     * The scaling function after decay start age
     * @param pCell a pointer to the cell to get scaling function value for
     */
    virtual double ScalingFunction(CellPtr pCell);

	/**
	 * Get method for member variable mDecayRate
	 * @return the value of mDecayRate
	 */
	double GetDecayRate() const;
	
	/**
	 * Set method for member variable mDecayRate
	 * @param decayRate the new value of mDecayRate
	 */
	void SetDecayRate(double decayRate);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ExponentialDecayPalssonAdhesionForce)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const ExponentialDecayPalssonAdhesionForce<DIM> * t, const unsigned int file_version)
		{
			double peak_force = t->GetPeakForce();
			double decay_start_age = t->GetDecayStartValue();
			double decay_rate = t->GetDecayRate();
			ar << peak_force;
			ar << decay_start_age;
			ar << decay_rate;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, ExponentialDecayPalssonAdhesionForce<DIM> * t, const unsigned int file_version)
		{
			double decay_rate, decay_start_value, peak_force;
			ar >> peak_force;
			ar >> decay_start_value;
			ar >> decay_rate;

			::new(t)ExponentialDecayPalssonAdhesionForce<DIM>(peak_force, decay_start_value, decay_rate);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*EXPONENTIALDECAYPALSSONADHESIONFORCE_HPP_*/