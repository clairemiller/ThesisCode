
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
	double mDecaySlope;

	friend class TestExponentialDecayPalssonAdhesionForce;
	// Add archiving functions
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<AbstractGradientPalssonAdhesionForce<DIM> >(*this);
		archive & mDecaySlope;
	}

public:
	/**
	 * Default constructor
	 */
	ExponentialDecayPalssonAdhesionForce(double peakForce, double decayStartAge, double decaySlope);

    /**
     * The scaling function after decay start age
     * @param pCell a pointer to the cell to get scaling function value for
     */
    virtual double ScalingFunction(CellPtr pCell);

	/**
	 * Get method for member variable mDecaySlope
	 * @return the value of mDecaySlope
	 */
	double GetDecaySlope() const;
	
	/**
	 * Set method for member variable mDecaySlope
	 * @param decaySlope the new value of mDecaySlope
	 */
	void SetDecaySlope(double decaySlope);
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
			double decay_slope = t->GetDecaySlope();
			ar << peak_force;
			ar << decay_start_age;
			ar << decay_slope;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, ExponentialDecayPalssonAdhesionForce<DIM> * t, const unsigned int file_version)
		{
			double decay_slope, decay_start_value, peak_force;
			ar >> peak_force;
			ar >> decay_start_value;
			ar >> decay_slope;

			::new(t)ExponentialDecayPalssonAdhesionForce<DIM>(peak_force, decay_start_value, decay_slope);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*EXPONENTIALDECAYPALSSONADHESIONFORCE_HPP_*/