
#ifndef STOCHASTICEXPONENTIALDECAYPALSSONADHESIONFORCE_HPP_
#define STOCHASTICEXPONENTIALDECAYPALSSONADHESIONFORCE_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractGradientPalssonAdhesionForce.hpp"


#include "FakePetscSetup.hpp"


template<unsigned DIM>
class StochasticExponentialDecayPalssonAdhesionForce : public AbstractGradientPalssonAdhesionForce<DIM>
{
private:
	// Member parameters
	double mMinDecayRate;
	double mMaxDecayRate;
	bool mIsInitialised;

	friend class TestStochasticExponentialDecayPalssonAdhesionForce;
	// Add archiving functions
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<AbstractGradientPalssonAdhesionForce<DIM> >(*this);
		archive & mMinDecayRate;
		archive & mMaxDecayRate;
		archive & mIsInitialised;
	}

public:
	/**
	 * Default constructor
	 */
	StochasticExponentialDecayPalssonAdhesionForce(double peakForce, double decayStartAge, double minDecayRate, double maxDecayRate);

	/**
	 * Special member method to assign a randomly generated decay rate to a cell
	 * @param pCell a cell pointer to the cell
	 */
	void AssignCellDecayRate(CellPtr pCell);

	/**
	 * Overridden method for calculate force between nodes to include initialisation
	 */
	virtual c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM, DIM>& rCellPopulation);
    
	/**
     * The scaling function after decay start age
     * @param pCell a pointer to the cell to get scaling function value for
     */
    virtual double ScalingFunction(CellPtr pCell);

	/**
	 * Get method for member variable mMinDecayRate
	 * @return the value of mMinDecayRate
	 */
	double GetMinDecayRate() const;
	
	/**
	 * Set method for member variable mMinDecayRate
	 * @param minDecayRate the new value of mMinDecayRate
	 */
	void SetMinDecayRate(double minDecayRate);

	/**
	 * Get method for member variable mMaxDecayRate
	 * @return the value of mMaxDecayRate
	 */
	double GetMaxDecayRate() const;
	
	/**
	 * Set method for member variable mMaxDecayRate
	 * @param maxDecayRate the new value of mMaxDecayRate
	 */
	void SetMaxDecayRate(double maxDecayRate);

	/**
	 * Get method to get whether the force is initialised or not
	 * @return the value of mIsInitialised
	 */
	bool IsInitialised() const;

	/**
	 * Set method for mIsInitialised
	 * @param isInitialised the new value of mIsInitialised
	 */
	void SetInitialised(bool isInitialised);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(StochasticExponentialDecayPalssonAdhesionForce)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const StochasticExponentialDecayPalssonAdhesionForce<DIM> * t, const unsigned int file_version)
		{
			double peak_force = t->GetPeakForce();
			double decay_start_age = t->GetDecayStartValue();
			double min_decay_rate = t->GetMinDecayRate();
			double max_decay_rate = t->GetMaxDecayRate();
			bool is_initialised = t->IsInitialised();
			ar << peak_force;
			ar << decay_start_age;
			ar << min_decay_rate;
			ar << max_decay_rate;
			ar << is_initialised;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, StochasticExponentialDecayPalssonAdhesionForce<DIM> * t, const unsigned int file_version)
		{
			double min_decay_rate, max_decay_rate, decay_start_value, peak_force;
			bool is_initialised;
			ar >> peak_force;
			ar >> decay_start_value;
			ar >> min_decay_rate;
			ar >> max_decay_rate;
			ar >> is_initialised;

			::new(t)StochasticExponentialDecayPalssonAdhesionForce<DIM>(peak_force, decay_start_value, min_decay_rate, max_decay_rate);
			t->SetInitialised(is_initialised);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*STOCHASTICEXPONENTIALDECAYPALSSONADHESIONFORCE_HPP_*/