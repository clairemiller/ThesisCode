#ifndef ABSTRACTGRADIENTPALSSONADHESIONFORCE_HPP_
#define ABSTRACTGRADIENTPALSSONADHESIONFORCE_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "PalssonAdhesionForce.hpp"

#include "FakePetscSetup.hpp"
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

template<unsigned DIM>
class AbstractGradientPalssonAdhesionForce : public PalssonAdhesionForce<DIM>
{
private:
	// Add archiving functions
	friend class boost::serialization::access;
	friend class TestAbstractGradientPalssonAdhesionForce;

	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<PalssonAdhesionForce<DIM> >(*this);
        archive & mPeakForce;
        archive & mDecayStartValue;
	}
protected:
	// Member variables
    double mPeakForce;
    double mDecayStartValue;

public:
	/**
	 * Default constructor
	 */
	AbstractGradientPalssonAdhesionForce(double peakForce=0.0, double decayStartValue=0.0);

	/**
	 * Default deconstructor
	*/
	virtual ~AbstractGradientPalssonAdhesionForce(){}

	/**
     * Calculates the force on each node.
     * @param rCellPopulation reference to the cell population
     */
	virtual c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * The scaling function after decay start age
     * @param pCell a pointer to the cell to get scaling function value for
     */
    virtual double ScalingFunction(CellPtr pCell) = 0;

    /**
     * Get method for the peak force
     * @return the value of mPeakForce
     */
    double GetPeakForce() const;

    /**
     * Get method for the decay start age
     * @return the value of mDecayStartValue
     */
    double GetDecayStartValue() const;

    /**
     * Set method for the peak force
     * @param peakForce the new value of mPeakForce
     */
    void SetPeakForce(double peakForce);

    /**
     * Set method for the decay start age
     * @param decayStartValue the new value of mDecayStartValue
     */
    void SetDecayStartValue(double decayStartValue);

    /**
    * Outputs force parameters to file.
    * @param rParamsFile the file stream to which the parameters are output
    */
    virtual void OutputForceParameters(out_stream& rParamsFile);
	 
	/**
	* Write any data necessary to a visualization setup file.
	* Used by AbstractCellBasedSimulation::WriteVisualizerSetupFile().
	* 
	* @param pVizSetupFile a visualization setup file
	*/
	virtual void WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile);
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractGradientPalssonAdhesionForce)


// namespace boost
// {
// namespace serialization
// {
// /**
//  * Serialize information required to construct a AbstractGradientPalssonAdhesionForce.
//  */
// template<class Archive, unsigned DIM>
// inline void save_construct_data(Archive & ar, const AbstractGradientPalssonAdhesionForce<DIM>* t, const unsigned int file_version)
// {
//     // Save data required to construct instance
//     double peakForce = t->GetPeakForce();
//     double decayStartValue = t->GetDecayStartValue();
    
//     ar & peakForce;
//     ar & decayStartValue;
// }

// /**
//  * De-serialize constructor parameters and initialise a AbstractGradientPalssonAdhesionForce.
//  */
// template<class Archive, unsigned DIM>
// inline void load_construct_data( Archive & ar, AbstractGradientPalssonAdhesionForce<DIM>* t, const unsigned int file_version)
// {
//     // Retrieve data from archive required to construct new instance
//     double peak_force, decay_start_value;
//     ar >> peak_force;
//     ar >> decay_start_value;

//     // Invoke inplace constructor to initialise instance
//     ::new(t)AbstractGradientPalssonAdhesionForce<DIM>(peak_force,decay_start_value);
// }
// }
// } // namespace ...

#endif /*ABSTRACTGRADIENTPALSSONADHESIONFORCE_HPP_*/
