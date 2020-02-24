
#ifndef LIREPULSIONFORCE_HPP_
#define LIREPULSIONFORCE_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"

#include "FakePetscSetup.hpp"


template<unsigned DIM>
class LiRepulsionForce : public RepulsionForce<DIM>
{
	/* Repulsive force taken from Li et al. 2013, Skin Stem Cell Hypotheses and Long Term 
	Clone Survival - Explored Using Agent-based Modelling */
private:
	// The contact stiffness
	double mContactStiffness;

	// Add archiving functions
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<RepulsionForce<DIM> >(*this);
		archive & mContactStiffness;
	}

public:
	/**
	 * Default constructor
	 */
	LiRepulsionForce();

	/**
	 * Overridden AddForceContribution method
	 * Checks that it is a node base cell population then runs the base method from AbstractTwoBody
	 * 
	 * @param rCellPopulation a reference to the cell population
	 */ 
	void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

	/**
	 * Overridden CalculateForceBetweenNodes method
	 *
	 * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     *
     * @return The force exerted on Node A by Node B.
	 */

	c_vector<double, DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM, DIM>& rCellPopulation);

	/**
	 * Set the contact stiffness
	 * @param contactStiffness the new value for the contact stiffness
	 */
	void SetContactStiffness(double contactStiffness);

	/**
	 * Get methods for the member functions
	 * @return the value of mContactStiffness
	 */
	double GetContactStiffness() const;

	/**
	 * Overridden OutputForceParameters method
	 */
	void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LiRepulsionForce)

#endif /*LIREPULSIONFORCE_HPP_*/