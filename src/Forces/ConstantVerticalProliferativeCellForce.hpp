
#ifndef CONSTANTVERTICALPROLIFERATIVECELLFORCE_HPP_
#define CONSTANTVERTICALPROLIFERATIVECELLFORCE_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"

#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

#include "FakePetscSetup.hpp"


template<unsigned DIM>
class ConstantVerticalProliferativeCellForce : public AbstractForce<DIM,DIM>
{
private:
	double mForceMagnitude;
	double mCutOffHeight;
	unsigned mVerticalDirection;

	// Add archiving functions
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version )
	{
		archive & boost::serialization::base_object<AbstractForce<DIM,DIM> >(*this);
		archive & mForceMagnitude;
		archive & mCutOffHeight;
		archive & mVerticalDirection;
	}

public:
	/**
	 * Default constructor
	 * @param strength the value of mStrength
	 * @param radius the value of mRadius
	 */
	ConstantVerticalProliferativeCellForce(double forceMagnitude, double cutOffHeight, unsigned verticalDirection = (DIM-1));

	/**
	 * Default deconstructor
	*/
	~ConstantVerticalProliferativeCellForce(){}

	/**
     * Calculates the force on each node.
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Outputs force parameters to file.
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

	/**
	 * Get methods for the member functions
	 * @return the member variable values
	 */
	double GetForceMagnitude() const;
	unsigned GetVerticalDirection() const;
	double GetCutOffHeight() const;

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ConstantVerticalProliferativeCellForce)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const ConstantVerticalProliferativeCellForce<DIM> * t, const unsigned int file_version)
		{
			double force_magnitude = t->GetForceMagnitude();
			ar << force_magnitude;
			double cutoff_height = t->GetCutOffHeight();            
			ar << cutoff_height;
			unsigned vertical_direction = t->GetVerticalDirection();
			ar << vertical_direction;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, ConstantVerticalProliferativeCellForce<DIM> * t, const unsigned int file_version)
		{
			double force_magnitude,cutoff_height;
			unsigned vertical_direction;

			ar >> force_magnitude;
			ar >> cutoff_height;
			ar >> vertical_direction;
			

		    ::new(t)ConstantVerticalProliferativeCellForce<DIM>(force_magnitude, cutoff_height, vertical_direction);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*CONSTANTVERTICALPROLIFERATIVECELLFORCE_HPP_*/
