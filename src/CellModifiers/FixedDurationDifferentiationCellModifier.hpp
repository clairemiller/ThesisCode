
#ifndef FIXEDDURATIONDIFFERENTIATIONCELLMODIFIER_HPP_
#define FIXEDDURATIONDIFFERENTIATIONCELLMODIFIER_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "CellsGenerator.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "RepulsionForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"


// Cell modifier that converts transit cells to differentiated cell after a specified duration of time
// Currently only implemented in 2D

template<unsigned DIM>
class FixedDurationDifferentiationCellModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{

	// Life span of a transit cell
	double mTransitTime;

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
		archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM> >(*this);

		// Archive member variable
		archive & mTransitTime;
	}

public:
	// Default constructor
	FixedDurationDifferentiationCellModifier(double transitTime);

	// Default destructor
	~FixedDurationDifferentiationCellModifier(){}

	// Override UpdateAtEndOfTImeStep method
	void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

	// Override output method
	void OutputSimulationModifierParameters(out_stream& rParamsFile);

	// Get method for transit time
	double GetTransitTime() const;

	// Set method for transit time
	void SetTransitTime(double transitTime);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixedDurationDifferentiationCellModifier)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const FixedDurationDifferentiationCellModifier<DIM> * t, const unsigned int file_version)
		{
		    double transit_time = t->GetTransitTime();
		    ar << transit_time;
		    
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, FixedDurationDifferentiationCellModifier<DIM> * t, const unsigned int file_version)
		{
            double transit_time;
		    ar >> transit_time;

		    ::new(t)FixedDurationDifferentiationCellModifier<DIM>(transit_time);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*FIXEDDURATIONDIFFERENTIATIONCELLMODIFIER_HPP_*/
