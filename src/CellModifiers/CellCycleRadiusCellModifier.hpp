
#ifndef CELLCYCLERADIUSCELLMODIFIER_HPP_
#define CELLCYCLERADIUSCELLMODIFIER_HPP_

#include "RandomNumberGenerator.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "CellCyclePhases.hpp"
#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"

template<unsigned DIM>
class CellCycleRadiusCellModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
		archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM> >(*this);
		
		// Archive member variables
		archive & mRadiusAtDivision;
		archive & mRadiusAfterDivision;
	}

protected:
	// The probability of division
	double mRadiusAtDivision;
	double mRadiusAfterDivision;

public:
	// Default constructor
	CellCycleRadiusCellModifier(double radiusAtDivision=1.0);

	// Default destructor
	~CellCycleRadiusCellModifier(){}

	// Override UpdatAtEndOfTimeStep method
	void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

	// Override output method
	void OutputSimulationModifierParameters(out_stream& rParamsFile);

	// Override initial setup method
	void SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory);

	/**
	 * Get mRadiusAtDivision
	 * @return value of mRadiusAtDivision
	 */
	double GetRadiusAtDivision() const;

	/**
	 * Set mProb_NonContactInhibited
	 * @param probNonContactInhibited the value of mProb_NonContactInhibited
	 */
	void SetRadiusAtDivision(double radiusAtDivision);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellCycleRadiusCellModifier)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const CellCycleRadiusCellModifier<DIM> * t, const unsigned int file_version)
		{
		     double radius_at_division = t->GetRadiusAtDivision();

		     ar << radius_at_division;	    
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, CellCycleRadiusCellModifier<DIM> * t, const unsigned int file_version)
		{
            double radius_at_division;
		    ar >> radius_at_division;

		    ::new(t)CellCycleRadiusCellModifier<DIM>(radius_at_division);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*CELLCYCLERADIUSCELLMODIFIER_HPP_*/
