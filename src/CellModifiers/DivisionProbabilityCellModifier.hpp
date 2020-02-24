
#ifndef DIVISIONPROBABILITYCELLMODIFIER_HPP_
#define DIVISIONPROBABILITYCELLMODIFIER_HPP_

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
class DivisionProbabilityCellModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
		archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM> >(*this);
		
		// Archive member variables
		archive & mProb_NonContactInhibited;
		archive & mProb_InhibitedStem;
		archive & mProb_InhibitedTransit;
	}

protected:
	// The probability of division
	double mProb_NonContactInhibited;
	double mProb_InhibitedStem;
	double mProb_InhibitedTransit;

public:
	// Default constructor
	DivisionProbabilityCellModifier(double probNonContactInhibited, double probInhibitedStem, double probInhibitedTransit);

	// Default destructor
	~DivisionProbabilityCellModifier(){}

	// Override UpdatAtEndOfTimeStep method
	void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

	// Override output method
	void OutputSimulationModifierParameters(out_stream& rParamsFile);

	/**
	 * Get mProb_NonContactInhibited
	 * @return value of mProb_NonContactInhibited
	 */
	double GetProbNonContactInhibited() const;

	/**
	 * Get mProb_InhibitedStem
	 * @return value of mProb_InhibitedStem
	 */
	double GetProbInhibitedStem() const;

		/**
	 * Get mProb_InhibitedTransit
	 * @return value of mProb_InhibitedTransit
	 */
	double GetProbInhibitedTransit() const;

	/**
	 * Set mProb_NonContactInhibited
	 * @param probNonContactInhibited the value of mProb_NonContactInhibited
	 */
	void SetProbNonContactInhibited(double probNonContactInhibited);

	/**
	 * Set mProb_InhibitedStem
	 * @param probInhibtedStem the value of mProb_InhibitedStem
	 */
	void SetProbInhibitedStem(double probInhibitedStem);

	/**
	 * Set mProb_InhibitedTransit
	 * @param probInhibitedTransit the value of mProb_InhibitedTransit
	 */
	void SetProbInhibitedTransit(double probInhibitedTransit);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DivisionProbabilityCellModifier)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const DivisionProbabilityCellModifier<DIM> * t, const unsigned int file_version)
		{
		     double prob_non_inhibited = t->GetProbNonContactInhibited();
		     double prob_inhibited_stem = t->GetProbInhibitedStem();
		     double prob_inhibited_transit = t->GetProbInhibitedTransit();

		     ar << prob_non_inhibited;
		     ar << prob_inhibited_stem;
		     ar << prob_inhibited_transit;
		    
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, DivisionProbabilityCellModifier<DIM> * t, const unsigned int file_version)
		{
            double probabilityNonContactInhibitedDivision, probabilityInibitedStemDivision, probabilityInhibitedTransitDivision;
		    ar >> probabilityNonContactInhibitedDivision;
		    ar >> probabilityInibitedStemDivision;
		    ar >> probabilityInhibitedTransitDivision;

		    ::new(t)DivisionProbabilityCellModifier<DIM>(probabilityNonContactInhibitedDivision, probabilityInibitedStemDivision, probabilityInhibitedTransitDivision);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*DIVISIONPROBABILITYCELLMODIFIER_HPP_*/
