
#ifndef HEIGHTDEPENDENTDIVISIONMODIFIER_HPP_
#define HEIGHTDEPENDENTDIVISIONMODIFIER_HPP_

#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "AbstractUndulatingBaseMembraneBoundaryCondition.hpp"

#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"

template<unsigned DIM>
class HeightDependentDivisionModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
		archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM> >(*this);
		
		// Archive member variables
		archive & mHeight;
		archive & mRevertAfterG1Phase;
		archive & mVertical;
	}

protected:
	double mHeight;
	bool mRevertAfterG1Phase;
	unsigned mVertical;
	boost::shared_ptr<AbstractUndulatingBaseMembraneBoundaryCondition<DIM> > mpMembrane;
	bool mUseMembrane;
	unsigned mMaxStemCells;

public:
	// Default constructor
	HeightDependentDivisionModifier(double height);

	// Default destructor
	~HeightDependentDivisionModifier(){}

	// Override UpdatAtEndOfTimeStep method
	void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

	// Override UpdateCellData
	void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

	// Override output method
	void OutputSimulationModifierParameters(out_stream& rParamsFile);

	// Override virtual function
	virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

	/**
	 * Get mHeight
	 * @return value of mHeight
	 */
	double GetHeight() const;

	bool GetRevertAfterG1Phase() const;

	unsigned GetVertical() const;

	bool GetUseMembrane() const;

	/**
	 * Set mHeight
	 * @param Height the value of mHeight
	 */
	void SetHeight(double height);

	void SetRevertAfterG1Phase(bool revertAfterG1Phase);

	void SetVertical(unsigned vertical);

	void SetMembrane(boost::shared_ptr<AbstractUndulatingBaseMembraneBoundaryCondition<DIM> > pMembrane);

	void SetMaxStemCells(unsigned maxCells);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(HeightDependentDivisionModifier)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const HeightDependentDivisionModifier<DIM> * t, const unsigned int file_version)
		{
		     double height = t->GetHeight();
		     bool revert_after_g1 = t->GetRevertAfterG1Phase();
		     unsigned vertical = t->GetVertical();

		     ar << height;
		     ar << revert_after_g1;
		     ar << vertical;
		    
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, HeightDependentDivisionModifier<DIM> * t, const unsigned int file_version)
		{
            double height;
			ar >> height;
			bool revert_after_g1;
			ar >> revert_after_g1;
			unsigned vertical;
		    ar >> vertical;

			::new(t)HeightDependentDivisionModifier<DIM>(height);
			t->SetRevertAfterG1Phase(revert_after_g1);
			t->SetVertical(vertical);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*HEIGHTDEPENDENTDIVISIONMODIFIER_HPP_*/
