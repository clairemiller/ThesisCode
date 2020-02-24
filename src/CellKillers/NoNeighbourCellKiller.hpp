
#ifndef NONEIGHBOURCELLKILLER_HPP_
#define NONEIGHBOURCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "LowNeighbourCellMutationState.hpp"

/**
 * Cell Killer to remove cells that have no neighbours within a certain radius
 */

template< unsigned DIM >

class NoNeighbourCellKiller : public AbstractCellKiller<DIM>
{
private:
	// Member variables
	double mNeighbourRadius;
	double mMinHeight;
	unsigned mVertical;

    friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
			archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
	}

public:
    /** 
     * Default constructor
     * @param pCellPopulation the cell population
     */
    NoNeighbourCellKiller(AbstractCellPopulation<DIM,DIM>* pCellPopulation, double neighbourRadius, double minHeight = 0.0, unsigned vertical = DIM-1);

    /**
     * Default deconstructor
     */
    ~NoNeighbourCellKiller(){}

    /**
     * Overridden virtual class
    */
    void CheckAndLabelCellsForApoptosisOrDeath();

	/**
	 * Get method for the neighbour radius
	 * @return mNeighbourRadius the neighbour radius
	 */
	double GetNeighbourRadius() const;

	/**
	 * Get method for the minimum height
	 * @return mMinHeight
	 */
	double GetMinimumHeight() const;

	/**
	 * Get method for the vertical direction
	 * @return mVertical
	 */
	unsigned GetVerticalDirection() const;

    /**
     * Overridden virtual class
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NoNeighbourCellKiller)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const NoNeighbourCellKiller<DIM> * t, const unsigned int file_version)
		{
			const AbstractCellPopulation<DIM,DIM>* const p_cell_population = t->GetCellPopulation();
		    ar << p_cell_population;
			double neighbour_radius = t->GetNeighbourRadius();
			ar << neighbour_radius;
			double min_height = t->GetMinimumHeight();
			ar << min_height;
			unsigned vertical = t->GetVerticalDirection();
			ar << vertical;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, NoNeighbourCellKiller<DIM> * t, const unsigned int file_version)
		{
			AbstractCellPopulation<DIM,DIM>* p_cell_population;
			ar >> p_cell_population;
			double neighbour_radius, min_height;
			ar >> neighbour_radius;
			ar >> min_height;
			unsigned vertical;
			ar >> vertical;

		    ::new(t)NoNeighbourCellKiller<DIM>(p_cell_population, neighbour_radius, min_height, vertical);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*NONEIGHBOURCELLKILLER_HPP_*/

