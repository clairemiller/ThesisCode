

// Kills a column of cells in 2D or 3D

#ifndef SCRATCHASSAYCELLKILLER_HPP_
#define SCRATCHASSAYCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCellPopulation.hpp"



/**
 * Cell Killer to replicate a Scratch Assay. Removes a column of cells at a specified location and time.
 * Currently implemented in 2D only.
 */

template< unsigned DIM >

class ScratchAssayCellKiller : public AbstractCellKiller<DIM>
{
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & archive, const unsigned int version)
	{
			archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
	}

protected:
	// Width of column to strip
	double mColWidth;

	// Centre of column to strip
	double mCentre;

	// Height to start column
	double mHeight;

    // Time to start the CellKiller 
    double mKillTime;

public:

	// Default constructor
	ScratchAssayCellKiller(AbstractCellPopulation<DIM,DIM>* pCellPopulation, double colWidth, double centre, double height, double killTime);

	// Return mColWidth
	double GetColWidth() const;

	// Return mCentre
	double GetCentre() const;

	// Return height
	double GetHeight() const;

    // Return mKillTime
    double GetKillTime() const;

	// Determine cells within column
	void CheckAndLabelCellsForApoptosisOrDeath();

	// Override OutputCellKillerParameters() method.
	void OutputCellKillerParameters(out_stream& rParamsFiles);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ScratchAssayCellKiller)

namespace boost
{
	namespace serialization
	{
		template<class Archive, unsigned DIM>
		inline void save_construct_data(Archive & ar, const ScratchAssayCellKiller<DIM> * t, const unsigned int file_version)
		{
			const AbstractCellPopulation<DIM,DIM>* const p_cell_population = t->GetCellPopulation();
		    ar << p_cell_population;
		    double col_width = t->GetColWidth();
		    ar << col_width;
		    double col_centre = t->GetCentre();
		    ar << col_centre;
		    double col_height = t->GetHeight();
		    ar << col_height;
            double kill_time = t->GetKillTime();
            ar << kill_time;
		}

		template<class Archive, unsigned DIM>
		inline void load_construct_data(Archive & ar, ScratchAssayCellKiller<DIM> * t, const unsigned int file_version)
		{
			AbstractCellPopulation<DIM,DIM>* p_cell_population;
			ar >> p_cell_population;
		    
            double col_width, col_centre, col_height, kill_time, duration;
		    ar >> col_width;
		    ar >> col_centre;
		    ar >> col_height;
            ar >> kill_time;

		    ::new(t)ScratchAssayCellKiller<DIM>(p_cell_population, col_width, col_height, col_centre, kill_time);
		}
	} // End namespace serialization
} // End namespace boost

#endif /*SCRATCHASSAYCELLKILLER_HPP_*/



