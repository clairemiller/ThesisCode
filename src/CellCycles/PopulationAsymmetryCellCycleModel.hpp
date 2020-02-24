
#ifndef POPULATIONASYMMETRYCELLCYCLEMODEL_HH_
#define POPULATIONASYMMETRYCELLCYCLEMODEL_HH_

#include "ChasteSerialization.hpp"
#include <cxxtest/TestSuite.h>

#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

class PopulationAsymmetryCellCycleModel : public AbstractSimplePhaseBasedCellCycleModel
{
private:
	/** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimplePhaseBasedCellCycleModel>(*this);

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;

        // Archive all the member variables
        archive & mPStemToSS;
        archive & mPStemToST; 
        archive & mPStemToTT; 
        archive & mPTransitToTT; 
        archive & mPTransitToTD; 
        archive & mPTransitToDD; 
    }

protected:
    // Probabilities of different daughter cell types
    double mPStemToSS;     // Probability of two stem from stem
    double mPStemToST;      // Probability of stem and transit from stem
    double mPStemToTT;      // Probability of two transit from stem
    double mPTransitToTT;   // Probability of two transit from transit
    double mPTransitToTD;   // Probability of transit and differentiated from transit
    double mPTransitToDD;   // Probability of two differentiated from transit

    /**
     * Whether the cell division was symmetric or asymmetric
     * Added for use with the division vector calculation in AbstractOffLatticeCellPopulation
     */
    bool mSymmetricDivision;

public:
    /* 
     * Default constructor - Probabilities are set manually using the 'Set' methods
     */
    PopulationAsymmetryCellCycleModel();

    // Copy constructor
    PopulationAsymmetryCellCycleModel(const PopulationAsymmetryCellCycleModel& rModel);

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /** 
     * Overriden ResetForDivision method 
     * - Determines proliferative type of the parent cell
     */
    void ResetForDivision();

    /** 
     * Overriden InitialiseDaughterCell method
     * - Determines proliferative type of daughter cell
     */
    void InitialiseDaughterCell();

    /**
     * Set mPStemToSS.
     *
     * @param pStemToSS the value of mPStemToSS
     */
    void SetProbStemToSS(double pStemToSS);

    /**
     * Set mPStemToST.
     *
     * @param pStemToST the value of mPStemToST
     */
    void SetProbStemToST(double pStemToST);

    /**
     * Set mPStemToTT.
     *
     * @param pStemToTT the value of mPStemToTT
     */
    void SetProbStemToTT(double pStemToTT);

    /**
     * Set mPTransitToTT
     *
     * @param pTransitToTT the value of mPTransitToTT
     */
    void SetProbTransitToTT(double pTransitToTT);

    /**
     * Set mPTransitToTD.
     *
     * @param pTransitToTD the value of mPTransitToTD
     */
    void SetProbTransitToTD(double pTransitToTD);

    /**
     * Set mPTransitToDD
     *
     * @param pTransitToDD the value of mPTransitToDD
     */
    void SetProbTransitToDD(double pTransitToDD);

    /**
     * Get mPStemToSS.
     *
     * @return the value of mPStemToSS
     */
    double GetProbStemToSS() const;

    /**
     * Get mPStemToST.
     *
     * @return the value of mPStemToST
     */
    double GetProbStemToST() const;

    /**
     * Get mPStemToTT.
     *
     * @return the value of mPStemToTT
     */
    double GetProbStemToTT() const;

    /**
     * Get mPTransitToTT
     *
     * @return the value of mPTransitToTT
     */
    double GetProbTransitToTT() const;

    /**
     * Get mPTransitToTD.
     *
     * @return the value of mPTransitToTD
     */
    double GetProbTransitToTD() const;

    /**
     * Get mPTransitToDD
     *
     * @return the value of mPTransitToDD
     */
    double GetProbTransitToDD() const;

    /**
     * Get mSymmetricDivision
     *
     * @return the value of mSymmetricDivision
     */
    bool GetIfSymmetricDivision() const;


    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(PopulationAsymmetryCellCycleModel)

#endif /*POPULATIONASYMMETRYCELLCYCLEMODEL_HH_*/
