

#include "PopulationAsymmetryCellCycleModel.hpp"

#include "SimulationTime.hpp"

// Default constructor
PopulationAsymmetryCellCycleModel::PopulationAsymmetryCellCycleModel()
: mPStemToSS(0.1), mPStemToST(0.8), mPStemToTT(0.1), mPTransitToTT(0.1), mPTransitToTD(0.8), mPTransitToDD(0.1), mSymmetricDivision(false)
{
    TS_ASSERT_DELTA(mPStemToSS+mPStemToST+mPStemToTT, 1.0, 1e-6);
    TS_ASSERT_DELTA(mPTransitToTT+mPTransitToTD+mPTransitToDD, 1.0, 1e-6);
}

// Copy constructor
PopulationAsymmetryCellCycleModel::PopulationAsymmetryCellCycleModel(const PopulationAsymmetryCellCycleModel& rModel)
 : mPStemToSS(rModel.mPStemToSS), mPStemToST(rModel.mPStemToST), mPStemToTT(rModel.mPStemToTT), 
    mPTransitToTT(rModel.mPTransitToTT), mPTransitToTD(rModel.mPTransitToTD), mPTransitToDD(rModel.mPTransitToDD),
    mSymmetricDivision(rModel.mSymmetricDivision)
{}

// Get function for use within the class only
bool PopulationAsymmetryCellCycleModel::GetIfSymmetricDivision() const
{
    return mSymmetricDivision;
}

// Copy constructor
AbstractCellCycleModel* PopulationAsymmetryCellCycleModel::CreateCellCycleModel()
{
    return new PopulationAsymmetryCellCycleModel(*this);
}

// Get functions
double PopulationAsymmetryCellCycleModel::GetProbStemToSS() const
{
    return mPStemToSS;
}

double PopulationAsymmetryCellCycleModel::GetProbStemToST() const
{
    return mPStemToST;
}

double PopulationAsymmetryCellCycleModel::GetProbStemToTT() const
{
    return mPStemToTT;
}

double PopulationAsymmetryCellCycleModel::GetProbTransitToTT() const
{
    return mPTransitToTT;
}

double PopulationAsymmetryCellCycleModel::GetProbTransitToTD() const
{
    return mPTransitToTD;
}

double PopulationAsymmetryCellCycleModel::GetProbTransitToDD() const
{
    return mPTransitToDD;
}

// Set functions
void PopulationAsymmetryCellCycleModel::SetProbStemToSS(double pStemToSS) 
{
    mPStemToSS = pStemToSS;
}

void PopulationAsymmetryCellCycleModel::SetProbStemToST(double pStemToST)
{
    mPStemToST = pStemToST;
}

void PopulationAsymmetryCellCycleModel::SetProbStemToTT(double pStemToTT)
{
    mPStemToTT = pStemToTT;
}

void PopulationAsymmetryCellCycleModel::SetProbTransitToTT(double pTransitToTT)
{
    mPTransitToTT = pTransitToTT;
}

void PopulationAsymmetryCellCycleModel::SetProbTransitToTD(double pTransitToTD)
{
    mPTransitToTD = pTransitToTD;
}

void PopulationAsymmetryCellCycleModel::SetProbTransitToDD(double pTransitToDD)
{
    mPTransitToDD = pTransitToDD;
}

// Override the ResetForDivision method
void PopulationAsymmetryCellCycleModel::ResetForDivision()
{
    // Generate a random number
    double v_rand = RandomNumberGenerator::Instance()->ranf();

	// Determine the original cell type
    // 1. Stem Cell:
	if ( mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
	{
		if ( v_rand <= mPStemToSS )
		{
            // Symmetric division with two stem cells
            mSymmetricDivision = true;
        }
		else if (v_rand <= (mPStemToSS + mPStemToST) )
        {
            // Asymmetric division: parent remains as stem cell and the daughter will be transit
            mSymmetricDivision = false;
        }
		else {
			// Symmetric division with two transit cells
            mSymmetricDivision = true;
            // Convert parent into transit
			boost::shared_ptr<AbstractCellProperty> p_transit_type =
            	mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        	mpCell->SetCellProliferativeType(p_transit_type);		
		}
	}
    // 2. Transit Cell:
    else if ( mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() )
    {
        if ( v_rand <= mPTransitToTT )
        {
            // Symmetric division staying as transit cells
            mSymmetricDivision = true;
        }
        else if ( v_rand <= (mPTransitToTT + mPTransitToTD) )
        { 
            // Asymmetric division to transit and differentiated
            mSymmetricDivision = false;
        }   
        else
        {
            // Symmetric division to two differentiated
            mSymmetricDivision = true;
            // Convert the parent cell to a differentiated cell
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff_type);          
        }
    }
    // Differentiated cell
    else
    {
        NEVER_REACHED;
    }

	AbstractSimplePhaseBasedCellCycleModel::ResetForDivision();
}

// Override the InitialiseDaughterCell
void PopulationAsymmetryCellCycleModel::InitialiseDaughterCell()
{

    // Get a random value
    double v_rand = RandomNumberGenerator::Instance()->ranf();

	/* We can determine the type of division occuring from the parent type 
    and the symmetric division flag
    If it is symmetric division then daughter remains the same as the parent
    Otherwise we need to modify the daughter type*/
    
    if ( !mSymmetricDivision )
    {
        // If the parent is a stem cell, the daughter becomes transit
        if ( mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>() )
        {
            boost::shared_ptr<AbstractCellProperty> p_transit_type =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_transit_type);
        }
        // If the parent is a transit cell, the daughter becomes differentiated
        else if (  mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() )
        {
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff_type);
        }
        // Should never reach this point
        else
        {
            NEVER_REACHED;
        }
    }

    // Call parent method
	AbstractSimplePhaseBasedCellCycleModel::InitialiseDaughterCell();
}

/**
 * Outputs cell cycle model parameters to file.
 *
 * @param rParamsFile the file stream to which the parameters are output
 */
void PopulationAsymmetryCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // Output member variables
    *rParamsFile << "\t\t\t<ProbabilityStemToSS>" << mPStemToSS << "</ProbabilityStemToSS>\n";  
    *rParamsFile << "\t\t\t<ProbabilityStemToST>" << mPStemToST << "</ProbabilityStemToST>\n";    
    *rParamsFile << "\t\t\t<ProbabilityStemToTT>" << mPStemToTT << "</ProbabilityStemToTT>\n";      
    *rParamsFile << "\t\t\t<ProbabilityTransitToTT>" << mPTransitToTT << "</ProbabilityTransitToTT>\n";  
    *rParamsFile << "\t\t\t<ProbabilityTransitToTD>" << mPTransitToTD << "</ProbabilityTransitToTD>\n";  
    *rParamsFile << "\t\t\t<ProbabilityTransitToDD>" << mPTransitToDD << "</ProbabilityTransitToDD>\n";   

    // Call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PopulationAsymmetryCellCycleModel)