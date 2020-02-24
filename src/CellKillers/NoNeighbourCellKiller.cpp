
#include "NoNeighbourCellKiller.hpp"

template<unsigned DIM>
NoNeighbourCellKiller<DIM>::NoNeighbourCellKiller(AbstractCellPopulation<DIM,DIM>* pCellPopulation, double neighbourRadius, double minHeight, unsigned vertical)
: AbstractCellKiller<DIM>(pCellPopulation), mNeighbourRadius(neighbourRadius), mMinHeight(minHeight), mVertical(vertical)
{
    NodeBasedCellPopulation<DIM>* pNodePopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation);
    if ( ! (bool) pNodePopulation )
    {
        EXCEPTION("No neighbour cell killer only implemented for node based cell populations.\n");
    }
    assert(vertical < DIM);
    assert(pNodePopulation->GetMechanicsCutOffLength() >= neighbourRadius);
}

template<unsigned DIM>
void NoNeighbourCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    // Loop over cells and label those with no neighbours as apoptotic
    for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
            cell_iter != this->mpCellPopulation->End();
            ++cell_iter)
    {
        // Check that it is at a sufficient height
        double cell_height = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[mVertical];
        if ( !(cell_iter->HasApoptosisBegun()) && cell_height  > mMinHeight )
        {
            //unsigned n_neighbours = cell_iter->GetCellData()->GetItem("n_neighbours");
            NodeBasedCellPopulation<DIM>* pNodePopulation = dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation);
            // Check if it is part of a pair and if so, kill both
            unsigned index = pNodePopulation->GetLocationIndexUsingCell(*cell_iter);
            std::set<unsigned> neighbours = pNodePopulation->GetNodesWithinNeighbourhoodRadius(index,mNeighbourRadius);
            // If it has no neighbours kill it
            if ( neighbours.size() < 2 )
            {
                // We actually want to look at a smaller radius than we did to calculate n_neighbours
                cell_iter -> StartApoptosis();
            }
            else if ( (cell_iter)->GetMutationState()->template IsType<LowNeighbourCellMutationState>() )
            {
                // Determine if it has no neighbours below
                bool no_neighbours_below = true;
                for( std::set<unsigned>::iterator nbr_it = neighbours.begin(); nbr_it != neighbours.end() && no_neighbours_below; ++nbr_it )
                {
                    CellPtr pNeighbourCell = this->mpCellPopulation->GetCellUsingLocationIndex(*nbr_it);
                    double nbr_height = this->mpCellPopulation->GetLocationOfCellCentre(pNeighbourCell)[mVertical];
                    if ( nbr_height < cell_height )
                    {
                        no_neighbours_below = false;
                    }
                }
                if ( no_neighbours_below )
                {
                    cell_iter ->StartApoptosis();
                }
            }
            
        }
    }
}

template<unsigned DIM>
double NoNeighbourCellKiller<DIM>::GetNeighbourRadius() const
{
    return mNeighbourRadius;
}

template<unsigned DIM>
double NoNeighbourCellKiller<DIM>::GetMinimumHeight() const
{
    return mMinHeight;
}

template<unsigned DIM>
unsigned NoNeighbourCellKiller<DIM>::GetVerticalDirection() const
{
    return mVertical;
}

template<unsigned DIM>
void NoNeighbourCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NeighbourRadius>" << mNeighbourRadius << "</NeighbourRadius>\n";
	*rParamsFile << "\t\t\t<MinHeight>" << mMinHeight << "</MinHeight>\n";
	*rParamsFile << "\t\t\t<Vertical>" << mVertical << "</Vertical>\n";
    // Call parent method
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class NoNeighbourCellKiller<1>;
template class NoNeighbourCellKiller<2>;
template class NoNeighbourCellKiller<3>;

// Serialization
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NoNeighbourCellKiller)