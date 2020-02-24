
#include "LiRepulsionForce.hpp"
#include "MeshBasedCellPopulation.hpp"

// Default constructor
template<unsigned DIM>
LiRepulsionForce<DIM>::LiRepulsionForce()
 : RepulsionForce<DIM>(), mContactStiffness(0.3)
{
}

// Overridden AddForceContribution method
template<unsigned DIM>
void LiRepulsionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
	// Check that it is a node based cell population
    if (!bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)))
    {
        EXCEPTION("LiRepulsionForce: currently only valid for node based population.\n");
    }
	
    // Call the base member class
	RepulsionForce<DIM>::AddForceContribution(rCellPopulation);
}

// Overridden CalculateForceBetweenNodes method
// NOTE: Currently only implemented for node based cell populations
template<unsigned DIM>
c_vector<double, DIM> LiRepulsionForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex, unsigned nodeBGlobalIndex, AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{   
	// We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, DIM> node_a_location = p_node_a->rGetLocation();
    c_vector<double, DIM> node_b_location = p_node_b->rGetLocation();

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, DIM> unit_difference;
    /*
     * We use the mesh method GetVectorFromAtoB() to compute the direction of the
     * unit vector along the line joining the two nodes, rather than simply subtract
     * their positions, because this method can be overloaded (e.g. to enforce a
     * periodic boundary in Cylindrical2dMesh).
     */
    unit_difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    // Calculate the spring length
    double rest_length = ( p_node_a->GetRadius() ) + ( p_node_b->GetRadius() );
    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(p_node_a->GetIndex());
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(p_node_b->GetIndex());
    double ageA = p_cell_A->GetAge();
    double ageB = p_cell_B->GetAge();
    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));
    /*
    * If the cells are both newly divided, then the rest length of the spring
    * connecting them grows linearly with time, until 1 hour after division.
    */
    if (ageA < (this->mMeinekeSpringGrowthDuration) && ageB < (this->mMeinekeSpringGrowthDuration))
    {
        AbstractCentreBasedCellPopulation<DIM,DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<DIM,DIM>*>(&rCellPopulation);
        std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

        if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
        {
            // Spring rest length increases from a small value to the normal rest length over 1 hour
            double lambda = this->mMeinekeDivisionRestingSpringLength;
            rest_length = lambda + (rest_length - lambda) * ageA/(this->mMeinekeSpringGrowthDuration);
        }
        // Spring already unmarked in RepulsionForce AddForceContribution
    }

    // Subtract the spring length and check they are overlapping
    distance_between_nodes = distance_between_nodes - rest_length;
    assert( distance_between_nodes <= 0.0 );

    // Calculate the force
    double Fij = mContactStiffness*distance_between_nodes*2.0/rest_length;
    
    // Convert to the vector and return
    return ( Fij*unit_difference );
}

template<unsigned DIM>
void LiRepulsionForce<DIM>::SetContactStiffness(double contactStiffness)
{
    assert(contactStiffness>0.0);
    mContactStiffness = contactStiffness;
}

template<unsigned DIM>
double LiRepulsionForce<DIM>::GetContactStiffness() const
{
	return mContactStiffness;
}

template<unsigned DIM>
void LiRepulsionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<ContactStiffness>" << mContactStiffness << "</ContactStiffness>\n";
	RepulsionForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class LiRepulsionForce<1>;
template class LiRepulsionForce<2>;
template class LiRepulsionForce<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LiRepulsionForce)