/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef VERTICALTOMEMBRANEDIVISIONRULE_HPP_
#define VERTICALTOMEMBRANEDIVISIONRULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCentreBasedDivisionRule.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "AbstractUndulatingBaseMembraneForce.hpp"

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule;

/**
 * A class to generate two daughter cell positions, located a distance
 * AbstractCentreBasedCellPopulation::mMeinekeDivisionSeparation apart,
 * along a random axis. The midpoint between the two daughter cell
 * positions corresponds to the parent cell's position.
 */
template<unsigned DIM>
class VerticalToMembraneDivisionRule : public AbstractCentreBasedDivisionRule<DIM, DIM>
{
private:
    /**
     * The specified division direction of the new daughter cell.
     * Initialized in the constructor.
     */
    boost::shared_ptr<AbstractUndulatingBaseMembraneForce<DIM> > mpMembrane;
    double mDivisionSpringLength;
    unsigned mVertical;
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<DIM, DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    VerticalToMembraneDivisionRule(boost::shared_ptr<AbstractUndulatingBaseMembraneForce<DIM> > pMembrane, double divisionSpringLength);

    /**
     * Empty destructor.
     */
    virtual ~VerticalToMembraneDivisionRule()
    {
    }

    /**
     * Get method for the membrane
     */
    const boost::shared_ptr<AbstractUndulatingBaseMembraneForce<DIM> > GetMembrane() const;

     /**
     * Get method for the spring length at division
     */
    double GetDivisionSpringLength() const;


    /**
     * Overridden CalculateCellDivisionVector() method.
     *
     * @param pParentCell  The cell to divide
     * @param rCellPopulation  The centre-based cell population
     *
     * @return the two daughter cell positions.
     */
    virtual std::pair<c_vector<double, DIM>, c_vector<double, DIM> > CalculateCellDivisionVector(CellPtr pParentCell,
        AbstractCentreBasedCellPopulation<DIM, DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VerticalToMembraneDivisionRule)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a FixedCentreBasedDivisionRule.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const VerticalToMembraneDivisionRule<DIM>* t, const unsigned int file_version)
{
    // Archive c_vector one component at a time
    boost::shared_ptr<AbstractUndulatingBaseMembraneForce<DIM> > membrane = t->GetMembrane();
    ar << membrane;

    double divisionSpringLength = t->GetDivisionSpringLength();
    ar << divisionSpringLength;
}

/**
 * De-serialize constructor parameters and initialize a FixedCentreBasedDivisionRule.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, VerticalToMembraneDivisionRule<DIM>* t, const unsigned int file_version)
{
    boost::shared_ptr<AbstractUndulatingBaseMembraneForce<DIM> > membrane;
    ar >> membrane;

    double divisionSpringLength;
    ar >> divisionSpringLength;
    
    // Invoke inplace constructor to initialise instance
    ::new(t)VerticalToMembraneDivisionRule<DIM>(membrane,divisionSpringLength);
}
}
} // namespace ...

#endif // VERTICALTOMEMBRANEDIVISIONRULE_HPP_
