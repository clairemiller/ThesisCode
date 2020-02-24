/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef LOCALCELLDENSITYCELLMODIFIER_HPP_
#define LOCALCELLDENSITYCELLMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include "AbstractCellPopulation.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NodesOnlyMesh.hpp"
#include "PeriodicNdNodesOnlyMesh.hpp"

/**
 * A modifier class to calculate domain density, from: 
 * Yamaguchi et al. (2017) Dynamical crossover in a stochastic model of cell fate decision.
 */
template<unsigned  DIM>
class LocalCellDensityCellModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
      archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM> >(*this);
      archive & mGridSize;
    }

    // Member parameters
    double mGridSize;
    c_vector<unsigned,DIM> mDomainSize;
    std::vector<unsigned> mLocalDensityGrid;
    c_vector<bool, DIM> mIsDimPeriodic;

    // Private member function: Access grid location
    unsigned& AccessGridCell(c_vector<unsigned, DIM> indices);

public:

    /**
     * Default constructor.
     */
    LocalCellDensityCellModifier(double gridSize = 1);

    /**
     * Destructor.
     */
    virtual ~LocalCellDensityCellModifier();

    /**
     * Specify what to do in the simulation at the end of each timestep.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Specify what to do in the simulation at the end of each output timestep,
     * after UpdateAtEndOfTimeStep() has been called.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Get the local cell density around specified location
     * 
     * @param Cell location
     * @param The radius of interest
     * @return the local cell density
     */
    double GetCellDensity(c_vector<double,DIM> cell_location);

    /**
     * Specify what to do in the simulation before the start of the time loop.
     * Required override
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Initialisation of cell data for cell creation
     * 
     * @param rCellPopulation the reference to the cell population
     * @param p the density to use for initialisation
     */
    void SetAllCellDataToValue(AbstractCellPopulation<DIM,DIM>& rCellPopulation, double p);

    /**
     * Output any simulation modifier parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);

    /**
     * Get method for the grid size variable
     * @return the value of mGridSize
     */
    double GetGridSize() const;

    /**
     * Set method for the grid size variable
     * @param the new value of mGridSize
     */
    void SetGridSize(double GridSize);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LocalCellDensityCellModifier)

#endif /*LOCALCELLDENSITYCELLMODIFIER_HPP_*/
