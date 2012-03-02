/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef OFFLATTICESIMULATIONWITHEQUILIBRIUMSTOPPINGEVENT_HPP_
#define OFFLATTICESIMULATIONWITHEQUILIBRIUMSTOPPINGEVENT_HPP_

#include <map>
#include "ChasteSerialization.hpp"

#include "CellBasedSimulation.hpp"

/**
 * An off-lattice cell-based simulation class with a stopping event.
 */
template<unsigned DIM>
class OffLatticeSimulationWithEquilibriumStoppingEvent : public OffLatticeSimulation<DIM>
{
private:

    /** Define a stopping event which says stop if mNumberOfStaticTimes > mNumberOfStaticTimesThreshold. */
    bool StoppingEventHasOccurred();

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<OffLatticeSimulation<DIM> >(*this);
        archive & mNumberOfStaticTimes;
        archive & mMovementThreshold;
        archive & mNumberOfStaticTimesThreshold;
    }

    /** Counter to see if in steady state **/
    unsigned mNumberOfStaticTimes;

public:

    /** Max squared movement of nodes in steady state. */
    double mMovementThreshold;

    /** Number of timesteps with movement below threshold before defined to be in steady state. */
    unsigned mNumberOfStaticTimesThreshold;

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population facade class (contains a mesh and cells)
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells Whether to initialise cells (defaults to true, set to false when loading
     *     from an archive)
     * @param movementThreshold Max squared movement of nodes in steady state (defaults to 1e-6)
     * @param numStaticTimesThreshold Number of timesteps with movement below threshold before defined
     *     to be in steady state (defaults to 1e3)
     */
    OffLatticeSimulationWithEquilibriumStoppingEvent(AbstractCellPopulation<DIM>& rCellPopulation,
												     bool deleteCellPopulationInDestructor=false,
												     bool initialiseCells=true,
                                                     double movementThreshold=1e-6,
                                                     unsigned numStaticTimesThreshold=1e3);

    /**
     * Overridden ApplyCellPopulationBoundaryConditions() method.
     *
     * Calculates the displacement of the cells and increments mNumberOfStaticTimes if neceasry
     *
     * @param rOldLocations the node locations at the previous time step
     */
    void ApplyCellPopulationBoundaryConditions(const std::vector< c_vector<double, DIM> >& rOldLocations);

    /**
     * @return mMovementThreshold
     */
    double GetMovementThreshold();

    /**
     * @return mNumberOfStaticTimesThreshold
     */
    unsigned GetNumberOfStaticTimesThreshold();

    /**
     * Overridden OutputSimulationParameters() method.
     * Outputs simulation parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OffLatticeSimulationWithEquilibriumStoppingEvent)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an OffLatticeSimulationWithEquilibriumStoppingEvent.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const OffLatticeSimulationWithEquilibriumStoppingEvent<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM> * p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise an OffLatticeSimulationWithEquilibriumStoppingEvent.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, OffLatticeSimulationWithEquilibriumStoppingEvent<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)OffLatticeSimulationWithEquilibriumStoppingEvent<DIM>(*p_cell_population, true, false);
}
}
} // namespace ...

#endif /*OFFLATTICESIMULATIONWITHEQUILIBRIUMSTOPPINGEVENT_HPP_*/
