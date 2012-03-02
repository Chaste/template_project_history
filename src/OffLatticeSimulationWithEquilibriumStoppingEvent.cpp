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

#include "OffLatticeSimulationWithEquilibriumStoppingEvent.hpp"

template<unsigned DIM>
bool OffLatticeSimulationWithEquilibriumStoppingEvent<DIM>::StoppingEventHasOccurred()
{
    return (mNumberOfStaticTimes > mNumberOfStaticTimesThreshold);
}

template<unsigned DIM>
OffLatticeSimulationWithEquilibriumStoppingEvent<DIM>::OffLatticeSimulationWithEquilibriumStoppingEvent(AbstractCellPopulation<DIM>& rCellPopulation,
    bool deleteCellPopulationInDestructor,
    bool initialiseCells,
    double movementThreshold,
    unsigned numStaticTimesThreshold)
    : OffLatticeSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells),
      mNumberOfStaticTimes(0),
      mMovementThreshold(movementThreshold),
      mNumberOfStaticTimesThreshold(numStaticTimesThreshold)
{
}

template<unsigned DIM>
void OffLatticeSimulationWithEquilibriumStoppingEvent<DIM>::ApplyCellPopulationBoundaryConditions(const std::vector< c_vector<double, DIM> >& rOldLocations)
{
    bool in_steady_state = true;

    // Iterate over all nodes and calculate the displacement.
    for (unsigned node_index=0; node_index<this->mrCellPopulation.GetNumNodes(); node_index++)
    {
        Node<DIM>* p_node = this->mrCellPopulation.GetNode(node_index);

        // Check if nodes have moved
        c_vector<double, DIM> displacement = p_node->rGetLocation() - rOldLocations[node_index];

        if (norm_2(displacement) > mMovementThreshold)
        {
            in_steady_state = false;
        }
    }

    /*
     * Check that we're not at the first timestep (as then rOldLocations are the same as
     * the current positions, except for the boundary nodes).
     */
    if ((in_steady_state) && (SimulationTime::Instance()->GetTimeStepsElapsed()>1))
    {
        mNumberOfStaticTimes++;
    }
    else
    {
        mNumberOfStaticTimes = 0;
    }
}

template<unsigned DIM>
double OffLatticeSimulationWithEquilibriumStoppingEvent<DIM>::GetMovementThreshold()
{
	return mMovementThreshold;
}

template<unsigned DIM>
unsigned OffLatticeSimulationWithEquilibriumStoppingEvent<DIM>::GetNumberOfStaticTimesThreshold()
{
	return mNumberOfStaticTimesThreshold;
}

template<unsigned DIM>
void OffLatticeSimulationWithEquilibriumStoppingEvent<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<MovementThreshold>"<< mMovementThreshold << "</MovementThreshold>\n";
    *rParamsFile << "\t\t<mNumberOfStaticTimesThreshold>"<< mNumberOfStaticTimesThreshold << "</mNumberOfStaticTimesThreshold>\n";

    // Call method on direct parent class
    OffLatticeSimulation<DIM>::OutputSimulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OffLatticeSimulationWithEquilibriumStoppingEvent<1>;
template class OffLatticeSimulationWithEquilibriumStoppingEvent<2>;
template class OffLatticeSimulationWithEquilibriumStoppingEvent<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OffLatticeSimulationWithEquilibriumStoppingEvent)
