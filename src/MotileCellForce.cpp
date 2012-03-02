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

#include "MotileCellForce.hpp"
#include "MotileMutationState.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned DIM>
MotileCellForce<DIM>::MotileCellForce()
    : AbstractForce<DIM>()
{
}

template<unsigned DIM>
MotileCellForce<DIM>::~MotileCellForce()
{
}

template<unsigned DIM>
void MotileCellForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                 AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("MotileCellForce is to be used with a VertexBasedCellPopulation only");
    }

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<p_cell_population->GetNumNodes(); node_index++)
    {
        bool is_motile_vertex = false;

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over the elements containing the vertex
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element and its index
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned element_index = p_element->GetIndex();

            if (p_cell_population->GetCellUsingLocationIndex(element_index)->template HasCellProperty<MotileMutationState>())
            {
                is_motile_vertex = true;
                break;
            }
        }

        // Only motile cells move
        if (is_motile_vertex)
        {
            c_vector<double,DIM> direction = zero_vector<double>(DIM);
            direction(0) = 1;
            double dt = SimulationTime::Instance()->GetTimeStep();

            double magnitude = 0.01/dt;

            rForces[node_index] += magnitude*direction;
        }
    }
}

template<unsigned DIM>
void MotileCellForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MotileCellForce<1>;
template class MotileCellForce<2>;
template class MotileCellForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MotileCellForce)
