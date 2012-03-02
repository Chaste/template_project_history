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

#include "ObstructionBoundaryCondition.hpp"
#include "VertexBasedCellPopulation.hpp"

template <unsigned DIM>
ObstructionBoundaryCondition<DIM>::ObstructionBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation)
{
    // This class is implemented only for 2D vertex-based cell populations
    assert(dynamic_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation));
    assert(DIM == 2);
}

template<unsigned DIM>
void ObstructionBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::vector< c_vector<double, DIM> >& rOldLocations)
{
    VertexBasedCellPopulation<DIM>* mpStaticCastCellPopulation = static_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation);

    // Iterate over all nodes and update their positions according to the boundary conditions
    unsigned num_nodes = mpStaticCastCellPopulation->GetNumNodes();
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_node = mpStaticCastCellPopulation->GetNode(node_index);
        c_vector<double, DIM> node_location = p_node->rGetLocation();

        double obsruction_radius = 1.5;

        c_vector<double, 2> obstruction_centre;
        obstruction_centre[0] = 3.0;
        obstruction_centre[1] = 5.0;

        double radius = norm_2(node_location-obstruction_centre);
        if (radius < obsruction_radius)
        {
            p_node->rGetModifiableLocation() = obstruction_centre + obsruction_radius/radius*(node_location-obstruction_centre);
        }
    }
}

template<unsigned DIM>
bool ObstructionBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    return true;
}

template<unsigned DIM>
void ObstructionBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class ObstructionBoundaryCondition<1>;
template class ObstructionBoundaryCondition<2>;
template class ObstructionBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ObstructionBoundaryCondition)
