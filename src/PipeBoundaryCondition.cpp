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

#include "PipeBoundaryCondition.hpp"
#include "VertexBasedCellPopulation.hpp"

template <unsigned DIM>
PipeBoundaryCondition<DIM>::PipeBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation)
{
    // This class is implemented only for 2D vertex-based cell populations
    assert(dynamic_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation));
    assert(DIM == 2);
}

template<unsigned DIM>
void PipeBoundaryCondition<DIM>::ImposeBoundaryCondition()
{
    VertexBasedCellPopulation<DIM>* mpStaticCastCellPopulation = static_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation);

    // Iterate over all nodes and update their positions according to the boundary conditions
    unsigned num_nodes = mpStaticCastCellPopulation->GetNumNodes();
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_node = mpStaticCastCellPopulation->GetNode(node_index);
        c_vector<double, DIM> node_location = p_node->rGetLocation();

        double x = node_location[0];
        double y = node_location[1];

        double inner_radius = 2.0;
        double outer_radius = 3*inner_radius;
        double height = 3.0;

        c_vector<double, 2> base_centre;
        c_vector<double, 2> top_centre;

        base_centre[0] = inner_radius;
        base_centre[1] = inner_radius;

        top_centre[0] = 3*inner_radius;
        top_centre[1] = inner_radius+height;

        if (y<inner_radius && x<3*inner_radius )
        {
            double radius = norm_2(node_location-base_centre);

            if (radius > inner_radius)
            {
                p_node->rGetModifiableLocation() = base_centre + inner_radius/radius*(node_location-base_centre);
            }
        }
        else if (y > inner_radius+height)
        {
            double radius = norm_2(node_location-top_centre);

            if (radius < inner_radius)
            {
                p_node->rGetModifiableLocation() = top_centre + inner_radius/radius*(node_location-top_centre);
            }
            if (radius > outer_radius)
            {
                p_node->rGetModifiableLocation() = top_centre + outer_radius/radius*(node_location-top_centre);
            }
        }
        else
        {
            if (x<3.0*inner_radius)
            {
                if (x < 0)
                {
                    p_node->rGetModifiableLocation()[0] = 0.0;
                }
                if (x > 2.0*inner_radius)
                {
                    p_node->rGetModifiableLocation()[0] = 2.0*inner_radius;
                }
            }
            else
            {
                if (x < 4.0*inner_radius)
                {
                    p_node->rGetModifiableLocation()[0] = 4.0*inner_radius;
                }
                if (x > 6.0*inner_radius)
                {
                    p_node->rGetModifiableLocation()[0] = 6.0*inner_radius;
                }
            }

        }
    }
}

template<unsigned DIM>
bool PipeBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    return true;
}

template<unsigned DIM>
void PipeBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class PipeBoundaryCondition<1>;
template class PipeBoundaryCondition<2>;
template class PipeBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PipeBoundaryCondition)
