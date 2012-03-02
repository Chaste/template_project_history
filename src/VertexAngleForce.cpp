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

#include "VertexAngleForce.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned DIM>
VertexAngleForce<DIM>::VertexAngleForce()
   : AbstractForce<DIM>(),
     mAngleRestrainingParameter(1.0)
{
}

template<unsigned DIM>
VertexAngleForce<DIM>::~VertexAngleForce()
{
}

template<unsigned DIM>
void VertexAngleForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                 AbstractCellPopulation<DIM>& rCellPopulation)
{
    assert(DIM == 2); // this method only works in 2D at present

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    /*
     * The force on each node is given by the interaction between the area and
     * the perimeter of the element containing the node.
     */

    // Iterate over elements in the cell population
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator element_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         element_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++element_iter)
    {
        unsigned element_index = element_iter->GetIndex();

        c_vector<double, DIM> centroid = p_cell_population->rGetMesh().GetCentroidOfElement(element_index);

        double dbl_num_nodes = (double) element_iter->GetNumNodes();

        unsigned num_nodes = element_iter->GetNumNodes();
        for (unsigned node_local_index = 0; node_local_index < num_nodes; node_local_index++)
        {
            unsigned node_global_index = element_iter->GetNodeGlobalIndex(node_local_index);

            c_vector<double, DIM> current_node = element_iter->GetNodeLocation(node_local_index);
            c_vector<double, DIM> next_node = element_iter->GetNodeLocation((node_local_index + 1)%(element_iter->GetNumNodes()));
            c_vector<double, DIM> previous_node = element_iter->GetNodeLocation((node_local_index + element_iter->GetNumNodes() - 1)%(element_iter->GetNumNodes()));

            c_vector<double, DIM> clockwise_unit_vector = p_cell_population->rGetMesh().GetVectorFromAtoB(current_node, previous_node);
            double clockwise_edge_length = norm_2(clockwise_unit_vector);
            clockwise_unit_vector /= clockwise_edge_length;
            c_vector<double, DIM> anti_clockwise_unit_vector = p_cell_population->rGetMesh().GetVectorFromAtoB(current_node, next_node);
            double anti_clockwise_edge_length = norm_2(anti_clockwise_unit_vector);
            anti_clockwise_unit_vector /= anti_clockwise_edge_length;

            // Calculate the angle subtended at the node.
            double determinant = 0.0;
            double angle =0.0;
            double target_angle =0.0;

            determinant = -clockwise_unit_vector(0)*anti_clockwise_unit_vector(1) + clockwise_unit_vector(1)*anti_clockwise_unit_vector(0);

            double dot_product =  inner_prod(-clockwise_unit_vector,anti_clockwise_unit_vector);

            if (dot_product < 1 - 1e-12) // as acos(1)=nan not 0
            {
                angle = acos(dot_product);
            }
            else
            {
                angle = 0.0;
            }

            target_angle = 2.0*M_PI/dbl_num_nodes;

            // Calculate the outward normal at the node
            c_vector<double, DIM> outward_normal = -0.5*clockwise_unit_vector - 0.5*anti_clockwise_unit_vector;
            double outward_normal_length = norm_2(outward_normal);
            if (outward_normal_length > 1e-12)
            {
                outward_normal /= outward_normal_length;
            }
            else
            {
                outward_normal(0) = -clockwise_unit_vector(1);
                outward_normal(1) = clockwise_unit_vector(0);
            }

            if (determinant < 0)
            {
                outward_normal = -outward_normal;
            }

            c_vector<double, DIM> angle_contribution = GetAngleRestrainingParameter() * pow((target_angle - angle),3) * outward_normal; // Force outward from centre for positive determinant and inward for negative determinant

            for (unsigned i=0; i<DIM; i++)
            {
                assert(!std::isnan(angle_contribution(i)));
            }

            rForces[node_global_index] += angle_contribution;
        }
    }
}

template<unsigned DIM>
double VertexAngleForce<DIM>::GetAngleRestrainingParameter()
{
    return mAngleRestrainingParameter;
}

template<unsigned DIM>
void VertexAngleForce<DIM>::SetAngleRestrainingParameter(double angleRestrainingParameter)
{
    mAngleRestrainingParameter = angleRestrainingParameter;
}

template<unsigned DIM>
void VertexAngleForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Output member variables then call method on direct parent class
	*rParamsFile <<  "\t\t\t<AngleRestrainingParameter>"<<  mAngleRestrainingParameter << "</AngleRestrainingParameter> \n";

	AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class VertexAngleForce<1>;
template class VertexAngleForce<2>;
template class VertexAngleForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexAngleForce)
