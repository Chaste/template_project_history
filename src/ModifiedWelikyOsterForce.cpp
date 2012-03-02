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

#include "ModifiedWelikyOsterForce.hpp"

template<unsigned DIM>
ModifiedWelikyOsterForce<DIM>::ModifiedWelikyOsterForce()
   : AbstractForce<DIM>(),
     mWelikyOsterAreaParameter(1.0),
     mWelikyOsterPerimeterParameter(1.0),
     mMatureCellTargetArea(1.0)
{
}

template<unsigned DIM>
ModifiedWelikyOsterForce<DIM>::~ModifiedWelikyOsterForce()
{
}

template<unsigned DIM>
void ModifiedWelikyOsterForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                 AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
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

        /******** Start of deformation force calculation ********/

        // Compute the area of this element
        double element_area = p_cell_population->rGetMesh().GetVolumeOfElement(element_index);

        double dbl_num_nodes = (double) element_iter->GetNumNodes();

        // Get the target area of the cell
        double cell_target_area = mMatureCellTargetArea; // GetTargetAreaOfCell(p_cell_population->GetCellUsingLocationIndex(element_index));

        // This is the modification scaling to allow target areas.
        double scaling = pow(cell_target_area,1.5)*4.0*GetWelikyOsterPerimeterParameter()*sin(M_PI/dbl_num_nodes)*sqrt(dbl_num_nodes*tan(M_PI/dbl_num_nodes));

        double deformation_coefficient = GetWelikyOsterAreaParameter()/element_area*scaling;

        /******** End of deformation force calculation *************/

        /******** Start of membrane force calculation ***********/

        // Compute the perimeter of the element
        double element_perimeter = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(element_index);

        double membrane_surface_tension_coefficient = GetWelikyOsterPerimeterParameter()*element_perimeter;

        /******** End of membrane force calculation **********/

        unsigned num_nodes = element_iter->GetNumNodes();
           for (unsigned node_local_index = 0; node_local_index < num_nodes; node_local_index++)
        {
            unsigned node_global_index = element_iter->GetNodeGlobalIndex(node_local_index);

            c_vector<double, DIM> current_node = element_iter->GetNodeLocation(node_local_index);
            c_vector<double, DIM> next_node = element_iter->GetNodeLocation((node_local_index + 1)%(element_iter->GetNumNodes()));
            c_vector<double, DIM> previous_node = element_iter->GetNodeLocation((node_local_index + element_iter->GetNumNodes() - 1)%(element_iter->GetNumNodes()));

            c_vector<double, DIM> clockwise_unit_vector = p_cell_population->rGetMesh().GetVectorFromAtoB(current_node, previous_node);
            clockwise_unit_vector /= norm_2(clockwise_unit_vector);
            c_vector<double, DIM> anti_clockwise_unit_vector = p_cell_population->rGetMesh().GetVectorFromAtoB(current_node, next_node);
            anti_clockwise_unit_vector /= norm_2(anti_clockwise_unit_vector);

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

            c_vector<double, DIM> deformation_contribution = deformation_coefficient * outward_normal;

            c_vector<double, DIM> membrane_surface_tension_contribution = membrane_surface_tension_coefficient * (clockwise_unit_vector + anti_clockwise_unit_vector);

            c_vector<double, DIM> force_on_node = deformation_contribution +
                                                  membrane_surface_tension_contribution;

            for (unsigned i=0; i<DIM; i++)
            {
                assert(!std::isnan(force_on_node(i)));
            }

            rForces[node_global_index] += force_on_node;
        }
    }
}

template<unsigned DIM>
double ModifiedWelikyOsterForce<DIM>::GetWelikyOsterAreaParameter()
{
    return mWelikyOsterAreaParameter;
}

template<unsigned DIM>
double ModifiedWelikyOsterForce<DIM>::GetWelikyOsterPerimeterParameter()
{
    return mWelikyOsterPerimeterParameter;
}

template<unsigned DIM>
void ModifiedWelikyOsterForce<DIM>::SetWelikyOsterAreaParameter(double welikyOsterAreaParameter)
{
    mWelikyOsterAreaParameter = welikyOsterAreaParameter;
}

template<unsigned DIM>
void ModifiedWelikyOsterForce<DIM>::SetWelikyOsterPerimeterParameter(double welikyOsterPerimeterParameter)
{
    mWelikyOsterPerimeterParameter = welikyOsterPerimeterParameter;
}

template<unsigned DIM>
void ModifiedWelikyOsterForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Output member variables
	*rParamsFile <<  "\t\t\t<ModifiedWelikyOsterAreaParameter>"<<  mWelikyOsterAreaParameter << "</WelikyOsterAreaParameter> \n";
	*rParamsFile <<  "\t\t\t<WelikyOsterPerimeterParameter>"<<  mWelikyOsterPerimeterParameter << "</WelikyOsterPerimeterParameter> \n";
	*rParamsFile <<  "\t\t\t<MatureCellTargetArea>"<<  mMatureCellTargetArea << "</MatureCellTargetArea> \n";

	// Call method on direct parent class
	AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class ModifiedWelikyOsterForce<1>;
template class ModifiedWelikyOsterForce<2>;
template class ModifiedWelikyOsterForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ModifiedWelikyOsterForce)
