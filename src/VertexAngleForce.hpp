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

#ifndef VERTEXANGLEFORCE_HPP_
#define VERTEXANGLEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"

/**
 * A force class for use in vertex dynamics simulations. Provides a force which aims to prevent
 * convex elements by adding a restraining force to try to make all elements regular.
 */
template<unsigned DIM>
class VertexAngleForce : public AbstractForce<DIM>
{
private:

    /** Coefficient of restraining force. Has units of kg s^-2 (cell size at equilibrium rest length)^-1. */
    double mAngleRestrainingParameter;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mAngleRestrainingParameter;
    }

public:

    /**
     * Constructor.
     */
    VertexAngleForce();

    /**
     * Destructor.
     */
    ~VertexAngleForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * @return mAngleRestrainingParameter.
     */
    double GetAngleRestrainingParameter();

    /**
     * Set mAngleRestrainingParameter.
     * 
     * @param angleRestrainingParameter the new value of mAngleRestrainingParameter
     */
    void SetAngleRestrainingParameter(double angleRestrainingParameter);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VertexAngleForce)

#endif /*VERTEXANGLEFORCE_HPP_*/
