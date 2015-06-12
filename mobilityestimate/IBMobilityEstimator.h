// Filename: IBMobilityEstimator.h
// Created on 29 Sept 2014 by Steven Delong
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBMobilityEstimator
#define included_IBMobilityEstimator

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LNodeSetData.h"
#include <ibamr/IBLagrangianForceStrategySet.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/AppInitializer.h>
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/app_namespaces.h>

// IBTK THIRD-PARTY INCLUDES
#include <muParser.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class IBMobilityEstimator is a subclass of IBHierarchyIntegrator
 * whose purpose is to estimate the mobility at given points.
 */
class IBMobilityEstimator
    : public IBHierarchyIntegrator
{
public:
    /*!
     * The constructor for class IBMobilityEstimator sets some default
     * values and reads in configuration information from input and restart
     * databases.
     *
     * \warning This simple example class does not support restarting.
     */
    IBMobilityEstimator(
        const std::string& object_name,
        Pointer<Database> input_db,
        Pointer<IBMethod> ib_method_ops,
        Pointer<INSHierarchyIntegrator> ins_hier_integrator,
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry);

    /*!
     * The destructor for class IBMobilityEstimator does
     * not do anything interesting.
     */
    ~IBMobilityEstimator();

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void
    preprocessIntegrateHierarchy(
        double current_time,
        double new_time,
        int num_cycles=1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void
    integrateHierarchy(
        double current_time,
        double new_time,
        int cycle_num=0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void
    postprocessIntegrateHierarchy(
        double current_time,
        double new_time,
        bool skip_synchronize_new_state_data,
        int num_cycles=1);

    /* !
     * Estimate mobility given the position of a particle and the setup of the
     * geometry from this object.
     */
    void
    estimateMobility(vector<double> position, double* mobility);
    
    /*!
     * Initialize any variables, communications algorithms, solvers, or other
     * data structures required by this time integrator object.
     */
    void
    initializeHierarchyIntegrator(
        Pointer<PatchHierarchy<NDIM> > hierarchy,
        Pointer<GriddingAlgorithm<NDIM> > gridding_alg);

    // function to normalize the eularian force by subtracting
    // totalForce/volume from each component
    void
    NormalizePointForces(void);

    // Overwrite the regrid function.
    void
    regridHierarchy();

    // Set up velocity BCs.
    void
    setVelocityBC(vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>*);

    // Spread forces.
    void
    spreadForce(const int f_data_idx,
                Pointer<LData> F_data,
                Pointer<LData> X_data,
                double time);

    // Interpolate velocity.
    void
    interpolateVelocity(const int u_data_idx,
                        Pointer<LData> U_data,
                        Pointer<LData> X_data,
                        double time);
    
    // Fill ghost cells. Used during force spreading.
    void
    fillGhostCells(int in, const double time);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBMobilityEstimator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBMobilityEstimator(
        const IBMobilityEstimator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMobilityEstimator&
    operator=(
        const IBMobilityEstimator& that);

    /* Fluid solver Variables
     * */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_u_star_var;
    
    /*
     * Pointers to Lagrangian data objects.
     */
    Pointer<LData> d_X_current_data, d_X_new_data, d_U_data,
        d_F_data;

    /*
     *  Mobility to be populated by "integrateHierarchy."
     */
    double* d_mobility;

    /* User-provided constants specified in the input file.
     */
    std::map<std::string,double> d_constants;

    int d_normalize_force_flag;
    
    // array of nonbonded interaction radii for different species.
    std::vector<std::vector<double> > d_interaction_radius;

    // array of nonbonded interaction radii for different species.
    std::vector<std::vector<std::vector<double> > > d_nonbonded_parameters;
    
    // output file for mobility estimate
    ofstream d_out_mob_file;

    // Vector of BC coeffs.  HACK: Make sure this is needed.
    vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> *d_u_bc_coefs;

};

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBMobilityEstimator
