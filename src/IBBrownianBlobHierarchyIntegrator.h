// Filename: IBBrownianBlobHierarchyIntegrator.h
// Created on 04 Apr 2013 by Boyce Griffith
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

#ifndef included_IBBrownianBlobHierarchyIntegrator
#define included_IBBrownianBlobHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include "NonbondedForceEvaluator.h"
#include "WallForceEvaluator.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LNodeSetData.h"
#include "ibamr/IBSpringForceSpec.h"
#include <ibamr/IBLagrangianForceStrategySet.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/AppInitializer.h>
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/app_namespaces.h>

#include "tbox/Pointer.h"
#include "tbox/Database.h"
#include "SAMRAIVectorReal.h"
#include <RobinBcCoefStrategy.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LData.h"
#include "string"

// IBTK THIRD-PARTY INCLUDES
#include <muParser.h>


/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class IBBrownianBlobHierarchyIntegrator is an implementation of a simple
 * first-order accurate, semi-implicit version of the immersed boundary method.
 */
class IBBrownianBlobHierarchyIntegrator
    : public IBHierarchyIntegrator
{
public:
    /*!
     * The constructor for class IBBrownianBlobHierarchyIntegrator sets some default
     * values and reads in configuration information from input and restart
     * databases.
     *
     * \warning This simple example class does not support restarting.
     */
    IBBrownianBlobHierarchyIntegrator(
        const std::string& object_name,
        Pointer<Database> input_db,
        Pointer<IBMethod> ib_method_ops,
        Pointer<INSHierarchyIntegrator> ins_hier_integrator,
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry);

    /*!
     * The destructor for class IBBrownianBlobHierarchyIntegrator does
     * not do anything interesting.
     */
    ~IBBrownianBlobHierarchyIntegrator();

    /*!
     * Generate gaussian N(0,1) random numbers at each node for each component
     * */

    void
    genrandn(Pointer<LData> Noise_data);
    
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

    /*!
     * Initialize any variables, communications algorithms, solvers, or other
     * data structures required by this time integrator object.
     */
    void
    initializeHierarchyIntegrator(
        Pointer<PatchHierarchy<NDIM> > hierarchy,
        Pointer<GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Function to do RFD spreading to approximate thermal drift.  This adds
     * to the eularian force variable corresponding to d_f_idx.
     */
    void
    spreadRFDForces(double time);

    /*  Spread Lagrangian force to the grid, not spreading forces from
     *  anchor points */
    void
    spreadForce(int d_f_idx,
                Pointer<LData> F_data,
                Pointer<LData> X_data,
                double time);

    /*  Interpolate velocities from the grid to the particles.  This won't
        interpolate at anchor points. */
    void
    interpolateVelocity(int d_u_idx,
                        Pointer<LData> U_data,
                        Pointer<LData> X_data,
                        double time);
    
    // register the ib_force_fcn directly with the hierarchy integrator 
    void
    registerIBLagrangianForceFunction(
        Pointer<IBLagrangianForceStrategySet> ib_force_fcn);
    
    /* function to set xi*/
    void
    SetXi(double inputXi);

    // function to set KT
    void
    SetkT(double inputkT);

    // function to set delta for random finite difference 
    void
    SetRfdelta(double inputrf, Pointer<CartesianGridGeometry<NDIM> > grid_geometry);

    // function to set dynamics, OVERDAMPED or RESOLVED
    void
    SetDynamics(string dynamics);

    // function to set the BG Flow function.
    void
    SetBGFlowFunction(const std::string& object_name,
                          Pointer<Database> input_db);

   // used in default case to set BG Flow to 0 everywhere.
    void
    SetBGFlowFunctionNull(const std::string& object_name);

    // set a lagrangian function
    void
    SetLagrangianFunction(const std::string& /*object_name*/,
                          Pointer<Database> input_db,
                          std::vector<mu::Parser>& parser);


    // set lagrangian function to void
    void
    SetLagrangianFunctionNull(const std::string& /*object_name*/,
                              std::vector<mu::Parser>& parser);
    
    
    // add the background flow
    void
    adBGFlow(Pointer<LData> d_U_data, Pointer<LData> d_X_data, double time);

    // helper functino to print out total force on Eularian Grid
    // in each direction.
    void
    CheckEularianForceTotals(void);

    // function to print total momentum and stochastic flux, as well as <v^2>
    void
    PrintDiagnostics(double current_time);

    //function to add forces to each point particle
    void
    addPointForces(Pointer<LData> d_F_data, Pointer<LData> d_X_current_data, double eval_time);

    
    // function to normalize the eularian force by subtracting
    // totalForce/volume from each component
    void
    NormalizePointForces(Pointer<LData> F_data);

    // print eularian side data to terminal.
    void
    printEularianData(int data_idx);
    
    // function to check that particles don't move more than one cellwidth
    // per timestep.  Causes an error if they do.
    void
    checkParticleMoveDistance(Pointer<LData> X_moved_data);

    // overwrite regridding criterion.  Regrid when any particle has moved more than
    // regrid_alpha*dx eularian distance since the last regrid
    
    bool
    atRegridPointSpecialized() const;

    // overwrite the regrid function to also update the d_X_last_regrid lagrangian
    // data.
    void
    regridHierarchy();

    // set the shear rate and origin for shear flow.  0.0 means no shear.
    // By default this is d/d_NDIM (v_x)
    // shear just adds to particle motion
    void
    setShear(double shear_rate_in, double shear_origin_in);

    
    // function to add a constant shear to the velocity
    void
    addShear(Pointer<LData> d_U_data,
             Pointer<LData> d_X_data,
             int axis,
             int direction);

    // function to add background flow velocity to particle.
    void
    addBGFlow(Pointer<LData> U_data, Pointer<LData> X_data, double eval_time);
    

    // function to read in parameters for nonbonded interactions
    void
    registerNonBondedInteractionParameters(Pointer<Database> input_db);
    

    // function to remove all nonbonded springs from the springspec objects
    // belonging to the Lagrangian data
    void
    removeNonBondedSprings(void);

    
    // function to remove all nonbonded springs from the springspec objects
    // belonging to the Lagrangian data
    void
    addNonBondedSprings(void);

    // function to remove the springs corresponding to indices in remove_idxs
    // from the springspec force_spec
    void
    removeSpringFromForceSpec(IBSpringForceSpec** force_spec_ptr,std::vector<int>& remove_idxs);
    
    // function to set nonbonded force index
    void
    setNonBondedForceFunctionIdx(int force_fcn_idx);

    // Function to initialize set of anchor indices based on node data.
    void
    initializeAnchorIndices(void);
    
    // Function to zero lagrangian velocities and forces for anchor points.
    void
    resetAnchorPointValues(Pointer<LData> U_data);


    // function to calculate nonbonded forces for periodic boundaries.
    void
    addNonbondedForces(Pointer<LData> X_data, Pointer<LData> F_data);
    
    
    // rebuild the lists of nodes near each wall
    void
    rebuildWalls();

    // add a wall.
    void
    addWall(Pointer<Database> wall_db,
            Pointer<CartesianGridGeometry<NDIM> > grid_geometry,
            double wall_ghost_dist);

    
    // add walls from input database.
    void
    registerWalls(Pointer<AppInitializer> app_initializer,
                 Pointer<CartesianGridGeometry<NDIM> > grid_geometry);
    
    // multiply quaternions to rotate particles.
    void
    multiplyQuaternion(Pointer<LData> Q_data, Pointer<LData> W_data);

    // fill ghost cells
    void
      fillGhostCells(int in, const double time);

    // set velocity boundary conditions
    void
    setVelocityBC(vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> *u_bc_coefs);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBBrownianBlobHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBBrownianBlobHierarchyIntegrator(
        const IBBrownianBlobHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBBrownianBlobHierarchyIntegrator&
    operator=(
        const IBBrownianBlobHierarchyIntegrator& that);


    /* Fluid solver Variables
     * */

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_u_star_var;

    
    int d_u_star_current_idx, d_u_star_new_idx, d_u_star_scratch_idx,
        d_nonbonded_force_idx;
    /*
     * Pointers to Lagrangian data objects.
     */
    Pointer<LData> d_X_current_data, d_X_new_data, d_U_data,
        d_F_data, d_X_star_data, d_F_new_data, d_U_star_data,
        d_Noise_data, d_Noise_2_data, d_RFD_Noise_data, d_RFD_loc_1_data,
        d_RFD_loc_2_data, d_X_last_regrid, d_Q_data, d_Q_star_data, d_W_data;

    // physical and computational parameters
    double d_kT, d_xi, d_rfdelta, d_physical_rho;

    // flag to set point/wall forces on particles or not
    int d_point_force_flag, d_wall_force_flag, d_bg_flow_flag, d_nonbonded_flag,
        d_nonbonded_by_cell_flag, d_normalize_force_flag;

    // variable to include RFD term or not
    bool d_use_rfd;

    // dynamics and time_stepping_type.
    string d_dynamics, d_time_stepping_type;

    // control how far a particle moves between regriddings,
    // regridding will happen when any particle has moved regrid_alpha*dx since the last
    // regrid
    double d_regrid_alpha;

    /* User-provided constants specified in the input file.
     */
    std::map<std::string,double> d_constants;

    /*!
     * The strings providing the data-setting functions which are evaluated by the
     * mu::Parser objects for pointForces.
     */
    std::vector<std::string> d_function_strings;

    /*!
     * The mu::Parser objects which evaluate the data-setting functions for pointforces.
     */
    std::vector<mu::Parser> d_parsers;


    /* Set of ints indicating which points are anchor points.  These points will not
     * interact with the fluid in any way.
     */
    std::set<int> d_anchor_point_local_idxs;

    /*!
     * The strings providing the data-setting functions which are evaluated by the
     * mu::Parser objects for background flow.
     */
    std::vector<std::string> d_bgflow_strings;

    /*!
     * The mu::Parser objects which evaluate the data-setting functions for background flow.
     */
    std::vector<mu::Parser> d_flow_parsers;

    /*!
     * The mu::Parser objects which evaluate the data-setting functions for nonbonded interactions.
     */

    // Force evaluator object for nonbonded interactions.
    NonbondedForceEvaluator* d_force_evaluator;

    // Force evaluator object for wall forces.
    WallForceEvaluator* d_wall_force_evaluator;
    
    
    // position and time for Lagrangian Functions
    double d_parser_time;
    Point d_parser_posn;
    
    // integers to determine when we print diagnostics.
    int d_diagnostic_counter;
    int d_diagnostic_interval;

    // lagrangian force generator
    SAMRAI::tbox::Pointer<IBLagrangianForceStrategySet> d_ib_force_fcn;

    // number of nonbonded species
    int d_num_species;
    
    // array of nonbonded interaction radii for different species.
    std::vector<std::vector<double> > d_interaction_radius;

    // array of nonbonded interaction radii for different species.
    std::vector<std::vector<std::vector<double> > > d_nonbonded_parameters;
    
    // output file for mobility estimate
    ofstream out_mob_file;

    // mobility component, changes each timestep.
    int mob_comp;

    // velocity boundary coefficients
    //vector<RobinBcCoefStrategy<NDIM>*> *d_u_bc_coefs;
    vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> *d_u_bc_coefs;

};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBBrownianBlobHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBBrownianBlobHierarchyIntegrator
