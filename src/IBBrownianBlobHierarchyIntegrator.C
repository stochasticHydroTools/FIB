// Filename: IBBrownianBlobHierarchyIntegrator.C
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


/////////////////////////////// INCLUDES /////////////////////////////////////

#include "IBBrownianBlobHierarchyIntegrator.h"
#include "NonbondedForceEvaluator.h"
#include "WallForceEvaluator.h"
#include "ibamr/IBAnchorPointSpec.h"
#include <ibamr/RNG.h>
#include <cmath>
#include <math.h> // ceil

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/ibamr_utilities.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBBrownianBlobHierarchyIntegrator::IBBrownianBlobHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<IBMethod> ib_method_ops,
    Pointer<INSHierarchyIntegrator> ins_hier_integrator,
    Pointer<CartesianGridGeometry<NDIM> > /*grid_geometry */)
    : IBHierarchyIntegrator(object_name, input_db, ib_method_ops, ins_hier_integrator, /*register_for_restart*/ false),
      d_xi(0.0),
      d_point_force_flag(0),
      d_wall_force_flag(0),
      d_bg_flow_flag(0),
      d_nonbonded_flag(0),
      d_normalize_force_flag(0),
      d_diagnostic_counter(0)
{
    
    // set time stepping type
    if(input_db->keyExists("time_stepping_type"))
    {
        d_time_stepping_type = input_db->getString("time_stepping_type");
    }
    else
    {
        pout << "No Time Stepping Type Specified.  Defaulting to Forward Euler." << std::endl;
        d_time_stepping_type = "FORWARD_EULER";
    }

    // set RFD flag
    if (input_db->keyExists("use_rfd"))
    {
        d_use_rfd = input_db->getBool("use_rfd");
    }
    else
    {
        d_use_rfd = true;
    }
    // Set up fraction of a cell width to use as a regrid distance.
    if(input_db->keyExists("regrid_alpha"))
    {
        d_regrid_alpha = input_db->getDouble("regrid_alpha");
    }
    else
    {
        pout << "Regrid alpha not in database.  Defaulting to 0.5 cell width "
            " regrid distance" << std::endl;
        d_regrid_alpha = 0.5;
    }
    
    // Set up diagnostic_interval, after how many timesteps do we print checks (momentum,
    // etc).
    // 0 indicates not to print checks.
    if(input_db->keyExists("diagnostic_interval"))
    {
        d_diagnostic_interval = input_db->getInteger("diagnostic_interval");
    }
    else
    {
        d_diagnostic_interval = 0;
    }
     // set up d_physical_rho, used only in momentum calculations for diagnostics (above).
    if(input_db->keyExists("physical_rho"))
    {
        d_physical_rho = input_db->getDouble("physical_rho");
    }
    else
    {
        d_physical_rho = 0.0;
        pout << "physical_rho missing from IBHierarchyIntegrator, setting it to 0.0.  Momentum calculations will be 0." << std::endl;
    }

    // do we normalize forces or not?  This should be true with a fully periodic domain.
    if (input_db->keyExists("normalize_force"))
    {
        d_normalize_force_flag = input_db->getInteger("normalize_force");
    }
    else
    {
        TBOX_ERROR("Must indicate if we should normalize forces or not with normalize_force (1 means normalize)");
    }

    // set lagrangian force functions for pointforces, nonbonded forces,
    // and background flow.

    // register IB point forces
    if (input_db->keyExists("PointForceFunction"))
    {
        cout << "found point force function." << std::endl;
        SetLagrangianFunction("F_Point_fcn",input_db->getDatabase("PointForceFunction"),d_parsers);
        d_point_force_flag = 1;
    }
    else
    { 
        SetLagrangianFunctionNull("F_Point_fcn",d_parsers);
        d_point_force_flag = 1;
    }


    // set nonbonded force switch.
    d_nonbonded_by_cell_flag = input_db->getIntegerWithDefault("nonbonded_by_cell",1);
    // TODO(delong): Check to make sure switch will work for springs v. cell search.

    // register background flow functions
    if (input_db->keyExists("BGFlowFunction"))
    {
        SetLagrangianFunction("BGFlow_fcn",input_db->getDatabase("BGFlowFunction"),d_flow_parsers);
        d_bg_flow_flag = 1;
    }
    else
    { 
        SetLagrangianFunctionNull("BGFlow_fcn",d_flow_parsers);
        d_bg_flow_flag = 0;
    }
    return;
}// IBBrownianBlobHierarchyIntegrator


IBBrownianBlobHierarchyIntegrator::~IBBrownianBlobHierarchyIntegrator()
{
    // Intentionally blank.
    return;
}// ~IBBrownianBlobHierarchyIntegrator

 void
IBBrownianBlobHierarchyIntegrator::genrandn(
    Pointer<LData> Noise_data)
{
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
 
    Vec Noise_data_vec = Noise_data->getVec();
    PetscScalar *noise;
    VecGetArray(Noise_data_vec,&noise);
     const Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_level_num);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
         cit != local_nodes.end(); ++cit)
    {
        
        LNode* const node_idx = *cit;
        int local_idx = node_idx->getLocalPETScIndex();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            RNG::genrandn(&(noise[local_idx*NDIM + d]));
        }
    }
    VecRestoreArray(Noise_data_vec,&noise);
    return;
}// genrandn

void
IBBrownianBlobHierarchyIntegrator::preprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const int num_cycles)
{
    IBHierarchyIntegrator::preprocessIntegrateHierarchy(current_time,
                                                        new_time,
                                                        num_cycles);
    const int coarsest_level_num = 0;
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();

    // Allocate Eulerian scratch and new data.
    for (int level_num = coarsest_level_num;
         level_num <= finest_level_num;
         ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level =
            d_hierarchy->getPatchLevel(level_num);
        level->allocatePatchData(d_u_idx, current_time);
        level->allocatePatchData(d_f_idx, current_time);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data    ,     new_time);
    }
     // Initialize the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    d_ins_hier_integrator->preprocessIntegrateHierarchy(current_time,
                                                        new_time,
                                                        ins_num_cycles);
    // Initialize IB data.
    d_ib_method_ops->preprocessIntegrateData(current_time, new_time,1);
    // NOTE: We assume here that all IB data are assigned to the finest level of
    // the AMR patch hierarchy.
    d_X_current_data  = l_data_manager->getLData(LDataManager::POSN_DATA_NAME,
                                                 finest_level_num);
    d_X_star_data     = l_data_manager->createLData("X_star", finest_level_num,
                                                    NDIM);
    d_X_new_data      = l_data_manager->createLData("X_new", finest_level_num,
                                                    NDIM);
    d_U_data          = l_data_manager->getLData(LDataManager::VEL_DATA_NAME,
                                                 finest_level_num);
    d_U_star_data     = l_data_manager->createLData("U_star", finest_level_num,
                                                    NDIM);
    d_Q_data          = l_data_manager->createLData("Q", finest_level_num,
                                                    ((NDIM-1)*NDIM)/2 + 1);
    d_Q_star_data     = l_data_manager->createLData("Q_star", finest_level_num,
                                                    ((NDIM-1)*NDIM)/2 + 1);
    d_W_data          = l_data_manager->createLData("W", finest_level_num,
                                                    ((NDIM-1)*NDIM)/2 + 1);
    d_F_data          = l_data_manager->createLData("F", finest_level_num,
                                                    NDIM);
    d_F_new_data      = l_data_manager->createLData("F_new", finest_level_num,
                                                    NDIM);
    d_Noise_data      = l_data_manager->createLData("Noise", finest_level_num,
                                                    NDIM);
    d_Noise_2_data    = l_data_manager->createLData("Noise_2", finest_level_num,
                                                    NDIM);
    d_RFD_Noise_data  = l_data_manager->createLData("RFD_Noise", finest_level_num,
                                                    NDIM);
    d_RFD_loc_1_data  = l_data_manager->createLData("RFD_loc_1", finest_level_num,
                                                    NDIM);
    d_RFD_loc_2_data  = l_data_manager->createLData("RFD_loc_2", finest_level_num,
                                                    NDIM);
    return;
}// preprocessIntegrateHierarchy

void
IBBrownianBlobHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time,
    const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);

    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    PetscErrorCode ierr;
    const double dt = new_time-current_time;
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();

    // WE HAVE NOT IMPLEMENTED ANYTHING BUT OVERDAMPED DYNAMICS YET
    if(d_dynamics == "RESOLVED")
    {
        TBOX_ERROR("RESOLVED dynamcis is not yet implemented. Please use OVERDAMPED dynamics.");
        // Here we implement a simple time integration scheme:
        //
        // (1) Compute Lagrangian forces and spread those forces to the grid: f =
        //     S[X^{n}] F^{n}.
        //
        // (2) Solve the fluid equations using the fluid solver registered with this
        //     class using the forcing computed in (1).
        //
        // (3) Interpolate the Eulerian velocity and update the positions of the
        //     Lagrangian points: X^{n+1} = X^{n} + dt*S^{*}[X^{n}] u^{n+1}.
        //
        // NOTE: We assume here that all IB data are assigned to the finest level of
        // the AMR patch hierarchy.
        
        
        ///////////////////////////////////////////////////////////////////////////
        // (1) Compute Lagrangian forces and spread those forces to the grid: f =
        //     S[X^{n}] F^{n}..
        
        // For simplicity, set the Lagrangian force density to equal zero.
        ierr = VecZeroEntries(d_F_data->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecZeroEntries(d_U_data->getVec());  IBTK_CHKERRQ(ierr);
     
        // Spread the forces to the grid.  We use the "current" Lagrangian position
        // data to define the locations from where the forces are spread.
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_hier_velocity_data_ops->setToScalar(d_u_idx, 0.0);

        spreadForce(d_f_idx, d_F_data, d_X_current_data, current_time);
	///spreadForce(d_f_big_gc_idx, d_F_data, d_X_current_data, current_time);

        // NOTE: Any additional Eulerian forcing should be computed here and added
        // to the data associated with d_f_idx.
        
        ///////////////////////////////////////////////////////////////////////////
        // (2) Solve the fluid equations using the fluid solver registered with this
        // class using the forcing computed in (1).
        const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
        for (int ins_cycle_num = 0; ins_cycle_num < ins_num_cycles; ++ins_cycle_num)
        {   
	    d_ins_hier_integrator->integrateHierarchy(current_time, new_time, ins_cycle_num);
        }
        
        ///////////////////////////////////////////////////////////////////////////
        // (3) Interpolate the Eulerian velocity and update the positions of the
        // Lagrangian points: X^{n+1} = X^{n} + dt*S^{*}[X^{n}] u^{n+1}.
        //
        // NOTE: We use the "new" velocity data (defined at time level n+1) to
        // determine the velocity of the Lagrangian nodes.  We could also use the
        // "current" data (defined at time level n) or some other velocity field
        // here.  We use the "current" Lagrangian position data to define the
        // locations to where the velocities are interpolated.
     
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int u_new_idx =
            var_db->mapVariableAndContextToIndex(
                d_ins_hier_integrator->getVelocityVariable(),
                d_ins_hier_integrator->getNewContext());
        d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
        interpolateVelocity(d_u_idx,
                            d_U_data,
                            d_X_current_data,
                            current_time);
        // generate noise
        genrandn(d_Noise_data);
        
        // basic FE scheme
        ierr = VecWAXPY(d_X_new_data->getVec(),
                        dt,
                        d_U_data->getVec(),
                        d_X_current_data->getVec());
        IBTK_CHKERRQ(ierr);
    } // if (dynamics == "RESOLVED")
    
    else if (d_dynamics == "OVERDAMPED")
    {
        // Delong: This is the first order method implemented in the
        // thermaldriftscheme.pdf and the second order (deterministic) method
        // implemented in thermal drift scheme
        // Test for rho, make sure it is 0 for stokes solve.
        const StokesSpecifications* d_problem_coefs =
            d_ins_hier_integrator->getStokesSpecifications();
        if (d_problem_coefs->getRho() != 0.0)
        {
            // For now enforce RHO = 0, use a separate variable for density calculations
            TBOX_ERROR("RHO must be equal to zero for steady solves.  Use Physical Rho"
                       "for the non-zero physical density");
        }
        
        if(d_time_stepping_type == "FORWARD_EULER" ||
           d_time_stepping_type == "TRAPEZOIDAL_RULE" ||
           d_time_stepping_type == "MIDPOINT_RULE") // allowed first order schemes
        {
            /***********
             *  FIRST ORDER METHOD
             * **********/
            // Note: Zeroing all of this may not be necessary
            // Zero lagrangian variables that will be added to later.
            ierr = VecZeroEntries(d_F_data->getVec());  IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(d_U_data->getVec());  IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(d_X_star_data->getVec()); IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(d_X_new_data->getVec());  IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(d_F_new_data->getVec());  IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(d_RFD_Noise_data->getVec());  IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(d_RFD_loc_1_data->getVec());  IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(d_RFD_loc_2_data->getVec());  IBTK_CHKERRQ(ierr);
            
            // Zero Eularian variables.
            d_hier_velocity_data_ops->setToScalar(d_f_idx,0.0);
            d_hier_velocity_data_ops->setToScalar(d_u_idx, 0.0);
            
            // Update position ghost nodes before calculating forces.
            d_X_current_data->beginGhostUpdate();
            d_X_current_data->endGhostUpdate();
            

            // Add RFD term for drift reproducing schemes.
            // d_use_rfd defaults to true.
            if (d_use_rfd)
            {
                spreadRFDForces(current_time);
            }

            // ADD LAGRANGIAN FORCES HERE
            if(d_point_force_flag)
            { // if we have point forces, then apply them.
                addPointForces(d_F_data, d_X_current_data, current_time);
            }
            // add spring forces. Note, these forces MUST NOT change the total force
            d_ib_force_fcn->computeLagrangianForce(d_F_data,
                                                   d_X_current_data,
                                                   d_U_data,
                                                   d_hierarchy,
                                                   finest_level_num,
                                                   current_time,
                                                   l_data_manager);
            // Spread the forces to the grid.  We use the "current" Lagrangian position
            // data to define the locations from where the forces are spread.
            // sync computed F before spreading to eularian grid
            d_F_data->beginGhostUpdate();
            d_F_data->endGhostUpdate();
            spreadForce(d_f_idx, d_F_data, d_X_current_data, current_time);

            // NOTE: Any additional Eulerian forcing should be computed here and added
            // to the data associated with d_f_idx.
            
            ///////////////////////////////////////////////////////////////////////////
            // (2) Solve the fluid equations using the fluid solver registered with this
            // class using the forcing computed in (1).
            const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();

            if(ins_num_cycles != 1)
            {
                TBOX_ERROR("For OVERDAMPED and EULER options, number of cycles must "
                           "equal 1. \n");
            }
            for (int ins_cycle_num = 0; ins_cycle_num < ins_num_cycles; ++ins_cycle_num)
            {
                d_ins_hier_integrator->integrateHierarchy(current_time,
                                                          new_time,
                                                          ins_cycle_num);
            }
            
            // Copy over eularian velocity to u_new.
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int u_new_idx = var_db->mapVariableAndContextToIndex(
                d_ins_hier_integrator->getVelocityVariable(),
                d_ins_hier_integrator->getNewContext());
            d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx, /*interior_only*/ true);

            // Interpolate eularian velocity back to particles
            interpolateVelocity(d_u_idx,
                                d_U_data,
                                d_X_current_data,
                                current_time);

            // Generate noise for particle motion.
            genrandn(d_Noise_data);

            // If we specify a background flow, add it to the particle velocities here.
            if(d_bg_flow_flag)
            {
                addBGFlow(d_U_data, d_X_current_data, current_time);
            }
            
            // sync lagrangian velocity
            d_U_data->beginGhostUpdate();
            d_U_data->endGhostUpdate();
            
            if(d_time_stepping_type == "FORWARD_EULER")
            {
                // just one forward step to X_new
                // trap/euler.  F holds total velocity temporarily.
                // NOTE: This scheme WILL NOT captpure the correct thermal
                // drift.
                ierr = VecAXPBYPCZ(d_F_data->getVec(),
                                   1.0,
                                   sqrt(2.0*d_kT*d_xi/dt),
                                   d_xi,d_U_data->getVec(),
                                   d_Noise_data->getVec());
                IBTK_CHKERRQ(ierr);
                // euler
                ierr = VecWAXPY(d_X_new_data->getVec(),
                                dt,d_F_data->getVec(),
                                d_X_current_data->getVec());
                IBTK_CHKERRQ(ierr); 
            }

            // Interpolate at the new position, TRAPEZOIDAL
            if (d_time_stepping_type == "TRAPEZOIDAL_RULE")
            {
                // predictor: d_X_star = d_X_current + dt*xi*F +
                //                       dt*d_U_data + sqrt(2kTxi dt)*Noise
                // F will store the result of xi*f + U + sqrt(2kT xi /dt)*noise

                // trap/euler
                ierr = VecAXPBYPCZ(d_F_data->getVec(),
                                   1.0,
                                   sqrt(2.0*d_kT*d_xi*dt),
                                   d_xi,d_U_data->getVec(),
                                   d_Noise_data->getVec());
                IBTK_CHKERRQ(ierr); 
                // trap
                ierr = VecWAXPY(d_X_star_data->getVec(),
                                dt,d_F_data->getVec(),
                                d_X_current_data->getVec());
                IBTK_CHKERRQ(ierr);
            
                d_X_star_data->beginGhostUpdate();
                d_X_star_data->endGhostUpdate();

                // check that we haven't moved past 2 regrid intervals
                // with the predictor step.
                checkParticleMoveDistance(d_X_star_data);
                
                // then we evaluate forces at X_star.
                // interp u at X_star to get U_star
                // then corrector: d_X_new = d_X_current + 0.5*dt*d_xi*(F + F_star) +
                //                           0.5*dt*(d_U_data + d_U_star) +
                //                           sqrt(2kT d_xi)*Noise
                
                // predictor corrector for J thermal drift piece.
                interpolateVelocity(d_u_idx,
                                    d_U_star_data,
                                    d_X_star_data,
                                    current_time);
                if(abs(d_xi)  > 0.0)
                {
                    //evaluate forces at X^* if d_xi != 0 so we need to
                    if(d_point_force_flag)
                    { // if we have pointforces, then apply them
                        // no need to normalize, no more fluid solves this step
                        addPointForces(d_F_new_data,d_X_star_data,new_time);
                    } // if(d_point_force_flag)

                    // spring forces at X^*
                    d_ib_force_fcn->computeLagrangianForce(d_F_new_data,
                                                           d_X_star_data,
                                                           d_U_data,
                                                           d_hierarchy,
                                                           finest_level_num,
                                                           new_time,
                                                           l_data_manager);
                    // sync forces
                    d_F_new_data->beginGhostUpdate();
                    d_F_new_data->endGhostUpdate();
                }
                // Add background flow to U_star if it exists.
                if(d_bg_flow_flag)
                {
                    addBGFlow(d_U_star_data,d_X_star_data,current_time + dt);
                }
                
                // add trapezoidal U piece from x^n
                ierr = VecAXPBYPCZ(d_X_new_data->getVec(),
                                   0.5*dt,
                                   1.0,
                                   1.0,
                                   d_U_data->getVec(),
                                   d_X_current_data->getVec());
                IBTK_CHKERRQ(ierr);
                // add trapezoidal d_xi F piece and noise from x^n
                ierr = VecAXPBYPCZ(d_X_new_data->getVec(),
                                   0.5*dt*d_xi,
                                   sqrt(0.5*d_kT*d_xi/dt),
                                   1.0,
                                   d_F_data->getVec(),
                                   d_Noise_data->getVec());
                IBTK_CHKERRQ(ierr);

                // add trapezoidal U piece from x^*
                ierr = VecAXPY(d_X_new_data->getVec(),
                               0.5*dt,
                               d_U_star_data->getVec());
                IBTK_CHKERRQ(ierr);
            
                // add trapezoidal d_xi F piece and noise from x^*
                ierr = VecAXPBYPCZ(d_X_new_data->getVec(),
                                   0.5*dt*d_xi,
                                   sqrt(0.5*d_kT*d_xi/dt),
                                   1.0,
                                   d_F_new_data->getVec(),
                                   d_Noise_data->getVec());
                IBTK_CHKERRQ(ierr);
            } // end trapezoidal corrector
            
            // interpolate at new position, MIDPOINT CORRECTOR
            if (d_time_stepping_type == "MIDPOINT_RULE")
            {
                // predictor: d_X_star = d_X_current + 0.5*dt*d_xi*F + 0.5*dt*d_U_data +
                //                      sqrt(d_kT d_xi dt)*Noise
                // F will store the result of d_xi*f + U + sqrt(4d_kT d_xi /dt)*noise
                ierr = VecAXPBYPCZ(d_F_data->getVec(),
                                   1.0,
                                   sqrt(4.0*d_kT*d_xi/dt),
                                   d_xi,
                                   d_U_data->getVec(),
                                   d_Noise_data->getVec());
                IBTK_CHKERRQ(ierr); // mid

                ierr = VecWAXPY(d_X_star_data->getVec(),
                                0.5*dt,
                                d_F_data->getVec(),
                                d_X_current_data->getVec());
                IBTK_CHKERRQ(ierr); // mid
                
                // sync at X_star
                d_X_star_data->beginGhostUpdate();
                d_X_star_data->endGhostUpdate();
                
                // check that we haven't moved past 2 regrid intervals with the predictor step
                checkParticleMoveDistance(d_X_star_data);
                
                // then we evaluate forces at X_star.
                // interp u at X_star to get U_star
                // then corrector: d_X_new = d_X_current + dt*d_xi F_star
                // + dt*d_U_star) + sqrt(d_kT d_xi)*(Noise_1 + Noise_2)
                // generate the second noise
                genrandn(d_Noise_2_data);
                interpolateVelocity(d_u_idx,
                                    d_U_star_data,
                                    d_X_star_data,
                                    current_time);


                // Apply forces at X^* if xi > 0
                if(abs(d_xi)  > 1.0e-9)
                {
                    if(d_point_force_flag)
                    { // if we have pointforces, then apply them
                        // no need to normalize, no more fluid solves this step
                        addPointForces(d_F_new_data,d_X_star_data,new_time);
                    } // if(d_point_force_flag)
                    // spring forces at X^*
                    d_ib_force_fcn->computeLagrangianForce(d_F_new_data,
                                                           d_X_star_data,
                                                           d_U_data,
                                                           d_hierarchy,
                                                           finest_level_num,
                                                           new_time,
                                                           l_data_manager);
                    // sync forces
                    d_F_new_data->beginGhostUpdate();
                    d_F_new_data->endGhostUpdate();
                }
                if(d_bg_flow_flag)
                {
                    addBGFlow(d_U_star_data,d_X_star_data,current_time + dt/2.0);
                }
                // sync at U_star
                d_U_star_data->beginGhostUpdate();
                d_U_star_data->endGhostUpdate();

                // add midpoint U piece from x^*
                ierr = VecWAXPY(d_X_new_data->getVec(),
                                dt,
                                d_U_star_data->getVec(),
                                d_X_current_data->getVec());
                IBTK_CHKERRQ(ierr);
                // add midpoint d_xi F piece and noise from x^*
                ierr = VecAXPBYPCZ(d_X_new_data->getVec(),
                                   dt*d_xi,
                                   sqrt(dt*1.0*d_kT*d_xi),
                                   1.0,
                                   d_F_new_data->getVec(),
                                   d_Noise_data->getVec());
                IBTK_CHKERRQ(ierr);
                
                // add second noise piece for midpoint
                ierr = VecAXPY(d_X_new_data->getVec(),
                               sqrt(dt*1.0*d_kT*d_xi),
                               d_Noise_2_data->getVec());
                IBTK_CHKERRQ(ierr);
            } // midpoint corrector

            // if it's time, print diagnostics
            PrintDiagnostics(current_time);
        } // end if d_time_stepping_type == "FORWARD_EULER" ,
          // "MIPOINT_RULE", or "TRAPEZOIDAL_RULE"
        else if (d_time_stepping_type == "SECOND_TRAPEZOIDAL_RULE" ||
                 d_time_stepping_type == "SECOND_MIDPOINT_RULE")
            // Allowed "second order" schemes.
        {
            /***************
             * SECOND ORDER (FOR ADDITIVE NOISE) METHOD
             * *************/
            // For simplicity, set the Lagrangian force density to equal zero.
            ierr = VecZeroEntries(d_F_data->getVec());  IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(d_U_data->getVec());  IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(d_F_new_data->getVec());  IBTK_CHKERRQ(ierr);
            ierr = VecZeroEntries(d_X_new_data->getVec());  IBTK_CHKERRQ(ierr);
                                                                                 
            // Spread the forces to the grid.  We use the "current" Lagrangian position
            // data to define the locations from where the forces are spread.
            d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
            d_hier_velocity_data_ops->setToScalar(d_u_idx,0.0);
            
            // sync X
            d_X_current_data->beginGhostUpdate();
            d_X_current_data->endGhostUpdate();

            if(d_point_force_flag)
            { // if we have pointforces, then apply them and
                // normalize the eularian forces.
                addPointForces(d_F_data,d_X_current_data,current_time);
            } // if(d_point_force_flag)

            // apply RFD term if it's trapezoidal.  midpoint RFD gets applied in
            // second solve
            if (d_use_rfd && d_time_stepping_type == "SECOND_TRAPEZOIDAL_RULE")
            {
                spreadRFDForces(current_time);
            }
            
            // add spring forces.
            d_ib_force_fcn->computeLagrangianForce(d_F_data,
                                                   d_X_current_data,
                                                   d_U_data,
                                                   d_hierarchy,
                                                   finest_level_num,
                                                   current_time,
                                                   l_data_manager);
            // sync F
            d_F_data->beginGhostUpdate();
            d_F_data->endGhostUpdate();

            // spread force
            spreadForce(d_f_idx,
                        d_F_data,
                        d_X_current_data,
                        current_time);

            // NOTE: Any additional Eulerian forcing should be computed here and added
            // to the data associated with d_f_idx.

            ///////////////////////////////////////////////////////////////////////////
            // (2) Solve the fluid equations using the fluid solver registered with this
            // class using the forcing computed in (1).
            const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
            if(ins_num_cycles !=2) {
                TBOX_ERROR("Number of INS Cycles must be 2 with OVERDAMPED dynamics"
                           " and SECOND_MIDPOINT_RULE or SECOND_TRAPEZOIDAL RULE TS_TYPE");
            }
            // first cycle
            d_ins_hier_integrator->integrateHierarchy(current_time, new_time,0);

            //copy over eularian velocity to u_new
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int u_star_idx= var_db->mapVariableAndContextToIndex(
                d_ins_hier_integrator->getVelocityVariable(),
                d_ins_hier_integrator->getNewContext());
            d_hier_velocity_data_ops->copyData(d_u_idx, u_star_idx);

            // interpolate eularian velocity back to particles
            interpolateVelocity(d_u_idx,
                                d_U_data,
                                d_X_current_data,
                                current_time);

            // (3) UPDATE LAGRANGIAN VARIABLES FOR PREDICTOR
            // generate noise
            genrandn(d_Noise_data);
            // predictor: d_X_star = d_X_current + dt*d_xi*F + dt*d_U_data
            //                       + sqrt(2kTd_xi dt)*Noise
            // F will store the result of d_xi*f + U + sqrt(2kT d_xi /dt)*noise

            // add BG FLow if it is present
            if(d_bg_flow_flag)
            {
                addBGFlow(d_U_data,d_X_current_data,current_time);
            }
            
            // sync U
            d_U_data->beginGhostUpdate();
            d_U_data->endGhostUpdate();
            
            if (d_time_stepping_type == "SECOND_MIDPOINT_RULE")
            {
                // calculate total move rate and store it in d_F_data
                ierr = VecAXPBYPCZ(d_F_data->getVec(),
                                   1.0,
                                   sqrt(4.0*d_kT*d_xi/dt),
                                   d_xi,
                                   d_U_data->getVec(),
                                   d_Noise_data->getVec());
                IBTK_CHKERRQ(ierr);
                
                // now move the particle with this rate from X^n to X^*, midpoint predictor
                ierr = VecWAXPY(d_X_star_data->getVec(),
                                0.5*dt,
                                d_F_data->getVec(),
                                d_X_current_data->getVec());
                IBTK_CHKERRQ(ierr);
            } // midpoint predictor

            if (d_time_stepping_type == "SECOND_TRAPEZOIDAL_RULE")
            {
                // calculate total move rate and store it in d_U_data
                ierr = VecAXPBYPCZ(d_U_data->getVec(),
                                   d_xi,
                                   sqrt(2.0*d_kT*d_xi/dt),
                                   1.0,
                                   d_F_data->getVec(),
                                   d_Noise_data->getVec());
                IBTK_CHKERRQ(ierr);
                // move the particle with this rate from X^n to X^*, trapezoid predictor
                ierr = VecWAXPY(d_X_star_data->getVec(),
                                dt,
                                d_U_data->getVec(),
                                d_X_current_data->getVec());
                IBTK_CHKERRQ(ierr);
            } // Trapezoidal predictor
            
            // Zero Eularian forces for use in next IB fluid solve,
            d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
            
            // sync X star for force computation
            d_X_star_data->beginGhostUpdate();
            d_X_star_data->endGhostUpdate();
            
            // check that we haven't moved past 2 regrid intervals with the predictor step
            checkParticleMoveDistance(d_X_star_data);
            
            //(4) SOLVE FLUID EQUATIONS AGAIN
            if(d_point_force_flag)
            {   // if we have pointforces, then apply them and
                // normalize the eularian forces.
                // now applying F* forces evaluated at X*
                addPointForces(d_F_new_data,d_X_star_data,new_time);
            }   // if(d_point_force_flag)
            if (d_use_rfd && d_time_stepping_type == "SECOND_MIDPOINT_RULE")
            {
                spreadRFDForces(current_time);
            }
            
            // add spring forces.
            d_ib_force_fcn->computeLagrangianForce(d_F_new_data,
                                                   d_X_star_data,
                                                   d_U_data,
                                                   d_hierarchy,
                                                   finest_level_num,
                                                   new_time,
                                                   l_data_manager);
            // sync F
            d_F_new_data->beginGhostUpdate();
            d_F_new_data->endGhostUpdate();            
            
            //spread this Force + Noise at X_star
            spreadForce(d_f_idx,
                        d_F_new_data,
                        d_X_star_data,
                        current_time);
            
            // Second stokes solve, cycle_num = 1
            d_ins_hier_integrator->integrateHierarchy(current_time, new_time,1);

            //copy over eularian velocity to u_new
            const int u_new_idx= var_db->mapVariableAndContextToIndex(
                d_ins_hier_integrator->getVelocityVariable(),
                d_ins_hier_integrator->getNewContext());
            d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);

            // interpolate eularian velocity back to particles
            interpolateVelocity(d_u_idx,
                                d_U_star_data,
                                d_X_star_data,
                                current_time);
            // (5) FINAL UPDATE OF X, CORRECTOR
            // add shear velocity
            if (d_time_stepping_type == "SECOND_TRAPEZOIDAL_RULE")
            {
                if(d_bg_flow_flag)
                {
                    addBGFlow(d_U_star_data,d_X_star_data,current_time + dt);
                }
                // use U_star to store the total move,
                // after this we are missing just one half of the noise term
                ierr = VecAXPBYPCZ(d_U_star_data->getVec(),
                                   0.5,
                                   0.5*d_xi,
                                   0.5,
                                   d_U_data->getVec(),
                                   d_F_new_data->getVec());
                IBTK_CHKERRQ(ierr);

                // add the correct half of the noise
                ierr = VecAXPY(d_U_star_data->getVec(),
                               sqrt(0.5*d_kT*d_xi/dt),
                               d_Noise_data->getVec());
                IBTK_CHKERRQ(ierr);

                // Final Corrector for X
                ierr = VecWAXPY(d_X_new_data->getVec(),
                                dt,
                                d_U_star_data->getVec(),
                                d_X_current_data->getVec());
                IBTK_CHKERRQ(ierr);
            } // Trapezoidal Final Corrector
            if (d_time_stepping_type == "SECOND_MIDPOINT_RULE")
            {
                if(d_bg_flow_flag)
                {
                    addBGFlow(d_U_star_data,d_X_star_data,current_time + dt/2.0);
                }
                // generate second noise
                genrandn(d_Noise_2_data);
                
                // use U_star to store the total move, missing just one half of the noise term
                ierr = VecAXPBYPCZ(d_U_star_data->getVec(),
                                   sqrt(d_kT*d_xi/dt),
                                   d_xi,
                                   1.0,
                                   d_Noise_data->getVec(),
                                   d_F_new_data->getVec());
                IBTK_CHKERRQ(ierr);

                // add half the noise, which for midpoint is a second noise
                ierr = VecAXPY(d_U_star_data->getVec(),
                               sqrt(d_kT*d_xi/dt),
                               d_Noise_2_data->getVec());
                IBTK_CHKERRQ(ierr);
                
                // Final Corrector for X
                ierr = VecWAXPY(d_X_new_data->getVec(),
                                dt,
                                d_U_star_data->getVec(),
                                d_X_current_data->getVec());
                IBTK_CHKERRQ(ierr);
            }
        } // end if d_time_stepping_type == SECOND_TRAPEZOIDAL_RULE, SECOND_MIDPOINT_RULE
        else
        {
            TBOX_ERROR("For overdamped dynamics, time_stepping_type must "
                       "be MIDPOINT_RULE, TRAPEZOIDAL_RULE, FORWARD_EULER, "
                       "SECOND_TRAPEZOIDAL_RULE, or SECOND_MIDPOINT_RULE."); 
        }
    } // end if dynamics == "OVERDAMPED"
    else
    {
        TBOX_ERROR("dynamics must be OVERDAMPED."); 
    }
    return;
}// integratehierarchy


void
IBBrownianBlobHierarchyIntegrator::postprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const bool skip_synchronize_new_state_data,
    const int num_cycles)
{
    IBHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time,
        new_time,
        skip_synchronize_new_state_data,
        num_cycles);
    
    const int coarsest_level_num = 0;
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    PetscErrorCode ierr;
    
    // Deallocate Eulerian scratch data.
    for (int level_num = coarsest_level_num;
         level_num <= finest_level_num;
         ++level_num)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        level->deallocatePatchData(d_u_idx);
        level->deallocatePatchData(d_f_idx);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }

    // Deallocate the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    d_ins_hier_integrator->
        postprocessIntegrateHierarchy(
            current_time,
            new_time,
            skip_synchronize_new_state_data,
            ins_num_cycles);

    // Reset and deallocate IB data.
    // NOTE: We assume here that all IB data are assigned to the finest level of
    // the AMR patch hierarchy.
    ierr = VecSwap(d_X_current_data->getVec(),
                   d_X_new_data->getVec());
    IBTK_CHKERRQ(ierr);
    
    d_X_current_data = NULL; 
    d_X_star_data    = NULL;
    d_X_new_data     = NULL;
    d_U_data         = NULL;
    d_U_star_data    = NULL;
    d_F_data         = NULL;
    d_RFD_Noise_data   = NULL;
    d_RFD_loc_1_data   = NULL;
    d_RFD_loc_2_data   = NULL;
    d_F_new_data    = NULL;
    d_Noise_data     = NULL;
    d_Noise_2_data     = NULL;
    return;
}// postprocessIntegrateHierarchy

void
IBBrownianBlobHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;
     // Setup the fluid solver for explicit coupling.
    // NOTE: This will use the data associated with d_f_idx to provide forcing
    // for the fluid equations.
    d_ins_hier_integrator->registerBodyForceFunction(new IBEulerianForceFunction(this));

    // NOTE: Any additional implementation-specific initialization should be
    // performed here.
    
    // Register intermediate eularian variables

    // steven: just made this 4 for now, needs to be support/2 + 1, so this works for 6 pt.
    // TODO(steven): set this up to be dynamic.
    const IntVector<NDIM> side_ghosts = 4;
    d_u_star_var = new SideVariable<NDIM,double>(d_object_name+"::u_star");
    registerVariable(d_u_star_scratch_idx, d_u_star_var, side_ghosts);

    // Finish initializing the hierarchy integrator.  This function call should
    // come at the end of this function.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);
    return;
}// initializeHierarchyIntegrator


void
IBBrownianBlobHierarchyIntegrator::spreadRFDForces(double time) {
    // random "forces" added for RFD.
    PetscErrorCode ierr;
    genrandn(d_RFD_Noise_data);
    // set locations where RFD noise will be spread, X +/- 0.5*delta*W
    ierr = VecWAXPY(d_RFD_loc_1_data->getVec(),0.5*d_rfdelta,d_RFD_Noise_data->getVec(),d_X_current_data->getVec()); IBTK_CHKERRQ(ierr);
    ierr = VecWAXPY(d_RFD_loc_2_data->getVec(),-0.5*d_rfdelta,d_RFD_Noise_data->getVec(),d_X_current_data->getVec()); IBTK_CHKERRQ(ierr);
    
    // ghost update here for spreading
    d_RFD_loc_1_data->beginGhostUpdate();
    d_RFD_loc_1_data->endGhostUpdate();
    d_RFD_loc_2_data->beginGhostUpdate();
    d_RFD_loc_2_data->endGhostUpdate();
    
    // scale noise to be kT/d_rfdelta * Noise for RFD term
    // scale by 2 for trapezoidal PC
    ierr = VecScale(d_RFD_Noise_data->getVec(),2.*d_kT/d_rfdelta);IBTK_CHKERRQ(ierr);
    
    // sync noise data before spread
    d_RFD_Noise_data->beginGhostUpdate();
    d_RFD_Noise_data->endGhostUpdate();
    
    // spread noise at first location
    spreadForce(d_f_idx, d_RFD_Noise_data, d_RFD_loc_1_data, time);
    // scale by negative 1 to subtract
    ierr = VecScale(d_RFD_Noise_data->getVec(),-1.0);IBTK_CHKERRQ(ierr);
    
    // sync noise data before second spread. Uneeded, but just to be sure.
    d_RFD_Noise_data->beginGhostUpdate();
    d_RFD_Noise_data->endGhostUpdate();
    
    // spread noise at second location.
    spreadForce(d_f_idx, d_RFD_Noise_data, d_RFD_loc_2_data, time);
}  // spreadRFDForce.
    
void
IBBrownianBlobHierarchyIntegrator::spreadForce(const int f_data_idx,
                                               Pointer<LData> F_data,
                                               Pointer<LData> X_data,
                                               double time) {
    // NOTE: Time may not matter at all here, but we include
    // the correct time anyway.
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    // Zero forces on anchor vertices manually.
    resetAnchorPointValues(F_data);

    // spread force
    l_data_manager->spread(f_data_idx,
                           F_data,
                           X_data,
                           d_u_phys_bdry_op,
                           finest_level_num,
                           getProlongRefineSchedules(d_object_name+"::f"),
                           time,
                           /*F_needs_ghost_fill*/ true,
                           /*X_needs_ghost_fill*/ false);
    if(d_normalize_force_flag)
    {
        NormalizePointForces(F_data);
    }
} // spreadForce

void
IBBrownianBlobHierarchyIntegrator::interpolateVelocity(const int u_data_idx,
                                                       Pointer<LData> U_data,
                                                       Pointer<LData> X_data,
                                                       double time) {
    // If walls -> fill ghost cells before interpolation.
    // Probably there is a better way to decide if the system is periodic.
    if(d_normalize_force_flag==0) fillGhostCells(u_data_idx, time); 

    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    // Interpolate Velocity.
    l_data_manager->interp(u_data_idx,
                           U_data,
                           X_data,
                           finest_level_num,
                           getCoarsenSchedules(
                               d_object_name+"::u::CONSERVATIVE_COARSEN"),
                           getGhostfillRefineSchedules(
                               d_object_name+"::u"),
                           time);

    // Zero velocity on anchor vertices.  resetAnchorPointValues takes a vector.
    resetAnchorPointValues(U_data);
} // interpolateVelocity

void
IBBrownianBlobHierarchyIntegrator::registerIBLagrangianForceFunction(
    Pointer<IBLagrangianForceStrategySet> ib_force_fcn)
{
    d_ib_force_fcn = ib_force_fcn;
    return;
}
    
void
IBBrownianBlobHierarchyIntegrator::SetXi(double xi)
{
    // set xi to the input value, bare diffusion.
    d_xi = xi;
    return;
}

void
IBBrownianBlobHierarchyIntegrator::SetkT(double kT)
{
    // set kT to the input value
    d_kT = kT;
    return;
}

void
IBBrownianBlobHierarchyIntegrator::SetRfdelta(double inputrf,
                                      Pointer<CartesianGridGeometry<NDIM> > grid_geometry)
{
    // set RFD delta to the input value multiplied by mesh width
    const double* dx = grid_geometry->getDx();
    double dx_min = dx[0];
    for (int k = 1; k < NDIM; ++k)
    {
        dx[k] < dx_min ? dx_min = dx[k] : dx_min = dx_min;
    }
    d_rfdelta = inputrf*dx_min;
    return;
}

void
IBBrownianBlobHierarchyIntegrator::SetDynamics(string inputDynamics)
{
    // set dynamics to the input value
    d_dynamics = inputDynamics;
    return;
}// SetDynamics

void
IBBrownianBlobHierarchyIntegrator::SetLagrangianFunction(
    const std::string& /* object_name */,
    Pointer<Database> input_db,
    std::vector<mu::Parser>& parser)
{
    // Set a general lagrangian function from the input file.
    // Used for Point Force function and Background Flow function.
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
#endif
    // Read in user-provided constants.
    SAMRAI::tbox::Array<std::string> db_key_names = input_db->getAllKeys();
    for (int k = 0; k < db_key_names.size(); ++k)
    {
        const std::string& name = db_key_names[k];
        if (input_db->isDouble(name))
        {
            d_constants[name] = input_db->getDouble(name);
        }
        else if (input_db->isFloat(name))
        {
            d_constants[name] = input_db->getFloat(name);
        }
        else if (input_db->isInteger(name))
        {
            d_constants[name] = input_db->getInteger(name);
        }
    }
    // Assume vector-valued function.
    
    int d = 0;
    std::string key_name = "function_0";
    while (input_db->isString(key_name))
    {
        d_function_strings.push_back(input_db->getString(key_name));
        parser.resize(parser.size()+1);
        parser.back().SetExpr(d_function_strings.back());
        ++d;
        std::ostringstream stream;
        stream << "function_" << d;
        key_name = stream.str();
    }
    
    // Define the default and user-provided constants.
    const double pi = 3.1415926535897932384626433832795;
    for (std::vector<mu::Parser>::iterator it = parser.begin(); it != parser.end(); ++it)
    {
        // Various names for pi.
        it->DefineConst("pi", pi);
        it->DefineConst("Pi", pi);
        it->DefineConst("PI", pi);
        
        // User-provided constants.
        for (std::map<std::string,double>::const_iterator map_cit = d_constants.begin();
             map_cit != d_constants.end();
             ++map_cit)
        {
            it->DefineConst(map_cit->first, map_cit->second);
        }
        
        // Variables.
        it->DefineVar("T", &d_parser_time);
        it->DefineVar("t", &d_parser_time);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream stream;
            stream << d;
            const std::string postfix = stream.str();
            it->DefineVar("X" + postfix, &(d_parser_posn[d]));
            it->DefineVar("x" + postfix, &(d_parser_posn[d]));
            it->DefineVar("X_" + postfix, &(d_parser_posn[d]));
            it->DefineVar("x_" + postfix, &(d_parser_posn[d]));
        }
    }
    return;
}//SetLagrangianFunction

void
IBBrownianBlobHierarchyIntegrator::SetLagrangianFunctionNull(
    const std::string& /* object_name */,
    std::vector<mu::Parser>& parser)
{
    // set Lagrangian function to null, add no point force functions
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
#endif
    // default to 0 force
    std::string null_fcn_str = "0.0";
    int d = 0;
    std::string key_name = "function_0";
    while (d < NDIM)
    {
        d_function_strings.push_back(null_fcn_str);
        parser.resize(parser.size()+1);
        parser.back().SetExpr(d_function_strings.back());
        ++d;
        std::ostringstream stream;
        stream << "function_" << d;
        key_name = stream.str();
    }
    return;
}//SetPointForceFunctionNull

// void
// IBBrownianBlobHierarchyIntegrator::CheckEularianForceTotals()
// {
//     // helper function to print out total forces in each direction on
//     // the eularian mesh
//     const int finest_level_num = d_hierarchy->getFinestLevelNumber();
//     double force_spreaded[NDIM];
//     for(int k = 0; k < NDIM; ++k)
//     {
//         force_spreaded[k] = 0.0;
//     }
    
//     int wgt_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
//     for (int level_num = 0; level_num <= finest_level_num; ++level_num)
//     { // iterate through levels
//         Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
//         for (PatchLevel<NDIM>::Iterator p(level); p; p++)
//         {//iterate through patches
//             Pointer<Patch<NDIM> > patch = level->getPatch(p());
//             // get data on this patch
//             Pointer<SideData<NDIM,double> > f_data =
//                 patch->getPatchData(d_f_idx);
//             Pointer<SideData<NDIM,double> > wgt_data =
//                 patch->getPatchData(wgt_idx);
//             const Box<NDIM>& patch_box = patch->getBox();
//             for (int axis = 0; axis < NDIM; ++axis)
//             {
//                 for (SideIterator<NDIM> it(patch_box,axis); it; it++)
//                 {
//                     SideIndex<NDIM> i_s = it();
//                     force_spreaded[axis] += (*f_data)(i_s)*(*wgt_data)(i_s);
//                 }
//             } //loop through axes
//         } // loop through patches
//     } // loop through levels
//     for (unsigned int d = 0; d < NDIM; ++d)
//     {
//         cout << "Total Force component " << d << ": "
//              << force_spreaded[d] << std::endl;
//     }
// } //CheckEularianForceTotals

void
IBBrownianBlobHierarchyIntegrator::PrintDiagnostics(double current_time)
{
    //function to print out total momentum, stochastic flux, and velocity magnitude
    // Must be called after the at least one INS integration step to allow
    // calculation of the stochastic forces.
    if (d_diagnostic_interval <= 0)
    {
        // no diagnostics, just return
        return;
    }
    else if (d_diagnostic_counter == (d_diagnostic_interval - 1))
    {
        const int finest_level_num = d_hierarchy->getFinestLevelNumber();    
    
        d_body_force_fcn->setDataOnPatchHierarchy(
            d_f_idx,
            d_ins_hier_integrator->getBodyForceVariable(),
            d_hierarchy,
            current_time);

        // varioables ot hold total force, total velocity, total velocity squared, and area
        // contributing to total velocity squared
        double ftotal[NDIM];
        double utotal[NDIM];
        
        // area averaged <v^2>
        double usquared = 0.0;
        double area = 0.0;
        for (int i = 0; i < NDIM; ++i)
        {
            ftotal[i] = 0.0;
            utotal[i] = 0.0;
        }
        
        int wgt_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
        for (int level_num = 0; level_num <= finest_level_num; ++level_num)
        { // iterate through levels
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {//iterate through patches
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                // get data on this patch
                Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(d_f_idx);
                Pointer<SideData<NDIM,double> > u_data = patch->getPatchData(d_u_idx);
                Pointer<SideData<NDIM,double> > wgt_data = patch->getPatchData(wgt_idx);
                const Box<NDIM>& patch_box = patch->getBox();
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIterator<NDIM> it(patch_box,axis); it; it++)
                    {
                        SideIndex<NDIM> i_s = it();
                        ftotal[axis] += (*f_data)(i_s)*(*wgt_data)(i_s);
                        utotal[axis] += (*u_data)(i_s)*(*wgt_data)(i_s);
                        usquared += (*u_data)(i_s)*(*u_data)(i_s)*(*wgt_data)(i_s)*
                            (*wgt_data)(i_s);
                        area += (*wgt_data)(i_s)*(*wgt_data)(i_s);
                    }
                } //loop through axes
            } // loop through patches
        } // loop through levels
        
        usquared = usquared/area;
        usquared = SAMRAI_MPI::sumReduction(usquared);
        
        cout << std::endl;
        cout << "Printing Diagnostics" << std::endl;
        for (int i = 0; i < NDIM; ++i)
        {
            // reduce across processors
            ftotal[i] = SAMRAI_MPI::sumReduction(ftotal[i]);
            utotal[i] = SAMRAI_MPI::sumReduction(utotal[i]);
            
            cout << "sum (f) stochastic forces component " << i << " is "
                 << ftotal[i] << std::endl;
            cout << "sum(\\rho u) in direction " << i << " is "
                 << d_physical_rho*utotal[i] << std::endl;
        }
        
        cout << "<v^2> = " << usquared << std::endl;
        cout << std::endl;
        d_diagnostic_counter = 0;
    }// else if (diagnostic_counter == d_diagnostic_interval - 1 )
    else{
        d_diagnostic_counter += 1;
    }
    return;
}//PrintDiagnostics

void
IBBrownianBlobHierarchyIntegrator::addPointForces(Pointer<LData> d_F_data, Pointer<LData> d_X_current_data, double eval_time)
{
    // set up l_data_manager and finest_level_num
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    
    // get access to petsc vector of forces.
    Vec forcePetscVec= d_F_data->getVec();
    PetscScalar *force;
    VecGetArray(forcePetscVec,&force);
            
    // get access to petsc vector of positions.
    Vec positionPetscVec= d_X_current_data->getVec();
    PetscScalar *lposition;
    VecGetArray(positionPetscVec,&lposition);
            
    // iterate through local nodes to add point forces to particles
    //delong: using finest level number, because we assume IB data are
    // defined on this level.
    
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_level_num);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
         cit != local_nodes.end(); ++cit)
    {
        LNode* const node_idx = *cit;
        int local_idx = node_idx->getLocalPETScIndex();
        const Vector periodic_disp = node_idx->getPeriodicDisplacement();
        d_parser_time = eval_time;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_parser_posn[d] = lposition[local_idx*NDIM + d] +
                periodic_disp[d];
        }
        for (unsigned int d = 0; d < NDIM; ++d)
        {
	    force[local_idx*NDIM + d] += d_parsers[d].Eval();
            // keep track of the total force in each direction so we can
            // average and subtract
        }
    }
    VecRestoreArray(d_F_data->getVec(),&force);
    return;
}//addPointForces

void
IBBrownianBlobHierarchyIntegrator::NormalizePointForces(Pointer<LData> F_data)
{
    // Normalize the applied point forces, so that the total force
    // in each component is 0.
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    double total_force[NDIM];
    for(int k = 0; k < NDIM; ++k)
    {
        total_force[k] = 0.0;
    }
    // Get Total force in each direction.
    PetscScalar *force_vec;
    VecGetArray(F_data->getVec(), &force_vec);
    for (int level_num = 0; level_num <= finest_level_num; ++level_num)
    { // iterate through levels
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {//iterate through patches
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            int node_data_idx = l_data_manager->getLNodePatchDescriptorIndex();
            Pointer<LNodeSetData> node_data = patch->getPatchData(node_data_idx);
            for (LNodeSetData::DataIterator it=node_data->data_begin(patch_box);
                 it!=node_data->data_end(); ++it)
            {
                LNode* node_idx = *it;
                int local_idx = node_idx->getLocalPETScIndex();
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    total_force[d] += force_vec[local_idx*NDIM + d];
                }
            }
        }
    }
    const double volume = this->d_hier_math_ops->getVolumeOfPhysicalDomain();
    // Now normalize eularian forces.
    for (int level_num = 0; level_num <= finest_level_num; ++level_num)
    { // iterate through levels
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {//iterate through patches
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            // get data on this patch
            Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(d_f_idx);
            const Box<NDIM>& patch_box = patch->getBox();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> it(patch_box,axis); it; it++)
                {
                    SideIndex<NDIM> i_s = it();
                    (*f_data)(i_s) -= 1.0*total_force[axis]/volume;
                }
            } //loop through axes
        } // loop through patches
    } // loop through levels
    return;
} // NormalizePointForces

void
IBBrownianBlobHierarchyIntegrator::printEularianData(int data_idx)
{
    // Display eularian data.

    // get finest level
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();    
    
    for (int level_num = 0; level_num <= finest_level_num; ++level_num)
    { // iterate through levels
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {//iterate through patches
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            // get data on this patch
            Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(data_idx);
            const Box<NDIM>& patch_box = patch->getBox();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                cout << "Axis = " << axis << std::endl;
                for (SideIterator<NDIM> it(patch_box,axis); it; it++)
                {
                    SideIndex<NDIM> i_s = it();
                      cout << (*f_data)(i_s) << ",";                  
                }
                cout << std::endl;
            } //loop through axes
        } // loop through patches
    } // loop through levels
    return;
} // printEularianData

void
IBBrownianBlobHierarchyIntegrator::checkParticleMoveDistance(Pointer<LData> X_moved_data)
{
    // Check how far each particle has moved.  If any of them have moved further than
    // 2 * regrid interval, throw an error and abort.

    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    // set up lagrangian data
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    Pointer<LData> current_data = l_data_manager->getLData(
        LDataManager::POSN_DATA_NAME,finest_level_num);
    Vec X_new_values = X_moved_data->getVec();
    Vec X_old_values = d_X_last_regrid->getVec();
    PetscScalar *new_position, *old_position;
    VecGetArray(X_new_values,&new_position);
    VecGetArray(X_old_values,&old_position);
    double distance_moved;
    double max_dist_moved[NDIM];
    for (int i = 0; i < NDIM;++i)
    {
        max_dist_moved[i] = 0.0;
    }
    
    for (int level_num = 0; level_num <= finest_level_num; ++level_num)
    { // iterate through levels
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {//iterate through patches
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const double * dx = pgeom->getDx();
            int node_data_idx =l_data_manager->getLNodePatchDescriptorIndex();
            Pointer<LNodeSetData> node_data = patch->getPatchData(node_data_idx);
            for (LNodeSetData::DataIterator it=node_data->data_begin(patch_box);
                 it!=node_data->data_end(); ++it)
            {
                LNode* node_idx = *it;
                int local_idx = node_idx->getLocalPETScIndex();
                distance_moved = 0.0;
                // check distance moved component wise
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    distance_moved = abs(new_position[local_idx*NDIM + d]  -
                                         old_position[local_idx*NDIM + d]);
                    if (distance_moved > 2.0*d_regrid_alpha*dx[d]){ 
                        cout << "Distance moved: " << distance_moved << std::endl;
                        TBOX_ERROR("Particle moved more than 2 regrid_intervals "
                                   "between regrids");
                    }
                    if (distance_moved > max_dist_moved[d])
                    {
                        max_dist_moved[d] = distance_moved;
                    }
                }
            }
        }
    }
    return;
} //CheckParticleMoveDistance

bool
IBBrownianBlobHierarchyIntegrator::atRegridPointSpecialized() const
{
    //Check how far the particle has moved, and if it's moved more than one box since the last
    //regrid, then regrid again.

    // make sure all processes return the same answer
  
    // always regrid at the first step
    if (d_integrator_step == 0){
        return true;
    }

    // get finest level
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    // set up lagrangian data
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    Pointer<LData> current_data = l_data_manager->getLData(
        LDataManager::POSN_DATA_NAME,finest_level_num);
    Vec X_new_values = current_data->getVec();
    Vec X_old_values = d_X_last_regrid->getVec();
    PetscScalar *new_position, *old_position;
    VecGetArray(X_new_values,&new_position);
    VecGetArray(X_old_values,&old_position);
    double distance_moved;
    double max_dist_moved[NDIM];
    for (int i = 0; i < NDIM;++i)
    {
        max_dist_moved[i] = 0.0;
    }
        
    for (int level_num = 0; level_num <= finest_level_num; ++level_num)
    { // iterate through levels
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {//iterate through patches
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const double * dx = pgeom->getDx();
            int node_data_idx =l_data_manager->getLNodePatchDescriptorIndex();
            Pointer<LNodeSetData> node_data = patch->getPatchData(node_data_idx);
            for (LNodeSetData::DataIterator it=node_data->data_begin(patch_box);
                 it!=node_data->data_end(); ++it)
            {
                LNode* node_idx = *it;
                int local_idx = node_idx->getLocalPETScIndex();
                distance_moved = 0.0;
                // check distance moved component wise
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    distance_moved = std::abs(new_position[local_idx*NDIM + d] -
                                              old_position[local_idx*NDIM + d]);
                    if (distance_moved > 2.0*d_regrid_alpha*dx[d]) { 
                        cout << "Distance moved: " << distance_moved << std::endl;
                        TBOX_ERROR("Particle moved more than 2 regrid_intervals "
                                   "between regrids");
                    }
                    if (distance_moved > max_dist_moved[d])
                    {
                        max_dist_moved[d] = distance_moved;
                    }
                }
            }
            // NOTE: This is slower than it could be.  We should be able to stop as soon as we
            // hit any particle that has moved more than the grid cell,
            // without checking the rest.
            for (int d = 0; d < NDIM; ++d)
            {
                max_dist_moved[d] = SAMRAI_MPI::maxReduction(max_dist_moved[d]);
                if (max_dist_moved[d] > dx[d]*d_regrid_alpha)
                {  // Particle has moved more than one regrid_alpha, regrid.
                    cout << "Regridding, Distance moved is : "
                         << max_dist_moved[d] << std::endl;
                    return true;
                }
            }
        }
    }
    return false;
} //atRegridPointSpecialized

void
IBBrownianBlobHierarchyIntegrator::regridHierarchy()
{
    IBHierarchyIntegrator::regridHierarchy();
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    // set up lagrangian data
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    PetscErrorCode ierr;
    d_X_last_regrid = l_data_manager->createLData("regrid_data",finest_level_num,
                                                  NDIM,
                                                  /*maintain_data*/ false);
    
    d_X_current_data = l_data_manager->getLData(LDataManager::POSN_DATA_NAME,
                                                finest_level_num);
    
    // initialize last regrid data
    ierr = VecCopy(d_X_current_data->getVec(),d_X_last_regrid->getVec());
    IBTK_CHKERRQ(ierr);

    // if we have nonbonded forces, manage their temporary springs
    // TODO(steven): Switch to change between these springs added on the fly and the
    // cell search for nonbonded forces.
    if(!d_nonbonded_by_cell_flag && d_nonbonded_flag)
    {
        // remove old springs
        removeNonBondedSprings();
        // add new springs
        addNonBondedSprings();
        // move spring data to vectors of springs in IBStandardForceGen
        d_ib_force_fcn->initializeLevelData(d_hierarchy,
                                            finest_level_num,
                                            0.0,
                                            0.0,
                                            l_data_manager);
    }
    d_X_current_data = NULL;
//    d_wall_force_evaluator->rebuildWalls(l_data_manager,d_hierarchy);

    initializeAnchorIndices();
    return;
} // regridHierarchy

void
IBBrownianBlobHierarchyIntegrator::addBGFlow(Pointer<LData> U_data, Pointer<LData> X_data, double eval_time)
{
    // set up l_data_manager and finest_level_num
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    
    // get access to petsc vector of forces.
    PetscScalar *velocity;
    VecGetArray(U_data->getVec(),&velocity);
            
    // get access to petsc vector of positions.
    PetscScalar *lposition;
    VecGetArray(X_data->getVec(),&lposition);
            
    // iterate through local nodes to add point forces to particles
    // delong: using finest level number, because we assume IB data are
    // defined on this level.
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_level_num);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
         cit != local_nodes.end(); ++cit)
    {
        LNode* const node_idx = *cit;
        int local_idx = node_idx->getLocalPETScIndex();
        const Vector periodic_disp = node_idx->getPeriodicDisplacement();
        d_parser_time = eval_time;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_parser_posn[d] = lposition[local_idx*NDIM + d] +
                periodic_disp[d];
        }
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            velocity[local_idx*NDIM + d] += d_flow_parsers[d].Eval();
        }
    }
    VecRestoreArray(U_data->getVec(),&velocity);
    // Zero anchor points.
    resetAnchorPointValues(U_data);
    return;
}

void
IBBrownianBlobHierarchyIntegrator::registerNonBondedInteractionParameters(Pointer<Database> input_db)
{
    // set up parameters and interaction radii arrays for use with nonbonded springs.
    // TODO: Make this multiple species, arbitrary number of parameters
    //       and make it apply to the by cell force additions

    if(input_db->keyExists("num_species"))
    {
        d_num_species = input_db->getDouble("num_species");
    }
    else
    {
        d_num_species = 0;
    }

    if(d_num_species > 0)
    {

        // resize vectors to be able to hold interaction radius and parameters
        // for num_species x num_species combinations.
        d_interaction_radius.resize(d_num_species);
        for (std::vector<std::vector<double> >::iterator cit =
                 d_interaction_radius.begin();
             cit != d_interaction_radius.end(); ++cit)
        {
            (*cit).resize(d_num_species);
        }
        d_nonbonded_parameters.resize(d_num_species);
        for (std::vector<std::vector<std::vector<double> > >::iterator cit =
                 d_nonbonded_parameters.begin();
             cit != d_nonbonded_parameters.end(); ++cit)
        {
            (*cit).resize(d_num_species);
        }
        if (input_db->keyExists("interaction_radius"))
        {
            // NOTE: Hard coded for 1 species right now
            d_interaction_radius[0][0] = input_db->
                getDouble("interaction_radius");
        }
        else
        {
            // NOTE: Hard coded for 1 species right now
            d_interaction_radius[0][0] = 0.0;
        }
        // first parameter is the interaction radius.
        d_nonbonded_parameters[0][0].push_back(d_interaction_radius[0][0]);

        // NOTE: for now jsut one species.
        if(input_db->keyExists("nonbonded_parameter"))
        {
            d_nonbonded_parameters[0][0].push_back(input_db->getDouble("nonbonded_parameter"));
        }
        else
        {
            cout << "WARNING: Second Parameter Defaulting to 0 \n" << std::endl;
            d_nonbonded_parameters[0][0].push_back(0.0);
        }
    }
    return;
}  //registerNonBondedInteractionParameters
                                                                    
void
IBBrownianBlobHierarchyIntegrator::removeNonBondedSprings(void)
{
    // remove all springs corresponding to spring type d_nonbonded_force_idx
    // delete empty springspecs afterwards.

    // set up l_data_manager and finest_level_num
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    
    // iterate through local nodes to check springspecs for non-bonded forces
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_level_num);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
         cit != local_nodes.end();
         ++cit)
    {
        std::vector<int> remove_idxs;
        // get springforcespec object, and inspect it to see if it has a nonbonded spring
        const LNode* node_idx = *cit;
        IBSpringForceSpec* force_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
        if(force_spec)
        {
            std::vector<int>& force_fcn_idxs = force_spec->getForceFunctionIndices();
            // instead of using iterators, just use an index.
            // We need this index to modify the other vectors
            // in springSpec
            for(int k = (force_fcn_idxs.size() - 1); k >= 0; --k)
            {
                // check if the force is of the non-bonded type.
                if (force_fcn_idxs[k] == d_nonbonded_force_idx)
                {
                    // keep track of all indices we want to remove
                    // WARNING: remove_idxs must be in decreasing order.
                    remove_idxs.push_back(k);
                }
            }
            // remove this spring from the springspec object.
            IBSpringForceSpec** force_spec_ptr = &force_spec;
            removeSpringFromForceSpec(force_spec_ptr,remove_idxs);
        }
    }
} // removeNonBondedSprings

void
IBBrownianBlobHierarchyIntegrator::addNonBondedSprings(void)
{
    // Add temporary non-bonded springs to the springspec objects belonging to LNodes.
    // When InitializeLevelSpringData is called, these springs will be recognized by
    // the IBStandardForceGen, which can then apply their effects.

    // assumes ONE level for now (dx is for the coarsest level,
    // but is used at the finest level.)

    // also relies on the domain being a single box!

   // set up l_data_manager and finest_level_num
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();

    // get grid geometry and relevant lower and upper limits.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    if (!grid_geom->getDomainIsSingleBox())
    {
        TBOX_ERROR("physical domain must be a single box...\n");
    }

    // These will only work if the domain is a single box.
    assert(grid_geom->getDomainIsSingleBox());
    const Index<NDIM>& lower = grid_geom->getPhysicalDomain()[0].lower();
    const Index<NDIM>& upper = grid_geom->getPhysicalDomain()[0].upper();
    const double* const x_lower = grid_geom->getXLower();
    const double* const x_upper = grid_geom->getXUpper();
    const double* const dx = grid_geom->getDx();

    // iterate through local nodes to check if nodes have neighbors close enough
    // for nonbonded forces 
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_level_num);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

    // get position petsc scalar array
    PetscScalar* position;
    VecGetArray(d_X_current_data->getVec(),&position);
        
    // TODO: need to iterate through species, then for each species iterate through
    // nodes of that species.
    // NOTE: species hard coded for now to be one species, called 0
    int species_num_1 = 0;
    int species_num_2 = 0;
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
         cit != local_nodes.end(); ++cit)
    {
        LNode* node_idx = *cit;
        IBSpringForceSpec* force_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
        // interaction distance of this type of particle.
        // TODO: Make interaction depend on species/input
        //      double interaction = d_interaction_radius[species_num_1][species_num_2];
        //get node location  cell index
        const int local_idx = node_idx->getLocalPETScIndex();
        
        // pointer to first component of position for this particle
        const double* const X = &(position[local_idx*NDIM]); 

        const CellIndex<NDIM> cell_idx =
            IndexUtilities::getCellIndex(X, x_lower, x_upper, dx, lower, upper);
        
        // grow a box at this index to include all relevant area.
        Box<NDIM> local_node_box(cell_idx,cell_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            local_node_box.grow(
                d,
                static_cast<int>((d_interaction_radius[species_num_1][species_num_2]
                                                  + 2.0*d_regrid_alpha)/dx[d])+1); 
        }
        // iterate through particles in the box, and add springs.
        for (int level_num = 0; level_num <= finest_level_num; ++level_num)
        {
            // iterate through levels, not actually necessary here.
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                //iterate through patches
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                
                // get box that intersects patch with particle neighborhood
                const Box<NDIM>& patch_box = patch->getBox();
                const Box<NDIM> intersect_box = patch_box * (local_node_box);
                int node_data_idx = l_data_manager->getLNodePatchDescriptorIndex();
                Pointer<LNodeSetData> node_data = patch->getPatchData(node_data_idx);
                for (LNodeSetData::DataIterator it=node_data->data_begin(intersect_box);
                     it!=node_data->data_end(); ++it)
                {
                    // iterate through nodes in this relevant area
                    // don't connect to yourself or previously viewed particles.
                    if((*it)->getLagrangianIndex() <= (*cit)->getLagrangianIndex()) continue;

                    // add spring connecting these nodes to current springspec
                     int slave_idx = (*it)->getLagrangianIndex();
                    // need to get parameters from input
                    std::vector<double> parameters = d_nonbonded_parameters[species_num_1][species_num_2];
                    
                    if (!force_spec)
                    {
                        // create vectors for spring spec
                        std::vector<int> slave_idxs(1,slave_idx);
                        std::vector<int> force_fcn_idxs(1,d_nonbonded_force_idx);
                        std::vector<std::vector<double> > parameters_vec(1,parameters);

                        std::vector<Pointer<Streamable> > data_items =
                            node_idx->getNodeData();
                        data_items.push_back(
                            new IBSpringForceSpec(node_idx->getLagrangianIndex(),
                                                  slave_idxs,
                                                  force_fcn_idxs,
                                                  parameters_vec));
                        node_idx->setNodeData(data_items);
                        force_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
                    }
                    else
                    {
                        //if the springspec object is there, add entries for the new spring
                        std::vector<int>& force_fcn_idxs =
                            force_spec->getForceFunctionIndices();
                        force_fcn_idxs.push_back(d_nonbonded_force_idx);

                        std::vector<int>& slave_idxs = force_spec->getSlaveNodeIndices();
                        slave_idxs.push_back(slave_idx);

                        std::vector<std::vector<double> >& parameters_vec =
                            force_spec->getParameters();
                        parameters_vec.push_back(parameters);
                    }
                } // loop through slave indices
            } // loop through patches
        } // loop through levels
    } // loop through master nodes (local nodes).
    return;
} //addNonBondedSprings

void
IBBrownianBlobHierarchyIntegrator::removeSpringFromForceSpec(
    IBSpringForceSpec** force_spec_ptr,std::vector<int>& remove_idxs)
{
    // iterate through tagged pieces, and remove them
    // the way we do this is by copying the last element over the one being
    // deleted, and then popping the last element.
    // THIS WILL REARRANGE THE ORDER OF THE ELEMENTS IN THESE VECTORS,
    // BUT ALL IN THE SAME WAY.
    // This is fine so long as the nodes to be removed are in
    // decreasing order.

    // get access to vectors of spring data from springspec object
    std::vector<int>& force_fcn_idxs = (*force_spec_ptr)->getForceFunctionIndices();
    std::vector<int>& slave_idxs = (*force_spec_ptr)->getSlaveNodeIndices();
    std::vector<std::vector<double> >& spring_parameters =
        (*force_spec_ptr)->getParameters();
    
    // go from the back, since we delete by copying and popping off the last element.
    for (std::vector<int>::const_iterator rit = remove_idxs.begin();
         rit != remove_idxs.end(); ++rit)
    {
        // force functions
        force_fcn_idxs[*rit] = force_fcn_idxs.back();
        force_fcn_idxs.pop_back();

        // slave indices
        slave_idxs[*rit] = slave_idxs.back();
        slave_idxs.pop_back();

        // parameters
        spring_parameters[*rit] = spring_parameters.back();
        spring_parameters.pop_back();
    }
}  //removeSpringFromForceSpec

void
IBBrownianBlobHierarchyIntegrator::setNonBondedForceFunctionIdx(int force_fcn_idx)
{
    d_nonbonded_force_idx = force_fcn_idx;
} //setNonBondedForceFunctionIdx

void
IBBrownianBlobHierarchyIntegrator::initializeAnchorIndices()
{
    // First clear the old anchor indices.
    d_anchor_point_local_idxs.clear();
    // Populate d_anchor_point_local_idxs from node data.
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    // Cast from IBStrategy to IBMethod.
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_level_num);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBAnchorPointSpec* const anchor_point_spec = node_idx->getNodeDataItem<IBAnchorPointSpec>();
        if (anchor_point_spec)
        {
            d_anchor_point_local_idxs.insert(node_idx->getLocalPETScIndex());
        }
    }
}

void
IBBrownianBlobHierarchyIntegrator::resetAnchorPointValues(
    Pointer<LData> U_data)
{
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    PetscErrorCode ierr;
    Vec U_vec = U_data->getVec();
    double* U_arr;
    ierr = VecGetArray(U_vec, &U_arr);
    IBTK_CHKERRQ(ierr);
    for (std::set<int>::const_iterator cit = d_anchor_point_local_idxs.begin();
         cit != d_anchor_point_local_idxs.end();
         ++cit)
    {
        const int& i = *cit;
        for (int d = 0; d < NDIM; ++d)
        {
            U_arr[NDIM * i + d] = 0.0;
        }
    }
    ierr = VecRestoreArray(U_vec, &U_arr);
    IBTK_CHKERRQ(ierr);
    return;
} // resetAnchorPointValues

void
IBBrownianBlobHierarchyIntegrator::multiplyQuaternion(Pointer<LData> Q_data, Pointer<LData> W_data)
{
    // code to multiply an increment d_W_data into the existing orientation of the
    // particles, d_Q_data
    //  inputs:
    //       Q_data - pointer to LData holding the current orientations of the particles
    //       W_data - increment to rotate the current particles, saved as quaternions.
    //       
    //  outputs:      
    //       Q_data - is updated to be the quaternion product of Q and W
    //
    //////////////////////////////////////////////////

    // set up l_data_manager and finest_level_num
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    
    // get access to petsc vector of forces.
    PetscScalar *orientation;
    VecGetArray(Q_data->getVec(),&orientation);
            
    // get access to petsc vector of positions.
    PetscScalar *qincrement;
    VecGetArray(W_data->getVec(),&qincrement);
            
    // iterate through local nodes to add point forces to particles.
    // Use finest level number, because we assume IB data are
    // defined on this level.
    
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_level_num);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        LNode* const node_idx = *cit;
        int local_idx = node_idx->getLocalPETScIndex();
        // 3d case for quaternion multiplication
#if (NDIM == 3)
        
        // temporarily hold orientation so we can use it for the products.
        double temp_orientation[4];
        // first entry is s1s2 - p1 dot p2
        for(int d = 0; d < 4; ++d)
        {

            temp_orientation[d] = orientation[local_idx*4 + d];
        }
        
        
        orientation[local_idx*4] = qincrement[local_idx*4]*temp_orientation[0] -
                                   qincrement[local_idx*4 + 1]*temp_orientation[1] -
                                   qincrement[local_idx*4 + 2]*temp_orientation[2] -
                                   qincrement[local_idx*4 + 3]*temp_orientation[3];

        
        //rest of the entries are s1p2 + s2p1 - p1 x p2
        for (unsigned int d = 0; d < 3; ++d)
        {
            int prev = (d-1 + 3) % 3;
            int next = (d+1) % 3;
            
            orientation[local_idx*4+d+1] = qincrement[local_idx*4]*temp_orientation[d+1] +
                temp_orientation[0]*qincrement[local_idx*4+d+1] +
                qincrement[local_idx*4+1+next]*temp_orientation[1+prev] -
                qincrement[local_idx*4+1+prev]*temp_orientation[1+next];
        }
#endif
        // 2D case for quaternion multiplication
#if (NDIM == 2)
        double temp_orientation[2];

        // hold onto current orientation
        for(int d = 0; d < 2; ++d)
        {
            temp_orientation[d] = orientation[local_idx*4 + d];
        }
        
        // entry one, new s1 = s1s2 + p1p2
        orientation[local_idx*2] = temp_orientation[0]*qincrement[local_idx*2] -
                            temp_orientation[1]*qincrement[local_idx*2 + 1];

        // entry two, new p1 = s1p2 + s2p1
        orientation[local_idx*2 + 1] = temp_orientation[0]*qincrement[local_idx*2 + 1] +
            temp_orientation[1]*qincrement[local_idx*2 + 0];
#endif
    }
    VecRestoreArray(Q_data->getVec(),&orientation);
    VecRestoreArray(W_data->getVec(),&qincrement);
    return;
} // multiplyQuaternion


    
void
IBBrownianBlobHierarchyIntegrator::fillGhostCells(int in, const double time){

  //Copied from spde/BoundaryJS/InterpOperator.C
  typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
  std::vector<InterpolationTransactionComponent> transaction_comps;
  InterpolationTransactionComponent x_component(in, 
						/*DATA_REFINE_TYPE*/ "NONE", //NONE
						/*USE_CF_INTERPOLATION*/ true, //true
						/*DATA_COARSEN_TYPE*/ "CUBIC_COARSEN", //Cubic_coarsen
						/*BDRY_EXTRAP_TYPE*/ "LINEAR", //Linear
						/*CONSISTENT_TYPE_2_BDRY*/ false, //false
						d_u_bc_coefs[0],
						/*d_fill_pattern*/ NULL); //Null
  transaction_comps.push_back(x_component);
  
  // Setup 
  int d_coarsest_ln = 0;
  int d_finest_ln   = d_hierarchy->getFinestLevelNumber();   

  const bool homogeneous_bc = false;
  SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_bdry_fill; 
  d_hier_bdry_fill = new IBTK::HierarchyGhostCellInterpolation();
  d_hier_bdry_fill->initializeOperatorState(transaction_comps, d_hierarchy, d_coarsest_ln, d_finest_ln);
  d_hier_bdry_fill->setHomogeneousBc(homogeneous_bc);

  // Fill ghost cells of x
  d_hier_bdry_fill->fillData(time);

  //Deallocate data
  d_hier_bdry_fill->deallocateOperatorState();
  d_hier_bdry_fill.setNull();
  transaction_comps.clear();
  
  //cout << "fillGhostCells ENDS" << endl;
  return;
}


// set velocity boundary conditions
void
IBBrownianBlobHierarchyIntegrator::setVelocityBC(vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> *u_bc_coefs){
    d_u_bc_coefs = u_bc_coefs;
    return;
}


