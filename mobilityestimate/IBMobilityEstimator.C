// Filename: IBMobilityEstimator.C
// Created on Sept 2014 by Steven Delong
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

#include "IBMobilityEstimator.h"

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

IBMobilityEstimator::IBMobilityEstimator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<IBMethod> ib_method_ops,
    Pointer<INSHierarchyIntegrator> ins_hier_integrator,
    Pointer<CartesianGridGeometry<NDIM> > /* grid_geometry */)
    : IBHierarchyIntegrator(object_name, input_db, ib_method_ops,
                            ins_hier_integrator, /*register_for_restart*/ false),
      d_normalize_force_flag(0)
{
    // do we normalize forces or not?  This should be true with a fully periodic domain.
    if (input_db->keyExists("normalize_force"))
    {
        d_normalize_force_flag = input_db->getInteger("normalize_force");
    }
    else
    {
        TBOX_ERROR("Must indicate if we should normalize forces or not with "
                   "normalize_force (1 means normalize)");
    }
    d_mobility = new double[NDIM*NDIM];
    
    d_out_mob_file.open("./mobility.dat",std::ofstream::out);
}

IBMobilityEstimator::~IBMobilityEstimator()
{
    delete[] d_mobility;
    return;
}// ~IBBrownianBlobHierarchyIntegrator

void
IBMobilityEstimator::estimateMobility(vector<double> particle_positions,
                                      double* mobility)
{
    regridHierarchy();
//    preprocessIntegrateHierarchy(0.0, 1.0, 1);
    int current_time = 0.0;
    const int coarsest_level_num = 0;
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    PetscErrorCode ierr;    
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
    }
    d_X_current_data  = l_data_manager->getLData(LDataManager::POSN_DATA_NAME,
                                                 finest_level_num);
    d_U_data          = l_data_manager->getLData(LDataManager::VEL_DATA_NAME,
                                                 finest_level_num);
    d_F_data          = l_data_manager->createLData("F", finest_level_num,
                                                    NDIM, false);    
     // no timestepping, just computing the mobility
     int ins_cycle_num = 0;
     for (int k = 0; k < NDIM*NDIM; ++k)
     {
         mobility[k] = 0.0;
     }
     // TODO: get number of particles here, and go through
     // one component of one particle at a time, to get the full
     // M = J L^-1 S
     for (int mob_comp = 0; mob_comp < NDIM; ++mob_comp) {
         // don't actually move the particle, just compute mobility
         // from an applied force in each direction.
         // Zero variables that we will add .
         ierr = VecZeroEntries(d_F_data->getVec());  IBTK_CHKERRQ(ierr);
         ierr = VecZeroEntries(d_U_data->getVec());  IBTK_CHKERRQ(ierr);
         d_hier_velocity_data_ops->setToScalar(d_f_idx,0.0);
         d_hier_velocity_data_ops->setToScalar(d_u_idx,0.0);
      
         // apply single force for mobility calc
         PetscScalar *force; 
         VecGetArray(d_F_data->getVec(),&force);
         // iterate through local nodes to add point forces to particles
         // using finest level number, because we assume IB data are
         // defined on this level.
         // Apply forces on all particles at mob_comp (which changes
         // from "time step" to "time step").
         const Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_level_num);
         const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
         for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
              cit != local_nodes.end();
              ++cit)
         {
             LNode* const node_idx = *cit;
             int local_idx = node_idx->getLocalPETScIndex();
             force[local_idx*NDIM + mob_comp] = 1.0;
         }
         VecRestoreArray(d_F_data->getVec(), &force);
         d_F_data->beginGhostUpdate();
         d_F_data->endGhostUpdate();

         // Set particle position.
         PetscScalar *position;
         VecGetArray(d_X_current_data->getVec(),&position);
         // iterate through local nodes to add point forces to particles
         // using finest level number, because we assume IB data are
         // defined on this level.
         // Apply forces on all particles at mob_comp (which changes
         // from "time step" to "time step").
         // There should just be one particle here.
         for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
              cit != local_nodes.end();
              ++cit)
         {
             LNode* const node_idx = *cit;
             int local_idx = node_idx->getLocalPETScIndex();
             for (int k = 0; k < NDIM; ++k) {
                 position[k] = particle_positions[k];
             }
         }
         VecRestoreArray(d_X_current_data->getVec(), &position);
         d_X_current_data->beginGhostUpdate();
         d_X_current_data->endGhostUpdate();
         // spread force
         spreadForce(d_f_idx,
                     d_F_data,
                     d_X_current_data,
                     current_time);
         // l_data_manager->spread(d_f_idx,
         //                        d_F_data,
         //                        d_X_current_data,
         //                        d_u_phys_bdry_op,
         //                        finest_level_num,
         //                        getProlongRefineSchedules(d_object_name+"::f"),
         //                        1.0 /* new time, unused */,
         //                        /*F_needs_ghost_fill*/ false,
         //                       /*X_needs_ghost_fill*/ false);

         // normalize force, but not when we have walls
         if(d_normalize_force_flag)
         {
             NormalizePointForces();
         }
         // Solve fluid equations for velocity, should be a stokes solve.
         d_ins_hier_integrator->preprocessIntegrateHierarchy(0.0, 1.0, 1);
         d_ins_hier_integrator->integrateHierarchy(0.0 /* current_time */,
                                                   1.0 /*new_time, unused*/,
                                                   1);
         // Copy over eularian velocity to u_new.
         VariableDatabase<NDIM>* var_db =
             VariableDatabase<NDIM>::getDatabase();
         const int u_new_idx= var_db->mapVariableAndContextToIndex(
             d_ins_hier_integrator->getVelocityVariable(),
             d_ins_hier_integrator->getNewContext());
         d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
      
         // interpolate eularian velocity back to particles
         interpolateVelocity(d_u_idx,
                             d_U_data,
                             d_X_current_data,
                             current_time);
         // l_data_manager->interp(d_u_idx,
         //                        d_U_data,
         //                        d_X_current_data,
         //                        finest_level_num,
         //                        getCoarsenSchedules(
         //                            d_object_name+"::u::CONSERVATIVE_COARSEN"),
         //                        getGhostfillRefineSchedules(
         //                            d_object_name+"::u"),
         //                        0.0 /* current_time */);
//         d_U_data->beginGhostUpdate();
//         d_U_data->endGhostUpdate();            
      
         // output velocity = column of mobility
         PetscScalar *velocity;
         VecGetArray(d_U_data->getVec(), &velocity);
         VecGetArray(d_X_current_data->getVec(), &position);
      
         // For the mobility estimate, we just print the mobility in
         // one line in column major order.  
         for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
              cit != local_nodes.end();
              ++cit)
         {
             LNode* const node_idx = *cit;
             int local_idx = node_idx->getLocalPETScIndex();
             for (int k = 0; k < NDIM; ++k)
             {
                 mobility[mob_comp*NDIM + k] = velocity[local_idx*NDIM + k];
             }
         }
      
         // keep particle where it is, and change the component along which force is added
         ierr = VecCopy(d_X_current_data->getVec(),
                        d_X_new_data->getVec());
         IBTK_CHKERRQ(ierr);
     }
    return;
}  // estimateMobility

void
IBMobilityEstimator::NormalizePointForces()
{
    // Normalize the applied point forces, so that the total force
    // in each component is 0.
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    double total_force[NDIM];
    for(int k = 0; k < NDIM; ++k)
    {
        total_force[k] = 0.0;
    }
    int wgt_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    for (int level_num = 0; level_num <= finest_level_num; ++level_num)
    { // iterate through levels
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {//iterate through patches.
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            // get data on this patch.
            Pointer<SideData<NDIM,double> > f_data = patch->getPatchData(d_f_idx);
            Pointer<SideData<NDIM,double> > wgt_data = patch->getPatchData(wgt_idx);
            const Box<NDIM>& patch_box = patch->getBox();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> it(patch_box,axis); it; it++)
                {
                    SideIndex<NDIM> i_s = it();
                    total_force[axis] += (*f_data)(i_s)*(*wgt_data)(i_s);
                }
            } //loop through axes
        } // loop through patches
    } // loop through levels

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
IBMobilityEstimator::preprocessIntegrateHierarchy(
    double current_time,
    double new_time,
    int num_cycles) {
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
    d_X_new_data      = l_data_manager->createLData("X_new", finest_level_num,
                                                    NDIM);
    d_U_data          = l_data_manager->getLData(LDataManager::VEL_DATA_NAME,
                                                 finest_level_num);
    d_F_data          = l_data_manager->createLData("F", finest_level_num,
                                                    NDIM);
}

void
IBMobilityEstimator::integrateHierarchy(
    double current_time,
    double new_time,
    int cycle_num) {
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    PetscErrorCode ierr;
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    // no timestepping, just computing the mobility
    int ins_cycle_num = 0;
    for (int k = 0; k < NDIM*NDIM; ++k)
    {
        d_mobility[k] = 0.0;
    }
    for (int mob_comp = 0; mob_comp < NDIM; ++mob_comp) {
        // don't actually move the particle, just compute mobility
        // from an applied force in each direction.
        // Zero variables that we will add .
        ierr = VecZeroEntries(d_F_data->getVec());  IBTK_CHKERRQ(ierr);
        ierr = VecZeroEntries(d_U_data->getVec());  IBTK_CHKERRQ(ierr);
        d_hier_velocity_data_ops->setToScalar(d_f_idx,0.0);
        d_hier_velocity_data_ops->setToScalar(d_u_idx,0.0);
        
        // apply single force for mobility calc
        PetscScalar *force, *position; 
        VecGetArray(d_F_data->getVec(),&force);
        // iterate through local nodes to add point forces to particles
        // using finest level number, because we assume IB data are
        // defined on this level.
        // Apply forces on all particles at mob_comp (which changes
        // from "time step" to "time step").
        const Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_level_num);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
             cit != local_nodes.end();
             ++cit)
        {
            LNode* const node_idx = *cit;
            int local_idx = node_idx->getLocalPETScIndex();
            force[local_idx*NDIM + mob_comp] = 1.0;
        }
        VecRestoreArray(d_F_data->getVec(), &force);
        d_F_data->beginGhostUpdate();
        d_F_data->endGhostUpdate();

        // spread force
        l_data_manager->spread(d_f_idx,
                               d_F_data,
                               d_X_current_data,
                               d_u_phys_bdry_op,
                               finest_level_num,
                               getProlongRefineSchedules(d_object_name+"::f"),
                               1.0 /* new time, unused */,
                               /*F_needs_ghost_fill*/ false,
                               /*X_needs_ghost_fill*/ false);

        // normalize force, but not when we have walls
        if(d_normalize_force_flag)
        {
            NormalizePointForces();
        }
        // Solve fluid equations for velocity, should be a stokes solve.
        d_ins_hier_integrator->integrateHierarchy(0.0 /* current_time */,
                                                  1.0 /*new_time, unused*/,
                                                  ins_cycle_num);
        //copy over eularian velocity to u_new
        VariableDatabase<NDIM>* var_db =
            VariableDatabase<NDIM>::getDatabase();
        const int u_new_idx= var_db->mapVariableAndContextToIndex(
            d_ins_hier_integrator->getVelocityVariable(),
            d_ins_hier_integrator->getNewContext());
        d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
        // interpolate eularian velocity back to particles
        l_data_manager->interp(d_u_idx,
                               d_U_data,
                               d_X_current_data,
                               finest_level_num,
                               getCoarsenSchedules(
                                   d_object_name+"::u::CONSERVATIVE_COARSEN"),
                               getGhostfillRefineSchedules(
                                   d_object_name+"::u"),
                               0.0 /* current_time */);
        d_U_data->beginGhostUpdate();
        d_U_data->endGhostUpdate();            

        // output velocity = column of mobility
        PetscScalar *velocity;
        VecGetArray(d_U_data->getVec(), &velocity);
        VecGetArray(d_X_current_data->getVec(), &position);
        
        // For the mobility estimate, we just print the mobility in
        // one line in column major order.
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
             cit != local_nodes.end();
             ++cit)
        {
            LNode* const node_idx = *cit;
            int local_idx = node_idx->getLocalPETScIndex();
            for (int k = 0; k < NDIM; ++k)
            {
                d_mobility[mob_comp*NDIM + k] = velocity[local_idx*NDIM + k];
                d_out_mob_file << velocity[local_idx*NDIM + k] << ",";
            }
            d_out_mob_file << std::endl;
        }
        // keep particle where it is, and change the component along which force is added
        ierr = VecCopy(d_X_current_data->getVec(),
                       d_X_new_data->getVec());
        IBTK_CHKERRQ(ierr);
    }
    return;
}

void
IBMobilityEstimator::postprocessIntegrateHierarchy(
    double current_time,
    double new_time,
    bool skip_synchronize_new_state_data,
    int num_cycles) {
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
    d_U_data         = NULL;
    d_F_data         = NULL;
    return;
}


void
IBMobilityEstimator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;
     // Setup the fluid solver for explicit coupling.
    // NOTE: This will use the data associated with d_f_idx to provide forcing
    // for the fluid equations.
    d_ins_hier_integrator->registerBodyForceFunction(new IBEulerianForceFunction(this));

    // Finish initializing the hierarchy integrator.  This function call should
    // come at the end of this function.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);
    return;
}// initializeHierarchyIntegrator


void
IBMobilityEstimator::regridHierarchy()
{
    IBHierarchyIntegrator::regridHierarchy();
    d_X_current_data = NULL;
    return;
}


// set velocity boundary conditions
void
IBMobilityEstimator::setVelocityBC(vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> *u_bc_coefs){
    d_u_bc_coefs = u_bc_coefs;
    return;
}


void
IBMobilityEstimator::spreadForce(const int f_data_idx,
                                 Pointer<LData> F_data,
                                 Pointer<LData> X_data,
                                 double time) {
    // NOTE: Time may not matter at all here, but we include
    // the correct time anyway.
    const int finest_level_num = d_hierarchy->getFinestLevelNumber();
    Pointer<IBMethod> p_ib_method_ops = d_ib_method_ops;
    LDataManager* l_data_manager = p_ib_method_ops->getLDataManager();
    
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
        NormalizePointForces();
    }
} // spreadForce

void
IBMobilityEstimator::interpolateVelocity(const int u_data_idx,
                                         Pointer<LData> U_data,
                                         Pointer<LData> X_data,
                                         double time) {
    //If walls -> fill ghost cells before interpolation
    if(d_normalize_force_flag==0) fillGhostCells(u_data_idx, time); //Probably there is a better way to decide if the system is periodic or not

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

} // interpolateVelocity

void
IBMobilityEstimator::fillGhostCells(int in, const double time){
    
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
  
  return;
}
