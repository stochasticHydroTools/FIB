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

// Config files

//#include <IBAMR_prefix_config.h>
//#include <IBTK_prefix_config.h> // these are only for debugging?
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include <ibamr/IBMethod.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBLagrangianForceStrategySet.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSStaggeredStochasticForcing.h>
#include <ibamr/RNG.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>


#include "IBMobilityEstimator.h"


// Function prototypes
void
output_data(
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
    Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
    LDataManager* l_data_manager,
    const int iteration_num,
    const double loop_time,
    const int output_level,
    const string& data_dump_dirname);

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(
    int argc,
    char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc,&argv,NULL,NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    {// cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool is_from_restart = app_initializer->isFromRestart();
        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        // set defaults and switches
        double regrid_alpha = 2.0;   // regrid parameter default
        string spread_fcn = "IB_4";  // default interpolation and spreading function
        bool is_vel_normalized = false;  // variable for enforcing velocity normalization
        
        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getStringWithDefault("solver_type", "STAGGERED");

        // use staggered solver with stochastic forces
        navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        SAMRAI::hier::IntVector<NDIM> gcw;
        double interaction = app_initializer->getComponentDatabase("IBHierarchyIntegrator")->getDoubleWithDefault("interaction_radius",0.0);
        
        for (int k = 0; k < NDIM; ++k)
        {
            gcw[k] = 4.0*regrid_alpha + interaction;
        }
        if(app_initializer->getComponentDatabase("IBMethod")->keyExists("delta_fcn"))
        {
            spread_fcn = app_initializer->getComponentDatabase("IBMethod")->getString("delta_fcn");
        }

        // Initialize l_data_manager
        LDataManager* l_data_manager = LDataManager::getManager("IBMethod::LDataManager",
                                                                spread_fcn,
                                                                spread_fcn,
                                                                gcw,
                                                                false);
        Pointer<IBMethod> ib_method_ops = new IBMethod(
            "IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<IBMobilityEstimator> time_integrator =
            new IBMobilityEstimator("IBMobilityEstimator",
                app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                ib_method_ops, navier_stokes_integrator,grid_geometry);
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
            "PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>(
                "StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
            "LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);

        // Create Eulerian initial condition specification objects.
        // NOTE: may not need these, test.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }
	time_integrator->setVelocityBC(&u_bc_coefs);
	        
        // Create stochastic forcing function specification object
        Pointer<INSStaggeredStochasticForcing> f_stoch_fcn = new INSStaggeredStochasticForcing("INSStaggeredStochasticForcing", app_initializer->getComponentDatabase("INSStaggeredStochasticForcing"), navier_stokes_integrator);
        time_integrator->registerBodyForceFunction(f_stoch_fcn);

        //////////////////////////////////////
        // check inputs for consistency:
        ////////////////////////////////////

        // check to make sure velocity is normalized
        if (app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator")->keyExists("normalize_velocity"))
        {
                
            is_vel_normalized = app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator")->getBool("normalize_velocity");
        }

        // Velocity must be normalized if all BCs are periodic, otherwise it must not be.
        // necessary for rho = 0 and periodic BCs only
        if (!is_vel_normalized)
        {
            
            if (periodic_shift.min() > 0)
            {
                TBOX_ERROR("NORMALIZE_VELOCITY must be TRUE with all periodic boundaries.");
            }
            
        }
        else
        {
            if (!(periodic_shift.min() > 0))
            {
                TBOX_ERROR("When boundaries are not all periodic, NORMALIZE_VELOCITY must be FALSE.");
            }
        }
        
        // Seed the random number generator.
        int seed = 0;
        if (input_db->keyExists("SEED"))
        {
            seed = input_db->getInteger("SEED");
        }
        else
        {
            TBOX_ERROR("Key data `seed' not found in input.");
        }
        RNG::parallel_seed(seed);

        // We are primarily interested in the Lagrangian data so we can skip writing Eulerian data:
        int output_level = 2; // -1=none, 0=txt, 1=silo+txt, 2=visit+silo+txt, 3=visit+SAMRAI+silo+txt, >3=verbose
        if (input_db->keyExists("OUTPUT_LEVEL"))
        {
            output_level = input_db->getInteger("OUTPUT_LEVEL");
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.	
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        double dt = time_integrator->getMaximumTimeStepSize();
        time_integrator->advanceHierarchy(dt);
        
        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        (void) l_data_manager; // silence unused warning.
    }// cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
}// main

void
output_data(
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
    Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
    LDataManager* l_data_manager,
    const int iteration_num,
    const double loop_time,
    const int output_level,
    const string& data_dump_dirname)
{
    if (output_level >= 0)
    {
        pout << "\nWriting output files at timestep # " <<  iteration_num << " t=" << loop_time << "\n";
    }
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;

    // Write Cartesian data.
    if (output_level >= 3)
    {
        Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
        hier_db->create(file_name);
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        ComponentSelector hier_data;
        hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(), navier_stokes_integrator->getCurrentContext()));
        hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(), navier_stokes_integrator->getCurrentContext()));
        patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
        hier_db->putDouble("loop_time", loop_time);
        hier_db->putInteger("iteration_num", iteration_num);
        hier_db->close();
    }

    // Write Lagrangian data.
    if (output_level >= 0)
    {
        const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
        Pointer<LData> X_data = l_data_manager->getLData("X", finest_hier_level);
        Vec X_petsc_vec = X_data->getVec();
        Vec X_lag_vec;
        Vec X_real_petsc_vec;
        PetscErrorCode ierr;
        VecDuplicate(X_petsc_vec, &X_lag_vec);
        VecDuplicate(X_petsc_vec,&X_real_petsc_vec);
        ierr = VecCopy(X_petsc_vec,X_real_petsc_vec); IBTK_CHKERRQ(ierr);
        PetscScalar *outposition;

        VecGetArray(X_real_petsc_vec,&outposition);

        // output virtual location, the "real" location in a periodic setting.
        const Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_hier_level);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
        
            LNode* const node_idx = *cit;
            int local_idx = node_idx->getLocalPETScIndex();
            const Vector& periodic_disp = node_idx->getPeriodicDisplacement();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                outposition[local_idx*NDIM + d] = outposition[local_idx*NDIM + d] + periodic_disp[d];
            }
        }
        VecRestoreArray(X_real_petsc_vec,&outposition);
        l_data_manager->scatterPETScToLagrangian(X_real_petsc_vec, X_lag_vec, finest_hier_level);
        file_name = data_dump_dirname + "/" + "X.";
        sprintf(temp_buf, "%05d", iteration_num);
        file_name += temp_buf;
        PetscViewer viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
        VecView(X_lag_vec, viewer);
        PetscViewerDestroy(&viewer);
        VecDestroy(&X_lag_vec);
    }
    return;
}// output_data


