/* Copyright (C) 2015-2017 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Header file for GPU solver functions
 *
*/
/** \file
 * \brief Header file for GPU solver functions
*/
#ifndef PHLEX_GPU_SOLVER_H_
#define PHLEX_GPU_SOLVER_H_
#include <cuda.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "phlex_solver.h"

#define CUDA_MAX_THREADS 1024

typedef struct {
  PMC_C_FLOAT * host_state;             // host pointer to the working state array
  PMC_C_FLOAT * dev_state;              // device pointer to the working state array
  PMC_C_FLOAT * host_env;               // host pointer to the working environmental state
  PMC_C_FLOAT * dev_env;                // device pointer to the working environmental state
  PMC_SOLVER_C_FLOAT * host_deriv;      // host pointer to the working deriv array
  PMC_SOLVER_C_FLOAT * dev_deriv;       // device pointer to the working deriv array
  PMC_SOLVER_C_FLOAT * host_jac;        // host pointer to the working Jacobian data
  PMC_SOLVER_C_FLOAT * dev_jac;         // device pointer to the working Jacobian data
  int deriv_size;                       // size of the derivative array
  int jac_size;                         // size of the Jacobian data array
  int deriv_block;                      // block to calculate derivative for this state
  int jac_block;                        // block to calculate Jacobian for this state
  int env_block;                        // block to update environmental state
  int deriv_start_id;                   // id of the first species in this state on the
                                        // shared derivative array
  int jac_start_id;                     // id of the first species in this state on the
                                        // shared Jacobian array
  void * host_rxn_dev_data;             // host pointer to reaction device data
  void * dev_rxn_dev_data;              // device pointer to reaction device data
} ModelDeviceData;

typedef struct {
  int deriv_blocks;                     // number of blocks to use during deriv calc
  int jac_blocks;                       // number of blocks to use during jac calc
  int env_blocks;                       // number of blocks to use during env update
  int deriv_threads;                    // number of threads to use during deriv calc
  int jac_threads;                      // number of threads to use during jac calc
  int env_threads;                      // number of threads to use during env update
  int n_states;                         // number of states to solve simultaneously
  ModelDeviceData * host_model_dev_data;// model device data array - one for each state
  ModelDeviceData * dev_model_dev_data; // model device data array - one for each state
} SolverDeviceData;

void phlex_gpu_solver_new( SolverData *solver_data );
void phlex_gpu_solver_update_env_state( SolverData *solver_data );
int phlex_gpu_solver_f( realtype t, N_Vector y, N_Vector deriv, void *solver_data );
int phlex_gpu_solver_Jac( realtype t, N_Vector y, N_Vector deriv, SUNMatrix J,
        void *solver_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3 );
void phlex_gpu_solver_solver_device_data_free_vp( void * solver_device_data );
void phlex_gpu_solver_solver_device_data_free( SolverDeviceData solver_device_data );
void phlex_gpu_solver_model_device_data_free( ModelDeviceData model_device_data );

#endif
