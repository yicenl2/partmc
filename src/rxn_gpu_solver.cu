/* Copyright (C) 2015-2018 Matthew Dawson
 * Licensed under the GNU General Public License version 2 or (at your
 * option) any later version. See the file COPYING for details.
 *
 * Reaction-specific functions for use by the solver
 *
*/
/** \file
 * \brief Reaction solver functions
*/
#include "cuda_util.h"
extern "C" {
#include "phlex_solver.h"
#include "rxn_gpu_solver.h"
#include "rxn_solver.h"
#include "rxns/rxn_gpu_arrhenius.h"
}

// Reaction types (Must match parameters defined in pmc_rxn_factory)
#define RXN_ARRHENIUS 1
#define RXN_TROE 2
#define RXN_CMAQ_H2O2 3
#define RXN_CMAQ_OH_HNO3 4
#define RXN_PHOTOLYSIS 5
#define RXN_HL_PHASE_TRANSFER 6
#define RXN_AQUEOUS_EQUILIBRIUM 7
#define RXN_ZSR_AEROSOL_WATER 8
#define RXN_PDFITE_ACTIVITY 9
#define RXN_SIMPOL_PHASE_TRANSFER 10
#define RXN_CONDENSED_PHASE_ARRHENIUS 11

/** \brief Assemble a set of indices for each reaction to solve with GPUs
 *
 * \param model_dev_data Model device data
 * \param host_rxn_data Pointer to the host reaction data
 */
void rxn_gpu_solver_new( ModelDeviceData * model_dev_data, void * orig_rxn_data )
{
  // Allocate space for a new RxnDeviceData object
  HANDLE_ERROR( cudaHostAlloc( (void**) &(model_dev_data->host_rxn_dev_data),
                               sizeof(RxnDeviceData),
                               cudaHostAllocWriteCombined |
                                  cudaHostAllocMapped
                             ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( (void**) &(model_dev_data->dev_rxn_dev_data),
                                          (void*) model_dev_data->host_rxn_dev_data,
                                          0 
                                        ) );
  RxnDeviceData *rxn_dev_data = (RxnDeviceData*) model_dev_data->host_rxn_dev_data;

  // Get the number of reactions
  int *rxn_data = (int*) (orig_rxn_data);
  int n_rxn = *(rxn_data++);

  // Count the number and size of reactions with GPU solver functions
  int n_gpu_rxn = 0;
  
  // save space for the number of reactions (as a double to maintain alignment)
  size_t size_gpu_rxn = sizeof(PMC_C_FLOAT);
  
  // Loop through the reactions and add up the space they require on the
  // rxn_data array
  for (int i_rxn = 0; i_rxn < n_rxn; i_rxn++ ) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);
    
    // Save the starting position of the current reaction's data
    char * first_datum = (char*) rxn_data;

      switch (rxn_type) {
        case RXN_AQUEOUS_EQUILIBRIUM :
          rxn_data = (int*) rxn_aqueous_equilibrium_skip(
                    (void*) rxn_data);
          break;
        case RXN_ARRHENIUS :
          n_gpu_rxn++;
          // save space for the reaction type (as a double to maintain alignment)
          size_gpu_rxn += sizeof(PMC_C_FLOAT);
          rxn_data = (int*) rxn_arrhenius_skip(
                    (void*) rxn_data);
          size_gpu_rxn += sizeof(char) * ((char*)rxn_data - first_datum);
          break;
        case RXN_CMAQ_H2O2 :
          rxn_data = (int*) rxn_CMAQ_H2O2_skip(
                    (void*) rxn_data);
          break;
        case RXN_CMAQ_OH_HNO3 :
          rxn_data = (int*) rxn_CMAQ_OH_HNO3_skip(
                    (void*) rxn_data);
          break;
        case RXN_CONDENSED_PHASE_ARRHENIUS :
          rxn_data = (int*) rxn_condensed_phase_arrhenius_skip(
                    (void*) rxn_data);
          break;
        case RXN_HL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_HL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_PDFITE_ACTIVITY :
          rxn_data = (int*) rxn_PDFiTE_activity_skip(
                    (void*) rxn_data);
          break;
        case RXN_PHOTOLYSIS :
          rxn_data = (int*) rxn_photolysis_skip(
                    (void*) rxn_data);
          break;
        case RXN_SIMPOL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_SIMPOL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_TROE :
          rxn_data = (int*) rxn_troe_skip(
                    (void*) rxn_data);
          break;
        case RXN_ZSR_AEROSOL_WATER :
          rxn_data = (int*) rxn_ZSR_aerosol_water_skip(
                    (void*) rxn_data);
          break;
      }
  }

  // Allocate space for the reaction parameters
  HANDLE_ERROR( cudaHostAlloc( (void**) &(rxn_dev_data->host_rxn_data),
                               size_gpu_rxn,
                               cudaHostAllocWriteCombined |
                                  cudaHostAllocMapped
                             ) );
  HANDLE_ERROR( cudaHostGetDevicePointer( (void**) &(rxn_dev_data->dev_rxn_data),
                                          (void*) rxn_dev_data->host_rxn_data,
                                          0 
                                        ) );

  if( n_gpu_rxn > 0 ) {
    // Allocate space for the starting index of each reaction
    HANDLE_ERROR( cudaHostAlloc( (void**) &(rxn_dev_data->host_rxn_data_start),
                                 n_gpu_rxn * sizeof(int),
                                 cudaHostAllocWriteCombined |
                                    cudaHostAllocMapped
                               ) );
    HANDLE_ERROR( cudaHostGetDevicePointer( (void**) &(rxn_dev_data->dev_rxn_data_start),
                                            (void*) rxn_dev_data->host_rxn_data_start,
                                            0
                                          ) );
  }

  int * host_rxn_data = (int*) rxn_dev_data->host_rxn_data;
  char * host_rxn_data_init = (char*) host_rxn_data;
  int * host_rxn_data_start = rxn_dev_data->host_rxn_data_start;
 
  // Set the number of reactions
  rxn_dev_data->n_rxn = *(host_rxn_data) = n_gpu_rxn;
  
  // Advance by size of a double to maintain alignment
  host_rxn_data += sizeof(PMC_C_FLOAT) / sizeof(int);

  // Reset the reaction data pointer
  rxn_data = (int*) (orig_rxn_data);
  n_rxn = *(rxn_data++);
  
  // Loop through the reactions, copying rxn data to the rxn data block
  for (int i_rxn = 0, i_gpu_rxn = 0; i_rxn < n_rxn && i_gpu_rxn < n_gpu_rxn; i_rxn++ ) {

    // Get the reaction type
    int rxn_type = *(rxn_data++);

      switch (rxn_type) {
        case RXN_AQUEOUS_EQUILIBRIUM :
          rxn_data = (int*) rxn_aqueous_equilibrium_skip(
                    (void*) rxn_data);
          break;
        case RXN_ARRHENIUS :
          // save the offset for the beginning of this reaction's data
          host_rxn_data_start[ i_gpu_rxn++ ] = 
                  ((char*) host_rxn_data) - host_rxn_data_init;
          
          // save the reaction type
          *(host_rxn_data) = rxn_type;
          
          // advance by the size of a double to maintain alignment
          host_rxn_data += sizeof(PMC_C_FLOAT) / sizeof(int);

          // copy the reaction data to the host array
          rxn_gpu_arrhenius_copy_data( (void*) rxn_data, (void*) host_rxn_data );

          // advance the pointers to the next reaction
          rxn_data = (int*) rxn_arrhenius_skip( (void*) rxn_data );
          host_rxn_data = (int*) rxn_arrhenius_skip( (void*) host_rxn_data );

          break;
        case RXN_CMAQ_H2O2 :
          rxn_data = (int*) rxn_CMAQ_H2O2_skip(
                    (void*) rxn_data);
          break;
        case RXN_CMAQ_OH_HNO3 :
          rxn_data = (int*) rxn_CMAQ_OH_HNO3_skip(
                    (void*) rxn_data);
          break;
        case RXN_CONDENSED_PHASE_ARRHENIUS :
          rxn_data = (int*) rxn_condensed_phase_arrhenius_skip(
                    (void*) rxn_data);
          break;
        case RXN_HL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_HL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_PDFITE_ACTIVITY :
          rxn_data = (int*) rxn_PDFiTE_activity_skip(
                    (void*) rxn_data);
          break;
        case RXN_PHOTOLYSIS :
          rxn_data = (int*) rxn_photolysis_skip(
                    (void*) rxn_data);
          break;
        case RXN_SIMPOL_PHASE_TRANSFER :
          rxn_data = (int*) rxn_SIMPOL_phase_transfer_skip(
                    (void*) rxn_data);
          break;
        case RXN_TROE :
          rxn_data = (int*) rxn_troe_skip(
                    (void*) rxn_data);
          break;
        case RXN_ZSR_AEROSOL_WATER :
          rxn_data = (int*) rxn_ZSR_aerosol_water_skip(
                    (void*) rxn_data);
          break;
      }
  }
}

/** \brief Update the environmental state for reactions with GPU solver functions
  *
  * \param sdd Solver device data
  * \param env New environmental state
  */
__global__ void rxn_gpu_update_env_state( SolverDeviceData sdd )
{
  
  // Get the state and reaction to solve for
  int i_state = 0;
  int i_rxn = threadIdx.x;
  int rxns_to_solve = 0;
  ModelDeviceData * mdd;
  RxnDeviceData * rdd;
  for( ; i_state < sdd.n_states; i_state++ ) {
    mdd = &( sdd.dev_model_dev_data[ i_state ] );
    if( mdd->env_block != blockIdx.x ) continue;
    rdd = ( RxnDeviceData* ) ( mdd->dev_rxn_dev_data );
    if( i_rxn < rdd->n_rxn ) { rxns_to_solve = 1; break; }
    i_rxn -= rdd->n_rxn;
  }
  if( rxns_to_solve == 0 ) return;

  // Get the number of reactions
  int * rxn_data = (int*) rdd->dev_rxn_data;
  int n_rxn = *(rxn_data);

  // Return if there are no reactions to update
  if( n_rxn == 0 ) return;

  // Advance the rxn data pointer to the reaction's data
  char *curr_rxn = ( (char*) rxn_data ) + 
                   rdd->dev_rxn_data_start[ i_rxn ];
  rxn_data = (int*) curr_rxn;

  // Get the reaction type
  int rxn_type = *(rxn_data);

  // Advance by the size of a double to maintain alignment
  rxn_data += sizeof(PMC_C_FLOAT) / sizeof(int);

  // Update environmental state for appropriate reaction type
  switch (rxn_type) {
    case RXN_ARRHENIUS :
      rxn_data = (int*) rxn_gpu_arrhenius_update_env_state( (void*) rxn_data, *mdd );
      break;
    default :
      printf("\nPartMC Internal Error: invalid rxn type in GPU update env state.\n"
             "block: %d thread: %d rxn: %d rxn_type: %d\n",
             blockIdx.x, threadIdx.x, i_rxn, rxn_type);
  }
  __syncthreads();

}

/** \brief Calculate the time derivative for reactions with GPU solver functions
 * 
 * \param sdd Device data needed for solving
 * \param time_step Current solver time step (s)
 */
__global__ void rxn_gpu_calc_deriv( SolverDeviceData sdd, PMC_C_FLOAT time_step)
{
  
  // Get the state and reaction to solve for
  int i_state = 0;
  int i_rxn = threadIdx.x;
  int rxns_to_solve = 0;
  ModelDeviceData * mdd;
  RxnDeviceData * rdd;
  for( ; i_state < sdd.n_states; i_state++ ) {
    mdd = &( sdd.dev_model_dev_data[ i_state ] );
    if( mdd->deriv_block != blockIdx.x ) continue;
    rdd = ( RxnDeviceData* ) ( mdd->dev_rxn_dev_data );
    if( i_rxn < rdd->n_rxn ) { rxns_to_solve = 1; break; }
    i_rxn -= rdd->n_rxn;
  }
  if( rxns_to_solve == 0 ) return;

  // Set up a shared derivative array
  __shared__ PMC_SOLVER_C_FLOAT shr_dev_deriv[ MAX_SHARED_ARRAY_SIZE_ ];
  PMC_SOLVER_C_FLOAT * state_deriv = &( shr_dev_deriv[ mdd->deriv_start_id ] );

  // Initialize the derivative array
  for( int i_spec = i_rxn; i_spec < mdd->deriv_size; i_spec += rdd->n_rxn )
    state_deriv[ i_spec ] = 0.0;
  __syncthreads();

  // Get the number of reactions
  int * rxn_data = (int*) rdd->dev_rxn_data;
  int n_rxn = *(rxn_data);

  // Return if there are no reactions to solve
  if( n_rxn == 0 ) return;

  if ( rxns_to_solve > 0 ) {
    
    // Advance the rxn data pointer to the reaction's data
    char *curr_rxn = ( (char*) rxn_data ) + 
                     rdd->dev_rxn_data_start[ i_rxn ];
    rxn_data = (int*) curr_rxn;

    // Get the reaction type
    int rxn_type = *(rxn_data);

    // Advance by the size of a double to maintain alignment
    rxn_data += sizeof(PMC_C_FLOAT) / sizeof(int);

    // Add derivative contribution from appropriate reaction type
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_gpu_arrhenius_calc_deriv_contrib( 
                            (void*) rxn_data, *mdd, state_deriv );
        break;
      default :
        printf("\nPartMC Internal Error: invalid rxn type in GPU solver.\n"
               "block: %d thread: %d rxn: %d rxn_type: %d\n",
               blockIdx.x, threadIdx.x, i_rxn, rxn_type);
    }
  }
  __syncthreads();

  // Add derivative contributions from this block to the primary deriv array
  for( int i_spec = i_rxn; i_spec < mdd->deriv_size; i_spec += rdd->n_rxn )
    atomicAdd( &( mdd->dev_deriv[ i_spec ] ), state_deriv[ i_spec ] );

}

/** \brief Calculate the Jacobian for reactions with GPU solver functions
 * 
 * \param sdd Device data needed for solving
 * \param time_step Current solver time step (s)
 */
__global__ void rxn_gpu_calc_jac( SolverDeviceData sdd, PMC_C_FLOAT time_step)
{
  
  // Get the state and reaction to solve for
  int i_state = 0;
  int i_rxn = threadIdx.x;
  int rxns_to_solve = 0;
  ModelDeviceData * mdd;
  RxnDeviceData * rdd;
  for( ; i_state < sdd.n_states; i_state++ ) {
    mdd = &( sdd.dev_model_dev_data[ i_state ] );
    if( mdd->jac_block != blockIdx.x ) continue;
    rdd = ( RxnDeviceData* ) ( mdd->dev_rxn_dev_data );
    if( i_rxn < rdd->n_rxn ) { rxns_to_solve = 1; break; }
    i_rxn -= rdd->n_rxn;
  }
  if( rxns_to_solve == 0 ) return;

  // Set up a shared Jacobian array
  __shared__ PMC_SOLVER_C_FLOAT shr_dev_jac[ MAX_SHARED_ARRAY_SIZE_ ];
  PMC_SOLVER_C_FLOAT * state_jac = &( shr_dev_jac[ mdd->jac_start_id ] );

  // Initialize the Jacobian array
  for( int i_elem = i_rxn; i_elem < mdd->jac_size; i_elem += rdd->n_rxn )
    state_jac[ i_elem ] = 0.0;
  __syncthreads();

  // Get the number of reactions
  int * rxn_data = (int*) rdd->dev_rxn_data;
  int n_rxn = *(rxn_data);

  // Return if there are no reactions to solve
  if( n_rxn == 0 ) return;

  if ( rxns_to_solve > 0 ) {
    
    // Advance the rxn data pointer to the reaction's data
    char *curr_rxn = ( (char*) rxn_data ) + 
                     rdd->dev_rxn_data_start[ i_rxn ];
    rxn_data = (int*) curr_rxn;

    // Get the reaction type
    int rxn_type = *(rxn_data);

    // Advance by the size of a double to maintain alignment
    rxn_data += sizeof(PMC_C_FLOAT) / sizeof(int);

    // Add Jacobian contribution from appropriate reaction type
    switch (rxn_type) {
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_gpu_arrhenius_calc_jac_contrib( 
                            (void*) rxn_data, *mdd, state_jac );
        break;
      default :
        printf("\nPartMC Internal Error: invalid rxn type in GPU Jac function.\n"
               "block: %d thread: %d rxn: %d rxn_type: %d\n",
               blockIdx.x, threadIdx.x, i_rxn, rxn_type);
    }
  }
  __syncthreads();

  // Add Jacobian contributions from this block to the primary Jacobian data array
  for( int i_elem = i_rxn; i_elem < mdd->jac_size; i_elem += rdd->n_rxn )
    atomicAdd( &( mdd->dev_jac[ i_elem ] ), state_jac[ i_elem ] );

}

/** \brief Print reaction data for reactions with GPU functions
  * 
  * \param host_rxn_data Pointer to the reaction data
  */
void rxn_gpu_solver_print( void * host_rxn_data )
{
  int * rxn_data = (int*) host_rxn_data;
  int n_rxn = *(rxn_data);
  rxn_data += sizeof(PMC_C_FLOAT) / sizeof(int);

  // Loop through the reactions, printing each
  for( int i_rxn = 0; i_rxn < n_rxn; i_rxn++ ) {

    printf( "\n Reaction %d position %d", i_rxn, 
        ((char*)rxn_data) - ((char*)host_rxn_data) );

    // Get the reaction type
    int rxn_type = *(rxn_data);
    rxn_data += sizeof(PMC_C_FLOAT) / sizeof(int);

    // Call the appropriate function
    switch (rxn_type) {
      case RXN_PDFITE_ACTIVITY :
        rxn_data = (int*) rxn_PDFiTE_activity_print(
                  (void*) rxn_data);
	break;
      case RXN_AQUEOUS_EQUILIBRIUM :
        rxn_data = (int*) rxn_aqueous_equilibrium_print(
                  (void*) rxn_data);
	break;
      case RXN_ARRHENIUS :
        rxn_data = (int*) rxn_arrhenius_print(
                  (void*) rxn_data);
	break;
      case RXN_CMAQ_H2O2 :
        rxn_data = (int*) rxn_CMAQ_H2O2_print(
                  (void*) rxn_data);
	break;
      case RXN_CMAQ_OH_HNO3 :
        rxn_data = (int*) rxn_CMAQ_OH_HNO3_print(
                  (void*) rxn_data);
	break;
      case RXN_CONDENSED_PHASE_ARRHENIUS :
        rxn_data = (int*) rxn_condensed_phase_arrhenius_print(
                  (void*) rxn_data);
	break;
      case RXN_HL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_HL_phase_transfer_print(
                  (void*) rxn_data);
	break;
      case RXN_PHOTOLYSIS :
        rxn_data = (int*) rxn_photolysis_print(
                  (void*) rxn_data);
	break;
      case RXN_SIMPOL_PHASE_TRANSFER :
        rxn_data = (int*) rxn_SIMPOL_phase_transfer_print(
                  (void*) rxn_data);
	break;
      case RXN_TROE :
        rxn_data = (int*) rxn_troe_print(
                  (void*) rxn_data);
	break;
      case RXN_ZSR_AEROSOL_WATER :
        rxn_data = (int*) rxn_ZSR_aerosol_water_print(
                  (void*) rxn_data);
	break;
    }
  }
}

/** \brief Free memory associated with a RxnDeviceData object
  *
  * \param rxn_dev_data Object to free
  */
void rxn_gpu_solver_free( void * rxn_dev_data )
{
  RxnDeviceData *rd = (RxnDeviceData*) rxn_dev_data;
  HANDLE_ERROR( cudaFreeHost( rd->host_rxn_data ) );
  HANDLE_ERROR( cudaFreeHost( rd->host_rxn_data_start ) );
  HANDLE_ERROR( cudaFreeHost( rd ) );
}
