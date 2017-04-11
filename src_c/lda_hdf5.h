/*
 * lda_hdf5.h
 *
 *  Created on: Apr 10, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_LDA_HDF5_H_
#define SRC_C_LDA_HDF5_H_

#include "spectral_density.h"

void hdf5_create_file( const char* hdf_filename );

size_t hdf5_get_dataset_array_length( const char *hdf_filename, const char* dataset_name );

double hdf5_get_sampling_frequency( const char* hdf_filename );

size_t hdf5_get_num_time_samples( const char* hdf_filename );

size_t hdf5_get_num_strains( const char* hdf_filename );

void hdf5_load_array( const char *hdf_filename, const char *dataset_name, double *data);

psd_t* hdf5_load_psd( const char *hdf_filename );

void hdf5_save_psd( const char *hdf_filename, psd_t *psd );



#endif /* SRC_C_LDA_HDF5_H_ */
