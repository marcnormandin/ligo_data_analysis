#ifndef SRC_C_LDA_HDF5_H_
#define SRC_C_LDA_HDF5_H_

#include <stddef.h>
#include <gsl/gsl_vector.h>

#if defined (__cplusplus)
extern "C" {
#endif

void hdf5_create_file( const char* hdf_filename );

size_t hdf5_get_dataset_array_length( const char *hdf_filename, const char* dataset_name );

double hdf5_get_sampling_frequency( const char* hdf_filename );

size_t hdf5_get_num_time_samples( const char* hdf_filename );

size_t hdf5_get_num_strains( const char* hdf_filename );

void hdf5_load_array( const char *hdf_filename, const char *dataset_name, double *data);

void hdf5_create_group(const char *hdf5_filename, const char* group_name);

void hdf5_save_array(const char *hdf5_filename, const char* group_name, const char *array_name, size_t len, double *array);

void hdf5_save_attribute_string( const char *hdf5_filename, const char *group_name, const char *attribute_name, const char *data);
void hdf5_save_attribute_double( const char *hdf5_filename, const char *group_name, const char *attribute_name, size_t len_array, const double *data );
void hdf5_save_attribute_ulong( const char *hdf5_filename, const char *group_name, const char *attribute_name, size_t len_array, const unsigned long *data );
void hdf5_save_attribute_gsl_vector( const char *hdf5_filename, const char *group_name, const char *attribute_name, const gsl_vector *data );

#if defined (__cplusplus)
}
#endif

#endif /* SRC_C_LDA_HDF5_H_ */
