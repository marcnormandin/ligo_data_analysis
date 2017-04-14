#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "hdf5_file.h"

void hdf5_create_file( const char* hdf_filename ) {
	assert(hdf_filename != NULL);

	hid_t file_id = H5Fcreate( hdf_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error. Unable to create the HDF5 file (%s). Aborting.\n",
				hdf_filename);
		exit(-1);
	}
	H5Fclose(file_id);
}

size_t hdf5_get_dataset_array_length( const char *hdf_filename, const char* dataset_name ) {
	assert(hdf_filename != NULL);
	assert(dataset_name != NULL);

	hid_t file_id, dataset_id, dspace_id;
	herr_t status;

	/* Open the file */
	file_id = H5Fopen( hdf_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error opening hdf5 file (%s) for reading. Aborting.\n", hdf_filename);
	}

	/* Read the dataset */
	dataset_id = H5Dopen2(file_id, dataset_name, H5P_DEFAULT);
	if (dataset_id < 0) {
		fprintf(stderr, "Error opening the dataset (%s) from the file (%s). Aborting.\n",
				dataset_name, hdf_filename);
		exit(-1);
	}

	/* Get the length of the dataset */
	dspace_id = H5Dget_space(dataset_id);
	int ndims = H5Sget_simple_extent_ndims(dspace_id);
	hssize_t len = H5Sget_simple_extent_npoints(dspace_id);

	H5Dclose(dataset_id);
	H5Fclose(file_id);

	return len;
}


size_t hdf5_get_num_strains( const char* hdf_filename ) {
	assert(hdf_filename != NULL);

	hid_t file_id, dataset_id, dspace_id;
	herr_t status;

	/* Open the file */
	file_id = H5Fopen( hdf_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error opening the hdf5 file (%s). Aborting.\n",
				hdf_filename);
		exit(-1);
	}

	/* Read the attribute */
	double num;
	const char *attribute_name = "num_strains";
	status = H5LTget_attribute_double( file_id, "/", attribute_name, &num);
	if (status < 0) {
		fprintf(stderr, "Error reading the attribute (%s) from the hdf5 file (%s). Aborting.\n",
				attribute_name, hdf_filename);
		exit(-1);
	}

	H5Fclose(file_id);

	return (int)num;
}

double hdf5_get_sampling_frequency( const char* hdf_filename )
{
	assert(hdf_filename != NULL);
	hid_t file_id, dataset_id, dspace_id;
	herr_t status;

	/* Open the file */
	file_id = H5Fopen( hdf_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error opening the hdf5 file (%s). Aborting.\n",
				hdf_filename);
		exit(-1);
	}

	/* Read the attribute */
	double fs;
	const char *attribute_name = "fs";
	status = H5LTget_attribute_double( file_id, "/", attribute_name, &fs);
	if (status < 0) {
		fprintf(stderr, "Error reading the attribute (%s) from the hdf5 file (%s). Aborting.\n",
				attribute_name, hdf_filename);
		exit(-1);
	}

	H5Fclose(file_id);

	return fs;
}

size_t hdf5_get_num_time_samples( const char* hdf_filename ) {
	assert(hdf_filename != NULL);

	hid_t file_id, dataset_id, dspace_id;
	herr_t status;

	/* Open the file */
	file_id = H5Fopen( hdf_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error opening the hdf5 file (%s). Aborting.\n",
				hdf_filename);
		exit(-1);
	}

	/* Read the attribute */
	double num;
	const char *attribute_name = "num_time_samples";
	status = H5LTget_attribute_double( file_id, "/", attribute_name, &num);
	if (status < 0) {
		fprintf(stderr, "Error reading the attribute (%s) from the hdf5 file (%s). Aborting.\n",
				attribute_name, hdf_filename);
		exit(-1);
	}

	H5Fclose(file_id);

	return (size_t)num;
}

void hdf5_load_array( const char *hdf_filename, const char *dataset_name, double *data) {
	assert(hdf_filename != NULL);
	assert(dataset_name != NULL);
	assert(data != NULL);

	hid_t file_id;
	herr_t status;

	file_id = H5Fopen( hdf_filename, H5F_ACC_RDWR, H5P_DEFAULT );
	if (file_id < 0) {
		fprintf(stderr, "Error: Unable to open the HDF5 file (%s). Aborting.\n", hdf_filename);
		exit(-1);
	}

	status = H5LTread_dataset_double( file_id, dataset_name, data );
	if (status < 0) {
		fprintf(stderr, "Error reading the dataset (%s) from the file (%s). Aborting.\n",
				dataset_name, hdf_filename);
		exit(-1);
	}

	H5Fclose( file_id );
}

void hdf5_create_group(const char *hdf5_filename, const char* group_name) {
	assert(hdf5_filename != NULL);
	assert(group_name != NULL);

	hid_t file_id, group_id;
	herr_t status;

	file_id = H5Fopen( hdf5_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error opening (%s) to create the group (%s). Aborting.\n",
				hdf5_filename, group_name);
		exit(-1);
	}

	group_id = H5Gcreate( file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (group_id < 0) {
		fprintf(stderr, "Error creating the group (%s) in the file (%s). Aborting.\n",
				group_name, hdf5_filename);
		exit(-1);
	}

	H5Gclose(group_id);
	H5Fclose(file_id);
}

void hdf5_save_array(const char *hdf5_filename, const char* group_name, const char *array_name, size_t len, double *array) {
	assert(hdf5_filename != NULL);
	assert(group_name != NULL);
	assert(array_name != NULL);
	assert(array != NULL);

	hid_t file_id, group_id;
	herr_t status;

	file_id = H5Fopen( hdf5_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error opening (%s) to create the group (%s). Aborting.\n",
				hdf5_filename, group_name);
		exit(-1);
	}

	group_id = H5Gopen2( file_id, group_name, H5P_DEFAULT);
	if (group_id < 0) {
		fprintf(stderr, "Error creating the group (%s) in the file (%s). Aborting.\n",
				group_name, hdf5_filename);
		exit(-1);
	}

	hsize_t dims[1];
	dims[0] = len;
	status = H5LTmake_dataset_double ( group_id, array_name, 1, dims, array );
	if (status < 0) {
		fprintf(stderr, "Error saving the dataset (//%s//%s) to the hdf5 file (%s). Aborting.\n",
				group_name, array_name, hdf5_filename);
		exit(-1);
	}

	H5Gclose(group_id);
	H5Fclose(file_id);
}
