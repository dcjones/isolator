
#ifndef ISOLATOR_HDF5_HPP
#define ISOLATOR_HDF5_HPP

/* Define a bunch of wrappers for HDF5 functions that check their output for
 * errors. */

#include <hdf5.h>
#include <hdf5_hl.h>
#include <logger.hpp>

// attributes
hid_t H5Aopen_checked(hid_t obj_id, const char* attr_name, hid_t aapl_id);
hid_t H5Dopen2_checked(hid_t loc_id, const char* name, hid_t dapl_id);


// datasets
void H5Dread_checked(hid_t dataset_id, hid_t mem_type_id,
                     hid_t mem_space_id, hid_t file_space_id,
                     hid_t xfer_plist_id, void* buf);

void H5Dwrite_checked(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id,
	                  hid_t file_space_id, hid_t xfer_plist_id,
	                  const void* buf);

hid_t H5Dcreate2_checked(hid_t loc_id, const char* name, hid_t dtype_id,
	                     hid_t space_id, hid_t lcpl_id, hid_t dcpl_id,
	                     hid_t dapl_id);


// dataspace
hid_t H5Screate_simple_checked(int rank, const hsize_t* current_dims,
	                           const hsize_t* maximum_dims);

void H5Sselect_hyperslab_checked(hid_t space_id, H5S_seloper_t op,
	                             const hsize_t* start, const hsize_t* stride,
	                             const hsize_t* count, const hsize_t* block);


#endif

