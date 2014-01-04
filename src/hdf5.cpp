
#include "hdf5.hpp"


hid_t H5Aopen_checked(hid_t obj_id, const char* attr_name, hid_t aapl_id)
{
    hid_t attr = H5Aopen(obj_id, attr_name, aapl_id);
    if (attr < 0) {
        Logger::abort("Failed to open HDF5 attribute \"%s\".", attr_name);
    }

    return attr;
}


hid_t H5Dopen2_checked(hid_t loc_id, const char* name, hid_t dapl_id)
{
    hid_t dataset = H5Dopen2(loc_id, name, dapl_id);
    if (dataset < 0) {
        Logger::abort("Failed to open the HDF5 dataset \"%s\".", name);
    }

    return dataset;
}


void H5Dread_checked(hid_t dataset_id, hid_t mem_type_id,
                            hid_t mem_space_id, hid_t file_space_id,
                            hid_t xfer_plist_id, void* buf)
{
    hid_t err = H5Dread(dataset_id, mem_type_id, mem_space_id, file_space_id,
                        xfer_plist_id, buf);
    if (err < 0) {
        Logger::abort("HDF5 read failed.");
    }
}


void H5Dwrite_checked(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id,
	                  hid_t file_space_id, hid_t xfer_plist_id,
	                  const void* buf)
{
	herr_t err = H5Dwrite(dataset_id, mem_type_id, mem_space_id,
		                  file_space_id, xfer_plist_id, buf);
	if (err < 0) {
		Logger::abort("HDF5 write failed.");
	}
}


hid_t H5Dcreate2_checked(hid_t loc_id, const char* name, hid_t dtype_id,
	                     hid_t space_id, hid_t lcpl_id, hid_t dcpl_id,
	                     hid_t dapl_id)
{
	hid_t dataset = H5Dcreate2(loc_id, name, dtype_id, space_id, lcpl_id,
		                       dcpl_id, dapl_id);
	if (dataset < 0) {
		Logger::abort("Failed to create the HDF5 dataset \"%s\".", name);
	}

	return dataset;
}


// dataspace
void H5Sselect_hyperslab_checked(hid_t space_id, H5S_seloper_t op,
	                             const hsize_t* start, const hsize_t* stride,
	                             const hsize_t* count, const hsize_t* block)
{
	herr_t err = H5Sselect_hyperslab(space_id, op, start, stride, count, block);
	if (err < 0) {
		Logger::abort("HDF5 hyperslab selection failed.");
	}
}


hid_t H5Screate_simple_checked(int rank, const hsize_t* current_dims,
	                           const hsize_t* maximum_dims)
{
	hid_t dataspace = H5Screate_simple(rank, current_dims, maximum_dims);
	if (dataspace < 0) {
		Logger::abort("HDF5 dataspace creation failed.");
	}

	return dataspace;
}

