#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "constants.h"
#include "hdf5.h"
#include "emfields.h"
#include "wlcs.h"

/******************************************************************************
 * Calculate electromagnetic fields for test case.
 * 
 * Input:
 * x, y, z, t: spatial positions and time.
 * Output:
 * emf: electromagnetic fields at (x,y,z,t) 
 ******************************************************************************/
void getemf_test(double x, double y, double z, double t, struct emfields *emf)
{
    emf->Bx = 0.0;
    emf->By = 0.0;
    emf->Bz = 1.0E0;
    emf->Ex = 0.0;
    emf->Ey = 0.0;
    emf->Ez = 0.0;
}

/******************************************************************************
 * Calculate or interpolate to get electromagnetic fields.
 * system_type = 0: test case.
 * system_type = 1: wire-loop current system case.
 * system_type = 2: force-free magnetic field case.
 * system_type = 3: MHD + test particle.
 *
 * Input:
 *  x, y, z, t: spatial positions and time.
 *  system_type: the system type.
 * Output:
 *  emf: electromagnetic fields at (x,y,z,t) 
 ******************************************************************************/
void get_emf(double x, double y, double z, double t, int system_type,
        struct emfields *emf)
{
    if (system_type == 0) {
        getemf_test(x, y, z, t, emf);
    }
    else if (system_type == 1) {
        getemf_wlcs(x, y, z, t, emf);
    }
    /* else if (system_type == 2){ */
    /*     getemf_ff(x, y, z, t, emf); */
    /* } */
    //printf("Bx, By, Bz, Ex, Ey, Ez: %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n",
    //        emf->Bx, emf->By, emf->Bz, emf->Ex, emf->Ey, emf->Ez);
}

/******************************************************************************
 * Create a compound data type containing the fields information for HDF5.
 *
 * Output:
 *  memtype: compound datatype for memory.
 *  filetype: compound datatype for the file.
 ******************************************************************************/
void create_fields_ctype(hid_t *memtype, hid_t *filetype)
{
    //herr_t status;
    /* Create the compound datatype in memory */
    *memtype = H5Tcreate(H5T_COMPOUND, sizeof(emfields));
    H5Tinsert(*memtype, "Bx", 
            HOFFSET(emfields, Bx), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "By", 
            HOFFSET(emfields, By), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "Bz", 
            HOFFSET(emfields, Bz), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "Ex", 
            HOFFSET(emfields, Ex), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "Ey", 
            HOFFSET(emfields, Ey), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "Ez", 
            HOFFSET(emfields, Ez), H5T_NATIVE_DOUBLE);
    
    /* Create the compound datatype for the file.  Because the standard */
    /* types we are using for the file may have different sizes than */
    /* the corresponding native types, we must manually calculate the */
    /* offset of each member. */
    
    *filetype = H5Tcreate(H5T_COMPOUND, 8*6);
    H5Tinsert (*filetype, "Bx", 0, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "By", 8, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "Bz", 8*2, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "Ex", 8*3, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "Ey", 8*4, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "Ez", 8*5, H5T_IEEE_F64BE);
}
