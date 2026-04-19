# cython: language_level=3
"""
Cython declarations for nifti_image and related types from
3rdparty/nifti/nifti1_io.h.
"""


cdef extern from "nifti1_io.h":
    ctypedef struct mat44:
        float m[4][4]

    ctypedef struct nifti1_extension:
        int    esize
        int    ecode
        char  *edata

    ctypedef int analyze_75_orient_code

    ctypedef struct nifti_image:
        int   ndim
        int   nx, ny, nz, nt, nu, nv, nw
        int   dim[8]
        size_t nvox
        int   nbyper
        int   datatype

        float dx, dy, dz, dt, du, dv, dw
        float pixdim[8]

        float scl_slope
        float scl_inter

        float cal_min
        float cal_max

        int   qform_code
        int   sform_code

        int   freq_dim
        int   phase_dim
        int   slice_dim

        int   slice_code
        int   slice_start
        int   slice_end
        float slice_duration

        float quatern_b, quatern_c, quatern_d
        float qoffset_x, qoffset_y, qoffset_z
        float qfac

        mat44 qto_xyz
        mat44 qto_ijk
        mat44 sto_xyz
        mat44 sto_ijk

        float toffset
        int   xyz_units
        int   time_units

        int   nifti_type
        int   intent_code
        float intent_p1, intent_p2, intent_p3
        char  intent_name[16]

        char  descrip[80]
        char  aux_file[24]

        char *fname
        char *iname
        int   iname_offset
        int   swapsize
        int   byteorder
        void *data

        int                 num_ext
        nifti1_extension   *ext_list
        analyze_75_orient_code analyze75_orient

    nifti_image *nifti_image_read(const char *hname, int read_data)
    void         nifti_image_free(nifti_image *nim)
    void         nifti_image_write(nifti_image *nim)
    nifti_image *nifti_simple_init_nim()
    nifti_image *nifti_copy_nim_info(const nifti_image *src)
    int          nifti_set_filenames(nifti_image *nim, const char *prefix,
                                     int check, int set_byte_order)

    # NIfTI datatype codes
    int DT_FLOAT32
    int NIFTI_TYPE_FLOAT32
