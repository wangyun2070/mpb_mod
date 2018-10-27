/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to disable sanity checks in code */
/* #undef CHECK_DISABLE */

/* Define to turn on debugging checks */
/* #undef DEBUG */

/* Define to use debugging malloc/free */
/* #undef DEBUG_MALLOC */

/* Define when using the profiler tool */
/* #undef ENABLE_PROF */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define if you have a BLAS library. */
#define HAVE_BLAS 1

/* Define if fenv.h declares this. */
#define HAVE_DECL_FEENABLEEXCEPT 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the `feenableexcept' function. */
#define HAVE_FEENABLEEXCEPT 1

/* Define to 1 if you have the `getopt' function. */
#define HAVE_GETOPT 1

/* Define to 1 if you have the <getopt.h> header file. */
#define HAVE_GETOPT_H 1

/* Define to 1 if you have the <guile/gh.h> header file. */
/* #undef HAVE_GUILE_GH_H */

/* Define to 1 if you have the `H5Pset_fapl_mpio' function. */
/* #undef HAVE_H5PSET_FAPL_MPIO */

/* Define to 1 if you have the `H5Pset_mpi' function. */
/* #undef HAVE_H5PSET_MPI */

/* Define if we have & link HDF5 */
/* #undef HAVE_HDF5 */

/* Define to 1 if you have the <hdf5.h> header file. */
/* #undef HAVE_HDF5_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if you have LAPACK library. */
#define HAVE_LAPACK 1

/* Define to 1 if you have the `ctl' library (-lctl). */
#define HAVE_LIBCTL 1

/* If we have the libctl_quiet variable */
#define HAVE_LIBCTL_QUIET 1

/* Define to 1 if you have the `dfftw' library (-ldfftw). */
/* #undef HAVE_LIBDFFTW */

/* Define to 1 if you have the `dfftw_mpi' library (-ldfftw_mpi). */
/* #undef HAVE_LIBDFFTW_MPI */

/* Define to 1 if you have the `dl' library (-ldl). */
/* #undef HAVE_LIBDL */

/* Define to 1 if you have the `drfftw' library (-ldrfftw). */
/* #undef HAVE_LIBDRFFTW */

/* Define to 1 if you have the `drfftw_mpi' library (-ldrfftw_mpi). */
/* #undef HAVE_LIBDRFFTW_MPI */

/* Define to 1 if you have the `efence' library (-lefence). */
/* #undef HAVE_LIBEFENCE */

/* Define to 1 if you have the `fftw' library (-lfftw). */
/* #undef HAVE_LIBFFTW */

/* Define to 1 if you have the `fftw3' library (-lfftw3). */
#define HAVE_LIBFFTW3 1

/* Define to 1 if you have the `fftw3f' library (-lfftw3f). */
/* #undef HAVE_LIBFFTW3F */

/* Define to 1 if you have the `fftw3f_mpi' library (-lfftw3f_mpi). */
/* #undef HAVE_LIBFFTW3F_MPI */

/* Define to 1 if you have the `fftw3f_omp' library (-lfftw3f_omp). */
/* #undef HAVE_LIBFFTW3F_OMP */

/* Define to 1 if you have the `fftw3l' library (-lfftw3l). */
/* #undef HAVE_LIBFFTW3L */

/* Define to 1 if you have the `fftw3l_mpi' library (-lfftw3l_mpi). */
/* #undef HAVE_LIBFFTW3L_MPI */

/* Define to 1 if you have the `fftw3l_omp' library (-lfftw3l_omp). */
/* #undef HAVE_LIBFFTW3L_OMP */

/* Define to 1 if you have the `fftw3_mpi' library (-lfftw3_mpi). */
/* #undef HAVE_LIBFFTW3_MPI */

/* Define to 1 if you have the `fftw3_omp' library (-lfftw3_omp). */
/* #undef HAVE_LIBFFTW3_OMP */

/* Define to 1 if you have the `fftw_mpi' library (-lfftw_mpi). */
/* #undef HAVE_LIBFFTW_MPI */

/* Define to 1 if you have the `guile' library (-lguile). */
/* #undef HAVE_LIBGUILE */

/* Define to 1 if you have the <libguile.h> header file. */
#define HAVE_LIBGUILE_H 1

/* Define to 1 if you have the `guile-ltdl' library (-lguile-ltdl). */
/* #undef HAVE_LIBGUILE_LTDL */

/* Define to 1 if you have the `ltdl' library (-lltdl). */
/* #undef HAVE_LIBLTDL */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `readline' library (-lreadline). */
/* #undef HAVE_LIBREADLINE */

/* Define to 1 if you have the `rfftw' library (-lrfftw). */
/* #undef HAVE_LIBRFFTW */

/* Define to 1 if you have the `rfftw_mpi' library (-lrfftw_mpi). */
/* #undef HAVE_LIBRFFTW_MPI */

/* Define to 1 if you have the `sfftw' library (-lsfftw). */
/* #undef HAVE_LIBSFFTW */

/* Define to 1 if you have the `sfftw_mpi' library (-lsfftw_mpi). */
/* #undef HAVE_LIBSFFTW_MPI */

/* Define to 1 if you have the `srfftw' library (-lsrfftw). */
/* #undef HAVE_LIBSRFFTW */

/* Define to 1 if you have the `srfftw_mpi' library (-lsrfftw_mpi). */
/* #undef HAVE_LIBSRFFTW_MPI */

/* Define to 1 if you have the `xfftw' library (-lxfftw). */
/* #undef HAVE_LIBXFFTW */

/* Define to 1 if you have the `xfftw_mpi' library (-lxfftw_mpi). */
/* #undef HAVE_LIBXFFTW_MPI */

/* Define to 1 if you have the `xrfftw' library (-lxrfftw). */
/* #undef HAVE_LIBXRFFTW */

/* Define to 1 if you have the `xrfftw_mpi' library (-lxrfftw_mpi). */
/* #undef HAVE_LIBXRFFTW_MPI */

/* Define to 1 if you have the `z' library (-lz). */
/* #undef HAVE_LIBZ */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define if you have & link an MPI library. */
/* #undef HAVE_MPI */

/* Have nlopt lib */
/* #undef HAVE_NLOPT */

/* Define to 1 if you have the <nlopt.h> header file. */
/* #undef HAVE_NLOPT_H */

/* Define if OpenMP is enabled */
/* #undef HAVE_OPENMP */

/* Define to 1 if you have the `scm_array_get_handle' function. */
#define HAVE_SCM_ARRAY_GET_HANDLE 1

/* Define to 1 if you have the `scm_is_array' function. */
#define HAVE_SCM_IS_ARRAY 1

/* Define to 1 if you have the `scm_make_smob_type' function. */
#define HAVE_SCM_MAKE_SMOB_TYPE 1

/* define if we have SCM_NEWSMOB */
#define HAVE_SCM_NEWSMOB 1

/* define if we have SCM_SMOB_DATA */
#define HAVE_SCM_SMOB_DATA 1

/* define if we have SCM_SMOB_PREDICATE */
#define HAVE_SCM_SMOB_PREDICATE 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strncmp' function. */
#define HAVE_STRNCMP 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to enable Kottke anisotropic smoothing */
#define KOTTKE 1

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* major version */
#define MPB_VERSION_MAJOR 1

/* minor version */
#define MPB_VERSION_MINOR 6

/* patch version */
#define MPB_VERSION_PATCH 0

/* Define if calling Fortran functions directly doesn't work. */
/* #undef NO_FORTRAN_FUNCTIONS */

/* Name of package */
#define PACKAGE "mpb"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "stevenj@alum.mit.edu"

/* Define to the full name of this package. */
#define PACKAGE_NAME "mpb"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "mpb 1.6.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "mpb"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.6.1"

/* Define to use complex fields and not to require inversion symmetry */
#define SCALAR_COMPLEX 1

/* Define to use long-double precision */
/* #undef SCALAR_LONG_DOUBLE_PREC */

/* Define to use single precision */
/* #undef SCALAR_SINGLE_PREC */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to use OpenMP threading. */
/* #undef USE_OPENMP */

/* Version number of package */
#define VERSION "1.6.1"

/* Define to support Hermitian/complex dielectric tensors. */
/* #undef WITH_HERMITIAN_EPSILON */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif
