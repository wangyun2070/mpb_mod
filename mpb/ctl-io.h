/* THIS FILE WAS AUTOMATICALLY GENERATED.  DO NOT MODIFY! */
/* generated from the file: mpb.scm */

#ifndef CTL_IO_H
#define CTL_IO_H

#include <ctl.h>

#ifdef __cplusplus
extern "C" {
#endif                          /* __cplusplus */

/******* Type declarations *******/

typedef struct material_type_struct {
enum { MATERIAL_TYPE_SELF, MATERIAL_GRID, MATERIAL_FUNCTION, MEDIUM_ANISOTROPIC, MEDIUM } which_subclass;
union {
struct material_grid_struct *material_grid_data;
struct material_function_struct *material_function_data;
struct medium_anisotropic_struct *medium_anisotropic_data;
struct medium_struct *medium_data;
} subclass;
} material_type;
#define MATERIAL_TYPE_ABSTRACT 1

typedef struct medium_struct {
number epsilon;
number mu;
} medium;

typedef struct medium_anisotropic_struct {
vector3 epsilon_diag;
cvector3 epsilon_offdiag;
vector3 epsilon_offdiag_imag;
vector3 mu_diag;
cvector3 mu_offdiag;
vector3 mu_offdiag_imag;
} medium_anisotropic;

typedef struct material_function_struct {
function material_func;
} material_function;

typedef struct material_grid_struct {
integer material_grid_kind;
number epsilon_min;
number epsilon_max;
number mu_min;
number mu_max;
vector3 size;
function matgrid_init;
SCM matgrid;
} material_grid;

typedef struct geometric_object_struct {
material_type material;
vector3 center;
enum { GEOMETRIC_OBJECT_SELF, PRISM, BLOCK, SPHERE, CYLINDER, COMPOUND_GEOMETRIC_OBJECT } which_subclass;
union {
struct prism_struct *prism_data;
struct block_struct *block_data;
struct sphere_struct *sphere_data;
struct cylinder_struct *cylinder_data;
struct compound_geometric_object_struct *compound_geometric_object_data;
} subclass;
} geometric_object;

typedef struct {
int num_items;
geometric_object *items;
} geometric_object_list;

typedef struct compound_geometric_object_struct {
geometric_object_list component_objects;
} compound_geometric_object;

typedef struct cylinder_struct {
vector3 axis;
number radius;
number height;
enum { CYLINDER_SELF, WEDGE, CONE } which_subclass;
union {
struct wedge_struct *wedge_data;
struct cone_struct *cone_data;
} subclass;
} cylinder;

typedef struct cone_struct {
number radius2;
} cone;

typedef struct wedge_struct {
number wedge_angle;
vector3 wedge_start;
vector3 e1;
vector3 e2;
} wedge;

typedef struct sphere_struct {
number radius;
} sphere;

typedef struct block_struct {
vector3 e1;
vector3 e2;
vector3 e3;
vector3 size;
matrix3x3 projection_matrix;
enum { BLOCK_SELF, ELLIPSOID } which_subclass;
union {
struct ellipsoid_struct *ellipsoid_data;
} subclass;
} block;

typedef struct {
int num_items;
vector3 *items;
} vector3_list;

typedef struct {
int num_items;
number *items;
} number_list;

typedef struct prism_struct {
vector3_list vertices;
vector3 centroid;
number height;
number_list workspace;
matrix3x3 m_c2p;
matrix3x3 m_p2c;
} prism;

typedef struct ellipsoid_struct {
vector3 inverse_semi_axes;
} ellipsoid;

typedef struct lattice_struct {
vector3 basis1;
vector3 basis2;
vector3 basis3;
vector3 size;
vector3 basis_size;
vector3 b1;
vector3 b2;
vector3 b3;
matrix3x3 basis;
matrix3x3 metric;
} lattice;

typedef struct {
int num_items;
cnumber *items;
} cnumber_list;

typedef struct {
int num_items;
SCM *items;
} SCM_list;

/******* Input variables *******/
extern integer dimensions;
extern material_type default_material;
extern vector3 geometry_center;
extern lattice geometry_lattice;
extern geometric_object_list geometry;
extern boolean ensure_periodicity;
extern vector3_list k_points;
extern integer num_bands;
extern number tolerance;
extern number target_freq;
extern integer mesh_size;
extern char* epsilon_input_file;
extern char* mu_input_file;
extern boolean force_mup;
extern boolean deterministicp;
extern boolean simple_preconditionerp;
extern integer eigensolver_flags;
extern integer eigensolver_block_size;
extern integer eigensolver_nwork;
extern boolean eigensolver_davidsonp;
extern number eigensolver_flops;
extern boolean negative_epsilon_okp;

/******* Output variables *******/
extern number eigensolver_flops;
extern number_list freqs;
extern integer iterations;
extern char* parity;

extern int num_read_input_vars;
extern int num_write_output_vars;

extern SCM read_input_vars(void);
extern SCM write_output_vars(void);
extern SCM destroy_input_vars(void);
extern SCM destroy_output_vars(void);

/******* external-functions *******/

extern vector3 compute_1_group_velocity_reciprocal(integer);
extern SCM compute_1_group_velocity_reciprocal_aux(SCM arg_scm_0);

extern vector3 compute_1_group_velocity(integer);
extern SCM compute_1_group_velocity_aux(SCM arg_scm_0);

extern number compute_1_group_velocity_component(vector3, integer);
extern SCM compute_1_group_velocity_component_aux(SCM arg_scm_0, SCM arg_scm_1);

extern number_list compute_group_velocity_component(vector3);
extern SCM compute_group_velocity_component_aux(SCM arg_scm_0);

extern number_list compute_yparities(void);
extern SCM compute_yparities_aux(void);

extern number_list compute_zparities(void);
extern SCM compute_zparities_aux(void);

extern number material_grids_min_tetm_gap(vector3, integer, number, number, integer, number);
extern SCM material_grids_min_tetm_gap_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2, SCM arg_scm_3, SCM arg_scm_4, SCM arg_scm_5);

extern number material_grids_mingap(vector3_list, integer, integer, number, number, integer, number);
extern SCM material_grids_mingap_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2, SCM arg_scm_3, SCM arg_scm_4, SCM arg_scm_5, SCM arg_scm_6);

extern number material_grids_maxgap(vector3_list, integer, integer, number, number, integer, number);
extern SCM material_grids_maxgap_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2, SCM arg_scm_3, SCM arg_scm_4, SCM arg_scm_5, SCM arg_scm_6);

extern number material_grids_approx_gradient(vector3, integer, integer, number);
extern SCM material_grids_approx_gradient_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2, SCM arg_scm_3);

extern void print_material_grids_deps_du_numeric(number);
extern SCM print_material_grids_deps_du_numeric_aux(SCM arg_scm_0);

extern void print_material_grids_deps_du(void);
extern SCM print_material_grids_deps_du_aux(void);

extern void print_material_grids_gradient(integer);
extern SCM print_material_grids_gradient_aux(SCM arg_scm_0);

extern void material_grids_match_epsilon_fileB(char*, number);
extern SCM material_grids_match_epsilon_fileB_aux(SCM arg_scm_0, SCM arg_scm_1);

extern void load_material_gridB(material_grid, char*, vector3);
extern SCM load_material_gridB_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2);

extern void save_material_grid(material_grid, char*);
extern SCM save_material_grid_aux(SCM arg_scm_0, SCM arg_scm_1);

extern void randomize_material_gridB(material_grid, number);
extern SCM randomize_material_gridB_aux(SCM arg_scm_0, SCM arg_scm_1);

extern cvector3 cvector_field_get_point_bloch(SCM, vector3);
extern SCM cvector_field_get_point_bloch_aux(SCM arg_scm_0, SCM arg_scm_1);

extern cvector3 cvector_field_get_point(SCM, vector3);
extern SCM cvector_field_get_point_aux(SCM arg_scm_0, SCM arg_scm_1);

extern cnumber cscalar_field_get_point(SCM, vector3);
extern SCM cscalar_field_get_point_aux(SCM arg_scm_0, SCM arg_scm_1);

extern number rscalar_field_get_point(SCM, vector3);
extern SCM rscalar_field_get_point_aux(SCM arg_scm_0, SCM arg_scm_1);

extern cnumber integrate_fieldL(function, SCM_list);
extern SCM integrate_fieldL_aux(SCM arg_scm_0, SCM arg_scm_1);

extern void field_mapLB(SCM, function, SCM_list);
extern SCM field_mapLB_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2);

extern void field_load(SCM);
extern SCM field_load_aux(SCM arg_scm_0);

extern void field_setB(SCM, SCM);
extern SCM field_setB_aux(SCM arg_scm_0, SCM arg_scm_1);

extern boolean fields_conformp(SCM, SCM);
extern SCM fields_conformp_aux(SCM arg_scm_0, SCM arg_scm_1);

extern SCM field_make(SCM);
extern SCM field_make_aux(SCM arg_scm_0);

extern void cvector_field_nonblochB(SCM);
extern SCM cvector_field_nonblochB_aux(SCM arg_scm_0);

extern SCM cvector_field_make(SCM);
extern SCM cvector_field_make_aux(SCM arg_scm_0);

extern SCM rscalar_field_make(SCM);
extern SCM rscalar_field_make_aux(SCM arg_scm_0);

extern boolean cur_fieldp(SCM);
extern SCM cur_fieldp_aux(SCM arg_scm_0);

extern void load_eigenvectors(char*);
extern SCM load_eigenvectors_aux(SCM arg_scm_0);

extern void save_eigenvectors(char*);
extern SCM save_eigenvectors_aux(SCM arg_scm_0);

extern SCM input_eigenvectors(char*, integer);
extern SCM input_eigenvectors_aux(SCM arg_scm_0, SCM arg_scm_1);

extern void output_eigenvectors(SCM, char*);
extern SCM output_eigenvectors_aux(SCM arg_scm_0, SCM arg_scm_1);

extern void scale_eigenvector(integer, cnumber);
extern SCM scale_eigenvector_aux(SCM arg_scm_0, SCM arg_scm_1);

extern SCM dot_eigenvectors(SCM, integer);
extern SCM dot_eigenvectors_aux(SCM arg_scm_0, SCM arg_scm_1);

extern void set_eigenvectors(SCM, integer);
extern SCM set_eigenvectors_aux(SCM arg_scm_0, SCM arg_scm_1);

extern SCM get_eigenvectors(integer, integer);
extern SCM get_eigenvectors_aux(SCM arg_scm_0, SCM arg_scm_1);

extern cnumber_list sqmatrix_eigvals(SCM);
extern SCM sqmatrix_eigvals_aux(SCM arg_scm_0);

extern SCM sqmatrix_diagm(cnumber_list);
extern SCM sqmatrix_diagm_aux(SCM arg_scm_0);

extern SCM sqmatrix_mult(SCM, SCM);
extern SCM sqmatrix_mult_aux(SCM arg_scm_0, SCM arg_scm_1);

extern cnumber sqmatrix_ref(SCM, integer, integer);
extern SCM sqmatrix_ref_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2);

extern integer sqmatrix_size(SCM);
extern SCM sqmatrix_size_aux(SCM arg_scm_0);

extern void set_kpoint_index(integer);
extern SCM set_kpoint_index_aux(SCM arg_scm_0);

extern integer get_kpoint_index(void);
extern SCM get_kpoint_index_aux(void);

extern boolean has_inversion_symp(void);
extern SCM has_inversion_symp_aux(void);

extern boolean has_hermitian_epsp(void);
extern SCM has_hermitian_epsp_aux(void);

extern number mpi_max(number);
extern SCM mpi_max_aux(SCM arg_scm_0);

extern integer mpi_proc_index(void);
extern SCM mpi_proc_index_aux(void);

extern integer mpi_num_procs(void);
extern SCM mpi_num_procs_aux(void);

extern boolean using_mpip(void);
extern SCM using_mpip_aux(void);

extern boolean mpi_is_masterp(void);
extern SCM mpi_is_masterp_aux(void);

extern void output_field_to_file(integer, char*);
extern SCM output_field_to_file_aux(SCM arg_scm_0, SCM arg_scm_1);

extern number compute_energy_in_object_list(geometric_object_list);
extern SCM compute_energy_in_object_list_aux(SCM arg_scm_0);

extern number compute_energy_integral(function);
extern SCM compute_energy_integral_aux(SCM arg_scm_0);

extern cnumber compute_field_integral(function);
extern SCM compute_field_integral_aux(SCM arg_scm_0);

extern number compute_energy_in_dielectric(number, number);
extern SCM compute_energy_in_dielectric_aux(SCM arg_scm_0, SCM arg_scm_1);

extern cnumber get_cscalar_point(vector3);
extern SCM get_cscalar_point_aux(SCM arg_scm_0);

extern cnumber get_bloch_cscalar_point(vector3);
extern SCM get_bloch_cscalar_point_aux(SCM arg_scm_0);

extern cvector3 get_field_point(vector3);
extern SCM get_field_point_aux(SCM arg_scm_0);

extern cvector3 get_bloch_field_point(vector3);
extern SCM get_bloch_field_point_aux(SCM arg_scm_0);

extern number get_energy_point(vector3);
extern SCM get_energy_point_aux(SCM arg_scm_0);

extern cmatrix3x3 get_epsilon_inverse_tensor_point(vector3);
extern SCM get_epsilon_inverse_tensor_point_aux(SCM arg_scm_0);

extern number get_epsilon_point(vector3);
extern SCM get_epsilon_point_aux(SCM arg_scm_0);

extern void compute_field_divergence(void);
extern SCM compute_field_divergence_aux(void);

extern number_list compute_field_energy(void);
extern SCM compute_field_energy_aux(void);

extern void fix_field_phase(void);
extern SCM fix_field_phase_aux(void);

extern void get_mu(void);
extern SCM get_mu_aux(void);

extern void get_epsilon(void);
extern SCM get_epsilon_aux(void);

extern void get_efield_from_dfield(void);
extern SCM get_efield_from_dfield_aux(void);

extern void get_bfield(integer);
extern SCM get_bfield_aux(SCM arg_scm_0);

extern void get_hfield(integer);
extern SCM get_hfield_aux(SCM arg_scm_0);

extern void get_dfield(integer);
extern SCM get_dfield_aux(SCM arg_scm_0);

extern void solve_kpoint(vector3);
extern SCM solve_kpoint_aux(SCM arg_scm_0);

extern void randomize_fields(void);
extern SCM randomize_fields_aux(void);

extern void set_parity(integer);
extern SCM set_parity_aux(SCM arg_scm_0);

extern boolean using_mup(void);
extern SCM using_mup_aux(void);

extern void init_params(integer, boolean);
extern SCM init_params_aux(SCM arg_scm_0, SCM arg_scm_1);

extern matrix3x3 square_basis(matrix3x3, vector3);
extern SCM square_basis_aux(SCM arg_scm_0, SCM arg_scm_1);

extern number range_overlap_with_object(vector3, vector3, geometric_object, number, integer);
extern SCM range_overlap_with_object_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2, SCM arg_scm_3, SCM arg_scm_4);

extern void display_geometric_object_info(integer, geometric_object);
extern SCM display_geometric_object_info_aux(SCM arg_scm_0, SCM arg_scm_1);

extern boolean point_in_periodic_objectp(vector3, geometric_object);
extern SCM point_in_periodic_objectp_aux(SCM arg_scm_0, SCM arg_scm_1);

extern vector3 normal_to_object(vector3, geometric_object);
extern SCM normal_to_object_aux(SCM arg_scm_0, SCM arg_scm_1);

extern boolean point_in_objectp(vector3, geometric_object);
extern SCM point_in_objectp_aux(SCM arg_scm_0, SCM arg_scm_1);


extern void export_external_functions(void);

/******* class input function prototypes *******/

extern void lattice_input(SCM so, lattice *o);
extern void ellipsoid_input(SCM so, ellipsoid *o);
extern void prism_input(SCM so, prism *o);
extern void block_input(SCM so, block *o);
extern void sphere_input(SCM so, sphere *o);
extern void wedge_input(SCM so, wedge *o);
extern void cone_input(SCM so, cone *o);
extern void cylinder_input(SCM so, cylinder *o);
extern void compound_geometric_object_input(SCM so, compound_geometric_object *o);
extern void geometric_object_input(SCM so, geometric_object *o);
extern void material_grid_input(SCM so, material_grid *o);
extern void material_function_input(SCM so, material_function *o);
extern void medium_anisotropic_input(SCM so, medium_anisotropic *o);
extern void medium_input(SCM so, medium *o);
extern void material_type_input(SCM so, material_type *o);

/******* class copy function prototypes *******/

extern void lattice_copy(const lattice *o0,lattice *o);
extern void ellipsoid_copy(const ellipsoid *o0,ellipsoid *o);
extern void prism_copy(const prism *o0,prism *o);
extern void block_copy(const block *o0,block *o);
extern void sphere_copy(const sphere *o0,sphere *o);
extern void wedge_copy(const wedge *o0,wedge *o);
extern void cone_copy(const cone *o0,cone *o);
extern void cylinder_copy(const cylinder *o0,cylinder *o);
extern void compound_geometric_object_copy(const compound_geometric_object *o0,compound_geometric_object *o);
extern void geometric_object_copy(const geometric_object *o0,geometric_object *o);
extern void material_grid_copy(const material_grid *o0,material_grid *o);
extern void material_function_copy(const material_function *o0,material_function *o);
extern void medium_anisotropic_copy(const medium_anisotropic *o0,medium_anisotropic *o);
extern void medium_copy(const medium *o0,medium *o);
extern void material_type_copy(const material_type *o0,material_type *o);

/******* class equal function prototypes *******/

extern boolean lattice_equal(const lattice *o0, const lattice *o);
extern boolean ellipsoid_equal(const ellipsoid *o0, const ellipsoid *o);
extern boolean prism_equal(const prism *o0, const prism *o);
extern boolean block_equal(const block *o0, const block *o);
extern boolean sphere_equal(const sphere *o0, const sphere *o);
extern boolean wedge_equal(const wedge *o0, const wedge *o);
extern boolean cone_equal(const cone *o0, const cone *o);
extern boolean cylinder_equal(const cylinder *o0, const cylinder *o);
extern boolean compound_geometric_object_equal(const compound_geometric_object *o0, const compound_geometric_object *o);
extern boolean geometric_object_equal(const geometric_object *o0, const geometric_object *o);
extern boolean material_grid_equal(const material_grid *o0, const material_grid *o);
extern boolean material_function_equal(const material_function *o0, const material_function *o);
extern boolean medium_anisotropic_equal(const medium_anisotropic *o0, const medium_anisotropic *o);
extern boolean medium_equal(const medium *o0, const medium *o);
extern boolean material_type_equal(const material_type *o0, const material_type *o);

/******* class destruction function prototypes *******/

extern void lattice_destroy(lattice o);
extern void ellipsoid_destroy(ellipsoid o);
extern void prism_destroy(prism o);
extern void block_destroy(block o);
extern void sphere_destroy(sphere o);
extern void wedge_destroy(wedge o);
extern void cone_destroy(cone o);
extern void cylinder_destroy(cylinder o);
extern void compound_geometric_object_destroy(compound_geometric_object o);
extern void geometric_object_destroy(geometric_object o);
extern void material_grid_destroy(material_grid o);
extern void material_function_destroy(material_function o);
extern void medium_anisotropic_destroy(medium_anisotropic o);
extern void medium_destroy(medium o);
extern void material_type_destroy(material_type o);


#ifdef __cplusplus
}                               /* extern "C" */
#endif                          /* __cplusplus */

#endif                          /* CTL_IO_H */

