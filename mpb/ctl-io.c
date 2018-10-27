/* THIS FILE WAS AUTOMATICALLY GENERATED.  DO NOT MODIFY! */
/* generated from the file: mpb.scm */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "ctl-io.h"

#ifdef CXX_CTL_IO
using namespace ctlio;
#endif

/******* Input variables *******/
integer dimensions;
material_type default_material;
vector3 geometry_center;
lattice geometry_lattice;
geometric_object_list geometry;
boolean ensure_periodicity;
vector3_list k_points;
integer num_bands;
number tolerance;
number target_freq;
integer mesh_size;
char* epsilon_input_file;
char* mu_input_file;
boolean force_mup;
boolean deterministicp;
boolean simple_preconditionerp;
integer eigensolver_flags;
integer eigensolver_block_size;
integer eigensolver_nwork;
boolean eigensolver_davidsonp;
number eigensolver_flops;
boolean negative_epsilon_okp;

/******* Output variables *******/
number eigensolver_flops;
number_list freqs;
integer iterations;
char* parity;

int num_read_input_vars = 0; /* # calls to read_input_vars */
int num_write_output_vars = 0; /* # calls to read_input_vars */

/******* class input functions *******/

void lattice_input(SCM so, lattice *o)
{
o->basis1 = vector3_object_property(so, "basis1");
o->basis2 = vector3_object_property(so, "basis2");
o->basis3 = vector3_object_property(so, "basis3");
o->size = vector3_object_property(so, "size");
o->basis_size = vector3_object_property(so, "basis-size");
o->b1 = vector3_object_property(so, "b1");
o->b2 = vector3_object_property(so, "b2");
o->b3 = vector3_object_property(so, "b3");
o->basis = matrix3x3_object_property(so, "basis");
o->metric = matrix3x3_object_property(so, "metric");
}

void ellipsoid_input(SCM so, ellipsoid *o)
{
o->inverse_semi_axes = vector3_object_property(so, "inverse-semi-axes");
}

void prism_input(SCM so, prism *o)
{
{
list lo_t = list_object_property(so, "vertices");
int i_t;
o->vertices.num_items = list_length(lo_t);
o->vertices.items = ((vector3 *) malloc(sizeof(vector3) * (o->vertices.num_items)));
for (i_t = 0; i_t < o->vertices.num_items; i_t++) {
o->vertices.items[i_t] = vector3_list_ref(lo_t, i_t);
}
}
o->centroid = vector3_object_property(so, "centroid");
o->height = number_object_property(so, "height");
{
list lo_t = list_object_property(so, "workspace");
int i_t;
o->workspace.num_items = list_length(lo_t);
o->workspace.items = ((number *) malloc(sizeof(number) * (o->workspace.num_items)));
for (i_t = 0; i_t < o->workspace.num_items; i_t++) {
o->workspace.items[i_t] = number_list_ref(lo_t, i_t);
}
}
o->m_c2p = matrix3x3_object_property(so, "m_c2p");
o->m_p2c = matrix3x3_object_property(so, "m_p2c");
}

void block_input(SCM so, block *o)
{
o->e1 = vector3_object_property(so, "e1");
o->e2 = vector3_object_property(so, "e2");
o->e3 = vector3_object_property(so, "e3");
o->size = vector3_object_property(so, "size");
o->projection_matrix = matrix3x3_object_property(so, "projection-matrix");
if (object_is_member("ellipsoid", so)) {
o->which_subclass = ELLIPSOID;
o->subclass.ellipsoid_data = ((ellipsoid *) malloc(sizeof(ellipsoid)));
ellipsoid_input(so, o->subclass.ellipsoid_data);
}
else 
o->which_subclass = BLOCK_SELF;
}

void sphere_input(SCM so, sphere *o)
{
o->radius = number_object_property(so, "radius");
}

void wedge_input(SCM so, wedge *o)
{
o->wedge_angle = number_object_property(so, "wedge-angle");
o->wedge_start = vector3_object_property(so, "wedge-start");
o->e1 = vector3_object_property(so, "e1");
o->e2 = vector3_object_property(so, "e2");
}

void cone_input(SCM so, cone *o)
{
o->radius2 = number_object_property(so, "radius2");
}

void cylinder_input(SCM so, cylinder *o)
{
o->axis = vector3_object_property(so, "axis");
o->radius = number_object_property(so, "radius");
o->height = number_object_property(so, "height");
if (object_is_member("wedge", so)) {
o->which_subclass = WEDGE;
o->subclass.wedge_data = ((wedge *) malloc(sizeof(wedge)));
wedge_input(so, o->subclass.wedge_data);
}
else if (object_is_member("cone", so)) {
o->which_subclass = CONE;
o->subclass.cone_data = ((cone *) malloc(sizeof(cone)));
cone_input(so, o->subclass.cone_data);
}
else 
o->which_subclass = CYLINDER_SELF;
}

void compound_geometric_object_input(SCM so, compound_geometric_object *o)
{
{
list lo_t = list_object_property(so, "component-objects");
int i_t;
o->component_objects.num_items = list_length(lo_t);
o->component_objects.items = ((geometric_object *) malloc(sizeof(geometric_object) * (o->component_objects.num_items)));
for (i_t = 0; i_t < o->component_objects.num_items; i_t++) {
geometric_object_input(object_list_ref(lo_t, i_t), &o->component_objects.items[i_t]);
}
}
}

void geometric_object_input(SCM so, geometric_object *o)
{
material_type_input(object_object_property(so, "material"), &o->material);
o->center = vector3_object_property(so, "center");
if (object_is_member("prism", so)) {
o->which_subclass = PRISM;
o->subclass.prism_data = ((prism *) malloc(sizeof(prism)));
prism_input(so, o->subclass.prism_data);
}
else if (object_is_member("block", so)) {
o->which_subclass = BLOCK;
o->subclass.block_data = ((block *) malloc(sizeof(block)));
block_input(so, o->subclass.block_data);
}
else if (object_is_member("sphere", so)) {
o->which_subclass = SPHERE;
o->subclass.sphere_data = ((sphere *) malloc(sizeof(sphere)));
sphere_input(so, o->subclass.sphere_data);
}
else if (object_is_member("cylinder", so)) {
o->which_subclass = CYLINDER;
o->subclass.cylinder_data = ((cylinder *) malloc(sizeof(cylinder)));
cylinder_input(so, o->subclass.cylinder_data);
}
else if (object_is_member("compound-geometric-object", so)) {
o->which_subclass = COMPOUND_GEOMETRIC_OBJECT;
o->subclass.compound_geometric_object_data = ((compound_geometric_object *) malloc(sizeof(compound_geometric_object)));
compound_geometric_object_input(so, o->subclass.compound_geometric_object_data);
}
else 
o->which_subclass = GEOMETRIC_OBJECT_SELF;
}

void material_grid_input(SCM so, material_grid *o)
{
o->material_grid_kind = integer_object_property(so, "material-grid-kind");
o->epsilon_min = number_object_property(so, "epsilon-min");
o->epsilon_max = number_object_property(so, "epsilon-max");
o->mu_min = number_object_property(so, "mu-min");
o->mu_max = number_object_property(so, "mu-max");
o->size = vector3_object_property(so, "size");
o->matgrid_init = function_object_property(so, "matgrid-init");
o->matgrid = SCM_object_property(so, "matgrid");
}

void material_function_input(SCM so, material_function *o)
{
o->material_func = function_object_property(so, "material-func");
}

void medium_anisotropic_input(SCM so, medium_anisotropic *o)
{
o->epsilon_diag = vector3_object_property(so, "epsilon-diag");
o->epsilon_offdiag = cvector3_object_property(so, "epsilon-offdiag");
o->epsilon_offdiag_imag = vector3_object_property(so, "epsilon-offdiag-imag");
o->mu_diag = vector3_object_property(so, "mu-diag");
o->mu_offdiag = cvector3_object_property(so, "mu-offdiag");
o->mu_offdiag_imag = vector3_object_property(so, "mu-offdiag-imag");
}

void medium_input(SCM so, medium *o)
{
o->epsilon = number_object_property(so, "epsilon");
o->mu = number_object_property(so, "mu");
}

void material_type_input(SCM so, material_type *o)
{
if (object_is_member("material-grid", so)) {
o->which_subclass = MATERIAL_GRID;
o->subclass.material_grid_data = ((material_grid *) malloc(sizeof(material_grid)));
material_grid_input(so, o->subclass.material_grid_data);
}
else if (object_is_member("material-function", so)) {
o->which_subclass = MATERIAL_FUNCTION;
o->subclass.material_function_data = ((material_function *) malloc(sizeof(material_function)));
material_function_input(so, o->subclass.material_function_data);
}
else if (object_is_member("medium-anisotropic", so)) {
o->which_subclass = MEDIUM_ANISOTROPIC;
o->subclass.medium_anisotropic_data = ((medium_anisotropic *) malloc(sizeof(medium_anisotropic)));
medium_anisotropic_input(so, o->subclass.medium_anisotropic_data);
}
else if (object_is_member("medium", so)) {
o->which_subclass = MEDIUM;
o->subclass.medium_data = ((medium *) malloc(sizeof(medium)));
medium_input(so, o->subclass.medium_data);
}
else 
o->which_subclass = MATERIAL_TYPE_SELF;
}

/******* class copy functions *******/

void lattice_copy(const lattice *o0,lattice *o)
{
o->basis1 = o0->basis1;
o->basis2 = o0->basis2;
o->basis3 = o0->basis3;
o->size = o0->size;
o->basis_size = o0->basis_size;
o->b1 = o0->b1;
o->b2 = o0->b2;
o->b3 = o0->b3;
o->basis = o0->basis;
o->metric = o0->metric;
}

void ellipsoid_copy(const ellipsoid *o0,ellipsoid *o)
{
o->inverse_semi_axes = o0->inverse_semi_axes;
}

void prism_copy(const prism *o0,prism *o)
{
{
int i_t;
o->vertices.num_items = o0->vertices.num_items;
o->vertices.items = ((vector3 *) malloc(sizeof(vector3) * (o->vertices.num_items)));
for (i_t = 0; i_t < o->vertices.num_items; i_t++) {
o->vertices.items[i_t] = o0->vertices.items[i_t];
}
}
o->centroid = o0->centroid;
o->height = o0->height;
{
int i_t;
o->workspace.num_items = o0->workspace.num_items;
o->workspace.items = ((number *) malloc(sizeof(number) * (o->workspace.num_items)));
for (i_t = 0; i_t < o->workspace.num_items; i_t++) {
o->workspace.items[i_t] = o0->workspace.items[i_t];
}
}
o->m_c2p = o0->m_c2p;
o->m_p2c = o0->m_p2c;
}

void block_copy(const block *o0,block *o)
{
o->e1 = o0->e1;
o->e2 = o0->e2;
o->e3 = o0->e3;
o->size = o0->size;
o->projection_matrix = o0->projection_matrix;
if (o0->which_subclass == ELLIPSOID) {
o->which_subclass = ELLIPSOID;
o->subclass.ellipsoid_data = ((ellipsoid *) malloc(sizeof(ellipsoid)));
ellipsoid_copy(o0->subclass.ellipsoid_data, o->subclass.ellipsoid_data);
}
else 
o->which_subclass = BLOCK_SELF;
}

void sphere_copy(const sphere *o0,sphere *o)
{
o->radius = o0->radius;
}

void wedge_copy(const wedge *o0,wedge *o)
{
o->wedge_angle = o0->wedge_angle;
o->wedge_start = o0->wedge_start;
o->e1 = o0->e1;
o->e2 = o0->e2;
}

void cone_copy(const cone *o0,cone *o)
{
o->radius2 = o0->radius2;
}

void cylinder_copy(const cylinder *o0,cylinder *o)
{
o->axis = o0->axis;
o->radius = o0->radius;
o->height = o0->height;
if (o0->which_subclass == WEDGE) {
o->which_subclass = WEDGE;
o->subclass.wedge_data = ((wedge *) malloc(sizeof(wedge)));
wedge_copy(o0->subclass.wedge_data, o->subclass.wedge_data);
}
else if (o0->which_subclass == CONE) {
o->which_subclass = CONE;
o->subclass.cone_data = ((cone *) malloc(sizeof(cone)));
cone_copy(o0->subclass.cone_data, o->subclass.cone_data);
}
else 
o->which_subclass = CYLINDER_SELF;
}

void compound_geometric_object_copy(const compound_geometric_object *o0,compound_geometric_object *o)
{
{
int i_t;
o->component_objects.num_items = o0->component_objects.num_items;
o->component_objects.items = ((geometric_object *) malloc(sizeof(geometric_object) * (o->component_objects.num_items)));
for (i_t = 0; i_t < o->component_objects.num_items; i_t++) {
geometric_object_copy(&o0->component_objects.items[i_t], &o->component_objects.items[i_t]);
}
}
}

void geometric_object_copy(const geometric_object *o0,geometric_object *o)
{
material_type_copy(&o0->material, &o->material);
o->center = o0->center;
if (o0->which_subclass == PRISM) {
o->which_subclass = PRISM;
o->subclass.prism_data = ((prism *) malloc(sizeof(prism)));
prism_copy(o0->subclass.prism_data, o->subclass.prism_data);
}
else if (o0->which_subclass == BLOCK) {
o->which_subclass = BLOCK;
o->subclass.block_data = ((block *) malloc(sizeof(block)));
block_copy(o0->subclass.block_data, o->subclass.block_data);
}
else if (o0->which_subclass == SPHERE) {
o->which_subclass = SPHERE;
o->subclass.sphere_data = ((sphere *) malloc(sizeof(sphere)));
sphere_copy(o0->subclass.sphere_data, o->subclass.sphere_data);
}
else if (o0->which_subclass == CYLINDER) {
o->which_subclass = CYLINDER;
o->subclass.cylinder_data = ((cylinder *) malloc(sizeof(cylinder)));
cylinder_copy(o0->subclass.cylinder_data, o->subclass.cylinder_data);
}
else if (o0->which_subclass == COMPOUND_GEOMETRIC_OBJECT) {
o->which_subclass = COMPOUND_GEOMETRIC_OBJECT;
o->subclass.compound_geometric_object_data = ((compound_geometric_object *) malloc(sizeof(compound_geometric_object)));
compound_geometric_object_copy(o0->subclass.compound_geometric_object_data, o->subclass.compound_geometric_object_data);
}
else 
o->which_subclass = GEOMETRIC_OBJECT_SELF;
}

void material_grid_copy(const material_grid *o0,material_grid *o)
{
o->material_grid_kind = o0->material_grid_kind;
o->epsilon_min = o0->epsilon_min;
o->epsilon_max = o0->epsilon_max;
o->mu_min = o0->mu_min;
o->mu_max = o0->mu_max;
o->size = o0->size;
o->matgrid_init = o0->matgrid_init;
o->matgrid = o0->matgrid;
}

void material_function_copy(const material_function *o0,material_function *o)
{
o->material_func = o0->material_func;
}

void medium_anisotropic_copy(const medium_anisotropic *o0,medium_anisotropic *o)
{
o->epsilon_diag = o0->epsilon_diag;
o->epsilon_offdiag = o0->epsilon_offdiag;
o->epsilon_offdiag_imag = o0->epsilon_offdiag_imag;
o->mu_diag = o0->mu_diag;
o->mu_offdiag = o0->mu_offdiag;
o->mu_offdiag_imag = o0->mu_offdiag_imag;
}

void medium_copy(const medium *o0,medium *o)
{
o->epsilon = o0->epsilon;
o->mu = o0->mu;
}

void material_type_copy(const material_type *o0,material_type *o)
{
if (o0->which_subclass == MATERIAL_GRID) {
o->which_subclass = MATERIAL_GRID;
o->subclass.material_grid_data = ((material_grid *) malloc(sizeof(material_grid)));
material_grid_copy(o0->subclass.material_grid_data, o->subclass.material_grid_data);
}
else if (o0->which_subclass == MATERIAL_FUNCTION) {
o->which_subclass = MATERIAL_FUNCTION;
o->subclass.material_function_data = ((material_function *) malloc(sizeof(material_function)));
material_function_copy(o0->subclass.material_function_data, o->subclass.material_function_data);
}
else if (o0->which_subclass == MEDIUM_ANISOTROPIC) {
o->which_subclass = MEDIUM_ANISOTROPIC;
o->subclass.medium_anisotropic_data = ((medium_anisotropic *) malloc(sizeof(medium_anisotropic)));
medium_anisotropic_copy(o0->subclass.medium_anisotropic_data, o->subclass.medium_anisotropic_data);
}
else if (o0->which_subclass == MEDIUM) {
o->which_subclass = MEDIUM;
o->subclass.medium_data = ((medium *) malloc(sizeof(medium)));
medium_copy(o0->subclass.medium_data, o->subclass.medium_data);
}
else 
o->which_subclass = MATERIAL_TYPE_SELF;
}

/******* class equal functions *******/

boolean lattice_equal(const lattice *o0, const lattice *o)
{
if (!vector3_equal(o->basis1, o0->basis1)) return 0;
if (!vector3_equal(o->basis2, o0->basis2)) return 0;
if (!vector3_equal(o->basis3, o0->basis3)) return 0;
if (!vector3_equal(o->size, o0->size)) return 0;
if (!vector3_equal(o->basis_size, o0->basis_size)) return 0;
if (!vector3_equal(o->b1, o0->b1)) return 0;
if (!vector3_equal(o->b2, o0->b2)) return 0;
if (!vector3_equal(o->b3, o0->b3)) return 0;
if (!matrix3x3_equal(o->basis, o0->basis)) return 0;
if (!matrix3x3_equal(o->metric, o0->metric)) return 0;
;
return 1;
}

boolean ellipsoid_equal(const ellipsoid *o0, const ellipsoid *o)
{
if (!vector3_equal(o->inverse_semi_axes, o0->inverse_semi_axes)) return 0;
;
return 1;
}

boolean prism_equal(const prism *o0, const prism *o)
{
{
int i_t;
if (o->vertices.num_items != o0->vertices.num_items) return 0;
for (i_t = 0; i_t < o->vertices.num_items; i_t++) {
if (!vector3_equal(o->vertices.items[i_t], o0->vertices.items[i_t])) return 0;
}
}
if (!vector3_equal(o->centroid, o0->centroid)) return 0;
if (o->height != o0->height) return 0;
{
int i_t;
if (o->workspace.num_items != o0->workspace.num_items) return 0;
for (i_t = 0; i_t < o->workspace.num_items; i_t++) {
if (o->workspace.items[i_t] != o0->workspace.items[i_t]) return 0;
}
}
if (!matrix3x3_equal(o->m_c2p, o0->m_c2p)) return 0;
if (!matrix3x3_equal(o->m_p2c, o0->m_p2c)) return 0;
;
return 1;
}

boolean block_equal(const block *o0, const block *o)
{
if (!vector3_equal(o->e1, o0->e1)) return 0;
if (!vector3_equal(o->e2, o0->e2)) return 0;
if (!vector3_equal(o->e3, o0->e3)) return 0;
if (!vector3_equal(o->size, o0->size)) return 0;
if (!matrix3x3_equal(o->projection_matrix, o0->projection_matrix)) return 0;
if (o0->which_subclass != o->which_subclass) return 0;
if (o0->which_subclass == ELLIPSOID) {
if (!ellipsoid_equal(o0->subclass.ellipsoid_data, o->subclass.ellipsoid_data)) return 0;
}
else ;
return 1;
}

boolean sphere_equal(const sphere *o0, const sphere *o)
{
if (o->radius != o0->radius) return 0;
;
return 1;
}

boolean wedge_equal(const wedge *o0, const wedge *o)
{
if (o->wedge_angle != o0->wedge_angle) return 0;
if (!vector3_equal(o->wedge_start, o0->wedge_start)) return 0;
if (!vector3_equal(o->e1, o0->e1)) return 0;
if (!vector3_equal(o->e2, o0->e2)) return 0;
;
return 1;
}

boolean cone_equal(const cone *o0, const cone *o)
{
if (o->radius2 != o0->radius2) return 0;
;
return 1;
}

boolean cylinder_equal(const cylinder *o0, const cylinder *o)
{
if (!vector3_equal(o->axis, o0->axis)) return 0;
if (o->radius != o0->radius) return 0;
if (o->height != o0->height) return 0;
if (o0->which_subclass != o->which_subclass) return 0;
if (o0->which_subclass == WEDGE) {
if (!wedge_equal(o0->subclass.wedge_data, o->subclass.wedge_data)) return 0;
}
else if (o0->which_subclass == CONE) {
if (!cone_equal(o0->subclass.cone_data, o->subclass.cone_data)) return 0;
}
else ;
return 1;
}

boolean compound_geometric_object_equal(const compound_geometric_object *o0, const compound_geometric_object *o)
{
{
int i_t;
if (o->component_objects.num_items != o0->component_objects.num_items) return 0;
for (i_t = 0; i_t < o->component_objects.num_items; i_t++) {
if (!geometric_object_equal(&o0->component_objects.items[i_t], &o->component_objects.items[i_t])) return 0;
}
}
;
return 1;
}

boolean geometric_object_equal(const geometric_object *o0, const geometric_object *o)
{
if (!material_type_equal(&o0->material, &o->material)) return 0;
if (!vector3_equal(o->center, o0->center)) return 0;
if (o0->which_subclass != o->which_subclass) return 0;
if (o0->which_subclass == PRISM) {
if (!prism_equal(o0->subclass.prism_data, o->subclass.prism_data)) return 0;
}
else if (o0->which_subclass == BLOCK) {
if (!block_equal(o0->subclass.block_data, o->subclass.block_data)) return 0;
}
else if (o0->which_subclass == SPHERE) {
if (!sphere_equal(o0->subclass.sphere_data, o->subclass.sphere_data)) return 0;
}
else if (o0->which_subclass == CYLINDER) {
if (!cylinder_equal(o0->subclass.cylinder_data, o->subclass.cylinder_data)) return 0;
}
else if (o0->which_subclass == COMPOUND_GEOMETRIC_OBJECT) {
if (!compound_geometric_object_equal(o0->subclass.compound_geometric_object_data, o->subclass.compound_geometric_object_data)) return 0;
}
else ;
return 1;
}

boolean material_grid_equal(const material_grid *o0, const material_grid *o)
{
if (o->material_grid_kind != o0->material_grid_kind) return 0;
if (o->epsilon_min != o0->epsilon_min) return 0;
if (o->epsilon_max != o0->epsilon_max) return 0;
if (o->mu_min != o0->mu_min) return 0;
if (o->mu_max != o0->mu_max) return 0;
if (!vector3_equal(o->size, o0->size)) return 0;
if (o->matgrid_init != o0->matgrid_init) return 0;
if (o->matgrid != o0->matgrid) return 0;
;
return 1;
}

boolean material_function_equal(const material_function *o0, const material_function *o)
{
if (o->material_func != o0->material_func) return 0;
;
return 1;
}

boolean medium_anisotropic_equal(const medium_anisotropic *o0, const medium_anisotropic *o)
{
if (!vector3_equal(o->epsilon_diag, o0->epsilon_diag)) return 0;
if (!cvector3_equal(o->epsilon_offdiag, o0->epsilon_offdiag)) return 0;
if (!vector3_equal(o->epsilon_offdiag_imag, o0->epsilon_offdiag_imag)) return 0;
if (!vector3_equal(o->mu_diag, o0->mu_diag)) return 0;
if (!cvector3_equal(o->mu_offdiag, o0->mu_offdiag)) return 0;
if (!vector3_equal(o->mu_offdiag_imag, o0->mu_offdiag_imag)) return 0;
;
return 1;
}

boolean medium_equal(const medium *o0, const medium *o)
{
if (o->epsilon != o0->epsilon) return 0;
if (o->mu != o0->mu) return 0;
;
return 1;
}

boolean material_type_equal(const material_type *o0, const material_type *o)
{
if (o0->which_subclass != o->which_subclass) return 0;
if (o0->which_subclass == MATERIAL_GRID) {
if (!material_grid_equal(o0->subclass.material_grid_data, o->subclass.material_grid_data)) return 0;
}
else if (o0->which_subclass == MATERIAL_FUNCTION) {
if (!material_function_equal(o0->subclass.material_function_data, o->subclass.material_function_data)) return 0;
}
else if (o0->which_subclass == MEDIUM_ANISOTROPIC) {
if (!medium_anisotropic_equal(o0->subclass.medium_anisotropic_data, o->subclass.medium_anisotropic_data)) return 0;
}
else if (o0->which_subclass == MEDIUM) {
if (!medium_equal(o0->subclass.medium_data, o->subclass.medium_data)) return 0;
}
else ;
return 1;
}

/******* class destruction functions *******/

void lattice_destroy(lattice o)
{
}

void ellipsoid_destroy(ellipsoid o)
{
}

void prism_destroy(prism o)
{
{
int index_t;
for (index_t = 0; index_t < o.vertices.num_items; index_t++) {
}
}
free(o.vertices.items);
{
int index_t;
for (index_t = 0; index_t < o.workspace.num_items; index_t++) {
}
}
free(o.workspace.items);
}

void block_destroy(block o)
{
if (o.which_subclass == ELLIPSOID) {
ellipsoid_destroy(*o.subclass.ellipsoid_data);
free(o.subclass.ellipsoid_data);
}
else { }
}

void sphere_destroy(sphere o)
{
}

void wedge_destroy(wedge o)
{
}

void cone_destroy(cone o)
{
}

void cylinder_destroy(cylinder o)
{
if (o.which_subclass == WEDGE) {
wedge_destroy(*o.subclass.wedge_data);
free(o.subclass.wedge_data);
}
else if (o.which_subclass == CONE) {
cone_destroy(*o.subclass.cone_data);
free(o.subclass.cone_data);
}
else { }
}

void compound_geometric_object_destroy(compound_geometric_object o)
{
{
int index_t;
for (index_t = 0; index_t < o.component_objects.num_items; index_t++) {
geometric_object_destroy(o.component_objects.items[index_t]);
}
}
free(o.component_objects.items);
}

void geometric_object_destroy(geometric_object o)
{
material_type_destroy(o.material);
if (o.which_subclass == PRISM) {
prism_destroy(*o.subclass.prism_data);
free(o.subclass.prism_data);
}
else if (o.which_subclass == BLOCK) {
block_destroy(*o.subclass.block_data);
free(o.subclass.block_data);
}
else if (o.which_subclass == SPHERE) {
sphere_destroy(*o.subclass.sphere_data);
free(o.subclass.sphere_data);
}
else if (o.which_subclass == CYLINDER) {
cylinder_destroy(*o.subclass.cylinder_data);
free(o.subclass.cylinder_data);
}
else if (o.which_subclass == COMPOUND_GEOMETRIC_OBJECT) {
compound_geometric_object_destroy(*o.subclass.compound_geometric_object_data);
free(o.subclass.compound_geometric_object_data);
}
else { }
}

void material_grid_destroy(material_grid o)
{
}

void material_function_destroy(material_function o)
{
}

void medium_anisotropic_destroy(medium_anisotropic o)
{
}

void medium_destroy(medium o)
{
}

void material_type_destroy(material_type o)
{
if (o.which_subclass == MATERIAL_GRID) {
material_grid_destroy(*o.subclass.material_grid_data);
free(o.subclass.material_grid_data);
}
else if (o.which_subclass == MATERIAL_FUNCTION) {
material_function_destroy(*o.subclass.material_function_data);
free(o.subclass.material_function_data);
}
else if (o.which_subclass == MEDIUM_ANISOTROPIC) {
medium_anisotropic_destroy(*o.subclass.medium_anisotropic_data);
free(o.subclass.medium_anisotropic_data);
}
else if (o.which_subclass == MEDIUM) {
medium_destroy(*o.subclass.medium_data);
free(o.subclass.medium_data);
}
else { }
}

/******* read input variables *******/

SCM read_input_vars(void)
{
if (num_read_input_vars++) destroy_input_vars();
dimensions = ctl_get_integer("dimensions");
material_type_input(ctl_get_object("default-material"), &default_material);
geometry_center = ctl_get_vector3("geometry-center");
lattice_input(ctl_get_object("geometry-lattice"), &geometry_lattice);
{
list lo_t = ctl_get_list("geometry");
int i_t;
geometry.num_items = list_length(lo_t);
geometry.items = ((geometric_object *) malloc(sizeof(geometric_object) * (geometry.num_items)));
for (i_t = 0; i_t < geometry.num_items; i_t++) {
geometric_object_input(object_list_ref(lo_t, i_t), &geometry.items[i_t]);
}
}
ensure_periodicity = ctl_get_boolean("ensure-periodicity");
{
list lo_t = ctl_get_list("k-points");
int i_t;
k_points.num_items = list_length(lo_t);
k_points.items = ((vector3 *) malloc(sizeof(vector3) * (k_points.num_items)));
for (i_t = 0; i_t < k_points.num_items; i_t++) {
k_points.items[i_t] = vector3_list_ref(lo_t, i_t);
}
}
num_bands = ctl_get_integer("num-bands");
tolerance = ctl_get_number("tolerance");
target_freq = ctl_get_number("target-freq");
mesh_size = ctl_get_integer("mesh-size");
epsilon_input_file = ctl_get_string("epsilon-input-file");
mu_input_file = ctl_get_string("mu-input-file");
force_mup = ctl_get_boolean("force-mu?");
deterministicp = ctl_get_boolean("deterministic?");
simple_preconditionerp = ctl_get_boolean("simple-preconditioner?");
eigensolver_flags = ctl_get_integer("eigensolver-flags");
eigensolver_block_size = ctl_get_integer("eigensolver-block-size");
eigensolver_nwork = ctl_get_integer("eigensolver-nwork");
eigensolver_davidsonp = ctl_get_boolean("eigensolver-davidson?");
eigensolver_flops = ctl_get_number("eigensolver-flops");
negative_epsilon_okp = ctl_get_boolean("negative-epsilon-ok?");
return SCM_UNSPECIFIED;
}

/******* write output variables *******/

SCM write_output_vars(void)
{
num_write_output_vars++;
ctl_set_number("eigensolver-flops", eigensolver_flops);
ctl_set_list("freqs", make_number_list(freqs.num_items, freqs.items));
ctl_set_integer("iterations", iterations);
ctl_set_string("parity", parity);
return SCM_UNSPECIFIED;
}

/******* destroy input variables *******/

SCM destroy_input_vars(void)
{
material_type_destroy(default_material);
lattice_destroy(geometry_lattice);
{
int index_t;
for (index_t = 0; index_t < geometry.num_items; index_t++) {
geometric_object_destroy(geometry.items[index_t]);
}
}
free(geometry.items);
{
int index_t;
for (index_t = 0; index_t < k_points.num_items; index_t++) {
}
}
free(k_points.items);
free(epsilon_input_file);
free(mu_input_file);
return SCM_UNSPECIFIED;
}

/******* destroy output variables *******/

SCM destroy_output_vars(void)
{
{
int index_t;
for (index_t = 0; index_t < freqs.num_items; index_t++) {
}
}
free(freqs.items);
free(parity);
return SCM_UNSPECIFIED;
}

/******* external-functions *******/

SCM compute_1_group_velocity_reciprocal_aux(SCM arg_scm_0)
{
SCM return_val_scm;
vector3 return_val_c;
integer arg_c_0;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = compute_1_group_velocity_reciprocal(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_vector3_to_scm(return_val_c);
return return_val_scm;
}

SCM compute_1_group_velocity_aux(SCM arg_scm_0)
{
SCM return_val_scm;
vector3 return_val_c;
integer arg_c_0;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = compute_1_group_velocity(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_vector3_to_scm(return_val_c);
return return_val_scm;
}

SCM compute_1_group_velocity_component_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
number return_val_c;
vector3 arg_c_0;
integer arg_c_1;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);
arg_c_1 = ctl_convert_integer_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = compute_1_group_velocity_component(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM compute_group_velocity_component_aux(SCM arg_scm_0)
{
SCM return_val_scm;
number_list return_val_c;
vector3 arg_c_0;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = compute_group_velocity_component(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_list_to_scm(make_number_list(return_val_c.num_items, return_val_c.items));
{
int index_t;
for (index_t = 0; index_t < return_val_c.num_items; index_t++) {
}
}
free(return_val_c.items);
return return_val_scm;
}

SCM compute_yparities_aux(void)
{
SCM return_val_scm;
number_list return_val_c;


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = compute_yparities();

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_list_to_scm(make_number_list(return_val_c.num_items, return_val_c.items));
{
int index_t;
for (index_t = 0; index_t < return_val_c.num_items; index_t++) {
}
}
free(return_val_c.items);
return return_val_scm;
}

SCM compute_zparities_aux(void)
{
SCM return_val_scm;
number_list return_val_c;


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = compute_zparities();

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_list_to_scm(make_number_list(return_val_c.num_items, return_val_c.items));
{
int index_t;
for (index_t = 0; index_t < return_val_c.num_items; index_t++) {
}
}
free(return_val_c.items);
return return_val_scm;
}

SCM material_grids_min_tetm_gap_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2, SCM arg_scm_3, SCM arg_scm_4, SCM arg_scm_5)
{
SCM return_val_scm;
number return_val_c;
vector3 arg_c_0;
integer arg_c_1;
number arg_c_2;
number arg_c_3;
integer arg_c_4;
number arg_c_5;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);
arg_c_1 = ctl_convert_integer_to_c(arg_scm_1);
arg_c_2 = ctl_convert_number_to_c(arg_scm_2);
arg_c_3 = ctl_convert_number_to_c(arg_scm_3);
arg_c_4 = ctl_convert_integer_to_c(arg_scm_4);
arg_c_5 = ctl_convert_number_to_c(arg_scm_5);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = material_grids_min_tetm_gap(arg_c_0, arg_c_1, arg_c_2, arg_c_3, arg_c_4, arg_c_5);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM material_grids_mingap_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2, SCM arg_scm_3, SCM arg_scm_4, SCM arg_scm_5, SCM arg_scm_6)
{
SCM return_val_scm;
number return_val_c;
vector3_list arg_c_0;
integer arg_c_1;
integer arg_c_2;
number arg_c_3;
number arg_c_4;
integer arg_c_5;
number arg_c_6;

{
list lo_t = ctl_convert_list_to_c(arg_scm_0);
int i_t;
arg_c_0.num_items = list_length(lo_t);
arg_c_0.items = ((vector3 *) malloc(sizeof(vector3) * (arg_c_0.num_items)));
for (i_t = 0; i_t < arg_c_0.num_items; i_t++) {
arg_c_0.items[i_t] = vector3_list_ref(lo_t, i_t);
}
}
arg_c_1 = ctl_convert_integer_to_c(arg_scm_1);
arg_c_2 = ctl_convert_integer_to_c(arg_scm_2);
arg_c_3 = ctl_convert_number_to_c(arg_scm_3);
arg_c_4 = ctl_convert_number_to_c(arg_scm_4);
arg_c_5 = ctl_convert_integer_to_c(arg_scm_5);
arg_c_6 = ctl_convert_number_to_c(arg_scm_6);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = material_grids_mingap(arg_c_0, arg_c_1, arg_c_2, arg_c_3, arg_c_4, arg_c_5, arg_c_6);

fflush(stdout); fflush(stderr);
{
int index_t;
for (index_t = 0; index_t < arg_c_0.num_items; index_t++) {
}
}
free(arg_c_0.items);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM material_grids_maxgap_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2, SCM arg_scm_3, SCM arg_scm_4, SCM arg_scm_5, SCM arg_scm_6)
{
SCM return_val_scm;
number return_val_c;
vector3_list arg_c_0;
integer arg_c_1;
integer arg_c_2;
number arg_c_3;
number arg_c_4;
integer arg_c_5;
number arg_c_6;

{
list lo_t = ctl_convert_list_to_c(arg_scm_0);
int i_t;
arg_c_0.num_items = list_length(lo_t);
arg_c_0.items = ((vector3 *) malloc(sizeof(vector3) * (arg_c_0.num_items)));
for (i_t = 0; i_t < arg_c_0.num_items; i_t++) {
arg_c_0.items[i_t] = vector3_list_ref(lo_t, i_t);
}
}
arg_c_1 = ctl_convert_integer_to_c(arg_scm_1);
arg_c_2 = ctl_convert_integer_to_c(arg_scm_2);
arg_c_3 = ctl_convert_number_to_c(arg_scm_3);
arg_c_4 = ctl_convert_number_to_c(arg_scm_4);
arg_c_5 = ctl_convert_integer_to_c(arg_scm_5);
arg_c_6 = ctl_convert_number_to_c(arg_scm_6);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = material_grids_maxgap(arg_c_0, arg_c_1, arg_c_2, arg_c_3, arg_c_4, arg_c_5, arg_c_6);

fflush(stdout); fflush(stderr);
{
int index_t;
for (index_t = 0; index_t < arg_c_0.num_items; index_t++) {
}
}
free(arg_c_0.items);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM material_grids_approx_gradient_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2, SCM arg_scm_3)
{
SCM return_val_scm;
number return_val_c;
vector3 arg_c_0;
integer arg_c_1;
integer arg_c_2;
number arg_c_3;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);
arg_c_1 = ctl_convert_integer_to_c(arg_scm_1);
arg_c_2 = ctl_convert_integer_to_c(arg_scm_2);
arg_c_3 = ctl_convert_number_to_c(arg_scm_3);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = material_grids_approx_gradient(arg_c_0, arg_c_1, arg_c_2, arg_c_3);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM print_material_grids_deps_du_numeric_aux(SCM arg_scm_0)
{
number arg_c_0;

arg_c_0 = ctl_convert_number_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
print_material_grids_deps_du_numeric(arg_c_0);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM print_material_grids_deps_du_aux(void)
{


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
print_material_grids_deps_du();

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM print_material_grids_gradient_aux(SCM arg_scm_0)
{
integer arg_c_0;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
print_material_grids_gradient(arg_c_0);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM material_grids_match_epsilon_fileB_aux(SCM arg_scm_0, SCM arg_scm_1)
{
char* arg_c_0;
number arg_c_1;

arg_c_0 = ctl_convert_string_to_c(arg_scm_0);
arg_c_1 = ctl_convert_number_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
material_grids_match_epsilon_fileB(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);
free(arg_c_0);

return SCM_UNSPECIFIED;
}

SCM load_material_gridB_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2)
{
material_grid arg_c_0;
char* arg_c_1;
vector3 arg_c_2;

material_grid_input(ctl_convert_object_to_c(arg_scm_0), &arg_c_0);
arg_c_1 = ctl_convert_string_to_c(arg_scm_1);
arg_c_2 = ctl_convert_vector3_to_c(arg_scm_2);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
load_material_gridB(arg_c_0, arg_c_1, arg_c_2);

fflush(stdout); fflush(stderr);
material_grid_destroy(arg_c_0);
free(arg_c_1);

return SCM_UNSPECIFIED;
}

SCM save_material_grid_aux(SCM arg_scm_0, SCM arg_scm_1)
{
material_grid arg_c_0;
char* arg_c_1;

material_grid_input(ctl_convert_object_to_c(arg_scm_0), &arg_c_0);
arg_c_1 = ctl_convert_string_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
save_material_grid(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);
material_grid_destroy(arg_c_0);
free(arg_c_1);

return SCM_UNSPECIFIED;
}

SCM randomize_material_gridB_aux(SCM arg_scm_0, SCM arg_scm_1)
{
material_grid arg_c_0;
number arg_c_1;

material_grid_input(ctl_convert_object_to_c(arg_scm_0), &arg_c_0);
arg_c_1 = ctl_convert_number_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
randomize_material_gridB(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);
material_grid_destroy(arg_c_0);

return SCM_UNSPECIFIED;
}

SCM cvector_field_get_point_bloch_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
cvector3 return_val_c;
SCM arg_c_0;
vector3 arg_c_1;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_vector3_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = cvector_field_get_point_bloch(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_cvector3_to_scm(return_val_c);
return return_val_scm;
}

SCM cvector_field_get_point_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
cvector3 return_val_c;
SCM arg_c_0;
vector3 arg_c_1;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_vector3_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = cvector_field_get_point(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_cvector3_to_scm(return_val_c);
return return_val_scm;
}

SCM cscalar_field_get_point_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
cnumber return_val_c;
SCM arg_c_0;
vector3 arg_c_1;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_vector3_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = cscalar_field_get_point(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_cnumber_to_scm(return_val_c);
return return_val_scm;
}

SCM rscalar_field_get_point_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
number return_val_c;
SCM arg_c_0;
vector3 arg_c_1;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_vector3_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = rscalar_field_get_point(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM integrate_fieldL_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
cnumber return_val_c;
function arg_c_0;
SCM_list arg_c_1;

arg_c_0 = ctl_convert_function_to_c(arg_scm_0);
{
list lo_t = ctl_convert_list_to_c(arg_scm_1);
int i_t;
arg_c_1.num_items = list_length(lo_t);
arg_c_1.items = ((SCM *) malloc(sizeof(SCM) * (arg_c_1.num_items)));
for (i_t = 0; i_t < arg_c_1.num_items; i_t++) {
arg_c_1.items[i_t] = SCM_list_ref(lo_t, i_t);
}
}

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = integrate_fieldL(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);
{
int index_t;
for (index_t = 0; index_t < arg_c_1.num_items; index_t++) {
}
}
free(arg_c_1.items);

return_val_scm = ctl_convert_cnumber_to_scm(return_val_c);
return return_val_scm;
}

SCM field_mapLB_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2)
{
SCM arg_c_0;
function arg_c_1;
SCM_list arg_c_2;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_function_to_c(arg_scm_1);
{
list lo_t = ctl_convert_list_to_c(arg_scm_2);
int i_t;
arg_c_2.num_items = list_length(lo_t);
arg_c_2.items = ((SCM *) malloc(sizeof(SCM) * (arg_c_2.num_items)));
for (i_t = 0; i_t < arg_c_2.num_items; i_t++) {
arg_c_2.items[i_t] = SCM_list_ref(lo_t, i_t);
}
}

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
field_mapLB(arg_c_0, arg_c_1, arg_c_2);

fflush(stdout); fflush(stderr);
{
int index_t;
for (index_t = 0; index_t < arg_c_2.num_items; index_t++) {
}
}
free(arg_c_2.items);

return SCM_UNSPECIFIED;
}

SCM field_load_aux(SCM arg_scm_0)
{
SCM arg_c_0;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
field_load(arg_c_0);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM field_setB_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM arg_c_0;
SCM arg_c_1;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_SCM_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
field_setB(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM fields_conformp_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
boolean return_val_c;
SCM arg_c_0;
SCM arg_c_1;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_SCM_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = fields_conformp(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_boolean_to_scm(return_val_c);
return return_val_scm;
}

SCM field_make_aux(SCM arg_scm_0)
{
SCM return_val_scm;
SCM return_val_c;
SCM arg_c_0;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = field_make(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_SCM_to_scm(return_val_c);
return return_val_scm;
}

SCM cvector_field_nonblochB_aux(SCM arg_scm_0)
{
SCM arg_c_0;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
cvector_field_nonblochB(arg_c_0);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM cvector_field_make_aux(SCM arg_scm_0)
{
SCM return_val_scm;
SCM return_val_c;
SCM arg_c_0;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = cvector_field_make(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_SCM_to_scm(return_val_c);
return return_val_scm;
}

SCM rscalar_field_make_aux(SCM arg_scm_0)
{
SCM return_val_scm;
SCM return_val_c;
SCM arg_c_0;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = rscalar_field_make(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_SCM_to_scm(return_val_c);
return return_val_scm;
}

SCM cur_fieldp_aux(SCM arg_scm_0)
{
SCM return_val_scm;
boolean return_val_c;
SCM arg_c_0;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = cur_fieldp(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_boolean_to_scm(return_val_c);
return return_val_scm;
}

SCM load_eigenvectors_aux(SCM arg_scm_0)
{
char* arg_c_0;

arg_c_0 = ctl_convert_string_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
load_eigenvectors(arg_c_0);

fflush(stdout); fflush(stderr);
free(arg_c_0);

return SCM_UNSPECIFIED;
}

SCM save_eigenvectors_aux(SCM arg_scm_0)
{
char* arg_c_0;

arg_c_0 = ctl_convert_string_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
save_eigenvectors(arg_c_0);

fflush(stdout); fflush(stderr);
free(arg_c_0);

return SCM_UNSPECIFIED;
}

SCM input_eigenvectors_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
SCM return_val_c;
char* arg_c_0;
integer arg_c_1;

arg_c_0 = ctl_convert_string_to_c(arg_scm_0);
arg_c_1 = ctl_convert_integer_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = input_eigenvectors(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);
free(arg_c_0);

return_val_scm = ctl_convert_SCM_to_scm(return_val_c);
return return_val_scm;
}

SCM output_eigenvectors_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM arg_c_0;
char* arg_c_1;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_string_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
output_eigenvectors(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);
free(arg_c_1);

return SCM_UNSPECIFIED;
}

SCM scale_eigenvector_aux(SCM arg_scm_0, SCM arg_scm_1)
{
integer arg_c_0;
cnumber arg_c_1;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);
arg_c_1 = ctl_convert_cnumber_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
scale_eigenvector(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM dot_eigenvectors_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
SCM return_val_c;
SCM arg_c_0;
integer arg_c_1;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_integer_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = dot_eigenvectors(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_SCM_to_scm(return_val_c);
return return_val_scm;
}

SCM set_eigenvectors_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM arg_c_0;
integer arg_c_1;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_integer_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
set_eigenvectors(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM get_eigenvectors_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
SCM return_val_c;
integer arg_c_0;
integer arg_c_1;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);
arg_c_1 = ctl_convert_integer_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = get_eigenvectors(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_SCM_to_scm(return_val_c);
return return_val_scm;
}

SCM sqmatrix_eigvals_aux(SCM arg_scm_0)
{
SCM return_val_scm;
cnumber_list return_val_c;
SCM arg_c_0;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = sqmatrix_eigvals(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_list_to_scm(make_cnumber_list(return_val_c.num_items, return_val_c.items));
{
int index_t;
for (index_t = 0; index_t < return_val_c.num_items; index_t++) {
}
}
free(return_val_c.items);
return return_val_scm;
}

SCM sqmatrix_diagm_aux(SCM arg_scm_0)
{
SCM return_val_scm;
SCM return_val_c;
cnumber_list arg_c_0;

{
list lo_t = ctl_convert_list_to_c(arg_scm_0);
int i_t;
arg_c_0.num_items = list_length(lo_t);
arg_c_0.items = ((cnumber *) malloc(sizeof(cnumber) * (arg_c_0.num_items)));
for (i_t = 0; i_t < arg_c_0.num_items; i_t++) {
arg_c_0.items[i_t] = cnumber_list_ref(lo_t, i_t);
}
}

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = sqmatrix_diagm(arg_c_0);

fflush(stdout); fflush(stderr);
{
int index_t;
for (index_t = 0; index_t < arg_c_0.num_items; index_t++) {
}
}
free(arg_c_0.items);

return_val_scm = ctl_convert_SCM_to_scm(return_val_c);
return return_val_scm;
}

SCM sqmatrix_mult_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
SCM return_val_c;
SCM arg_c_0;
SCM arg_c_1;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_SCM_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = sqmatrix_mult(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_SCM_to_scm(return_val_c);
return return_val_scm;
}

SCM sqmatrix_ref_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2)
{
SCM return_val_scm;
cnumber return_val_c;
SCM arg_c_0;
integer arg_c_1;
integer arg_c_2;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);
arg_c_1 = ctl_convert_integer_to_c(arg_scm_1);
arg_c_2 = ctl_convert_integer_to_c(arg_scm_2);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = sqmatrix_ref(arg_c_0, arg_c_1, arg_c_2);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_cnumber_to_scm(return_val_c);
return return_val_scm;
}

SCM sqmatrix_size_aux(SCM arg_scm_0)
{
SCM return_val_scm;
integer return_val_c;
SCM arg_c_0;

arg_c_0 = ctl_convert_SCM_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = sqmatrix_size(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_integer_to_scm(return_val_c);
return return_val_scm;
}

SCM set_kpoint_index_aux(SCM arg_scm_0)
{
integer arg_c_0;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
set_kpoint_index(arg_c_0);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM get_kpoint_index_aux(void)
{
SCM return_val_scm;
integer return_val_c;


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = get_kpoint_index();

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_integer_to_scm(return_val_c);
return return_val_scm;
}

SCM has_inversion_symp_aux(void)
{
SCM return_val_scm;
boolean return_val_c;


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = has_inversion_symp();

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_boolean_to_scm(return_val_c);
return return_val_scm;
}

SCM has_hermitian_epsp_aux(void)
{
SCM return_val_scm;
boolean return_val_c;


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = has_hermitian_epsp();

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_boolean_to_scm(return_val_c);
return return_val_scm;
}

SCM mpi_max_aux(SCM arg_scm_0)
{
SCM return_val_scm;
number return_val_c;
number arg_c_0;

arg_c_0 = ctl_convert_number_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = mpi_max(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM mpi_proc_index_aux(void)
{
SCM return_val_scm;
integer return_val_c;


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = mpi_proc_index();

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_integer_to_scm(return_val_c);
return return_val_scm;
}

SCM mpi_num_procs_aux(void)
{
SCM return_val_scm;
integer return_val_c;


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = mpi_num_procs();

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_integer_to_scm(return_val_c);
return return_val_scm;
}

SCM using_mpip_aux(void)
{
SCM return_val_scm;
boolean return_val_c;


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = using_mpip();

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_boolean_to_scm(return_val_c);
return return_val_scm;
}

SCM mpi_is_masterp_aux(void)
{
SCM return_val_scm;
boolean return_val_c;


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = mpi_is_masterp();

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_boolean_to_scm(return_val_c);
return return_val_scm;
}

SCM output_field_to_file_aux(SCM arg_scm_0, SCM arg_scm_1)
{
integer arg_c_0;
char* arg_c_1;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);
arg_c_1 = ctl_convert_string_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
output_field_to_file(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);
free(arg_c_1);

return SCM_UNSPECIFIED;
}

SCM compute_energy_in_object_list_aux(SCM arg_scm_0)
{
SCM return_val_scm;
number return_val_c;
geometric_object_list arg_c_0;

{
list lo_t = ctl_convert_list_to_c(arg_scm_0);
int i_t;
arg_c_0.num_items = list_length(lo_t);
arg_c_0.items = ((geometric_object *) malloc(sizeof(geometric_object) * (arg_c_0.num_items)));
for (i_t = 0; i_t < arg_c_0.num_items; i_t++) {
geometric_object_input(object_list_ref(lo_t, i_t), &arg_c_0.items[i_t]);
}
}

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = compute_energy_in_object_list(arg_c_0);

fflush(stdout); fflush(stderr);
{
int index_t;
for (index_t = 0; index_t < arg_c_0.num_items; index_t++) {
geometric_object_destroy(arg_c_0.items[index_t]);
}
}
free(arg_c_0.items);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM compute_energy_integral_aux(SCM arg_scm_0)
{
SCM return_val_scm;
number return_val_c;
function arg_c_0;

arg_c_0 = ctl_convert_function_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = compute_energy_integral(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM compute_field_integral_aux(SCM arg_scm_0)
{
SCM return_val_scm;
cnumber return_val_c;
function arg_c_0;

arg_c_0 = ctl_convert_function_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = compute_field_integral(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_cnumber_to_scm(return_val_c);
return return_val_scm;
}

SCM compute_energy_in_dielectric_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
number return_val_c;
number arg_c_0;
number arg_c_1;

arg_c_0 = ctl_convert_number_to_c(arg_scm_0);
arg_c_1 = ctl_convert_number_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = compute_energy_in_dielectric(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM get_cscalar_point_aux(SCM arg_scm_0)
{
SCM return_val_scm;
cnumber return_val_c;
vector3 arg_c_0;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = get_cscalar_point(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_cnumber_to_scm(return_val_c);
return return_val_scm;
}

SCM get_bloch_cscalar_point_aux(SCM arg_scm_0)
{
SCM return_val_scm;
cnumber return_val_c;
vector3 arg_c_0;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = get_bloch_cscalar_point(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_cnumber_to_scm(return_val_c);
return return_val_scm;
}

SCM get_field_point_aux(SCM arg_scm_0)
{
SCM return_val_scm;
cvector3 return_val_c;
vector3 arg_c_0;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = get_field_point(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_cvector3_to_scm(return_val_c);
return return_val_scm;
}

SCM get_bloch_field_point_aux(SCM arg_scm_0)
{
SCM return_val_scm;
cvector3 return_val_c;
vector3 arg_c_0;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = get_bloch_field_point(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_cvector3_to_scm(return_val_c);
return return_val_scm;
}

SCM get_energy_point_aux(SCM arg_scm_0)
{
SCM return_val_scm;
number return_val_c;
vector3 arg_c_0;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = get_energy_point(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM get_epsilon_inverse_tensor_point_aux(SCM arg_scm_0)
{
SCM return_val_scm;
cmatrix3x3 return_val_c;
vector3 arg_c_0;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = get_epsilon_inverse_tensor_point(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_cmatrix3x3_to_scm(return_val_c);
return return_val_scm;
}

SCM get_epsilon_point_aux(SCM arg_scm_0)
{
SCM return_val_scm;
number return_val_c;
vector3 arg_c_0;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = get_epsilon_point(arg_c_0);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM compute_field_divergence_aux(void)
{


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
compute_field_divergence();

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM compute_field_energy_aux(void)
{
SCM return_val_scm;
number_list return_val_c;


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = compute_field_energy();

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_list_to_scm(make_number_list(return_val_c.num_items, return_val_c.items));
{
int index_t;
for (index_t = 0; index_t < return_val_c.num_items; index_t++) {
}
}
free(return_val_c.items);
return return_val_scm;
}

SCM fix_field_phase_aux(void)
{


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
fix_field_phase();

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM get_mu_aux(void)
{


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
get_mu();

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM get_epsilon_aux(void)
{


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
get_epsilon();

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM get_efield_from_dfield_aux(void)
{


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
get_efield_from_dfield();

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM get_bfield_aux(SCM arg_scm_0)
{
integer arg_c_0;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
get_bfield(arg_c_0);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM get_hfield_aux(SCM arg_scm_0)
{
integer arg_c_0;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
get_hfield(arg_c_0);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM get_dfield_aux(SCM arg_scm_0)
{
integer arg_c_0;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
get_dfield(arg_c_0);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM solve_kpoint_aux(SCM arg_scm_0)
{
vector3 arg_c_0;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
solve_kpoint(arg_c_0);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM randomize_fields_aux(void)
{


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
randomize_fields();

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM set_parity_aux(SCM arg_scm_0)
{
integer arg_c_0;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
set_parity(arg_c_0);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM using_mup_aux(void)
{
SCM return_val_scm;
boolean return_val_c;


#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = using_mup();

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_boolean_to_scm(return_val_c);
return return_val_scm;
}

SCM init_params_aux(SCM arg_scm_0, SCM arg_scm_1)
{
integer arg_c_0;
boolean arg_c_1;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);
arg_c_1 = ctl_convert_boolean_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
init_params(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return SCM_UNSPECIFIED;
}

SCM square_basis_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
matrix3x3 return_val_c;
matrix3x3 arg_c_0;
vector3 arg_c_1;

arg_c_0 = ctl_convert_matrix3x3_to_c(arg_scm_0);
arg_c_1 = ctl_convert_vector3_to_c(arg_scm_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = square_basis(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);

return_val_scm = ctl_convert_matrix3x3_to_scm(return_val_c);
return return_val_scm;
}

SCM range_overlap_with_object_aux(SCM arg_scm_0, SCM arg_scm_1, SCM arg_scm_2, SCM arg_scm_3, SCM arg_scm_4)
{
SCM return_val_scm;
number return_val_c;
vector3 arg_c_0;
vector3 arg_c_1;
geometric_object arg_c_2;
number arg_c_3;
integer arg_c_4;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);
arg_c_1 = ctl_convert_vector3_to_c(arg_scm_1);
geometric_object_input(ctl_convert_object_to_c(arg_scm_2), &arg_c_2);
arg_c_3 = ctl_convert_number_to_c(arg_scm_3);
arg_c_4 = ctl_convert_integer_to_c(arg_scm_4);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = range_overlap_with_object(arg_c_0, arg_c_1, arg_c_2, arg_c_3, arg_c_4);

fflush(stdout); fflush(stderr);
geometric_object_destroy(arg_c_2);

return_val_scm = ctl_convert_number_to_scm(return_val_c);
return return_val_scm;
}

SCM display_geometric_object_info_aux(SCM arg_scm_0, SCM arg_scm_1)
{
integer arg_c_0;
geometric_object arg_c_1;

arg_c_0 = ctl_convert_integer_to_c(arg_scm_0);
geometric_object_input(ctl_convert_object_to_c(arg_scm_1), &arg_c_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
display_geometric_object_info(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);
geometric_object_destroy(arg_c_1);

return SCM_UNSPECIFIED;
}

SCM point_in_periodic_objectp_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
boolean return_val_c;
vector3 arg_c_0;
geometric_object arg_c_1;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);
geometric_object_input(ctl_convert_object_to_c(arg_scm_1), &arg_c_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = point_in_periodic_objectp(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);
geometric_object_destroy(arg_c_1);

return_val_scm = ctl_convert_boolean_to_scm(return_val_c);
return return_val_scm;
}

SCM normal_to_object_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
vector3 return_val_c;
vector3 arg_c_0;
geometric_object arg_c_1;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);
geometric_object_input(ctl_convert_object_to_c(arg_scm_1), &arg_c_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = normal_to_object(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);
geometric_object_destroy(arg_c_1);

return_val_scm = ctl_convert_vector3_to_scm(return_val_c);
return return_val_scm;
}

SCM point_in_objectp_aux(SCM arg_scm_0, SCM arg_scm_1)
{
SCM return_val_scm;
boolean return_val_c;
vector3 arg_c_0;
geometric_object arg_c_1;

arg_c_0 = ctl_convert_vector3_to_c(arg_scm_0);
geometric_object_input(ctl_convert_object_to_c(arg_scm_1), &arg_c_1);

#ifdef HAVE_SCM_FLUSH_ALL_PORTS
scm_flush_all_ports();
#endif
return_val_c = point_in_objectp(arg_c_0, arg_c_1);

fflush(stdout); fflush(stderr);
geometric_object_destroy(arg_c_1);

return_val_scm = ctl_convert_boolean_to_scm(return_val_c);
return return_val_scm;
}

void export_external_functions(void)
{
gh_new_procedure("compute-1-group-velocity-reciprocal-aux", (SCM (*)()) compute_1_group_velocity_reciprocal_aux, 1, 0, 0);
gh_new_procedure("compute-1-group-velocity-aux", (SCM (*)()) compute_1_group_velocity_aux, 1, 0, 0);
gh_new_procedure("compute-1-group-velocity-component-aux", (SCM (*)()) compute_1_group_velocity_component_aux, 2, 0, 0);
gh_new_procedure("compute-group-velocity-component-aux", (SCM (*)()) compute_group_velocity_component_aux, 1, 0, 0);
gh_new_procedure("compute-yparities-aux", (SCM (*)()) compute_yparities_aux, 0, 0, 0);
gh_new_procedure("compute-zparities-aux", (SCM (*)()) compute_zparities_aux, 0, 0, 0);
gh_new_procedure("material-grids-min-tetm-gap-aux", (SCM (*)()) material_grids_min_tetm_gap_aux, 6, 0, 0);
gh_new_procedure("material-grids-mingap-aux", (SCM (*)()) material_grids_mingap_aux, 7, 0, 0);
gh_new_procedure("material-grids-maxgap-aux", (SCM (*)()) material_grids_maxgap_aux, 7, 0, 0);
gh_new_procedure("material-grids-approx-gradient-aux", (SCM (*)()) material_grids_approx_gradient_aux, 4, 0, 0);
gh_new_procedure("print-material-grids-deps-du-numeric-aux", (SCM (*)()) print_material_grids_deps_du_numeric_aux, 1, 0, 0);
gh_new_procedure("print-material-grids-deps-du-aux", (SCM (*)()) print_material_grids_deps_du_aux, 0, 0, 0);
gh_new_procedure("print-material-grids-gradient-aux", (SCM (*)()) print_material_grids_gradient_aux, 1, 0, 0);
gh_new_procedure("material-grids-match-epsilon-file!-aux", (SCM (*)()) material_grids_match_epsilon_fileB_aux, 2, 0, 0);
gh_new_procedure("load-material-grid!-aux", (SCM (*)()) load_material_gridB_aux, 3, 0, 0);
gh_new_procedure("save-material-grid-aux", (SCM (*)()) save_material_grid_aux, 2, 0, 0);
gh_new_procedure("randomize-material-grid!-aux", (SCM (*)()) randomize_material_gridB_aux, 2, 0, 0);
gh_new_procedure("cvector-field-get-point-bloch-aux", (SCM (*)()) cvector_field_get_point_bloch_aux, 2, 0, 0);
gh_new_procedure("cvector-field-get-point-aux", (SCM (*)()) cvector_field_get_point_aux, 2, 0, 0);
gh_new_procedure("cscalar-field-get-point-aux", (SCM (*)()) cscalar_field_get_point_aux, 2, 0, 0);
gh_new_procedure("rscalar-field-get-point-aux", (SCM (*)()) rscalar_field_get_point_aux, 2, 0, 0);
gh_new_procedure("integrate-fieldL-aux", (SCM (*)()) integrate_fieldL_aux, 2, 0, 0);
gh_new_procedure("field-mapL!-aux", (SCM (*)()) field_mapLB_aux, 3, 0, 0);
gh_new_procedure("field-load-aux", (SCM (*)()) field_load_aux, 1, 0, 0);
gh_new_procedure("field-set!-aux", (SCM (*)()) field_setB_aux, 2, 0, 0);
gh_new_procedure("fields-conform?-aux", (SCM (*)()) fields_conformp_aux, 2, 0, 0);
gh_new_procedure("field-make-aux", (SCM (*)()) field_make_aux, 1, 0, 0);
gh_new_procedure("cvector-field-nonbloch!-aux", (SCM (*)()) cvector_field_nonblochB_aux, 1, 0, 0);
gh_new_procedure("cvector-field-make-aux", (SCM (*)()) cvector_field_make_aux, 1, 0, 0);
gh_new_procedure("rscalar-field-make-aux", (SCM (*)()) rscalar_field_make_aux, 1, 0, 0);
gh_new_procedure("cur-field?-aux", (SCM (*)()) cur_fieldp_aux, 1, 0, 0);
gh_new_procedure("load-eigenvectors-aux", (SCM (*)()) load_eigenvectors_aux, 1, 0, 0);
gh_new_procedure("save-eigenvectors-aux", (SCM (*)()) save_eigenvectors_aux, 1, 0, 0);
gh_new_procedure("input-eigenvectors-aux", (SCM (*)()) input_eigenvectors_aux, 2, 0, 0);
gh_new_procedure("output-eigenvectors-aux", (SCM (*)()) output_eigenvectors_aux, 2, 0, 0);
gh_new_procedure("scale-eigenvector-aux", (SCM (*)()) scale_eigenvector_aux, 2, 0, 0);
gh_new_procedure("dot-eigenvectors-aux", (SCM (*)()) dot_eigenvectors_aux, 2, 0, 0);
gh_new_procedure("set-eigenvectors-aux", (SCM (*)()) set_eigenvectors_aux, 2, 0, 0);
gh_new_procedure("get-eigenvectors-aux", (SCM (*)()) get_eigenvectors_aux, 2, 0, 0);
gh_new_procedure("sqmatrix-eigvals-aux", (SCM (*)()) sqmatrix_eigvals_aux, 1, 0, 0);
gh_new_procedure("sqmatrix-diagm-aux", (SCM (*)()) sqmatrix_diagm_aux, 1, 0, 0);
gh_new_procedure("sqmatrix-mult-aux", (SCM (*)()) sqmatrix_mult_aux, 2, 0, 0);
gh_new_procedure("sqmatrix-ref-aux", (SCM (*)()) sqmatrix_ref_aux, 3, 0, 0);
gh_new_procedure("sqmatrix-size-aux", (SCM (*)()) sqmatrix_size_aux, 1, 0, 0);
gh_new_procedure("set-kpoint-index-aux", (SCM (*)()) set_kpoint_index_aux, 1, 0, 0);
gh_new_procedure("get-kpoint-index-aux", (SCM (*)()) get_kpoint_index_aux, 0, 0, 0);
gh_new_procedure("has-inversion-sym?-aux", (SCM (*)()) has_inversion_symp_aux, 0, 0, 0);
gh_new_procedure("has-hermitian-eps?-aux", (SCM (*)()) has_hermitian_epsp_aux, 0, 0, 0);
gh_new_procedure("mpi-max-aux", (SCM (*)()) mpi_max_aux, 1, 0, 0);
gh_new_procedure("mpi-proc-index-aux", (SCM (*)()) mpi_proc_index_aux, 0, 0, 0);
gh_new_procedure("mpi-num-procs-aux", (SCM (*)()) mpi_num_procs_aux, 0, 0, 0);
gh_new_procedure("using-mpi?-aux", (SCM (*)()) using_mpip_aux, 0, 0, 0);
gh_new_procedure("mpi-is-master?-aux", (SCM (*)()) mpi_is_masterp_aux, 0, 0, 0);
gh_new_procedure("output-field-to-file-aux", (SCM (*)()) output_field_to_file_aux, 2, 0, 0);
gh_new_procedure("compute-energy-in-object-list-aux", (SCM (*)()) compute_energy_in_object_list_aux, 1, 0, 0);
gh_new_procedure("compute-energy-integral-aux", (SCM (*)()) compute_energy_integral_aux, 1, 0, 0);
gh_new_procedure("compute-field-integral-aux", (SCM (*)()) compute_field_integral_aux, 1, 0, 0);
gh_new_procedure("compute-energy-in-dielectric-aux", (SCM (*)()) compute_energy_in_dielectric_aux, 2, 0, 0);
gh_new_procedure("get-cscalar-point-aux", (SCM (*)()) get_cscalar_point_aux, 1, 0, 0);
gh_new_procedure("get-bloch-cscalar-point-aux", (SCM (*)()) get_bloch_cscalar_point_aux, 1, 0, 0);
gh_new_procedure("get-field-point-aux", (SCM (*)()) get_field_point_aux, 1, 0, 0);
gh_new_procedure("get-bloch-field-point-aux", (SCM (*)()) get_bloch_field_point_aux, 1, 0, 0);
gh_new_procedure("get-energy-point-aux", (SCM (*)()) get_energy_point_aux, 1, 0, 0);
gh_new_procedure("get-epsilon-inverse-tensor-point-aux", (SCM (*)()) get_epsilon_inverse_tensor_point_aux, 1, 0, 0);
gh_new_procedure("get-epsilon-point-aux", (SCM (*)()) get_epsilon_point_aux, 1, 0, 0);
gh_new_procedure("compute-field-divergence-aux", (SCM (*)()) compute_field_divergence_aux, 0, 0, 0);
gh_new_procedure("compute-field-energy-aux", (SCM (*)()) compute_field_energy_aux, 0, 0, 0);
gh_new_procedure("fix-field-phase-aux", (SCM (*)()) fix_field_phase_aux, 0, 0, 0);
gh_new_procedure("get-mu-aux", (SCM (*)()) get_mu_aux, 0, 0, 0);
gh_new_procedure("get-epsilon-aux", (SCM (*)()) get_epsilon_aux, 0, 0, 0);
gh_new_procedure("get-efield-from-dfield-aux", (SCM (*)()) get_efield_from_dfield_aux, 0, 0, 0);
gh_new_procedure("get-bfield-aux", (SCM (*)()) get_bfield_aux, 1, 0, 0);
gh_new_procedure("get-hfield-aux", (SCM (*)()) get_hfield_aux, 1, 0, 0);
gh_new_procedure("get-dfield-aux", (SCM (*)()) get_dfield_aux, 1, 0, 0);
gh_new_procedure("solve-kpoint-aux", (SCM (*)()) solve_kpoint_aux, 1, 0, 0);
gh_new_procedure("randomize-fields-aux", (SCM (*)()) randomize_fields_aux, 0, 0, 0);
gh_new_procedure("set-parity-aux", (SCM (*)()) set_parity_aux, 1, 0, 0);
gh_new_procedure("using-mu?-aux", (SCM (*)()) using_mup_aux, 0, 0, 0);
gh_new_procedure("init-params-aux", (SCM (*)()) init_params_aux, 2, 0, 0);
gh_new_procedure("square-basis-aux", (SCM (*)()) square_basis_aux, 2, 0, 0);
gh_new_procedure("range-overlap-with-object-aux", (SCM (*)()) range_overlap_with_object_aux, 5, 0, 0);
gh_new_procedure("display-geometric-object-info-aux", (SCM (*)()) display_geometric_object_info_aux, 2, 0, 0);
gh_new_procedure("point-in-periodic-object?-aux", (SCM (*)()) point_in_periodic_objectp_aux, 2, 0, 0);
gh_new_procedure("normal-to-object-aux", (SCM (*)()) normal_to_object_aux, 2, 0, 0);
gh_new_procedure("point-in-object?-aux", (SCM (*)()) point_in_objectp_aux, 2, 0, 0);
}

