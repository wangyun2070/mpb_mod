/* This file was automatically generated --- DO NOT EDIT */

/* For 1d, 2d, and 3d, quadrature points and weights on a unit sphere.
   There are num_sphere_quad[dim-1] points i, with the i-th point at
   (x,y,z) = (sphere_quad[dim-1][i][ 0, 1, 2 ]), and with a quadrature
   weight sphere_quad[dim-1][i][3]. */

#define NQUAD 50
static const int num_sphere_quad[3] = { 2, 12, 50 };

static const double sphere_quad[3][NQUAD][4] = {
    { {1,0,0,0.5}, {-1,0,0,0.5} },
    {
        { 1, 0, 0, 0.083333333333333328707 },
        { -1, 0, 0, 0.083333333333333328707 },
        { 0, 1, 0, 0.083333333333333328707 },
        { 0, -1, 0, 0.083333333333333328707 },
        { -0.49999999999999977796, 0.86602540378443881863, 0, 0.083333333333333328707 },
        { 0.5, -0.86602540378443859659, 0, 0.083333333333333328707 },
        { -0.86602540378443870761, 0.49999999999999994449, 0, 0.083333333333333328707 },
        { 0.86602540378443837454, -0.50000000000000044409, 0, 0.083333333333333328707 },
        { -0.50000000000000044409, -0.86602540378443837454, 0, 0.083333333333333328707 },
        { 0.50000000000000011102, 0.86602540378443859659, 0, 0.083333333333333328707 },
        { 0.86602540378443870761, 0.49999999999999994449, 0, 0.083333333333333328707 },
        { -0.86602540378443881863, -0.49999999999999977796, 0, 0.083333333333333328707 },
    },
    {
        { 0, -1, 0, 0.012698412698412698402 },
        { 0, 1, 0, 0.012698412698412698402 },
        { -1, 0, 0, 0.012698412698412698402 },
        { 1, 0, 0, 0.012698412698412698402 },
        { 0, 0, 1, 0.012698412698412698402 },
        { 0, 0, -1, 0.012698412698412698402 },
        { -0.57735026918962573106, -0.57735026918962573106, -0.57735026918962573106, 0.021093750000000001388 },
        { 0.57735026918962573106, 0.57735026918962573106, 0.57735026918962573106, 0.021093750000000001388 },
        { -0.57735026918962573106, -0.57735026918962573106, 0.57735026918962573106, 0.021093750000000001388 },
        { 0.57735026918962573106, 0.57735026918962573106, -0.57735026918962573106, 0.021093750000000001388 },
        { -0.57735026918962573106, 0.57735026918962573106, -0.57735026918962573106, 0.021093750000000001388 },
        { 0.57735026918962573106, -0.57735026918962573106, 0.57735026918962573106, 0.021093750000000001388 },
        { -0.57735026918962573106, 0.57735026918962573106, 0.57735026918962573106, 0.021093750000000001388 },
        { 0.57735026918962573106, -0.57735026918962573106, -0.57735026918962573106, 0.021093750000000001388 },
        { 0.70710678118654757274, -0.70710678118654757274, 0, 0.022574955908289243145 },
        { -0.70710678118654757274, 0.70710678118654757274, 0, 0.022574955908289243145 },
        { 0.70710678118654757274, 0, 0.70710678118654757274, 0.022574955908289243145 },
        { -0.70710678118654757274, 0, -0.70710678118654757274, 0.022574955908289243145 },
        { 0, -0.70710678118654757274, -0.70710678118654757274, 0.022574955908289243145 },
        { 0, 0.70710678118654757274, 0.70710678118654757274, 0.022574955908289243145 },
        { 0.70710678118654757274, 0, -0.70710678118654757274, 0.022574955908289243145 },
        { -0.70710678118654757274, 0, 0.70710678118654757274, 0.022574955908289243145 },
        { 0, 0.70710678118654757274, -0.70710678118654757274, 0.022574955908289243145 },
        { 0, -0.70710678118654757274, 0.70710678118654757274, 0.022574955908289243145 },
        { -0.70710678118654757274, -0.70710678118654757274, 0, 0.022574955908289243145 },
        { 0.70710678118654757274, 0.70710678118654757274, 0, 0.022574955908289243145 },
        { -0.90453403373329088755, -0.30151134457776362918, -0.30151134457776362918, 0.020173335537918869742 },
        { 0.90453403373329088755, 0.30151134457776362918, 0.30151134457776362918, 0.020173335537918869742 },
        { -0.30151134457776362918, 0.90453403373329088755, -0.30151134457776362918, 0.020173335537918869742 },
        { 0.30151134457776362918, -0.90453403373329088755, 0.30151134457776362918, 0.020173335537918869742 },
        { -0.30151134457776362918, -0.30151134457776362918, 0.90453403373329088755, 0.020173335537918869742 },
        { 0.30151134457776362918, 0.30151134457776362918, -0.90453403373329088755, 0.020173335537918869742 },
        { 0.30151134457776362918, -0.90453403373329088755, -0.30151134457776362918, 0.020173335537918869742 },
        { -0.30151134457776362918, 0.90453403373329088755, 0.30151134457776362918, 0.020173335537918869742 },
        { -0.30151134457776362918, 0.30151134457776362918, -0.90453403373329088755, 0.020173335537918869742 },
        { 0.30151134457776362918, -0.30151134457776362918, 0.90453403373329088755, 0.020173335537918869742 },
        { 0.90453403373329088755, 0.30151134457776362918, -0.30151134457776362918, 0.020173335537918869742 },
        { -0.90453403373329088755, -0.30151134457776362918, 0.30151134457776362918, 0.020173335537918869742 },
        { 0.30151134457776362918, -0.30151134457776362918, -0.90453403373329088755, 0.020173335537918869742 },
        { -0.30151134457776362918, 0.30151134457776362918, 0.90453403373329088755, 0.020173335537918869742 },
        { -0.30151134457776362918, -0.30151134457776362918, -0.90453403373329088755, 0.020173335537918869742 },
        { 0.30151134457776362918, 0.30151134457776362918, 0.90453403373329088755, 0.020173335537918869742 },
        { -0.30151134457776362918, -0.90453403373329088755, 0.30151134457776362918, 0.020173335537918869742 },
        { 0.30151134457776362918, 0.90453403373329088755, -0.30151134457776362918, 0.020173335537918869742 },
        { -0.90453403373329088755, 0.30151134457776362918, 0.30151134457776362918, 0.020173335537918869742 },
        { 0.90453403373329088755, -0.30151134457776362918, -0.30151134457776362918, 0.020173335537918869742 },
        { -0.90453403373329088755, 0.30151134457776362918, -0.30151134457776362918, 0.020173335537918869742 },
        { 0.90453403373329088755, -0.30151134457776362918, 0.30151134457776362918, 0.020173335537918869742 },
        { 0.30151134457776362918, 0.90453403373329088755, 0.30151134457776362918, 0.020173335537918869742 },
        { -0.30151134457776362918, -0.90453403373329088755, -0.30151134457776362918, 0.020173335537918869742 },
    }
};