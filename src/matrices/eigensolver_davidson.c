/* Copyright (C) 1999-2014 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* This file contains an alternative eigensolver, currently experimental,
   based on the Davidson method (a preconditioned variant of Lanczos):

   M. Crouzeix, B. Philippe, and M. Sadkane, "The Davidson Method,"
   SIAM J. Sci. Comput. 15, no. 1, pp. 62-76 (January 1994). */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "config.h"
#include <mpiglue.h>
#include <mpi_utils.h>
#include <check.h>
#include <scalar.h>
#include <matrices.h>
#include <blasglue.h>

#include "eigensolver.h"

extern void eigensolver_get_eigenvals_aux(evectmatrix Y, real *eigenvals,
                                          evectoperator A, void *Adata,
                                          evectmatrix Work1, evectmatrix Work2,
                                          sqmatrix U, sqmatrix Usqrt,
                                          sqmatrix Uwork);

#define STRINGIZEx(x) #x /* a hack so that we can stringize macro values */
#define STRINGIZE(x) STRINGIZEx(x)

/**************************************************************************/

#define EIGENSOLVER_MAX_ITERATIONS 100000
#define FEEDBACK_TIME 4.0 /* elapsed time before we print progress feedback */

/**************************************************************************/

void eigensolver_davidson(evectmatrix Y, real *eigenvals,
			  evectoperator A, void *Adata,
			  evectpreconditioner K, void *Kdata,
			  evectconstraint constraint, void *constraint_data,
			  evectmatrix Work[], int nWork,
			  real tolerance, int *num_iterations,
			  int flags,
			  real target)
{
     int nbasis, q;
     evectmatrix *AV, *V;
	 //evectmatrix *AV_z, *V_z;  //Work matrix for D field corresponding to AV and V
     sqmatrix VAV, S, Swork, U, S2, S3, I;
	 //sqmatrix S2_z, S3_z;
     mpiglue_clock_t prev_feedback_time;
     int iteration = 0, ibasis = 0;
     real *eigenvals2, prev_E = 0;

     prev_feedback_time = MPIGLUE_CLOCK;
     
#ifdef DEBUG
     flags |= EIGS_VERBOSE;
#endif

     CHECK(nWork >= 4, "not enough workspace");

	 /*Work are the input matrix array. Each quarters are used to store different data.
	 ndata keep invariant while the size of Work needs to be doubled*/
     nbasis = nWork / 2;
     V = Work;                      //field
     AV = Work + nbasis;			//The result after applying the maxwell operator on field
	 //V_z = Work + 2 * nbasis;		//D field
	 //AV_z = Work + 3 * nbasis;		//The result after applying the maxwell operator on D field

     q = Y.p * nbasis;           //Y.p should be doubled itself. While calculating VAV, S, Swork, we should consider D and H fields together
     VAV = create_sqmatrix(q);
     S = create_sqmatrix(q);
     Swork = create_sqmatrix(q);

     sqmatrix_resize(&VAV, 0, 0);
     sqmatrix_resize(&S, 0, 0);
     sqmatrix_resize(&Swork, 0, 0);

     CHK_MALLOC(eigenvals2, real, q);

     U = create_sqmatrix(Y.p);
     S2 = create_sqmatrix(Y.p);
     S3 = create_sqmatrix(Y.p);

     //U_z = create_sqmatrix(Z.p);
     //S2_z = create_sqmatrix(Z.p);
     //S3_z = create_sqmatrix(Z.p);

     I = create_sqmatrix(0);

     if (constraint)
	  constraint(Y, constraint_data); //Apply each constraints on the H field, which are included in constraint_data
	 //if (constraint)
	  //constraint(Z, constraint_data);  //Assume D field is constrained in the same way as H field

     evectmatrix_XtX(U, Y, S3);   //U=YtY, S3 is a scratch matrix
     CHECK(sqmatrix_invert(U, 1, S3), "singular YtY at start");
     sqmatrix_sqrt(S2, U, S3); /* S2 = 1/sqrt(Yt*Y) */ //in the sense that the eigenvalues of U are square rooted. It seems U has been overwritten as the eigenvectors of the original U.
     evectmatrix_XeYS(V[0], Y, S2, 1); /* V[0] = orthonormalize Y */ //V[0]=YS2, and "1" means S2 is assumed to be Hermitian

	 //evectmatrix_XtX(U_z, Z, S3_z);   //U_z=ZtZ, S3_z is a scratch matrix
     //CHECK(sqmatrix_invert(U_z, 1, S3_z), "singular ZtZ at start");
     //sqmatrix_sqrt(S2_z, U_z, S3_z); /* S2 = 1/sqrt(Zt*Z) */ //in the sense that the eigenvalues of U_z are square rooted. It seems U_z has been overwritten as the eigenvectors of the original U_z.
     //evectmatrix_XeYS(V_z[0], Z, S2_z, 1); /* V_z[0] = orthonormalize Z */ //V[0]=YS2, and "1" means S2 is assumed to be Hermitian

     do {
	  real E;
	  int itarget, i;

	  //A(V[ibasis], AV[ibasis], Adata, 0, Y;)
	  //Apply the maxwell operator on D and H field
	  //Adata is maxwell_data type variable, 0 and Y are not used
	  A(V[ibasis], AV[ibasis], Adata, 0, Y);  

	  q = Y.p * (ibasis + 1);    //Consider D and H fields together
	  sqmatrix_resize(&VAV, q, 1);
	  sqmatrix_resize(&S, q, 0);
	  sqmatrix_resize(&Swork, q, 0);

	  for (i = 0; i <= ibasis; ++i) {
		  //V[i] * AV[ibasis] is stored in VAV at an offset Y.p * (q * i + ibasis) (? offset is the index where it begins)
	       evectmatrixXtY_sub(VAV, Y.p * (q * i + ibasis),
				  V[i], AV[ibasis], S3); 
	  }
	  //Take the upper right half of VAV to S and fill the lower left half with the adjoint of the upper
	  sqmatrix_copy_upper2full(S, VAV); 

	  //Find the eigenvals2 of S, Swork is a work matrix, (? S is overwritten by its eigenvectors)
	  sqmatrix_eigensolve(S, eigenvals2, Swork);

	  /* find index itarget of start of "window" around the
	     target frequency : */
	  if (target == 0.0) /* not attempting targeted eigensolver */
	       itarget = 0;
	  else {
	       /* note that this technique seems to have convergence trouble */
	       for (itarget = 0; itarget + Y.p < q &&
			 fabs(target - eigenvals2[itarget]) >
			 fabs(target - eigenvals2[itarget + Y.p]); ++itarget)
		    ;
	  }

	  for (E = 0.0, i = 0; i < Y.p; ++i) {
	       E += (eigenvals[i] = eigenvals2[itarget + i]);
	  }
	  mpi_assert_equal(E);

	  /* compute Y = best eigenvectors */
	  //V[0] is the normalized input Y, (? S seems to be the eigenvector matrix of the original S gotten from VAV)
	  /*If i is nonzero, Y=1*Y+1*V[i]*S, use only V[i].p x V[i].p submatrix, 
	  begining at element indexed (itarget * q + Y.p * i)*/
	  /*If i is zero, Y=1*V[i]*S, use only V[i].p x V[i].p submatrix, 
	  begining at element indexed (itarget * q + Y.p * i)*/
	  for (i = 0; i <= ibasis; ++i) {
	       evectmatrix_aXpbYS_sub(i ? 1.0 : 0.0, Y,
				      1.0, V[i],
				      S, itarget * q + Y.p * i, 1);
	  }
	  
	  if (iteration > 0 && mpi_is_master() &&
	      ((flags & EIGS_VERBOSE) ||
	       MPIGLUE_CLOCK_DIFF(MPIGLUE_CLOCK, prev_feedback_time)
	       > FEEDBACK_TIME)) {
               printf("    iteration %4d: "
                      "trace = %0.16g (%g%% change)\n", iteration, E,
		      200.0 * fabs(E - prev_E) / (fabs(E) + fabs(prev_E)));
               fflush(stdout); /* make sure output appears */
               prev_feedback_time = MPIGLUE_CLOCK; /* reset feedback clock */
          }

	  if (iteration > 0 &&
              fabs(E - prev_E) < tolerance * 0.5 * (fabs(E) + 
						    fabs(prev_E) + 1e-7))
               break; /* convergence!  hooray! */
	  
	  /* compute new directions from residual & update basis: */
	  {
	       int ibasis2 = (ibasis + 1) % nbasis;

	       /* compute V[ibasis2] = AY */
#if 1
	  /*If i is nonzero, V[ibasis2]=1*V[ibasis2]+1*AV[i]*S, use only AV[i].p x AV[i].p submatrix, 
	  begining at element indexed (itarget * q + Y.p * i)*/
	  /*If i is zero, V[ibasis2]=1*AV[i]*S, use only AV[i].p x AV[i].p submatrix, 
	  begining at element indexed (itarget * q + Y.p * i)*/
	       for (i = 0; i <= ibasis; ++i) {
		    evectmatrix_aXpbYS_sub(i ? 1.0 : 0.0, V[ibasis2],
					   1.0, AV[i],
					   S, itarget * q + Y.p * i, 1);
	       }
#else
	       A(Y, V[ibasis2], Adata, 1, Y);
#endif
	  
	       /* handle restart case: */
	       if (ibasis2 == 0) {
		    evectmatrix_copy(AV[0], V[0]); //AV[0]=V[0]
		    evectmatrix_copy(V[0], Y);	   //V[0]=Y
		    sqmatrix_resize(&VAV, Y.p, 0);
		    evectmatrix_XtY(VAV, V[0], AV[0], S3); //VAV=V[0]*AV[0]
		    ibasis2 = 1;
		    evectmatrix_copy(V[ibasis2], AV[0]);  //V[ibasis2]=AV[0]
	       }

	       /* V[ibasis2] = residual = AY - Y * eigenvals */
		   //V[ibasis2].data+=-1 * Y.data *diag(eigenvals), where V[ibasis2].data,Y.data are Y.n x Y.p and eigenvals are real
	       matrix_XpaY_diag_real(V[ibasis2].data, 
				     -1.0, Y.data,
				     eigenvals, Y.n, Y.p);

	       /* AV[ibasis2] = precondition V[ibasis2]: */
	       if (K != NULL)
		    K(V[ibasis2], AV[ibasis2], Kdata, Y, eigenvals, I);
	       else
		    evectmatrix_copy(AV[ibasis2], V[ibasis2]);

	       /* project by the constraints, if any: */
	       if (constraint)
		    constraint(AV[ibasis2], constraint_data);

	       /* orthogonalize against previous V: */
	       for (i = 0; i < ibasis2; ++i) {
		    evectmatrix_XtY(U, V[i], AV[ibasis2], S3);//U=adjoint(V[i])*AV[ibasis2], S3 is a scratch matrix
		    evectmatrix_XpaYS(AV[ibasis2], -1.0, V[i], U, 0);//AV[ibasis2]+=-1.0*V[i]*U
	       }

	       /* orthonormalize within itself: */
	       evectmatrix_XtX(U, AV[ibasis2], S3);//U=adjoint(AV[ibasis2]) * AV[ibasis2]
	       CHECK(sqmatrix_invert(U, 1, S3), "non-independent AV subspace");
	       sqmatrix_sqrt(S2, U, S3);//S2=sqrt(U), S3 is a work matrix
	       evectmatrix_XeYS(V[ibasis2], AV[ibasis2], S2, 1);//V[ibasis2]=AV[ibasis2]*S2

	       ibasis = ibasis2;
	  }

	  prev_E = E;
     } while (++iteration < EIGENSOLVER_MAX_ITERATIONS);

     CHECK(iteration < EIGENSOLVER_MAX_ITERATIONS,
           "failure to converge after "
           STRINGIZE(EIGENSOLVER_MAX_ITERATIONS)
           " iterations");

     evectmatrix_XtX(U, Y, S3);
     CHECK(sqmatrix_invert(U, 1, S3), "singular YtY at end");
     eigensolver_get_eigenvals_aux(Y, eigenvals, A, Adata,
				   V[0], AV[0], U, S3, S2);

     free(eigenvals2);

     destroy_sqmatrix(VAV);
     destroy_sqmatrix(S);
     destroy_sqmatrix(Swork);
     destroy_sqmatrix(U);
     destroy_sqmatrix(S2);
     destroy_sqmatrix(S3);
     destroy_sqmatrix(I);

     *num_iterations = iteration;     
}
