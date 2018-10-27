; Photonic crystal slab consisting of a triangular lattice of air
; holes in a finite-thickness dielectric slab, optionally with a
; substrate on one side of the slab.  See the paper: S. G. Johnson,
; S. Fan, P. R. Villeneuve, J. D. Joannopoulos, L. A. Kolodziejski,
; "Guided modes in photonic crystal slabs," PRB 60, 5751 (August
; 1999).

; Note that this structure has mirror symmetry throught the z=0 plane,
; and we are looking at k-vectors in the xy plane only.  Thus, we can
; break up the modes into even and odd (analogous to TE and TM), using
; the run-zeven and run-zodd functions.

(define-param h 0.5) ; the thickness of the slab
(define-param eps 12.0) ; the dielectric constant of the slab
(define-param loweps 1.0) ; the dielectric constant of the substrate
(define-param r 0.3) ; the radius of the holes
(define-param supercell-h 4) ; height of the supercell

; triangular lattice with vertical supercell:
(set! geometry-lattice (make lattice (size 1 1 supercell-h)
                         (basis1 (/ (sqrt 3) 2) 0.5)
                         (basis2 (/ (sqrt 3) 2) -0.5)))

(set! geometry
      (list (make block (material (make dielectric (epsilon loweps)))
		  (center 0 0 (* 0.25 supercell-h))
		  (size infinity infinity (* 0.5 supercell-h)))
	    (make block (material (make dielectric (epsilon eps)))
		  (center 0) (size infinity infinity h))
	    (make cylinder (material air)
		  (center 0) (radius r) (height supercell-h))))

; 1st Brillouin zone of a triangular lattice:
(define Gamma (vector3 0 0 0))
(define M (vector3 0 0.5 0))
(define K (vector3 (/ -3) (/ 3) 0))

(define-param only-K false) ; run with only-K=true to only do this k-point
(define-param k-interp 4)   ; the number of k points to interpolate
(if only-K
    (set! k-points (list K))
    (set! k-points (interpolate k-interp (list Gamma M K Gamma))))

(set-param! resolution (vector3 32 32 16))
(set-param! num-bands 9)

; Run even and odd bands, outputting fields only at the K point:
(if (= loweps 1.0)
    (begin ; we only have even/odd classification for symmetric structure
      (run-zeven (output-at-kpoint K output-hfield-z))
      (run-zodd (output-at-kpoint K output-dfield-z)))
    (run (output-at-kpoint K output-hfield-z) display-zparities))

(display-eigensolver-stats)
