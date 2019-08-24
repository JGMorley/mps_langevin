# mps-miscellany
Misc MPS-related functions written in MATLAB.

---
Convention for storing an MPS
---
When working with non-translationally invariant or Vidal form MPS we need at 
least two tensors to specify our MPS. In this case we store an infinite, uniform 
MPS like

  mps = {FORM, T1, T2, ...,TN}

FORM indicates whether we're in Vidal, left-canonical or right-canonical form, 
taking for these cases the values 0, 1 and 2 respectively. The number of 
elements of mps then determines the size of the unit cell.


If FORM = 3 then the tensors T1,...,TN should be
  Gamma1, lambda1, ..., GammaM, lambdaM
for a uniform infinite MPS with periodicity over M sites.

If FORM = 1 then the tensors T1,...,TN should be
  A1,...,AM,Lambda
for a uniform infinite MPS with periodicity over M sites, where 
  Lambda = lambda.'*conj(lambda) (Lambda indexed from top to bottom of drawing)

If FORM = 2 then the tensors T1,...,TN should be 
  A1,...,AM,Lambda
for a uniform infinite MPS with periodicity over M sites, where
  Lambda = lambda*lambda' (Lambda indexed from top to bottom of drawing)

If we are using translationally invariant, non-Vidal form MPS then we can just
use a single tensor A, and the normalizeMPS function to find the left- and 
right- environment.
---

---
Normalization convention
---
> When in canonical form:
I. Schmidt spectrum at each bond index should be normalized
     ie trace(lambda*lambda') = 1
   This is achieved by setting
     lambda = lambda/sqrt(trace(lambda*lambda'))

> Otherwise:
I. Left and right environments should have a mutual overlap of unity

> In either case:
II. Eigenvalues in canonical form conditions should be unity.
    This is achieved here by setting 
      A -> A/sqrt(trace(EI)/D) for right canonical form A,
      A -> A/sqrt(trace(IE)/D) for left canonical form A, and
      G -> G/sqrt(trace(RI)/D) for Vidal form G(amma), 
    where E, R are the usual transfer matrices and I the 'identity' matrix
    used in the canonical form conditions in eg. DOI: 10.1103/PhysRevB.78.155117

Taken together these imply that the state is normalized in the conventional 
sense, ie. it has norm equal to 1.
---