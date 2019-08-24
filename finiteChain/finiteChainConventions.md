Conventions for finite chain mps code
-------------------------------------

MPS cell-arrays
---------------

MPS should be stored of cell arrays of the form

{L0, {A1, L1}, {A2, L2}, ...,{AN, LN}},

where the Ak are contracted to form the state and should satisfy

           | [d,Dright]         for  k = 1
size(Ak) = | [Dleft, Dright, d] for  1 < k < N
           | [Dleft, d]         for  k = N,

where Dleft/Dright is the bond dimension to the left/right respectively and d
is the physical dimension for the kth site.

If the MPS is in left canonical form the right environment of the nth site is
given by Ln. If the MPS is in right canonical form the left environment of the 
nth site is given by L(n-1). In both cases Ln is the square of the nth 'lambda'
matrix used in Vidal form.

This convention will also be able to describe Vidal form, by storing the nth
Gamma tensor as An. NB in this case the state is formed by contracting the A
tensors _and_ the L matrices, schematically: psi = GLGLGL...LG.

Hamiltonian tensors
-------------------

For an N-site system the Hamiltonian tensor should have 2N indices such that

H(s'1,...,s'N,s1,...,sN) = <s'1...s'N|\hat{H}|s1...sN>
