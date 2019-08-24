function output = kronlist( matrix_carray )
%KRONLIST Return tensor product of elemnts in matrix_carray
% eg. kronlist{A,B,C}   = kron(A, kron(B, C))
% eg. kronlist{A,B,C,D} = kron(A, kron(B, kron(C, D))) 

number_matrices = size(matrix_carray,2);

output = 1.;
for n=number_matrices:-1:1
    output = kron(matrix_carray{n}, output);
end

end

