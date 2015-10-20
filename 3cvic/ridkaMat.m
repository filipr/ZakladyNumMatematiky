function ridkaMat
% porovnani zaplneni QR rozkladu (pomoci mG-S) a LU rozkladu (tedy maticoveho zapisu Gaussovy eliminace)

load('west0479.mat')
A = west0479;

[Q,R] = mGS(A); 
[L,U] = lu(A); 

velikost_A = size(A)

pct = 100 / numel(A);


clf; 
figure(1);
spy(A), title('A - ridka symetricka matice')
nz = nnz(A);
xlabel(sprintf('nonzeros=%d (%.3f%%)',nz,nz*pct));

pause;


figure(2);
spy(L), title('Dolni trojuhelnikova matice z LU rozkladu')
nz = nnz(L);
xlabel(sprintf('nonzeros=%d (%.3f%%)',nz,nz*pct));

pause;

figure(3);
spy(Q), title('OG matice z QR rozkladu')
nz = nnz(Q);
xlabel(sprintf('nonzeros=%d (%.3f%%)',nz,nz*pct));

pause;

figure(4);
spy(U), title('Horni trojuhelnikova matice z LU rozkladu')
nz = nnz(U);
xlabel(sprintf('nonzeros=%d (%.3f%%)',nz,nz*pct));

pause;

figure(5);
spy(R), title('Horni trojuhelnikova matice z QR rozkladu')
nz = nnz(R);
xlabel(sprintf('nonzeros=%d (%.3f%%)',nz,nz*pct));

end

% modifikovaný Gram-Schmidt
function [Q, R] = mGS(A)
m = size(A,2); % pocet sloupcù matice A

Q = zeros(size(A));
R = zeros(size(A,2));

R(1,1) = norm(A(:,1));
Q(:,1) = A(:,1)/R(1,1);

for k = 2:m
    temp = A(:,k); % k-tý sloupec matice A, který bude ortogonalizován
    for j = 1:k-1 % v MGS se tomuto for-cyklu nelze, narozdíl od cCG, vyhnout
        R(j,k) = Q(:,j)'*temp;
        temp = temp - R(j,k)*Q(:,j);
    end
    R(k,k) = norm(temp);
    Q(:,k) = temp/R(k,k);
end
end
