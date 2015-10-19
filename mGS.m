
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
