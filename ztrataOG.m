%% Ztráta ortogonality v rùzných variantách Gram-Schmidtovy ortogonalizace
function GramSchmidt

n = 100; % velikost matice = pocet vektoru k ortogonalizaci
% vyzkoušejte si i jinou volbu, napríklad n = 50

% výbìr matice pro ilustraci ztráty ortogonality; viz skripta se stromem (3.24)
A = gallery('lauchli',n,1e-7);
cislo_podminenosti = cond(A)

% setrojení "ortogonální" matice pomocí rùzných variant GS; viz níže
Qc  =  cGS(A);
Qm  =  mGS(A);
Qic = icGS(A);
Qhh = householder(A);

% ztráta ortogonality merena Frobeniovou normou
I = eye(n); % jednotkova matice, teoreticky Q'*Q = I
ztrata_cCG = norm(I - Qc'*Qc, 'fro');
ztrata_mCG = norm(I - Qm'*Qm, 'fro');
ztrata_icCG = norm(I - Qic'*Qic, 'fro');

ztrata_hh = norm(I - Qhh'*Qhh, 'fro');

% vykreslení obrázkù zobrazujících ztrátu ortogonality mezi jednotlivými
%    sloupci matice Q
figure(1);
imagesc(abs(I - Qc'*Qc));
colorbar;
title(['Ztráta ortogonality CG-S || I - Q^{T}Q || = ', num2str(ztrata_cCG)]);

figure(2);
imagesc(abs(I - Qm'*Qm));
colorbar;
title(['Ztráta ortogonality MG-S || I - Q^{T}Q || = ', num2str(ztrata_mCG)]);

figure(3);
imagesc(abs(I - Qic'*Qic));
colorbar;
title(['Ztráta ortogonality ICG-S || I - Q^{T}Q || = ', num2str(ztrata_icCG)]);


figure(4);
imagesc(abs(I - Qhh'*Qhh));
colorbar;
title(['Ztráta ortogonality QR pomoci hh odrazu || I - Q^{T}Q || = ', num2str(ztrata_hh)]);

end




% klasický Gram-Schmidt
function [Q, R] = cGS(A)
m = size(A,2); % pocet sloupcù matice A

Q = zeros(size(A));
R = zeros(size(A,2));

% normalizace prvního vektoru
R(1,1) = norm(A(:,1));
Q(:,1) = A(:,1)/R(1,1);

for k = 2:m
    temp = A(:,k); % k-tý sloupec matice A, který bude ortogonalizován
    % následující dva radky jsou efektivním zápisem (v MATLABu)
    %   ortogonalizace vektoru  temp  vuci již ortogonalizovaným vektorum
    coef = Q(:,1:(k-1))'*temp; % vypocet koeficientu q_i^T * a_k
    temp = temp - Q(:,1:(k-1))*coef;    % ortogonalizace vuci predchozim vektorum q_i
    R(1:k-1,k) = coef'; 
    R(k,k) = norm(temp); 
    Q(:,k) = temp/R(k,k);
end
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


% iterovaný (klasický) Gram-Schmidt
function [Q, R] = icGS(A)
m = size(A,2); % pocet sloupcu matice A

iter = 2; % pocet iteraci GS procesu - mely by stacit 2

Q = zeros(size(A));
R = zeros(size(A,2));

R(1,1) = norm(A(:,1));
Q(:,1) = A(:,1)/R(1,1);

for k = 2:m
    temp = A(:,k); % k-tý sloupce matice A, který bude ortogonalizován
    r = 0;
    for ell = 1:iter % proto se ríká tomuto GS iterovaný
        rtemp = Q(:,1:(k-1))'*temp; 
        temp = temp - Q(:,1:(k-1))*rtemp;
        r = r + rtemp;
    end
    R(1:k-1, k) = r';
    R(k,k) = norm(temp);
    Q(:,k) = temp/R(k,k);
end
end

% QR pomoci hh reflexi
function [Q,R] = householder(A) 


[m,n] = size(A);
Q = eye(m,m);
R = A;
%H = eye(n,m);
% muzeme eliminovat prvnich (m-1) sloupcu
for i = 1:m-1
    x = R(i:m,i); % not A :-)
    e = zeros(m-i+1,1);
    e(1,1)=1;
    q = x + sign(x(1))*norm(x).*e; 
    q = q ./ norm(q); 
    %dimenze musí odpovidat
    H = eye(m-i+1,m-i+1) - 2.0.*q*q'; % matice HH reflexe

    HH = eye(m,m); 
    %staci oriznuta matice o posledni radky (stejne se nasobi nulami)
    % NESTACI matice HH reflexi mezi sebou nasobime
    HH(i:m,i:m) = H(1:m-i+1,1:m-i+1);

    Q = Q* HH';   %update

    R = HH*R; %update
    %chyba = norm(Q*R-A);
end

% nemusíme ukládat posledni radek R a posledni sloupec Q 
if m>n  
    R = R(1:n,1:n);
    Q = Q(1:m,1:n);
end

end 
