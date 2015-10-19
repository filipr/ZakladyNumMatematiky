% Algoritmus QR rozkladu - presne tak jak se nema delat :-) 
% V praxi se vzdy vyhybame sestavovani a NASOBENI prislusnych matic 
% H*A -> nuluje sloupce matice A - tedy Q = H_1' H_2' ... H_{m-1}'
%pridano v GIT
function [Q,R] = hh(A) 


[m,n] = size(A);
Q = eye(m,m);
R = A;

% muzeme eliminovat prvnich (m-1) sloupcu
for i = 1:m-1
    x = R(i:m,i); % not A :-)
    e = zeros(m-i+1,1);
    e(1,1)=1;
    % stabilni volba normaloveho vektoru
    q = x + sign(x(1))*norm(x).*e; 
    q = q ./ norm(q); 
    %dimenze musi odpovidat
    H = eye(m-i+1,m-i+1) - 2.0.*q*q'; % matice HH reflexe

    HH = eye(m,m); 
    % NESTACI oriznuta matice HH reflexe mezi sebou nasobime
    HH(i:m,i:m) = H(1:m-i+1,1:m-i+1);

    Q = Q* HH';   %update

    R = HH*R; %update
    
end

% nemusíme ukládat posledni radek R a posledni sloupec Q 
if m>n  
    R = R(1:n,1:n);
    Q = Q(1:m,1:n);
end

chyba = norm(Q*R-A)

figure(2);
title('Matice R');
surf(R);

end 
