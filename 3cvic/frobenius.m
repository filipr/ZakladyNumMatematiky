function y = frobenius(A)  

y = 0;
[m,n] = size(A);

for i = 1:m 
    for j = 1:n 
        y =  y + A(i,j).^2; 
    end

end

y = sqrt(y);
end