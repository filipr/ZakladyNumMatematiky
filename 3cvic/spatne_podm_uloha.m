%% Priklad spatne podminene ulohy

A = [1 2; 2 3.999];
b = [4; 7.999];

x_reseni = A\b % reseni soustavy vestavenym resicem

chyba1 = norm(b - A*x_reseni) % kontrola, je to skutecne reseni

%% perturbace prave strany
b_perturbed = [4.001; 7.998];

x2 = A\b_perturbed

chyba2 = norm(b_perturbed - A*x2) % kontrola, je to skutecne reseni

%% perturbace matice

A_perturbed = [1 2.001; 2 3.998];

x3 = A_perturbed\b

chyba3 = norm(b - A_perturbed*x3) % kontrola, je to skutecne reseni

%% cislo podminenosti matice A
podminenost_matice_A = cond(A)