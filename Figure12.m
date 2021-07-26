clear;

sizep = 100;

phi = linspace(0,2*pi,sizep)';
E = zeros(4,sizep);

for iter = 1:1:sizep
  t = exp(i*phi(iter)/4);
  H = zeros(4,4);
  H(1,2) = t; H(1,4) = t';  % atom 1 is connected to 2 and 4. t' according to Peierls subsitution. H is still Hermitian
  H(2,3) = t; H(2,1) = t';
  H(3,4) = t; H(3,2) = t';
  H(4,1) = t; H(4,3) = t';
  
E(:,iter) = eig(H); % exctracting the energy for that particular value of v and w
end

figure;
hold on;

for iter = 1:1:4
  plot(phi,E(iter,:));
end

xlabel('\Phi_0','Fontsize',16);
ylabel('Energy Eigenvalues','Fontsize',16);
