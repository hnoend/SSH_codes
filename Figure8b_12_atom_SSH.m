clear;

n = 12;  %number of atoms

v = 0.1;
w = 1;

  H = zeros(n,n);  

  for m = 1:2:n   %obtaining the 12x12 Hamiltonian
    H(m,m+1) = v;
    H(m+1,m) = v';
    if(m<n-1)
      H(m+1,m+2) = w;
    end
    if(m>1)
      H(m,m-1) = w';
    end
  end
  %uncomment this for periodic boundary condition
  %H(1,n) = w';
  %H(n,1) = w;
  
[V,E] = eig(H);   %obtaining the eigenvalues
E = eig(H);   % check what these energies are

figure;
hold
plot(V(:,1).*V(:,1),'r');
plot(V(:,6).*V(:,6),'b');
xlabel('Site number','Fontsize',16);
ylabel('|\Psi|^2','Fontsize',16);