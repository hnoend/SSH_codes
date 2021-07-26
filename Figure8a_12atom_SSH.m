clear;
n = 12;  %number of atoms

sizep = 100;

vp = linspace(0,4,sizep)';
E = zeros(n,sizep);
w = 1;

for iter = 1:1:sizep
  v = vp(iter);
  H = zeros(n,n);

  for m = 1:2:n
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
  %[V,E] = eig(H);
E(:,iter) = eig(H); % exctracting the energy for that particular value of v and w
end

figure;
hold on;

for iter = 1:1:n
  plot(vp,E(iter,:));
end

xlabel('v','Fontsize',16);
ylabel('Energy','Fontsize',16);
