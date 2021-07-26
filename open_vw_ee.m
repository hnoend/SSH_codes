function entropy = open_vw_ee(v,w,n,meow,filling)
% v,w are the hopping parameters
% n is the number of atoms
% meow is the site index where one wants to partition, meow = n/2 for halfway partition
% filling is n/2 for half-filled system

H = zeros(n,n);
Corrs2 = zeros(n,n);
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
 % uncomment for periodic system
 % H(1,n) = w';
 % H(n,1) = w;
  
[V,E] = eig(H);

Corrs2 = V(1:n,1:filling)*V(1:n,1:filling)';

%truncating the correlator matrix
subcorrs = Corrs2(1:meow,1:meow);
ck = eig(subcorrs);


entropy = 0;

for m = 1:1:meow
  if(ck(m) !=0 && ck(m)!=1)
      entropy = entropy - ck(m)*log(ck(m)) - (1-ck(m))*log(1-ck(m));
  endif
endfor

end