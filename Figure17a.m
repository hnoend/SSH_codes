clear;
% 4 atoms 2^4 = 16
vec=[zeros(1,15) 1]';kBT=0.000156;beta=0;
H=zeros(16);N=zeros(16);
c1=zeros(16);c2=zeros(16);c3=zeros(16);c4=zeros(16);
%annihilation operators c1,c2,c3,c4 according to Appendix 2
for k1=1:2
for k2=1:2
for k3=1:2
         kp=[1 k3-1 k2-1 k1-1]*[8;4;2;1];
         k=[0 k3-1 k2-1 k1-1]*[8;4;2;1];
         c4(k+1,kp+1)=(-1)^(k1-1+k2-1+k3-1); %k3-1+1
         kp=[k3-1 1 k2-1 k1-1]*[8;4;2;1];
         k=[k3-1 0 k2-1 k1-1]*[8;4;2;1];
         c3(k+1,kp+1)=(-1)^(k1-1+k2-1);
         kp=[k3-1 k2-1 1 k1-1]*[8;4;2;1];
         k=[k3-1 k2-1 0 k1-1]*[8;4;2;1];
         c2(k+1,kp+1)=(-1)^(k1-1);  % k1-1+1
         kp=[k3-1 k2-1 k1-1 1]*[8;4;2;1];
         k=[k3-1 k2-1 k1-1 0]*[8;4;2;1];
         c1(k+1,kp+1)=1;
end
end
end
d1=c1';d2=c2';d3=c3';d4=c4';%creation operators

sizep = 100;
delta = linspace(0,1,sizep)';
entE = zeros(sizep,1);

for iter = 1:1:sizep
v = 1-delta(iter);
w = delta(iter);

%the ssh hamiltonian
H = v*d1*c2+w*d2*c3+v*d3*c4+w*d4*c1 + v*d2*c1+w*d3*c2+v*d4*c3+w*d1*c4; %periodic
%H = v*d1*c2+w*d2*c3+v*d3*c4 + v*d2*c1+w*d3*c2+v*d4*c3; %open condition



[V,E] = eig(H);
E = eig(H);

den1 = V(:,1)*V(:,1)'; %the ground state

entropy = 0;
creo = TrX(den1,2,[4,4]); %taking partial trace over the first half
cka = eig(creo);  %calculating the eigenvalues
for m = 1:1:4
  if(cka(m) !=0 && cka(m)!=1)
      entropy = entropy - cka(m)*log(cka(m));
  end
end


entE(iter) = entropy;

end

figure;
plot(delta,entE);
xlabel('\delta','Fontsize', 16);
ylabel('Entanglement Entropy','Fontsize', 16);

