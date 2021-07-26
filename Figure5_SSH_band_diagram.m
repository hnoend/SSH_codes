clear;
%parameters v and w
v = 0.48;
w = 0.52;

sizep = 100;
E = zeros(sizep,1);
k = linspace(-pi,pi,sizep)';

for iter = 1:1:sizep
  E(iter) = sqrt(v*v+w*w+2*v*w*cos(k(iter))); %obtaining the energy for a periodic SSH chain
end

%plotting
figure;
hold;
plot(k,-E,'b');
plot(k,E,'b--');
xlabel('k (wavenumber)','Fontsize',16);
ylabel('Energy','Fontsize',16);

