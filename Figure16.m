clear;

sizep = 100;

k = linspace(-pi,pi,sizep);

E = zeros(sizep,2);
lambda = 0.3; % change value of lambda here

for iter = 1:1:sizep
  H1 = [0 exp(-i*k(iter));exp(i*k(iter)) 0];
  H2 = [0 exp(-2*i*k(iter));exp(2*i*k(iter)) 0];
  
  H = lambda*H1 + (1-lambda)*H2;
  
  E(iter,:) = eig(H);
end

figure;
hold on;
plot(k,E(:,1),'r');
plot(k,E(:,2),'b');
xlabel('ka \in (-\pi,\pi]','Fontsize', 16);
ylabel('E','Fontsize', 16);
