clear;

delta = linspace(0,1,100);

n = 100;
en = zeros(100,1);
delta = delta';
for p = 1:1:100
  en(p) = open_vw_ee(1-delta(p),delta(p),n,n/2,n/2-1); % calculating EE, with one less than half filling
 % using correlator matrix
endfor
en = real(en);

figure;
plot(delta(2:100),en(2:100));
xlabel('\delta','Fontsize', 16);
ylabel('Entanglement Entropy','Fontsize', 16);