clear;

v = 99; % setting v/w = 0.99
w = 100;

L = 100 ; % number of unit cells, number of atoms is 2L

kc = [];

angc = [];
berry2 = 1;
for m = 1:1:L
  k = 2*pi*m/L-pi;
  knext = 2*pi*(m-1)/L - pi;
  r = sqrt(v*v + w*w + 2*v*w*cos(k));
  rnext = sqrt(v*v + w*w + 2*v*w*cos(knext));
  
  t = v + w*exp(i*k);
  tnext = v + w*exp(-i*knext);
  overlap = 0.5*(r*rnext + t*tnext);  %obtaining the overlap angle for each k
  angc = [angc, imag(log(overlap))];%*180/pi];
  berry2 = berry2*overlap;
  kc = [kc,k];
end 
 

figure;
  h1 = plot(kc,angc,'r');
  xlabel('ka \in (-\pi,\pi]','Fontsize', 16);
  ylabel('Overlap angle','Fontsize', 16);
  