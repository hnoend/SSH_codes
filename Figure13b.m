clear;

v = 2.5; % setting v/w = 0.25, change this for plots in Figure13c as well
w = 10;

L = 8 ; % number of unit cells, number of atoms is 2L

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
  
  h1 = plot(kc,angc,'r');
  set(h1,'linewidth',[1.5]);
  set(h1,'marker','o');
  set(h1,'markersize',7);
  
  hold on;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
v = 5;
w = 10;
 
%berry_ac = [];
%w_ac = [];
%for v = 1:0.12:10

L = 8;

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
  overlap = 0.5*(r*rnext + t*tnext);
  angc = [angc, imag(log(overlap))];%*180/pi];
  berry2 = berry2*overlap;
  kc = [kc,k];
  end
  
  h2 = plot(kc,angc,'g');
  set(h2,'linewidth',[1.5]);
  set(h2,'marker','o');
  set(h2,'markersize',7);
  
  hold on;


  
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 v = 7.5;
w = 10;

L = 8;

kc = [];

angc = [];
berry2 = 1;
for m = 1:1:L
  k = 2*pi*m/L -pi;
  knext = 2*pi*(m-1)/L -pi;
  r = sqrt(v*v + w*w + 2*v*w*cos(k));
  rnext = sqrt(v*v + w*w + 2*v*w*cos(knext));
  
  t = v + w*exp(i*k);
  tnext = v + w*exp(-i*knext);
  overlap = 0.5*(r*rnext + t*tnext);
  angc = [angc, imag(log(overlap))];%*180/pi];
  berry2 = berry2*overlap;
  kc = [kc,k];
  end
  
  h3 = plot(kc,angc,'b');
  set(h3,'linewidth',[1.5]);
  set(h3,'marker','o');
  set(h3,'markersize',7);
  
  hold on;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  xlabel('ka \in (-\pi,\pi]','Fontsize', 16);
  ylabel('Overlap angle','Fontsize', 16);
  
  
 
  bruh = legend('v/w 0.25','v/w = 0.5','v/w = 0.75');
  set(bruh,'Fontsize',14);
  set(bruh,'Location','best');
  