clear;

sizet = 100;
delta = linspace(0,1,sizet); % parameter delta
berryp = zeros(sizet,1);  % geometric phase
L = 3;  %number of unit cells, number of atoms is 2L
for iter = 1:1:sizet
  
v = 1-delta(iter); % parameters v and w according to delta
w = delta(iter); 
angc = [];
for m = 1:1:L
                k = 2*pi*m/L;
                knext = 2*pi*(m-1)/L;
                r = sqrt(v*v + w*w + 2*v*w*cos(k));
                rnext = sqrt(v*v + w*w + 2*v*w*cos(knext));
  
                t = v + w*exp(i*k);
                tnext = v + w*exp(-i*knext);
                overlap = 0.5*(r*rnext + t*tnext); % getting overlap angel for each k
                angc = [angc, imag(log(overlap))*180/pi];
            end 
            berry = sum(angc);
          if(round(berry) != 180)
            berry = 0;
          end
berryp(iter) = berry/180; % getting the geometric phase modulo pi

end

figure;
plot(delta,berryp);
xlabel('\delta','Fontsize',16);
ylabel('Berry Phase modulo \pi','Fontsize',16);