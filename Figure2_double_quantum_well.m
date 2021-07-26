clear;
hcut = (6.64e-34)/(2*pi); %defining constants
me = 9.1e-31;qe = 1.6e-19; beta = 1e-9;
a = 1;    % length of well in nm
b = 0.1;  % distance between wells in nm
E1 = hcut*hcut*pi*pi/(2*me*qe*beta*beta*a*a); % getting energy of ground state of a 1D infinite box of length a
V = 3*E1;   % defining potential of finite well with respect to that energy

g = 2*me*qe*V*beta*beta/(hcut*hcut); % defining k^2 + p^2 = 2mV/(hbar)^2 = g


pt = fe_solve(a,b,g); % solving for the p that gives ground state energy
kt = sqrt(g-pt*pt); %obtaining the corresponding k

%plotting the wavefunction

%getting the coefficients first assuming B is 100 and then normalizing later
B = 100;
term5 = kt*cos(pt*b/2)*(exp(kt*b/2)-exp(-kt*b/2))+pt*sin(pt*b/2)*(exp(kt*b/2)+exp(-kt*b/2));
term6 = kt*sin(pt*b/2)*(exp(kt*b/2)-exp(-kt*b/2))-pt*cos(pt*b/2)*(exp(kt*b/2)+exp(-kt*b/2));
C = -B*term5/term6;
A = (B*cos(pt*(a+b/2)) + C*sin(pt*(a+b/2)))/exp(-kt*(a+b/2));
D = (B*cos(pt*b/2) + C*sin(pt*b/2))/(exp(+kt*b/2)+exp(-kt*b/2));

sizp = 100*(a+b)/max(min(a,b),0.1);
sizp = 1000;
x = linspace(-5,5,sizp);x = x'; %plotting over a distane of -5 to 5nm
dx = x(3)-x(2);
psi = zeros(sizp,1);
for loop = 1:1:sizp
  if(x(loop)<-a-b/2)
    psi(loop) = A*exp(kt*x(loop));
  elseif(x(loop)>=-a-b/2 && x(loop)<-b/2)
    psi(loop) = B*cos(pt*x(loop)) -C*sin(pt*x(loop));
  elseif(x(loop)>=-b/2 && x(loop)<b/2)
    psi(loop) = D*(exp(kt*x(loop))+exp(-kt*x(loop)));
  elseif(x(loop)>=b/2 && x(loop)<a+b/2)
    psi(loop) = B*cos(pt*x(loop)) +C*sin(pt*x(loop));
  else
    psi(loop) = A*exp(-kt*x(loop));
  end
end

norm = psi'*psi*dx;
psi = psi/sqrt(norm); %normalizing the wavefunction

figure;
plot(x,psi);
xlabel('x in nm','Fontsize',16);
ylabel('\Psi_0','Fontsize',16);
