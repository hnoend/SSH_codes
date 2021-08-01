clear all;
close all;
clc;

%m1 is a loop to vary v and w from 0 to 1
m1 = 0.1:0.02:0.9;
%setting fundamental constants
e = 1.6E-19;
h = 6.6E-34;
temp = e*e/h;

%the loops are meant to plot different types of linear conductance 
% v and w are varied via loops
% we can change length of the chain keeping coupling constant to observe
% variation of linear conductance with length 
% we can also change coupling to contacts keeping length constant to observe
% variation of linear conductance with coupling strength 

for loop2 = 1:length(m1) % varies w from 0.1 to 0.9 in steps of 0.02
    loop2
    G1(loop2) = linear_response_G(m1(loop2),1-(m1(loop2)),4,1); % unit size is 10
end
G1 = G1./temp;

for loop2 = 1:length(m1) % varies w from 0.1 to 0.9 in steps of 0.02
    loop2
    G2(loop2) = linear_response_G(m1(loop2),1-(m1(loop2)),8,1); % unit size is 10
end
G2 = G2./temp;
% 
for loop2 = 1:length(m1)% varies w from 0.1 to 0.9 in steps of 0.02
    loop2
    G3(loop2) = linear_response_G(m1(loop2),1-(m1(loop2)),16,1); % unit size is 10
end
G3 = G3./temp;



%the following function calculates linear conductance for given v,w,unit
%and gamma

function [G] = linear_response_G(v,w,unit,gamma)
%energy varies from -Em to +Em, E = sqrt(v^2 + w^2 + 2vwcos(k)) for an SSH
%chain so E_max = -E_min = (v+w). The multiplication factor of 1.1 gives
%some buffer around it
Em = (v+w)*1.1;
%chemical potential is set at 0
mu=0;

%setting other physical constants
hbar=1.06e-34;q=1.6e-19;m=.25*9.1e-31;kT=.025;IE=(q*q)/(2*pi*hbar);
h = 6.62E-34;a=1;

%setting up the Hamiltonian
H = kron(eye(unit),v.*[0,1;1,0]);
c1 =3;
Np =2*unit;
for j = 2:2:(unit*2 -1)
    H(j,c1) = w;
    H(j+1,c1-1)=w;
    c1= c1+2;
end


%Bias
NV=1; % points in IV
%small  variation in voltage to ensure that we are in linear region of
%conductance
VV=linspace(0,.0001,NV); % variation of voltage


zplus=1i*1e-8;
Ef=0; %fermi energy 

for iV=1:NV
V=VV(iV);mu1=Ef+(V/2);mu2=Ef-(V/2);

%Energy grid for Green’s function method
NE = 10000;
E=linspace(-Em,Em,NE);dE=E(2)-E(1);

%Transmission
I=0;%Current
for k=1:NE
    
%     gamma = 0.5;
    sig1 = zeros(Np, Np); sig2 = zeros(Np, Np);
    sig1(1,1) = -1i*0.5*gamma;
    sig2(Np,Np)=-1i*0.5*gamma;
    Gamma_1 = 1i*(sig1 - sig1');
    Gamma_2 = 1i*(sig2 - sig2');
    G=inv(((E(k)+zplus)*eye(Np))-H-sig1-sig2);
    T12=real(trace(Gamma_1*G*Gamma_2*G'));
    f1=1/(1+exp((E(k)-mu1)/kT));
    f2=1/(1+exp((E(k)-mu2)/kT));

    TM(k)=T12;
    I= I+(dE*IE*TM(k)*(f1-f2));
    end
    II(iV)=I; G = I/V;
end

end

%G values can now be plotted against v or w
