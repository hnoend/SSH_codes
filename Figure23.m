clear all;
clc;
close all;

% parameters of SSH chain, change w and v for different topological phases
w=0.5;v=0.5;
%change units to change length of chain, grid number corresponds to energy
%matrix for NEGF computation, the higher the finer are computations
unit = 16; grid = 100000;
%other parameters initialized
hbar=1;q=1;kT=0.002;

a=1;
%define SSH hamiltonian
H = kron(eye(unit),v.*[0,1;1,0]);
c1 =3;
Np =2*unit;
for j = 2:2:(unit*2 -1)
    H(j,c1) = w;
    H(j+1,c1-1)=w;
    c1= c1+2;
end
energy_H = eig(H);

%the range of energy depends on the values of w and v, we would want to
%cover the entire spectrum possible
E=linspace(-1.1,1.1,grid);
dE=E(2)-E(1);
zplus=1i*1e-8;

%the negf 
for loop1=1:length(E)
    
    gamma = 0.05; % set it to w/10   
    sig1 = zeros(Np, Np); sig2 = zeros(Np, Np);
    sig1(1,1) = -1i*0.5*gamma;
    sig2(Np,Np)=-1i*0.5*gamma;
    Gamma_1 = 1i*(sig1 - sig1');
    Gamma_2 = 1i*(sig2 - sig2');
    %G 
    G=inv(((E(loop1)+zplus)*eye(Np))-H-sig1 - sig2);
    Ai=1j*(G-G');
    DOS = trace(Ai);
    ALDOS1(:,loop1)=diag(Ai)./(2*pi);
end
for j = 1:unit*2
    ALDOS1(j,:) = ALDOS1(j,:)./sum(ALDOS1(j,:));
end

%plotting LDOS

%setting floor values to 1E-7 to avoid computational discrepancies
indices = find(abs(ALDOS1)<1E-7);
ALDOS1(indices) = 1E-7;

ALDOS2 = log(ALDOS1);
[x,y] = meshgrid(E,1:unit*2);
pcolor(y,x,ALDOS2);
colormap(hot);
colorbar;
shading interp
box on
set(gca,'FontSize',14);
set(gca,'Fontsize',[14]);
ylabel('E','Fontsize',14);
xlabel('Site Number','Fontsize',14);
title('Local DOS','Fontsize',14);
