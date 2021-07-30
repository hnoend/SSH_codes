clear all;clc;
% parameters of SSH chain, change w and v for different topological phases
unit = 5;
v=1;w=0.5;

%setting up the Hamiltonian of the SSH chain
H1 = kron(eye(unit),v.*[0,1;1,0]);
c1 =3;
Np =2*unit;
for j = 2:2:(unit*2 -1)
    H1(j,c1) = w;
    H1(j+1,c1-1)=w;
    c1= c1+2;
end

%set another set of v and w
v=0.5;w=1;
H2 = kron(eye(unit),v.*[0,1;1,0]);
c1 =3;
Np =2*unit;
for j = 2:2:(unit*2 -1)
    H2(j,c1) = w;
    H2(j+1,c1-1)=w;
    c1= c1+2;
end

%merge of PA-alpha and PA-beta
H_merge1 = zeros(20,20);
H_merge1(1:10,1:10)=H1;
H_merge1(11:20,11:20) = H2;
H_merge1(10,11)=1;H_merge1(11,10)=1; % H_merge1(11,12)=1;H_merge1(12,11)=1;
% H_merge1(22,21)=1;H_merge1(21,22)=1;
H_merge1(1,20)=0.5;H_merge1(20,1)=0.5;

H_1_eig = eig(H_merge1)


%merge of PA-alpha and PA-alpha
H_merge2 = zeros(20,20);
H_merge2(1:10,1:10)=H1;
H_merge2(11:20,11:20) = H1;
H_merge2(10,11)=0.5;H_merge2(11,10)=0.5;%H_merge2(11,12)=0.5;H_merge2(12,11)=0.5;
% H_merge2(22,21)=0.5;H_merge2(21,22)=0.5;H_merge2(1,22)=0.5;H_merge2(22,1)=0.5;
H_merge2(1,20)=0.5;H_merge2(20,1)=0.5;
H_2_eig = eig(H_merge2)

%eigenvalues in H_1_eig and H_2_eig can now be plotted

