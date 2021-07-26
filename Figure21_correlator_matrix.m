%run Hubbard_initialize.m before running this code!
sizep = 100;

U = 10;
delta = linspace(0,0.99,sizep);



openee =zeros(sizep,1);
    for iter = 1:1:sizep
      v = 1-delta(iter);
      w = delta(iter);
      
        %Hssh_up_spin   = v*d1*c3+w*d3*c5+v*d5*c7+w*d7*c1 + v*d3*c1+w*d5*c3+v*d7*c5+w*d1*c7; % periodic
        %Hssh_down_spin = v*d2*c4+w*d4*c6+v*d6*c8+w*d8*c2 + v*d4*c2+w*d6*c4+v*d8*c6+w*d2*c8; %periodic
        Hssh_up_spin   = v*d1*c3+w*d3*c5+v*d5*c7+ v*d3*c1+w*d5*c3+v*d7*c5; % open
        Hssh_down_spin = v*d2*c4+w*d4*c6+v*d6*c8+ v*d4*c2+w*d6*c4+v*d8*c6; % open

        Hu = U*d1*c1*d2*c2 + U*d3*c3*d4*c4 + U*d5*c5*d6*c6 + U*d7*c7*d8*c8;

        H = Hssh_up_spin+Hssh_down_spin+Hu;
        H4=H([i4],[i4]);
        [V4,D4]=eig(H4);
        
        gs = zeros(dim,1);
        mat = V4(:,1);
        p = size(mat,1);
        for m = 1:1:p
            ind = i4(m);
            gs(ind) = mat(m);
        end

%calculating all the correlations in the 2^N space        
cor = zeros(n,n);

cor(1,1) = gs'*d1*c1*gs;cor(1,2) = gs'*d1*c2*gs;cor(1,3) = gs'*d1*c3*gs;cor(1,4) = gs'*d1*c4*gs;
cor(1,5) = gs'*d1*c5*gs;cor(1,6) = gs'*d1*c6*gs;cor(1,7) = gs'*d1*c7*gs;cor(1,8) = gs'*d1*c8*gs;

cor(2,1) = gs'*d2*c1*gs;cor(2,2) = gs'*d2*c2*gs;cor(2,3) = gs'*d2*c3*gs;cor(2,4) = gs'*d2*c4*gs;
cor(2,5) = gs'*d2*c5*gs;cor(2,6) = gs'*d2*c6*gs;cor(2,7) = gs'*d2*c7*gs;cor(2,8) = gs'*d2*c8*gs;

cor(3,1) = gs'*d3*c1*gs;cor(3,2) = gs'*d3*c2*gs;cor(3,3) = gs'*d3*c3*gs;cor(3,4) = gs'*d3*c4*gs;
cor(3,5) = gs'*d3*c5*gs;cor(3,6) = gs'*d3*c6*gs;cor(3,7) = gs'*d3*c7*gs;cor(3,8) = gs'*d3*c8*gs;

cor(4,1) = gs'*d4*c1*gs;cor(4,2) = gs'*d4*c2*gs;cor(4,3) = gs'*d4*c3*gs;cor(4,4) = gs'*d4*c4*gs;
cor(4,5) = gs'*d4*c5*gs;cor(4,6) = gs'*d4*c6*gs;cor(4,7) = gs'*d4*c7*gs;cor(4,8) = gs'*d4*c8*gs;

cor(5,1) = gs'*d5*c1*gs;cor(5,2) = gs'*d5*c2*gs;cor(5,3) = gs'*d5*c3*gs;cor(5,4) = gs'*d5*c4*gs;
cor(5,5) = gs'*d5*c5*gs;cor(5,6) = gs'*d5*c6*gs;cor(5,7) = gs'*d5*c7*gs;cor(5,8) = gs'*d5*c8*gs;

cor(6,1) = gs'*d6*c1*gs;cor(6,2) = gs'*d6*c2*gs;cor(6,3) = gs'*d6*c3*gs;cor(6,4) = gs'*d6*c4*gs;
cor(6,5) = gs'*d6*c5*gs;cor(6,6) = gs'*d6*c6*gs;cor(6,7) = gs'*d6*c7*gs;cor(6,8) = gs'*d6*c8*gs;

cor(7,1) = gs'*d7*c1*gs;cor(7,2) = gs'*d7*c2*gs;cor(7,3) = gs'*d7*c3*gs;cor(7,4) = gs'*d7*c4*gs;
cor(7,5) = gs'*d7*c5*gs;cor(7,6) = gs'*d7*c6*gs;cor(7,7) = gs'*d7*c7*gs;cor(7,8) = gs'*d7*c8*gs;

cor(8,1) = gs'*d8*c1*gs;cor(8,2) = gs'*d8*c2*gs;cor(8,3) = gs'*d8*c3*gs;cor(8,4) = gs'*d8*c4*gs;
cor(8,5) = gs'*d8*c5*gs;cor(8,6) = gs'*d8*c6*gs;cor(8,7) = gs'*d8*c7*gs;cor(8,8) = gs'*d8*c8*gs;

subcor = cor(1:4,1:4);


ck = eig(subcor);


centropy = 0;

for m = 1:1:4
  if(ck(m) !=0 && ck(m)!=1)
      centropy = centropy - ck(m)*log(ck(m)) - (1-ck(m))*log(1-ck(m));
  endif
endfor
              
        openee(iter) = real(centropy);
        
end

figure;
plot(delta,openee);
xlabel('\delta','Fontsize', 16);
ylabel('Entanglement Entropy','Fontsize', 16);     