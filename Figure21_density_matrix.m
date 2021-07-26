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
        den1 = gs*gs';

        entropy = 0;
        creo = TrX(den1,2,[16,16]);
        cka = eig(creo);
        for m = 1:1:16
              if(cka(m) !=0 && cka(m)!=1)
                    entropy = entropy - cka(m)*log(cka(m));
              endif
        endfor
              
        openee(iter) = real(entropy);
end

%figure;
plot(delta,openee);
xlabel('\delta','Fontsize', 16);
ylabel('Entanglement Entropy','Fontsize', 16);     