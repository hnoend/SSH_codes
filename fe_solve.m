function root2 = fe_solve(a,b,g)
  %using bisection method to obtain the ground state energy as the root
  p1 = 0.1; p2 = p1; %observation that ground state energy is always more than 0.1
  while(fev(a,b,g,p2) >= 0) % obtaining right end starting point for bisection method
      p2 = p2+0.2;
  end
  while(abs(p1-p2)>0.001)  % error tolerance to stop
      pn = (p1+p2)/2;
      if(fev(a,b,g,pn) >= 0)
          p1 = pn;
      else
          p2 = pn;
      end
  end
  root2 = (p1+p2)/2;   
end
