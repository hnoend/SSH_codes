function fe = fev(a,b,g,p)
  % the equation for even solutions of the well.
k = sqrt(g-p*p);
fe = k*k*(1-exp(-k*b))*sin(p*a) - p*p*(1+exp(-k*b))*sin(p*a) + 2*k*p*cos(p*a);  
end
