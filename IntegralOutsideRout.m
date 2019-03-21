function IntegralOutsideRout = IntegralOutsideRout(N_max, Const, x0, y0, p)
n = 10000;
a = 0;  
b = 100;  
h = (b-a) / n;  
x=a:h:b;
m  = [-N_max:N_max];
if(N_max == 0) m = 0; end
q = sqrt(1 - p.^2);
q = q.* (2*(imag(q) <= 0)-1);
Q = k0.* dx * q;
    
    H2m  = besselh(m, 2, Q);
    dH2m = (H2m.* m)./ Q  - besselh(m + 1, 2, Q); 
    
    FieldOutsideEPHI =-1i*Const(3)*dH2m-1i*Const(4)*p*m*H2m/Q;
    FieldOutsideERHO = p*Const(4)*dH2m+Const(3)*m*H2m/Q;
    FieldOutsideHPHI =-Const(4)*dH2m+Const(3)*p*m*H2m/Q;
    FieldOutsideHRHO = 1i*p*Const(3)*dH2m+1i*Const(4)*m*H2m/Q;
    
IntegralOutsideRout = trapz(x, x.*(FieldOutsideERHO.*FieldOutsideHPHI-FieldOutsideEPHI.*FieldOutsideHRHO));


end

