function FieldOutsideEX = FieldOutsideEX(N_max, Const,k_0, p, x0, y0, x, y)
m  = [-N_max:N_max];
if(N_max == 0) m = 0; end
q = sqrt(1 - p.^2);
q = q.* (2*(imag(q) <= 0)-1);

dx = x - x0;
dy = y - y0;
rho = sqrt(dx.^2 + dy.^2);
phi1 = atan(abs(dy)/abs(dx));
phi    = (phi1).* (dx>0).* (dy>0); 
phi    = phi + (pi - phi1).* (dx<=0).* (dy>0); 
phi    = phi + (pi + phi1).* (dx<=0).* (dy<=0); 
phi    = phi + (2*pi - phi1).* (dx>0).* (dy<=0);

    Q = k_0.* rho * q;
    
    H2m  = besselh(N_max, 2, Q);
    dH2m = (H2m.* N_max)./ Q  - besselh(N_max + 1, 2, Q); 
    
    
    FieldOutsideEPHI =-1i*Const(3)*dH2m-1i*Const(4)*p*N_max*H2m/Q;
    FieldOutsideERHO = p*Const(4)*dH2m+Const(3)*N_max*H2m/Q;
    FieldOutsideEX = FieldOutsideERHO.*cos(phi)- FieldOutsideEPHI.*sin(phi);
end



