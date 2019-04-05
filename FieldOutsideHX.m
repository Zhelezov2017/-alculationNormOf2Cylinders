function FieldOutsideHX = FieldOutsideHX(N_max, Const,k_0, p, x0, y0, x, y)
m  = [-N_max:N_max];
if(N_max == 0) m = 0; end
q = sqrt(1 - p.^2);
q = q.* (2*(imag(q) <= 0)-1);

FieldOutsideHPHI = zeros(1001,1001);
FieldOutsideHRHO = zeros(1001,1001);
for ii = 1:1001
    for ij = 1:1001
        dx = x(ii)-x0;
        dy = y(ij)-y0;
        rho = sqrt(dx.^2 + dy.^2);
        phi1 = atan(abs(dy)/abs(dx));
        phi    = (phi1).* (dx>0).* (dy>0); 
        phi    = phi + (pi - phi1).* (dx<=0).* (dy>0); 
        phi    = phi + (pi + phi1).* (dx<=0).* (dy<=0); 
        phi    = phi + (2*pi - phi1).* (dx>0).* (dy<=0);

        Q = k_0.* rho * q;
    
        H2m  = besselh(N_max, 2, Q);
        dH2m = (H2m.* N_max)./ Q  - besselh(N_max + 1, 2, Q); 
    
        FieldOutsideHPHI(ii,ij) = -Const(4)*dH2m+Const(3)*p*N_max*H2m/Q;
        FieldOutsideHRHO(ii,ij) = 1i*p*Const(3)*dH2m+1i*Const(4)*N_max*H2m/Q;
        FieldOutsideHX(ii,ij) = cos(phi)*FieldOutsideHRHO(ii,ij) - sin(phi)*FieldOutsideHPHI(ii,ij);
    end
end


