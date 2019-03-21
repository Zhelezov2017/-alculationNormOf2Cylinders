function MagneticField = MagneticField(N_max, Const, k_0, p, rho, EE, GG, HH)
m  = [-N_max:N_max];
if(N_max == 0) m = 0; end
q = sqrt(1 - p.^2);
q = q.* (2*(imag(q) <= 0)-1);

    mainq = EE.^2 - GG.^2 + EE.* HH - (HH + EE).* p.^2;
    radq = sqrt((HH - EE).^2 * p.^4 + 2 * ((GG.^2).* (HH + EE) - EE.* (HH - EE).^2) * p.^2 +...
        (EE.^2 - GG.^2 - EE.* HH).^2);
    
    q1 = sqrt(0.5 * (mainq - radq)./ EE);
    q2 = sqrt(0.5 * (mainq + radq)./ EE);

    n1 = -(EE.* (p.* GG).^(-1)).* (p.^2 + q1.^2 + (GG.^2)./ EE - EE);
    n2 = -(EE.* (p * GG).^(-1)).* (p.^2 + q2.^2 + (GG.^2)./ EE - EE);
    
    alp1 = -1 + (p.^2 + q1^2 - EE)./ GG;
    alp2 = -1 + (p.^2 + q2^2 - EE)./ GG;
    
    bet1 = 1 + p./ n1;
    bet2 = 1 + p./ n2;
    
    D1 = q1* k_0;
    D2 = q2* k_0;
    Jm0Q1  = besselj(N_max, D1*rho);
    Jm1Q1  = besselj(N_max+1, D1*rho);
    
    Jm0Q2  = besselj(N_max, D2*rho);
    Jm1Q2  = besselj(N_max+1, D2*rho);
    
    
    MagneticField = -(Const(1)*n1*(Jm1Q1- bet1*N_max*Jm0Q1/(D1*rho))+Const(2)*n2*(Jm1Q2- bet2*N_max* Jm0Q2/(D2*rho)));
   
end