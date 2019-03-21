function [H_z_out, dHz_dx_out,  dHz_dy_out, ...
     E_z_out, dEz_dx_out, dEz_dy_out ] = fieldsGyrotropicCrystal2D_eigenwave(x, y, m, a_0, ee, k_0, p, w_0, BD)

global EE GG HH c
%%%% переходим от декатровых координат к цилиндрическим
rho   = sqrt(x.^2 + y.^2);
phi1   = atan(abs(y)./ abs(x));
phi    = (phi1).* (x>0).* (y>0); 
phi    = phi + (pi - phi1).* (x<=0).* (y>0); 
phi    = phi + (pi + phi1).* (x<=0).* (y<=0); 
phi    = phi + (2*pi - phi1).* (x>0).* (y<=0);

ka = k_0 * a_0;

q = sqrt(1 - p.^2);
q = q.* (2*(imag(q) <= 0)-1);
Q = k_0.* a_0 * q;

%%%%%% вычисляем величини для вычисления внутренних коэффициентов
    mainq = EE^2 - GG^2 + EE * HH - (HH + EE) * p.^2;
    radq = sqrt((HH - EE)^2 * p.^4 + 2 * (GG^2 * (HH + EE) - EE * (HH - EE)^2) * p.^2 +...
        (EE^2 - GG^2 - EE * HH)^2);
    q1 = sqrt(0.5 * (mainq - radq) / EE);
    q2 = sqrt(0.5 * (mainq + radq) / EE);

    n1 = -(EE * (p * GG).^(-1)).* (p.^2 + q1.^2 + GG^2 / EE - EE);
    n2 = -(EE * (p * GG).^(-1)).* (p.^2 + q2.^2 + GG^2 / EE - EE);
    
    alp1 = -1 + (p.^2 + q1^2 - EE)./ GG;
    alp2 = -1 + (p.^2 + q2^2 - EE)./ GG;
    
    bet1 = 1 + p./ n1;
    bet2 = 1 + p./ n2;

    Q1 = q1 * ka;
    Q2 = q2 * ka;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_z_in  = 0;
H_z_out = 0;
dHz_dx_out = 0;
dHz_dx_in  = 0;
dHz_dy_out = 0;
dHz_dy_in  = 0;

E_z_in  = 0;
E_z_out = 0;
dEz_dx_out = 0;
dEz_dx_in  = 0;
dEz_dy_out = 0;
dEz_dy_in  = 0;

mMax = size(m,2);
%%% цикл по азимутальному индексу
for jm = 1 : size(m,2)
    
% %%%%%%%%%%%%%%%%%% вычисляем коэффициенты B_1 и B_2 %%%%%%%%%%%%%%%%%%%%%%%
%     JM1    = besselj(m(jm)+1, Q1);
%     JM2    = besselj(m(jm)+1, Q2);
%     Jm1    = besselj(m(jm), Q1);
%     Jm2    = besselj(m(jm), Q2);
%     Jm1_Q1 = Jm1./ Q1;
%     Jm2_Q2 = Jm2./ Q2;
%     
%     Ez1    = 1i/HH * n1 * q1 * Jm1;
%     Ez2    = 1i/HH * n2 * q2 * Jm2;
%     Ephi1  = 1i * (JM1 + alp1 * m(jm) * Jm1_Q1);
%     Ephi2  = 1i * (JM2 + alp2 * m(jm) * Jm2_Q2);
%     Hz1    = - q1 * Jm1;
%     Hz2    = - q2 * Jm2;
%     Hphi1  = - n1 * (JM1 - bet1 * m(jm) * Jm1_Q1);
%     Hphi2  = - n2 * (JM2 - bet2 * m(jm) * Jm2_Q2);
% 
%     
%     
%     
%     H2m  = besselh(m(jm), 2, Q);
%     dH2m = (H2m * m(jm))./ Q  - besselh(m(jm) + 1, 2, Q);    
%     
%     Jm_out  = besselj(m(jm), Q);
%     dJm_out = Jm_out * m(jm) / Q  - besselj(m(jm) + 1, Q);
%     
%     Ez_inc  =     Jm_out;
%     Ez_sct  = q.* H2m;    
%     Hz_inc  =     Jm_out;
%     Hz_sct  = q.* H2m; 
%     
%     A = 1./ (k_0 * (1 - p.^2));
%     Ephi_incE  = -    (A.* ((p.* (m(jm)./a_0)).* Jm_out));
%     Ephi_sctE  = -q.* (A.* ((p.* (m(jm)./a_0)).* H2m));
%     Ephi_incH  = -    (A.* (- 1i * k_0 * q.* dJm_out));
%     Ephi_sctH  = -q.* (A.* (- 1i * k_0 * q.* dH2m));
%     
%     Hphi_incH  = -    (A.* ((p.* (m(jm)./a_0)).* Jm_out));
%     Hphi_sctH  = -q.* (A.* ((p.* (m(jm)./a_0)).* H2m));
%     Hphi_incE  = -    (A.* (  1i * k_0 * q.* dJm_out));
%     Hphi_sctE  = -q.* (A.* (  1i * k_0 * q.* dH2m));
% 
% %%%% граничные условия
% mat_boundConditions = [ Ez1      Ez2     -Ez_sct         0;                 %%% Ez
%                         Hz1      Hz2     0              -Hz_sct;            %%% Hz
%                         Ephi1    Ephi2   -Ephi_sctE     -Ephi_sctH;         %%% Ephi
%                         Hphi1    Hphi2   -Hphi_sctE     -Hphi_sctH];        %%% Ephi
% 
% [S_EE, S_HE, S_EH, S_HH] = funSm_TE_freeIncidence(m(jm), w_0, a_0, ee, p, EE, GG, HH, c);
% CE = S_EE * D(jm) + S_HE * D(jm + mMax);
% CH = S_EH * D(jm) + S_HH * D(jm + mMax);
% rightSide = [CE * Ez_inc;
%              CH * Hz_inc; 
%              CE * Ephi_incE + CH * Ephi_incH; 
%              CE * Hphi_incE + CH * Hphi_incH];
% 
% S = mat_boundConditions \ rightSide;
B1 = BD((jm - 1) * 4 + 1);
B2 = BD((jm - 1) * 4 + 2);
DE = BD((jm - 1) * 4 + 3);
DH = BD((jm - 1) * 4 + 4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jm1     = besselj(m(jm), q1.* k_0 * rho.* (abs(rho) < a_0));
Jm2     = besselj(m(jm), q2.* k_0 * rho.* (abs(rho) < a_0));
H2m    = besselh(m(jm), 2, k_0 * q * rho);
H2M    = besselh(m(jm) + 1,2, k_0 * q * rho);
J2m_in1 = besselj(m(jm), (q1 * k_0) * rho.* (abs(rho) < a_0));
J2m_in2 = besselj(m(jm), (q2 * k_0) * rho.* (abs(rho) < a_0));
J2M_in1 = besselj(m(jm) + 1, (q1 * k_0) * rho.* (abs(rho) < a_0));
J2M_in2 = besselj(m(jm) + 1, (q2 * k_0) * rho.* (abs(rho) < a_0));

%%% вычисление поля снаружи в внутри одновременно 

    %%%% вычисление электрического поля и его производных 
    E_z_in    = E_z_in + 1i/HH * (n1 * q1 * Jm1 * B1 + n2 * q2 * Jm2 * B2).*...
        exp(-1i * m(jm) * phi).* (abs(rho) < a_0);
    E_z_out = E_z_out + q.* (DE.* H2m).* exp(-1i * m(jm) * phi);
    
    dEz_dx_out = dEz_dx_out + DE.* (((H2m.* m(jm)./ (rho.^2)).* (x + 1i * y) -...
        k_0 * H2M.* x./ rho)).* exp(-1i * m(jm) * phi);
    dEz_dx_in =  dEz_dx_in  + 1i/HH * (B1 * n1 * q1 * (((J2m_in1.* m(jm)./ (rho.^2)).* (x + 1i * y) -...
         q1.* k_0 * J2M_in1.* x./ rho)) +...
         B2 * n2 * q2 * (((J2m_in2.* m(jm)./ (rho.^2)).* (x + 1i * y) -...
         q2.* k_0 * J2M_in2.* x./ rho))).* exp(-1i * m(jm) * phi).* (abs(rho) < a_0);

    dEz_dy_out = dEz_dy_out + DE.* ((H2m.* m(jm)./ (rho.^2)).* (y - 1i * x) -...
        k_0 * H2M.* y./ rho).* exp(-1i * m(jm) * phi);
    dEz_dy_in =  dEz_dy_in  + 1i/HH * (B1 * n1 * q1 * ((J2m_in1.* m(jm)./ (rho.^2)).* (y - 1i * x) -...
         q1.* k_0 * J2M_in1.* y./ rho) + ...
         B2 * n2 * q2 * ((J2m_in2.* m(jm)./ (rho.^2)).* (y - 1i * x) -...
         q2.* k_0 * J2M_in2.* y./ rho)).* exp(-1i * m(jm) * phi).* (abs(rho) < a_0);   
    
    %%%% вычисление магнитного поля и его производных
    H_z_in  = H_z_in  + (- q1 * B1 * Jm1 - q2 * B2 * Jm2).*...
        exp(-1i * m(jm) * phi).* (abs(rho) < a_0);    
    H_z_out = H_z_out + q.* (DH.* H2m).* exp(-1i * m(jm) * phi);
    
    dHz_dx_out = dHz_dx_out + DH.* (((H2m.* m(jm)./ (rho.^2)).* (x + 1i * y) -...
        k_0 * H2M.* x./ rho)).* exp(-1i * m(jm) * phi);
    dHz_dx_in =  dHz_dx_in  + (-B1 * q1 * (((J2m_in1.* m(jm)./ (rho.^2)).* (x + 1i * y) -...
         q1.* k_0 * J2M_in1.* x./ rho)) -...
         B2 * q2 * (((J2m_in2.* m(jm)./ (rho.^2)).* (x + 1i * y) -...
         q2.* k_0 * J2M_in2.* x./ rho))).* exp(-1i * m(jm) * phi).* (abs(rho) < a_0);

    dHz_dy_out = dHz_dy_out + DH.* ((H2m.* m(jm)./ (rho.^2)).* (y - 1i * x) -...
        k_0 * H2M.* y./ rho).* exp(-1i * m(jm) * phi);
    dHz_dy_in =  dHz_dy_in  + (-B1 * q1 * ((J2m_in1.* m(jm)./ (rho.^2)).* (y - 1i * x) -...
         q1.* k_0 * J2M_in1.* y./ rho) - ...
         B2 * q2 * ((J2m_in2.* m(jm)./ (rho.^2)).* (y - 1i * x) -...
         q2.* k_0 * J2M_in2.* y./ rho)).* exp(-1i * m(jm) * phi).* (abs(rho) < a_0); 
     
%%%% поле и его производные по координатам необходимо домножать на
%%%% (abs(rho) < a_0) для того, что бы внутренние поля цилиндров не
%%%% суммировались друг с другом (что не правильно)
end


