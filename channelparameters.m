% % plasma chanel parameters:

%%%%%%%% for ionosphere %%%%%%%%%%%%%%%%
% plasma chanel parameters:
global R EE GG HH a_0 c w w_H w_p typeOfmedia

if(isempty(typeOfmedia))
    typeOfmedia = 'Plasma';
end

if(strcmp(typeOfmedia, 'Plasma'))
n_e = 1e6;
e_0 = 4.8e-10;
m_e = 9.1e-28;
c = 3e10;   % velocity of light

n_e = 1e13; % концентрация электронов
H_0 = 800;
a_0 = 2.5;

n_e = 1e6; % концентрация электронов
H_0 = 0.5;
a_0 = 100;


w_H = (e_0 * H_0) / (m_e * c);
w_p = sqrt(4 * pi * e_0^2 * n_e / m_e);
% w_p = 5.5e4*(n_e)^(1/2);


% e_0 = 4.8e-10;
% T = 2.1494;
% T_e = T * 1.6e-12;
% r_D = (T_e / (4 * pi * e_0^2 * n_e))^(1/2);
M_Ar=39.95*1.66e-24;
O_p = sqrt(4 * pi * e_0^2 * n_e / M_Ar);
O_H = w_H * m_e / M_Ar;
% Par = 1;
c = 3e10;   % velocity of light

Nu_e = 0 * w_H;
Nu_i = 0;
w_LH = w_H * sqrt((O_p^2 + O_H^2) / (w_p^2 + w_H^2));
w_UH = sqrt(w_p^2 + w_H^2);

% w_0 = 0.01 * w_H;

% w_0 = w;
% % w_0 = [0.1:0.001:3]*w_H;
% EE = (w_0.^2 - w_UH.^2).* (w_0.^2 - w_LH.^2).*...
%             ((w_0.^2  - w_H^2).* (w_0.^2 - O_H^2)).^(-1) -...
%             w_p.^2 * 1i*Nu_e./ ((w_H^2 - w_0.^2).* w_0);
        
%%%%% EE при учете потерь
EE = 1 + w_p.^2 * (w_0 - 1i*Nu_e)./ ((w_H.^2 - (w_0 - 1i*Nu_e).^2).* w_0) +...
         O_p.^2 * (w_0 - 1i*Nu_i)./ ((O_H.^2 - (w_0 - 1i*Nu_i).^2).* w_0)
        
GG =  -w_p.^2 * w_H./ ((w_H.^2 - (w_0 - 1i*Nu_e).^2).* w_0) +...
              O_p.^2 * O_H./ ((O_H.^2 - (w_0).^2).* w_0);
          
HH = 1 - w_p.^2./ ((w_0 - 1i*Nu_e).* w_0) - O_p.^2./ ((w_0).* w_0);

Pcin  = sqrt(EE - ((EE + HH).* GG.^2).* (EE - HH).^(-2)  +...
         2 * ((EE - HH).^(-2)).* sqrt((EE.* HH.* GG.^2).* (GG.^2 - (EE - HH).^2)))
    %prop.cnst for conic refr. waves
Pin=sqrt(EE-GG)

% figure 
% plot(w_0/w_H,real(EE), w_0/w_H, imag(EE),'r')
% figure 
% plot(w_0/w_H,real(GG), w_0/w_H, imag(GG),'r')

% EE = 1 - ((I_VW.* v).* (I_VW.^2 - sqrt(u.* U))).* ((I_VW.^2 - u).* (I_VW.^2 - U)).^(-1);
% GG = (((I_VW.^2).* v).* sqrt(u)).* ((I_VW.^2 - u).* (I_VW.^2 - U)).^(-1);
% HH = 1 - v.* I_VW.^(-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EE = 1 - (w_p.^2)./ w.^2;
% GG = 0;
% HH = EE;


% EE = -1.14 - 0.1 * 1i;
% % EE = -1.14;
% % EE = -1 - 0.5 * 1i;
% % % EE = -1;
% % % EE = 4;
% % % EE = 2;
% % 
% % EE = 10000;
% GG = 0;
% HH = EE;

elseif(strcmp(typeOfmedia, 'Graphen'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Расчёт диэлектрической проницаемсоти графена-диэлектрической
%%%%%%%%%%%% слоистой среды

%%%%% поверхностная проводимость графена. Взято из стать
% @ARTICLE{Gusynin2007,
%   author = {V P Gusynin and S G Sharapov and J P Carbotte},
%   title = {Magneto-optical conductivity in graphene},
%   journal = {Journal of Physics: Condensed Matter},
%   year = {2007},
%   volume = {19},
%   pages = {026222},
%   number = {2},
%   abstract = {Landau level quantization in graphene reflects the Dirac nature of
% 	its quasiparticles and has been found to exhibit an unusual integer
% 	quantum Hall effect. In particular, the lowest Landau level can be
% 	thought of as shared equally by electrons and holes, and this leads
% 	to characteristic behaviour of the magneto-optical conductivity as
% 	a function of frequency О© for various values of the chemical potential
% 	Ој. Particular attention is paid to the optical spectral weight under
% 	various absorption peaks and its redistribution as Ој is varied.
% 	We also provide results for magnetic field B as well as chemical
% 	potential sweeps at selected fixed frequencies, which can be particularly
% 	useful for possible measurements in graphene. Both diagonal and Hall
% 	conductivities are considered.},
%   url = {http://stacks.iop.org/0953-8984/19/i=2/a=026222}
% }

c        = 3e10;                         %%% скоротсть света в см/с
Delta    = 0;
vF       = 10^8;                         %%% в см/c
e_0      = 4.8e-10;                      %%% в СГСЭ
B_tesla  = 1;                            %%% внешнее магнитное поле в теслах
B_0      = B_tesla * 1e4;                %%% внешнее магнитное поле в эрстедах

mu_eV    = 58.67 * 1e-3;                 %%% химический потенциал в эВ
mu       = mu_eV * (1.6 * 1e-12);        %%% химический потенциал в эрг

% mu_K     = 680;                          %%% химический потенциал в кельвинах
% mu       = mu_K * 1.3805 * 1e-16;        %%% химический потенциал в эрг

T_K      = 10;                           %%% температура в кельвинах
T        = T_K * 1.3805 * 1e-16;         %%% температура в эрг
h        = 6.6748 * 1e-27;               %%% постоянная Планка в  эрг с
h_dirac  = 1.0545 * 1e-27;               %%% дираковская постоянная Планка в  эрг с

% h = 6.582119 * 1e-16;     %%% постоянная Планка  h = 6.582119 * 1e-16 эВ с

GaM_n_eV = 1.3 * 1e-3;                   %%% энергия рассения электрона в эВ
GaM_n    = GaM_n_eV * (1.6 * 1e-12);     %%% энергия рассения электрона в эрг
% GaM_n_K  = 15;                           %%% энергия рассения электрона в кельвинах
% GaM_n    = GaM_n_K  * 1.3805 * 1e-16;    %%% энергия рассения электрона в эрг
GaM_n1   = GaM_n;

w0_eV    = 1e-2;                         %%% энергия падающего излучения в эВ
w0       = w0_eV * (1.6 * 1e-12);        %%% энергия падающего излучения в эрг

% w_H      = 1.7e7 * T_K * 1e4;            %%% циклотронная частот
% E_land   = h_dirac * w_H;                %%% первый уровень Ландау (полуклассический)

%%% функция вичисления энергии уровня Ландау
Mn     = @(n) sign(n) * sqrt(Delta.^2 + 2 *  abs(n) * h_dirac * vF.^2 * e_0 * B_0 / c);
%%% функция вичисления распределения Ферми
nF     = @(Mnn) 1./ (exp((Mnn - mu)/T) + 1);

w_0 = w0;

% w_0 = [0.1:0.1:30]*w0;
% 
for iw = 1:size(w_0,2)
    w0 = -w_0(iw);
n_max        = 25; %%% максимальный рассматриваемый уровень Ландау
sigmaDial    = 0;
sigmaNondiag = 0;
for n = 0:n_max
    
M_n  = Mn(n) + 1e-30; %%% прибавляем эту малую величину для избавления деления на ноль
M_n1 = Mn(n + 1);
    
sigmaDial = sigmaDial - e_0^2 * vF.^2 * e_0 * B_0 * h_dirac / (h * 1i * c ).*...
      ((1-Delta^2/(1e-20+M_n.* M_n1)).* ((nF(M_n) - nF(M_n1)) + (nF(-M_n1)-nF(-M_n))).*...
       (1./(M_n - M_n1 + w0 + 1i*(GaM_n + GaM_n1)) - 1./(M_n - M_n1 - w0 -1i*(GaM_n + GaM_n1)))./(M_n1 - M_n) +...
       (1+Delta^2/(1e-20+M_n.* M_n1)).* ((nF(-M_n)- nF(M_n1)) + (nF(-M_n1)-nF(M_n))).*...
       (1./(M_n + M_n1 + w0 + 1i*(GaM_n + GaM_n1)) - 1./(M_n + M_n1 - w0 -1i*(GaM_n + GaM_n1)))./(M_n1 + M_n));
   
sigmaNondiag = sigmaNondiag + e_0^2 * vF^2 * e_0* B_0  * h_dirac /(h * c).*...
    (((nF(M_n)-nF(M_n1))-(nF(-M_n1)-nF(-M_n))).*...
    ((1-Delta.^2/(M_n.* M_n1)).*(1/(M_n - M_n1 + w0 + 1i * (GaM_n + GaM_n1)) + 1/(M_n - M_n1 - w0 - 1i*(GaM_n+GaM_n1)))/(M_n1-M_n)-...
     (1+Delta.^2/(M_n.* M_n1)).*(1/(M_n + M_n1 + w0 + 1i * (GaM_n + GaM_n1)) + 1/(M_n + M_n1 - w0 - 1i*(GaM_n+GaM_n1)))/(M_n1+M_n)));


end
sigm(iw)        = sigmaDial;
sigmaNon(iw)    = sigmaNondiag;
end


% norm_w     = sqrt(2 *  h_dirac * vF.^2 * e_0 * B_0 / c);
% % norm_w     = (1.6 * 1e-12);

% norm_sigma = e_0^2 / h;
% norm_sigma = e_0^2 / h / (GaM_n/sqrt(h_dirac * vF.^2 * e_0 * B_0 / c));
% % norm_sigma = c;

EE_0 = 4.5;              %%% диэлектрическа проницамость диэлектрика между графеновыми дисками
                         %%% В нашем случае SiO2
d    = 1e-6;             %%% толщина диэлектрического слоя в см
w_0   = w_0 / h_dirac;     %%% частота падающего излучения в радиан/c
       
EE = EE_0 - 1i * 4 * pi * sigm./ (w_0 * d);
GG = 4 * pi * sigmaNon./ (w_0 * d);
HH = EE_0* ones(size(w_0));


end

% ww0 = 4 * pi / d / (1 + EE_0) * (e_0^2/(4 * h_dirac))
% 
% ww0 = ww0 * h_dirac /  (1.6 * 1e-12)



% figure 
% plot(h_dirac*w_0/(1.6 * 1e-12),real(EE), h_dirac*w_0/(1.6 * 1e-12), imag(EE),'r')
% figure 
% plot(h_dirac*w_0/(1.6 * 1e-12),real(GG), h_dirac*w_0/(1.6 * 1e-12), imag(GG),'r')



% figure
% plot(w_0/norm_w,real(sigm)/norm_sigma,w_0/norm_w,imag(sigm)/norm_sigma,'r')
% 
% figure
% plot(w_0/norm_w,real(sigmaNon)/norm_sigma,w_0/norm_w,imag(sigmaNon)/norm_sigma,'r')

% sigm     = -1i*sigm;
% sigmaNon = sigmaNon;
% 
% figure
% plot(w_0/norm_w,real(sigm)/norm_sigma,w_0/norm_w,imag(sigm)/norm_sigma,'r')
% 
% figure
% plot(w_0/norm_w,real(sigmaNon)/norm_sigma,w_0/norm_w,imag(sigmaNon)/norm_sigma,'r')

% clear all
 

















