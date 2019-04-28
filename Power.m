%%%%% Эта программа считает мощность
%%%%% 

clear all
tic

%%%% задаём параметры системы
global a_0 ee k_0 Theta N_max cylXY L E0z H0z p EE GG HH N_cylinders 
GPC_systemParameters

N_max = 0;
[BD] = crystalCoefficients2D_eigenwave(N_max, w_0, a_0, cylXY, k_0, EE, GG, HH, c, p);

E_usually = ElectricField(N_max, BD, k_0, p, a_0, EE, GG, HH);% + ElectricFieldEphi(N_max, BD, k_0, p, a_0, EE, GG, HH);
IntegralFull = -1271246357504.24 - 967098839412.943i;


%Задаем  J
d = 2;
I_0 = 4;
j = I_0 / d;

E_unussualy = ElectricField(N_max, BD, k_0, -p, a_0, EE, -GG, HH); %ElectricFieldEphi(N_max, BD, k_0, -p, a_0, EE, -GG, HH);




%2 формулы для P
a_p =  a_0 * 2 *pi *I_0/(IntegralFull*1i*k_0*p)*(exp(d)-exp(-d));
P_1 = 1/4 * IntegralFull * a_p^2;%[796.877305450313 - 606.223264834755i]
P_2 = -a_0 *pi *j* E_usually*(exp(d)-exp(-d))/(-1i*k_0*p);%[-33.8572446978129 - 2006644.89641296i]
