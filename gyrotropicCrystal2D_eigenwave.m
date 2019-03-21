%%%%% Эта программа позовляет вычислять распределения поля при рассеяние
%%%%% плоской электомагнитной волны TE типа на нескольких цилиндрах

clear all
tic

%%%% задаём параметры системы
global a_0 ee k_0 Theta N_max cylXY L E0z H0z p EE GG HH N_cylinders 
GPC_systemParameters

n = 10000;
a = -20000;  
b = 20000;
cr = -20000;
d = 20000;
h = (b-a) / n; 
x=a:h:b;
y=cr:h:d;

N_max = 0;
p = 8.298; 

N_cylinders = 2;
L = 3 * a_0;
Outside = 0;
[BD] = crystalCoefficients2D_eigenwave(N_max, w_0, a_0, cylXY, k_0, EE, GG, HH, c, p);

RateCalculationMatrix = RateCalculationMatrix(N_max, k_0, p, a_0, EE, GG, HH , BD, N_cylinders );

RateCalculationMatrixOutInt = RateCalculationMatrixOutInt(N_max, k_0, p, a_0, EE, GG, HH , BD, N_cylinders );



Integral = Integral(RateCalculationMatrix);
BeenFieldEX =  BeenFieldEX(N_cylinders, N_max, cylXY, BD,k_0, p, x, y);
BeenFieldEY =  BeenFieldEY(N_cylinders, N_max, cylXY, BD,k_0, p, x, y);
BeenFieldHX =  BeenFieldHX(N_cylinders, N_max, cylXY, BD,k_0, p, x, y);
BeenFieldHY =  BeenFieldHY(N_cylinders, N_max, cylXY, BD,k_0, p, x, y);
VectorProduct = VectorProduct( BeenFieldHX, BeenFieldHY, BeenFieldEX, BeenFieldEY );
IntegralOut = IntegralOutside(VectorProduct, x, y);


IntegralOutIntVent = IntegralNew(RateCalculationMatrixOutInt);


IntegralFull = Integral + IntegralOut - IntegralOutIntVent;
valueSet = zeros(4,(2*N_max+1)*N_cylinders);

%for jj = 0: N_cylinders-1
%     for jm = 1:2*N_max+1
%         DC = zeros(1,4);
%         for jh = 1:4
%             DC(1,jh)=BD(4*(2*N_max+1)*jj+jh+4*jm-4);
%             valueSet(jh,(2*N_max+1)*jj-jm+1)=DC(1,jh);
%         end 
%         R(1:4)=valueSet(1:4,jj);     
%     end
%end     





   