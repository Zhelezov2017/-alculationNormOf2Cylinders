%%%%% Эта программа позовляет вычислять распределения поля при рассеяние
%%%%% плоской электомагнитной волны TE типа на нескольких цилиндрах

clear all
tic

%%%% задаём параметры системы
global a_0 ee k_0 Theta N_max cylXY L E0z H0z p EE GG HH N_cylinders 
GPC_systemParameters

n = 2000;
a = -2000;  
b = 2000;
cr = -2000;
d = 2000;
h = (b-a) / n; 
x=a:h:b;
y=cr:h:d;

N_max = 3;
p = 4; 

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
Sum = 0;

Man = sum(BeenFieldEX)+sum(BeenFieldEY)+sum(BeenFieldHX)+sum(BeenFieldHY);
Man(1,1001) = 0;
%IntegralOut = IntegralOutside(Sum , x, y);
IntOut = IntegralOutside(Man, x, y);

IntegralOutIntVent = IntegralNew(RateCalculationMatrixOutInt);


IntegralFull = Integral + IntOut - IntegralOutIntVent;

   

  

    




   