function IntegralNew = IntegralNew(RateCalculationMatrixOutInt)
n = 10000;
a = 0;  
b = 100;  
h = (b-a) / n;  
x=a:h:b; 
IntegralNew = trapz(x, x*RateCalculationMatrixOutInt);


end

