function Integral = Integral(RateCalculationMatrix)
n = 10000;
a = 0;  
b = 100;  
h = (b-a) / n;  
x=a:h:b; 
Integral = trapz(x, x*(RateCalculationMatrix));


end

