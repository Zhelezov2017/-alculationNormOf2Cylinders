function IntegralOutside = IntegralOutside(VectorProduct, x, y)


IntegralOutsideDX = trapz(x,x.*VectorProduct);
IntegralOutside = IntegralOutsideDX;

end
