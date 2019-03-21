function IntegralOutside = IntegralOutside(VectorProduct, x, y)

IntegralOutsideDX = trapz(sqrt(x.^2+y.^2),sqrt(x.^2+y.^2).* VectorProduct);

IntegralOutside = IntegralOutsideDX;

end
