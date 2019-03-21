function RateCalculationMatrix = RateCalculationMatrix(N_max, k_0, p, rho  , EE, GG, HH, BD, N_cylinders )
m = [-N_max: N_max];  
if(N_max == 0) m = 0; end
    %Функции внутри численно
    valueSet = zeros(4,(2*N_max+1)*N_cylinders);
    for jj = 1: N_cylinders
        for jm = 1:2*N_max+1
            DC = zeros(1,4);
            for jh = 1:4
                for jn = 1:4
                    DC(1,jn)=BD(4*(2*N_max+1)*jj-4*jm-jn+5);
                end
                valueSet(jh,(2*N_max+1)*jj-jm+1)=DC(1,5-jh);
            end 
            
        end
    end  
   
    
    IntegralForm1 = 0;
    IntegralForm2 = 0;
    for jj = 1:2*N_max+1
        %R= values(M,{2});
        R(1:4)=valueSet(1:4,jj);
        M(1:4)=valueSet(1:4,jj+2*N_max+1);
        MagneticField1 = MagneticField(m(jj), R, k_0, p, rho, EE, GG, HH);
        ElectricField1 = ElectricField(m(jj), R, k_0, p, rho, EE, GG, HH);
        MagneticField2 = MagneticField(m(jj), M, k_0, p, rho, EE, GG, HH);
        ElectricField2 = ElectricField(m(jj), M, k_0, p, rho, EE, GG, HH);
        MagnetOnElect1 = MagnetOnElect( ElectricField1, MagneticField1); 
        MagnetOnElect2 = MagnetOnElect( ElectricField2, MagneticField2);
        InternalForm1 = IntegralForm1 + MagnetOnElect1;     
        InternalForm2 = IntegralForm2 + MagnetOnElect2;
        RateCalculationMatrix =+ InternalForm1 + InternalForm2;
    end
    
   
end