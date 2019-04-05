function BeenFieldHY = BeenFieldHY(N_cylinders, N_max, cylXY, Const,k_0, p, x, y )
m = [-N_max: N_max];

valueSet = zeros(4,(2*N_max+1)*N_cylinders);
Summ = 0;   
 for jj = 1: N_cylinders
      for jm = 1:2*N_max+1
          DC = zeros(1,4);
          for jh = 1:4
              for jn = 1:4
                  DC(1,jn)=Const(4*(2*N_max+1)*jj-4*jm-jn+5);
              end
              valueSet(jh,(2*N_max+1)*jj-jm+1)=DC(1,5-jh);
          end
      end
 end  
  
for jj = 0: N_cylinders-1
      for jm = 1:2*N_max+1
          R(1:4)=valueSet(1:4,jj*(2*N_max+1) + jm );
          Summ =+ FieldOutsideHY(m(jm), R, k_0, p, cylXY(jj+1,1), cylXY(jj+1,2), x, y);
      end
end


BeenFieldHY = Summ;


end



