% calcul du chargement en fonction du temps
% 
% 
% 
% 

function [Ftps1] = F_ricker(eff,x)

       ts=3;
       tp=3;
      if (x<=3*ts) ;
        Ftps1=(2*pi^2*(x-ts).^2./tp^2-1).*exp(-pi^2*(x-ts).^2./tp^2).*eff;
      else
        Ftps1 = 0 ;  
      end
