function [X, T, G] = global_topology(L1,h,Nex1,Ney1,Lef1,L2,Nex2,Ney2,Lef2,param_med,param_pml)
% This function define the global topology and for each node it defines if
% it is part of SDA or SDB and associates its parameters. If it is an
% interface node, the node will be marked as an interface node and will
% have a special treatment

% Start coordinates
% down left
cx1 = 0; 
cy1 = 0;
% Up right
cx2 = L1+L2; 
cy2 = h;
% Vector of coordinates in x and y direction
xp = cx1:(cx2-cx1)/(Nex1+Nex2):cx2;
yp = cy1:(cy2-cy1)/Ney1:cy2;
% Coordinate matrix
ind = 1;
for ii = 1:(Ney1+1)
    for jj = 1:(Nex1+Nex2+1)
        X(ind,:) = [xp(1,jj) yp(1,ii)];
        G{ind} = param_med;
        ind  = ind+1;
    end
end
% Connectivity
ip = 1 ; ll = 1 ;
for ii = 1:Ney1
    for jj = 0:((Nex1+Nex2)-1)
        T(ll,:) = [(ip+jj) (ip+jj+1) (ip+jj+(Nex1+Nex2)+2) (ip+jj+(Nex1+Nex2)+1) 1] ;
        ll = ll + 1 ;
    end
    ip = ip + (Nex1+1) ;
end
%% Parameter for each element (PML or Med)
for ii = 1:size(T,1)
   Xe = X(T(ii,1:size(T,2)-1));
   % medium
   if Xe(1,1) < L1
       G{ii} = param_med;
   % PML
   else
       G{ii} = param_pml;
   end
end


end
