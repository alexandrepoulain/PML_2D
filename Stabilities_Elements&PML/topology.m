function [X,T] = topology(xstart,ystart,xend,yend,Nex,Ney)
% This function return the coordinates of the nodes and the connectivity
% matrix.
% Inputs: xstart: x coordinate of the starting point
%         ystart: y coordinate of the starting point
%         xend: x coordinate of the end point
%         yend: y coordinate of the end point
%         Nex: number of elements in x direction
%         Ney: number of elements in y direction
% Start coordinates
% down left
cx1 = xstart; 
cy1 = ystart;
% Up right
cx2 = xend; 
cy2 = yend;
% Vector of coordinates in x and y direction
xp = cx1:(cx2-cx1)/Nex:cx2;
yp = cy1:(cy2-cy1)/Ney:cy2;
% Coordinate matrix
ind = 1;
for ii = 1:(Ney+1)
    for jj = 1:(Nex+1)
        X(ind,:) = [xp(1,jj) yp(1,ii)];
        ind  = ind+1;
    end
end
% Connectivity
ip = 1 ; ll = 1 ;
for ii = 1:Ney
    for jj = 0:(Nex-1)
        T(ll,:) = [(ip+jj) (ip+jj+1) (ip+jj+Nex+2) (ip+jj+Nex+1) 1] ;
        ll = ll + 1 ;
    end
    ip = ip + (Nex+1) ;
end
end

