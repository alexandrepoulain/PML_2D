function [globX, globT, paramX] = global_topology(XA,TA,XB,TB,interfaceA...
    , interfaceB,param_Med,param_PML,param_interface)
% This function define the global topology and for each node it defines if
% it is part of SDA or SDB and associates its parameters. If it is an
% interface node, the node will be marked as an interface node and will
% have a special treatment

nn_tot = size(XA,1) + size(XB,1) - length(interfaceA); % total nb nodes 
globX = zeros(nn_tot,size(XA,2));
globX(1:size(XA,1),:) = XA(:,:);
paramX(1:size(XA,1),:) = param_Med; 
% for each node in SDB
for ii=1:size(XB,1)
    if ~ ismember(ii,interfaceB)
        globX(size(XA,1)+ii,1) = (XB(ii,1));
        globX(size(XA,1)+ii,2) = (XB(ii,2));
        paramX(size(XA,1)+ii,:) = param_PML;
    end
        
end
% Interface node parameter
paramX(interfaceA,:) = param_interface; 
% global connectivity matrix
globT = zeros(size(TA,1)+size(TB,1),5); 
globT(1:size(TA,1),:) = TA(:,:);
ind =1;
for ii=1:size(TB,1)
    Te = TB(ii,1:4);
    if ismember(interfaceB, Te)==[true,true]
        for jj=1:size(Te,2)
            if Te(1,jj) == interfaceB(1)  
                globT(size(TA,1)+ii,jj) = interfaceA(1);
            elseif Te(1,jj) == interfaceB(2)  
                globT(size(TA,1)+ii,jj) = interfaceA(2);
            else
                globT(size(TA,1)+ii,jj) = Te(1,jj)+size(XA,1);
            end
        end
        index_interface_elements(size(TA,1)+ind) = ii;
        ind = ind+1;
    else
        globT(size(TA,1)+ii,:) = TB(ii,:)  + ...
        [size(XA,1) size(XA,1) size(XA,1) size(XA,1) 0];
    end
end
end

