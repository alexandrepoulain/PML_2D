function [M_global, K_global, C_global] = construct_global(XA, TA, XB, TB,...
    interfaceA, interfaceB, MA,KA,CA,MB,KB,CB,KKB,CCB, dof)
% This function construct the global matrices to bind both sub-domain
% together
% Inputs: XA, TA, XB, TB: Positions and connectivity matrices for both
%                         sub-domains.
%         interfaceA, interfaceB: indices of the interface node for both 
%                                 subdomains 
%         MA,KA,CA: mass, stiffness, damping matrices for SDA
%         MB,KB,CB: mass, stiffness, damping matrices for SDB
%         KKB, CCB: effective stiffness and damping for SDB
%         dof: number degree of freedom
% Initialization
nn_tot = size(XA,1) + size(XB,1) - length(interfaceA); % total nb nodes 
global_X = zeros(nn_tot,size(XA,2));
global_X(1:size(XA,1),:) = XA(:,:);
% for each node in SDB
for ii=1:size(XB,1)
    if ~ ismember(ii,interfaceB)
        global_X(size(XA,1)+ii,1) = (XB(ii,1));
        global_X(size(XA,1)+ii,2) = (XB(ii,2));
    end
end
% global connectivity matrix
global_T = zeros(size(TA,1)+size(TB,1),5); 
global_T(1:size(TA,1),:) = TA(:,:);
ind =1;
for ii=1:size(TB,1)
    Te = TB(ii,1:4);
    if ismember(interfaceB, Te)==[true,true]
        for jj=1:size(Te,2)
            if Te(1,jj) == interfaceB(1)  
                global_T(size(TA,1)+ii,jj) = interfaceA(1);
            elseif Te(1,jj) == interfaceB(2)  
                global_T(size(TA,1)+ii,jj) = interfaceA(2);
            else
                global_T(size(TA,1)+ii,jj) = Te(1,jj)+size(XA,1);
            end
        end
        index_interface_elements(size(TA,1)+ind) = ii;
        ind = ind+1;
    else
        global_T(size(TA,1)+ii,:) = TB(ii,:)  + ...
        [size(XA,1) size(XA,1) size(XA,1) size(XA,1) 0];
    end
end
M_global = zeros(nn_tot*dof,nn_tot*dofnn_tot*dof);     % mass matrix
K_global = zeros(nn_tot*dof,nn_tot*dofnn_tot*dof);     % stiffness matrix
C_global = zeros(nn_tot*dof,nn_tot*dofnn_tot*dof);     % damping matrix
nn = size(global_T,2)-1;
% for SDA
M_global(1:size(XA,1)*dof,1:size(XA,1)*dof) = MA;
K_global(1:size(XA,1)*dof,1:size(XA,1)*dof) = KA;
C_global(1:size(XA,1)*dof,1:size(XA,1)*dof) = CA;
% for SDB: we go through each node
for ii=size(TA,1)+1:size(global_T,1)
    Te = global_T(ii,:);
    me = zeros(dof*nn,dof*nn);
    ke = zeros(dof*nn,dof*nn);
    ce = zeros(dof*nn,dof*nn);
    for k=1:nn
        % find the corresponding element in the matrices
        me = me +
    end
    
end


end

