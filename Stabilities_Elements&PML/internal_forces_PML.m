function [vect_past,P_past] = internal_forces_PML(dof,ng,Mat_Sig_n,vect_Sig_n,Mat_epsi_n,...
        vect_epsi_n,Mat_Epsi_n,vect_Epsi_n,P_past,TB,Nef_2)
% This function calculate the vector of internal forces in the PML.
% Inputs: sig: Matrix of the stress at each Gauss point 
%         Sig: Matrix of the integral of the stress at each gauss point
%         XB: Coordinates of each nodes in the PML
%         TB: Connectivity matrix
%         dof: Number of degree of freedom

switch ng
    case 2
        vect_past = Mat_Sig_n{1}*vect_Sig_n{1} + Mat_Sig_n{2}*vect_Sig_n{2} + ...
                     Mat_Sig_n{3}*vect_Sig_n{3} +  Mat_Sig_n{4}*vect_Sig_n{4} ;
        vect_past = vect_past + Mat_epsi_n{1}*vect_epsi_n{1} +  Mat_epsi_n{2}*vect_epsi_n{2} + ...
                     Mat_epsi_n{3}*vect_epsi_n{3} +  Mat_epsi_n{4}*vect_epsi_n{4} ;
        vect_past = vect_past - Mat_Epsi_n{1}*vect_Epsi_n{1} - Mat_Epsi_n{2}*vect_Epsi_n{2} - ...
                     Mat_Epsi_n{3}*vect_Epsi_n{3} - Mat_Epsi_n{4}*vect_Epsi_n{4} ;
    case 3
        vect_past = Mat_Sig_n{1}*vect_Sig_n{1} + Mat_Sig_n{2}*vect_Sig_n{2} + ...
                 Mat_Sig_n{3}*vect_Sig_n{3} +  Mat_Sig_n{4}*vect_Sig_n{4} + ...
                 Mat_Sig_n{5}*vect_Sig_n{5} + Mat_Sig_n{6}*vect_Sig_n{6} + ...
                 Mat_Sig_n{7}*vect_Sig_n{7} +  Mat_Sig_n{8}*vect_Sig_n{8} + ...
                 Mat_Sig_n{9}*vect_Sig_n{9};
        vect_past = vect_past + Mat_epsi_n{1}*vect_epsi_n{1} +  Mat_epsi_n{2}*vect_epsi_n{2} + ...
                 Mat_epsi_n{3}*vect_epsi_n{3} +  Mat_epsi_n{4}*vect_epsi_n{4} + ...
                 Mat_epsi_n{5}*vect_epsi_n{5} +  Mat_epsi_n{6}*vect_epsi_n{6} + ...
                 Mat_epsi_n{7}*vect_epsi_n{7} +  Mat_epsi_n{8}*vect_epsi_n{8} + ...
                 Mat_epsi_n{9}*vect_epsi_n{9};
        vect_past = vect_past - Mat_Epsi_n{1}*vect_Epsi_n{1} - Mat_Epsi_n{2}*vect_Epsi_n{2} - ...
                 Mat_Epsi_n{3}*vect_Epsi_n{3} - Mat_Epsi_n{4}*vect_Epsi_n{4} - ...
                 Mat_Epsi_n{5}*vect_Epsi_n{5} - Mat_Epsi_n{6}*vect_Epsi_n{6} - ...
                 Mat_Epsi_n{7}*vect_Epsi_n{7} - Mat_Epsi_n{8}*vect_Epsi_n{8} - ...
                 Mat_Epsi_n{9}*vect_Epsi_n{9};
end
% Put the force at the nodes
ne = size(TB,2) - 1;
for i=1:Nef_2
    for s = 1:ne
        P_past((TB(i,s)-1)*dof+1)=P_past((TB(i,s)-1)*dof+1) + vect_past((i-1)*8+(s-1)*dof+1);
        P_past((TB(i,s)-1)*dof+2)=P_past((TB(i,s)-1)*dof+2) + vect_past((i-1)*8+(s-1)*dof+2);         
    end
end


end
