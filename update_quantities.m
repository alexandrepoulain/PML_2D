function [P_past,vect_epsi_n_2,vect_epsi_n,vect_sig_n,vect_Epsi_n,...
    vect_Sig_n] = update_quantities(dt,nnB,dof,ng,Nef_2,global_pos_elemB,u2...
    ,v2,vect_epsi_n_2,vect_epsi_n,Mat_Bepsi_n,Mat_BQ_n,Mat_Fepsi_n,...
    Mat_FQ_n,vect_sig_n,Mat_Dsig,vect_Epsi_n,vect_Sig_n)
% This function update the physical quantities at the end of a time loop
% Inputs: global_pos_elem: in this matrix each column corresponds to an
%                          element. Each element contains the global 
%                          position of each of its node. For each node,
%                          2 indices correspond to the 2 directions 


P_past = zeros(nnB*dof,1); % reset P_past to 0

% if we are in the case of 2 gauss points in each direction
switch ng
    case 2
        for i=1:Nef_2
            % Get back the effective displacement for the different nodes
            disp_ef = u2(global_pos_elemB(:,i)); 
            % and velocity
            vit_ef = v2(global_pos_elemB(:,i));
            % calculation of the first part of the strain
            vect_epsi_n_2{1}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{1}(:,:,i)*vit_ef) + (Mat_BQ_n{1}(:,:,i)*disp_ef) ;
            vect_epsi_n_2{2}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{2}(:,:,i)*vit_ef) + (Mat_BQ_n{2}(:,:,i)*disp_ef) ;
            vect_epsi_n_2{3}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{3}(:,:,i)*vit_ef) + (Mat_BQ_n{3}(:,:,i)*disp_ef) ;
            vect_epsi_n_2{4}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{4}(:,:,i)*vit_ef) + (Mat_BQ_n{4}(:,:,i)*disp_ef) ;
        end
        % Calculation of the second part of the strain at time n+1 (in Voigt notation)
        vect_epsi_n{1} = vect_epsi_n_2{1} + Mat_Fepsi_n{1}*vect_epsi_n{1}-Mat_FQ_n{1}*vect_Epsi_n{1};
        vect_epsi_n{2} = vect_epsi_n_2{2} + Mat_Fepsi_n{2}*vect_epsi_n{2}-Mat_FQ_n{2}*vect_Epsi_n{2};
        vect_epsi_n{3} = vect_epsi_n_2{3} + Mat_Fepsi_n{3}*vect_epsi_n{3}-Mat_FQ_n{3}*vect_Epsi_n{3};
        vect_epsi_n{4} = vect_epsi_n_2{4} + Mat_Fepsi_n{4}*vect_epsi_n{4}-Mat_FQ_n{4}*vect_Epsi_n{4};
        % Calculation of the stress tensor at time n+1 (In Voigt notation) 
        vect_sig_n{1} = Mat_Dsig{1}*vect_epsi_n{1};
        vect_sig_n{2} = Mat_Dsig{2}*vect_epsi_n{2};
        vect_sig_n{3} = Mat_Dsig{3}*vect_epsi_n{3};
        vect_sig_n{4} = Mat_Dsig{4}*vect_epsi_n{4};
        % update of the integral of the strain 
        vect_Epsi_n{1} = vect_Epsi_n{1} + dt*vect_epsi_n{1};
        vect_Epsi_n{2} = vect_Epsi_n{2} + dt*vect_epsi_n{2};
        vect_Epsi_n{3} = vect_Epsi_n{3} + dt*vect_epsi_n{3};
        vect_Epsi_n{4} = vect_Epsi_n{4} + dt*vect_epsi_n{4};
        % update of the integral of the stress
        vect_Sig_n{1} = vect_Sig_n{1} + dt*vect_sig_n{1};
        vect_Sig_n{2} = vect_Sig_n{2} + dt*vect_sig_n{2};
        vect_Sig_n{3} = vect_Sig_n{3} + dt*vect_sig_n{3};
        vect_Sig_n{4} = vect_Sig_n{4} + dt*vect_sig_n{4};
    case 3
        for i=1:Nef_2
            % Get back the effective displacement for the different nodes
            disp_ef = u2(global_pos_elemB(:,i)) ; 
            % and velocity
            vit_ef = v2(global_pos_elemB(:,i)) ;
            % calculation of the first part of the strain
            vect_epsi_n_2{1}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{1}(:,:,i)*vit_ef) + (Mat_BQ_n{1}(:,:,i)*disp_ef);
            vect_epsi_n_2{2}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{2}(:,:,i)*vit_ef) + (Mat_BQ_n{2}(:,:,i)*disp_ef);
            vect_epsi_n_2{3}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{3}(:,:,i)*vit_ef) + (Mat_BQ_n{3}(:,:,i)*disp_ef);
            vect_epsi_n_2{4}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{4}(:,:,i)*vit_ef) + (Mat_BQ_n{4}(:,:,i)*disp_ef);
            vect_epsi_n_2{5}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{5}(:,:,i)*vit_ef) + (Mat_BQ_n{5}(:,:,i)*disp_ef);
            vect_epsi_n_2{6}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{6}(:,:,i)*vit_ef) + (Mat_BQ_n{6}(:,:,i)*disp_ef);
            vect_epsi_n_2{7}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{7}(:,:,i)*vit_ef) + (Mat_BQ_n{7}(:,:,i)*disp_ef);
            vect_epsi_n_2{8}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{8}(:,:,i)*vit_ef) + (Mat_BQ_n{8}(:,:,i)*disp_ef);
            vect_epsi_n_2{9}((i-1)*3+1:(i-1)*3+3,1) = (Mat_Bepsi_n{9}(:,:,i)*vit_ef) + (Mat_BQ_n{9}(:,:,i)*disp_ef);
        end  
        % Calculation of the second part of the strain at time n+1 (in Voigt notation)
        vect_epsi_n{1} = vect_epsi_n_2{1} + Mat_Fepsi_n{1}*vect_epsi_n{1}-Mat_FQ_n{1}*vect_Epsi_n{1} ;
        vect_epsi_n{2} = vect_epsi_n_2{2} + Mat_Fepsi_n{2}*vect_epsi_n{2}-Mat_FQ_n{2}*vect_Epsi_n{2} ;
        vect_epsi_n{3} = vect_epsi_n_2{3} + Mat_Fepsi_n{3}*vect_epsi_n{3}-Mat_FQ_n{3}*vect_Epsi_n{3} ;
        vect_epsi_n{4} = vect_epsi_n_2{4} + Mat_Fepsi_n{4}*vect_epsi_n{4}-Mat_FQ_n{4}*vect_Epsi_n{4} ;
        vect_epsi_n{5} = vect_epsi_n_2{5} + Mat_Fepsi_n{5}*vect_epsi_n{5}-Mat_FQ_n{5}*vect_Epsi_n{5} ;
        vect_epsi_n{6} = vect_epsi_n_2{6} + Mat_Fepsi_n{6}*vect_epsi_n{6}-Mat_FQ_n{6}*vect_Epsi_n{6} ;
        vect_epsi_n{7} = vect_epsi_n_2{7} + Mat_Fepsi_n{7}*vect_epsi_n{7}-Mat_FQ_n{7}*vect_Epsi_n{7} ;
        vect_epsi_n{8} = vect_epsi_n_2{8} + Mat_Fepsi_n{8}*vect_epsi_n{8}-Mat_FQ_n{8}*vect_Epsi_n{8} ;
        vect_epsi_n{9} = vect_epsi_n_2{9} + Mat_Fepsi_n{9}*vect_epsi_n{9}-Mat_FQ_n{9}*vect_Epsi_n{9} ;
        % Calculation of the stress tensor at time n+1 (In Voigt notation) 
        vect_sig_n{1} = Mat_Dsig{1}*vect_epsi_n{1};
        vect_sig_n{2} = Mat_Dsig{2}*vect_epsi_n{2};
        vect_sig_n{3} = Mat_Dsig{3}*vect_epsi_n{3};
        vect_sig_n{4} = Mat_Dsig{4}*vect_epsi_n{4};
        vect_sig_n{5} = Mat_Dsig{5}*vect_epsi_n{5};
        vect_sig_n{6} = Mat_Dsig{6}*vect_epsi_n{6};
        vect_sig_n{7} = Mat_Dsig{7}*vect_epsi_n{7};
        vect_sig_n{8} = Mat_Dsig{8}*vect_epsi_n{8};
        vect_sig_n{9} = Mat_Dsig{9}*vect_epsi_n{9};
        % update of the integral of the strain 
        vect_Epsi_n{1} = vect_Epsi_n{1} + DT*vect_epsi_n{1};
        vect_Epsi_n{2} = vect_Epsi_n{2} + DT*vect_epsi_n{2};
        vect_Epsi_n{3} = vect_Epsi_n{3} + DT*vect_epsi_n{3};
        vect_Epsi_n{4} = vect_Epsi_n{4} + DT*vect_epsi_n{4};
        vect_Epsi_n{5} = vect_Epsi_n{5} + DT*vect_epsi_n{5};
        vect_Epsi_n{6} = vect_Epsi_n{6} + DT*vect_epsi_n{6};
        vect_Epsi_n{7} = vect_Epsi_n{7} + DT*vect_epsi_n{7};
        vect_Epsi_n{8} = vect_Epsi_n{8} + DT*vect_epsi_n{8};
        vect_Epsi_n{9} = vect_Epsi_n{9} + DT*vect_epsi_n{9};
        % update of the integral of the stress
        vect_Sig_n{1} = vect_Sig_n{1} + DT*vect_sig_n{1};
        vect_Sig_n{2} = vect_Sig_n{2} + DT*vect_sig_n{2};
        vect_Sig_n{3} = vect_Sig_n{3} + DT*vect_sig_n{3};
        vect_Sig_n{4} = vect_Sig_n{4} + DT*vect_sig_n{4};
        vect_Sig_n{5} = vect_Sig_n{5} + DT*vect_sig_n{5};
        vect_Sig_n{6} = vect_Sig_n{6} + DT*vect_sig_n{6};
        vect_Sig_n{7} = vect_Sig_n{7} + DT*vect_sig_n{7};
        vect_Sig_n{8} = vect_Sig_n{8} + DT*vect_sig_n{8};
        vect_Sig_n{9} = vect_Sig_n{9} + DT*vect_sig_n{9};
end

end

