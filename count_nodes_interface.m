function [interface_A,interface_B] = count_nodes_interface(L1,XA,XB)
% This functions count the number of nodes at the interface in A and B
% Inputs: L1: length of the soil
%         XA: Position of the nodes in A
%         XB: position of the nodes in B
% Number of nodes in A
disp("Number of nodes in A:");
nodes_A = length(XA(:,1))
disp("Number of nodes in B:");
nodes_B = length(XB(:,1))

ind_node = 1 ;
for i=1:nodes_A
  if (XA(i,1) == L1)
      interface_A(ind_node) = i;
      ind_node = ind_node + 1;
  end
end
ind_node = 1 ;
for i=1:nodes_B
  if (XB(i,1) == L1)
      interface_B(ind_node) = i;
      ind_node = ind_node + 1;
  end
end
end

