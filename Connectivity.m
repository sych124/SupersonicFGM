%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Element Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:x_ele
    for j=1:y_ele
        ele_indx=i+x_ele*(j-1);                                 % element number
        node_indx(ele_indx,1)=x_node*2*(j-1)+i*2-1;             % 1st node number of the element
        node_indx(ele_indx,5)=node_indx(ele_indx,1)+1;          % 5th node number of the element
        node_indx(ele_indx,2)=node_indx(ele_indx,5)+1;          % 2nd node number of the element
        node_indx(ele_indx,8)=node_indx(ele_indx,1)+x_node;     % 8th node number of the element
        node_indx(ele_indx,9)=node_indx(ele_indx,8)+1;     % 9th node number of the element
        node_indx(ele_indx,6)=node_indx(ele_indx,9)+1;     % 6th node number of the element
        node_indx(ele_indx,4)=node_indx(ele_indx,8)+x_node;     % 4th node number of the element
        node_indx(ele_indx,7)=node_indx(ele_indx,4)+1;     % 7th node number of the element
        node_indx(ele_indx,3)=node_indx(ele_indx,7)+1;     % 3rd node number of the element
             
    end
end
