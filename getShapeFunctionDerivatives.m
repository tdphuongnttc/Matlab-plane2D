function [deriX,deriY] = getShapeFunctionDerivatives(nnel,deriXi,deriEta,invJ)
%invJ : inverse Jacobian
%nnel: number of side of element
 for i=1:nnel
    deriX(i)=invJ(1,1)*deriXi(i)+invJ(1,2)*deriEta(i);
    deriY(i)=invJ(2,1)*deriXi(i)+invJ(2,2)*deriEta(i);
 end



