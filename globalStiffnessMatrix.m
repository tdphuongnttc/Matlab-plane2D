function [K]= globalStiffnessMatrix(nelem,elemdata,nnel,nodeCoor,ndof,edof,ngp,K)
%nelem :number of element

%calculate x and w matrix
[x,w]=Gauss(ngp);
A=elemdata(:,1:4); 
 %create x coordinate and y coordinate matrix from A
 [xx,yy]=coordinateMatrix(A,nelem,nodeCoor);
 %create global index base on local index
 [index]= getGlobalIndex(ndof,nelem,nnel,A);

for i=1:nelem     
     k = zeros(edof,edof); %k cua tung element 8x8
     for igx=1:ngp
        % xi and wtx number for outer intergral
        ptxi = x(igx);
        wtx = w(igx); 
        
        for igy=1:ngp 
             % eta and wty number for inner intergral
            pteta = x(igy);
            wty = w(igy);
            [deriXi]=createShapeFunctionDerivativeXi(pteta);
            [deriEta] = createShapeFunctionDerivativeEta(ptxi);
            [detJ,invJ,J] = Jacobian(deriXi,deriEta,xx(i,:),yy(i,:),nnel);%,,,,,
            [deriX,deriY] = getShapeFunctionDerivatives(nnel,deriXi,deriEta,invJ);
            [B] =  getOperators(edof,nnel,deriX,deriY);
            [C,E,nuy,h]=materialMatrix(i,elemdata);
            k = k+B'*C*B*wtx*wty*detJ*h;
        end      
     end
     K(index(i,:),index(i,:))=K(index(i,:),index(i,:))+k;
    
end