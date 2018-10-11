function [r_p]=forceVector(nelem,elemdata,elembody,nnel,nodeCoor,ndof,edof,sdof,ngp,K)
% r_p : body force
%nelem : number of element
%nnel: number of side of element
%create weight and xi of gauss point
[x,w]=Gauss(ngp);
%create matrix with ordinate of node per element
A=elemdata(:,1:4);
%body force of element
r_p_e=zeros(edof,1);
%body force of all element
r_p=zeros(sdof,1);
for i=1:nelem     
  [xx,yy]=coordinateMatrix(A,nelem,nodeCoor);
       for igx=1:ngp
        ptxi = x(igx); %exact coordinates gauss point
        wtx = w(igx); %exact weight number
        for igy=1:ngp            
            pteta = x(igy);
            wty = w(igy);
            [shape]=createShapeFunction(ptxi,pteta);
            [deriXi]=createShapeFunctionDerivativeXi(pteta);
            [deriEta] = createShapeFunctionDerivativeEta(ptxi);
            [detJ,invJ,J] = Jacobian(deriXi,deriEta,xx(i,:),yy(i,:),nnel);
            [N] = ansatzFuntionMatrix(shape);
            b=N*elembody(i,:)';
            r_p_e=r_p_e+wtx*wty*N'*b*elemdata(i,5)*detJ;            
        end      
       end
    %create global index
    [index]= getGlobalIndex(ndof,nelem,nnel,A);
    %assembly body force
    r_p(index(i,:))=r_p(index(i,:))+r_p_e;    
end

