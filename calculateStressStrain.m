function [strain,stress]=calculateStressStrain(nelem,edof,U,nnel,ngp,elemdata,nodeCoor,ndof)
   
%calculate x and w matrix
A=elemdata(:,1:4);
[xx,yy]=coordinateMatrix(A,nelem,nodeCoor);
[index]= getGlobalIndex(ndof,nelem,nnel,A);  

[x,w]=Gauss(ngp);%..........
 strain=zeros(nelem*ngp*ngp*3,1); 
  stress=zeros(nelem*ngp*ngp*3,1); 
for i=1:nelem
    %calculate strain at node
    %for j=1:length(nodes(1,1))
      strainelement=zeros(ngp*ngp*3,1);
      stresselement=zeros(ngp*ngp*3,1);
      
       [C,E,nuy,h]=materialMatrix(i,elemdata);
        
       disp(['DISPLACEMENT, STRESS, STRAIN IN ELEMENT ',num2str(i)])
        U_element=U(index(i,:))
       for igx=1:ngp
        % xi and wtx number for outer intergral
            ptxi = x(igx);
           for igy=1:ngp 
             % eta and wty number for inner intergral
                eta = x(igy);          
                [deriXi]=createShapeFunctionDerivativeXi(eta);
                [deriEta] = createShapeFunctionDerivativeEta(ptxi);
                [detJ,invJ,J] = Jacobian(deriXi,deriEta,xx(i,:),yy(i,:),nnel);
                [deriX,deriY] = getShapeFunctionDerivatives(nnel,deriXi,deriEta,invJ);
                B =  getOperators(edof,nnel,deriX,deriY);
                disp(['stress and strain in elememt  ',num2str(i),'  xi=  ',num2str(ptxi),'  eta=  ',num2str(eta)]);
                strain1=B*U_element
                stress1=C*strain1
                strainelement((ngp*(igx-1)+igy)*3-2:(ngp*(igx-1)+igy)*3)= strainelement((ngp*(igx-1)+igy)*3-2:(ngp*(igx-1)+igy)*3)+strain1;
                stresselement((ngp*(igx-1)+igy)*3-2:(ngp*(igx-1)+igy)*3)=stresselement((ngp*(igx-1)+igy)*3-2:(ngp*(igx-1)+igy)*3)+stress1;
           end
       end
       
       strainelement=  strainelement
       stresselement= stresselement
       strain((i*3*ngp*ngp-(3*ngp*ngp-1)):(i*3*ngp*ngp))=strain((i*3*ngp*ngp-(3*ngp*ngp-1)):(i*3*ngp*ngp)) + strainelement;
       stress((i*3*ngp*ngp-(3*ngp*ngp-1)):(i*3*ngp*ngp))=stress((i*3*ngp*ngp-(3*ngp*ngp-1)):(i*3*ngp*ngp)) + stresselement;
end
