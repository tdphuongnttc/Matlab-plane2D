function [r_n]= loadVector(nelem,elemdata,elemtrac,elemtracID,nnel,nodeCoor,ndof,edof,ngp,r_n,sdof)
    %create traction vector
    %nelem : number of element
    %nnel: number of side of element
    %ndof: degree of node
    % degree of freedom of element
    
    A=elemdata(:,1:4); 
    [indexForce]= getGlobalIndex(ndof,nelem,nnel,A)
    [x,w]=Gauss(ngp);
    [xx,yy]=coordinateMatrix(A,nelem,nodeCoor);
for i=1:nelem
    r_ele=zeros(edof,1);
    r_ele_reduce=zeros(nnel,1);
    for side=1:nnel
       for ig=1:ngp %calculate integral
         pt = x(ig);%exact coordinates gauss point
         wt = w(ig); %exact weight number
        if(elemtracID(i,side)==1) % identify side have traction
            if (side==1)
            [shape]=createShapeFunction(pt,-1);
            [deriXi]=createShapeFunctionDerivativeXi(-1);
            [deriEta] = createShapeFunctionDerivativeEta(pt);
            [detJ,invJ,J] = Jacobian(deriXi,deriEta,xx(i,:),yy(i,:),nnel);
            [N] = ansatzFuntionMatrix(shape);
            [C,E,nuy,h]=materialMatrix(i,elemdata);
             r_ele_reduce=r_ele_reduce+wt*N(:,1:4)'*N(:,1:4)*elemtrac(i,(side*4-3):side*4)'*detJ*h;          
             r_ele(1:4,1)=r_ele_reduce(:);       
               elseif(side==2)
            [shape]=createShapeFunction(1,pt);
            [deriXi]=createShapeFunctionDerivativeXi(pt);
            [deriEta] = createShapeFunctionDerivativeEta(1);                   
            [detJ,invJ,J] = Jacobian(deriXi,deriEta,xx(i,:),yy(i,:),nnel);
            [N] = ansatzFuntionMatrix(shape);
            [C,E,nuy,h]=materialMatrix(i,elemdata);
             r_ele_reduce=r_ele_reduce+wt*N(:,3:6)'*N(:,3:6)*elemtrac(i,(side*4-3):side*4)'*detJ*h;          
             r_ele(3:6,1)=r_ele_reduce(:);         
              elseif(side==3)                      
             [shape]=createShapeFunction(pt,1);
            [deriXi]=createShapeFunctionDerivativeXi(1);
            [deriEta] = createShapeFunctionDerivativeEta(pt);
            [N] = ansatzFuntionMatrix(shape);
            [C,E,nuy,h]=materialMatrix(i,elemdata);
             r_ele_reduce=r_ele_reduce+wt*N(:,5:8)'*N(:,5:8)*elemtrac(i,(side*4-3):side*4)'*detJ*h;          
             r_ele(5:8,1)=r_ele_reduce(:)        
                else(side==4)
            [shape]=createShapeFunction(pt,-1);
            [deriXi]=createShapeFunctionDerivativeXi(pt);
            [deriEta] = createShapeFunctionDerivativeEta(-1);                      
            [detJ,invJ,J] = Jacobian(deriXi,deriEta,xx(i,:),yy(i,:),nnel);
            [N] = ansatzFuntionMatrix(shape);
            [C,E,nuy,h]=materialMatrix(i,elemdata);
             r_ele_reduce(1:2,1)=r_ele_reduce+wt*N(:,7:8)'*N(:,7:8)*elemtrac(i,(side*4-3):side*4)'*detJ*h;      
             r_ele_reduce(3:4,1)=r_ele_reduce+wt*N(:,1:2)'*N(:,1:2)*elemtrac(i,(side*4-3):side*4)'*detJ*h;          
             r_ele(7:8,1)=r_ele_reduce(1:2)
             r_ele(1:2,1)=r_ele_reduce(3:4)
          end
        end
      end        
    end
    %assembly force
    r_n(indexForce(i,:))=r_n(indexForce(i,:))+r_ele;     
end
