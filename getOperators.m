function [B] =  getOperators(edof,nnel,deriX,deriY)
%create Operators Matrix
%edof=nnel*ndof; B=(3,8)
%nnel number side of element
%create operator matrix 3x8
%deriXi ; derivative X follow xi at element i
%deriEta ; derivative X follow eta at element i

%create matrix operator
B=zeros(3,edof);
for i=1:nnel
B(1,2*i-1)=deriX(i);
B(2,2*i-1)=0;
B(3,2*i-1)=deriY(i);
B(1,2*i)=0;
B(2,2*i)=deriY(i);
B(3,2*i)=deriX(i);
end

