function [detJ,invJ,J] = Jacobian(deriXi,deriEta,xx,yy,nnel)
%matrix jacobian for each element
%nnel number side of element
%xx X coordinate at element i
%yy Y coordinate at element i
%deriXi ; derivative X follow xi at element i
%deriEta ; derivative X follow eta at element i

J = zeros(2,2);
for i=1:nnel 
    J(1,1)=J(1,1)+deriXi(i)*xx(i);
    J(1,2)=J(1,2)+deriXi(i)*yy(i);
    J(2,1)=J(2,1)+deriEta(i)*xx(i);
    J(2,2)=J(2,2)+deriEta(i)*yy(i);
end
%determinate J
detJ=det(J);               
%inverse J
invJ=inv(J);


