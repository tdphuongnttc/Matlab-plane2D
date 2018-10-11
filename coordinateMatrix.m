function [xx,yy]=coordinateMatrix(A,nelem,nodeCoor)
% create  coordinate xx and yy matrix with element
for i=1:nelem
xx(i,:)=nodeCoor(A(i,:),1);
yy(i,:)=nodeCoor(A(i,:),2);
end