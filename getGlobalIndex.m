function [index]= getGlobalIndex(ndof,nelem,nnel,A)
%create index matrix nelem x (nnelxndof);
%nelem : number of element
%nnel: number of side of element
index=zeros(nelem,nnel*ndof);
for i=1:nelem
    indicies=A(i,:);
           index(i,:)=[indicies(1)*2-1 indicies(1)*2 indicies(2)*2-1 indicies(2)*2 indicies(3)*2-1 indicies(3)*2 indicies(4)*2-1 indicies(4)*2]; 
 end
