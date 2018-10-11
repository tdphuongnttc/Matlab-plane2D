% FEM code for 2D plane stress problems using quadrilateral elements
% clear memory
clear;
% clear screen
clc;
%==========================================================================
% Open input file for reading data
finput = fopen('input_prob911.txt');

% read total number of nodes
nnode = fscanf(finput,'%d',[1]); % tinh ra  nnode=9

% read nodal coordinates
for i=1:nnode
    nodeCoor(i,:) = fscanf(finput,'%f %f',[2])% the hien toa do cac nut
end

% read total number of elements
nelem = fscanf(finput,'%d',[1]); % nelem =4 the hien 04 element
for i=1:nelem
    % read element connectivity, thickness,E,nu
    elemdata(i,:) = fscanf(finput,'%d %d %d %d %f %f %f',[7]); %[1,2,5,4] [2,3,6,5] [4,5,8,7] [5,6,9,8]
    % read nodal values of body force
    elembody(i,:) = fscanf(finput,'%f %f %f %f %f %f %f %f',[8]);
    % read traction id of element sides (1: traction applied, 0:no traction)
    elemtracID(i,:) = fscanf(finput,'%d %d %d %d',[4]);
    % read nodal tractions (t1i,t2i,t1j,t2j) applied on 4 sides of element
    elemtrac(i,:) = fscanf(finput,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',[16]);
end

% read number of dofs under Dirichlet BC
nUc = fscanf(finput,'%d',[1]);  % 6 bac tu do
% read node number, component (1 or 2), prescribed displacement
for i=1:nUc
    constraint(i,:) = fscanf(finput,'%d %d %f',[3]); % [1:1,2] [4:1,2] [7:1,2] the hien cac rang buoc tai nut 1,4,7
end

%Close the input file
fclose(finput);
%%
nesscessaryInput;
%%
%calculate K
[K]= globalStiffnessMatrix(nelem,elemdata,nnel,nodeCoor,ndof,edof,ngp,K);
%arrange matrix
[Kff,Kcc,Kucf]= arrangeMatrix(K,nUc,nnode,constraint,ndof);
%calculate body force
[r_p]=forceVector(nelem,elemdata,elembody,nnel,nodeCoor,ndof,edof,sdof,ngp,K);
%loadVector
[r_n]= loadVector(nelem,elemdata,elemtrac,elemtracID,nnel,nodeCoor,ndof,edof,ngp,r_n,sdof)
%sum load vector 
r=r_p+r_n;

%%
%calculate U

bcDof=zeros(length(constraint(:,1)),1)
for i=1:length(constraint(:,1))
   bcDof(i) = (constraint(i,1)-1)*2+constraint(i,2)
end
[rc,rf]=arrangeForce(gdof,bcDof,r)


    U = zeros(gdof,1)
    Uc=U(bcDof)
    Uf=inv(Kff)*(rf-Kucf'*Uc)
    
 for i=1:length(bcDof)
    for j= 1:gdof
        if j == bcDof(i)
        U(j)=Uc(i);
        end
   end    
 end
 activeDof=setdiff([1:gdof]',[bcDof])
 for i=1:length(activeDof)
    for j= 1:gdof
        if j ==activeDof(i)
        U(j)=Uf(i);
        end
    end    
 end
U

%% calculate stress and strain
   [strain,stress]=calculateStressStrain(nelem,edof,U,nnel,ngp,elemdata,nodeCoor,ndof)