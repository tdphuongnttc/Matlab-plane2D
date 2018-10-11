% nelem : number of elements
% nnode : number of nodes
bcval = constraint(:,3);
ngp =2; % Number of Gauss points
nodes = elemdata(:,1:4); %node per element
ndof = 2; % degree of freedom of node
gdof=ndof*nnode;    %    gdof=ndof*nnode=2*9=18 total number of degrees of freedom % 2 nhan 9 = 18
nnel =4; %number of side
edof = nnel*ndof; %4*2=8 degree of freedom of element
sdof = nnode*ndof; %9*2=18 degree of freedeom of system
K = zeros(sdof,sdof); %global stiffness matrix
r_n = zeros(sdof,1); % %force vector for traction