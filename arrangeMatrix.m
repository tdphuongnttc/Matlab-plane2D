function [Kff,Kcc,Kucf]=...
    arrangeMatrix(K,nUc,nnode,constraint,ndof)
% arrange K, the same 2D truss
Uc=zeros(nUc,1); % ndisp =3
nUf=nnode*ndof-nUc; % nnode=5*2-3=7
Uf=zeros(nUf,1); %[ 0 0 0 0 0 0 0] bien cho bac tu do

for i=1:nUc % i chay tu 1 den 3
   index=ndof*(constraint(i,1)-1)+constraint(i,2);
   % i=1; index=2*(1-1)+1=1
   % i=2; index=2*(1-1)+2=2
   %i=3 , index=2*(3-1)+2=6
   Uc(i)=index; %Uc(1)=1; Uc(2)=2;Uc(3)=6;
end
 
k=1; %luc dau k =1
for i=1:nnode % i chay tu 1 den 5
 for j=1:ndof % j chay tu 1 den 2, doi lai 2 bac tu do va 3 bac tu do
 index =ndof*(i-1)+j; % i=1 ; index =2(1-1)+1=1;
 iUcc=0; % gan bien iUcc=0
       for m=1:nUc %k =1 , if index=Uc(1), ndisp=3
        if index == Uc(m)
         iUcc=1;
         break
        end    
       end
       
   if iUcc==0
    Uf(k)=index;
    k=k+1;
   end
  end
 end


 for i=1:nUf
     for j=1:nUf
         Kff(i,j)=K(Uf(i),Uf(j));
     end
 end
 
 for i=1:nUc
     for j=1:nUc
         Kcc(i,j)=K(Uc(i),Uc(j));
     end
 end
 
 for i=1:nUc
     for j=1:nUf
         Kucf(i,j)=K(Uc(i),Uf(j));
     end
 end
 