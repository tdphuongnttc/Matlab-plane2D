function [C,E,nuy,h]=materialMatrix(i,elemdata)
%create C,E matrix
%get h and nuy for every element
%C material matrix for plain stress
%E modul
%nuy: poisson ration
%h;thickness of plane

h = elemdata(i,5);
E = elemdata(i,6);
nuy = elemdata(i,7);
C = (E/(1-nuy^2))*[1 nuy 0 ; nuy 1 0; 0 0 (1-nuy)/2];