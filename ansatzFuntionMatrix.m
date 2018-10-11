function [N] = ansatzFuntionMatrix(shape)
%create matrix with ansatzFunction
N=[shape(1) 0 shape(2) 0 shape(3) 0 shape(4) 0;0 shape(1) 0 shape(2) 0 shape(3) 0 shape(4)];