function [deriXi]=createShapeFunctionDerivativeXi(eta)
%derivative shape function follow xi
deriXi=[(-1/4)*(1-eta) (1/4)*(1-eta) (1/4)*(1+eta) (-1/4)*(1+eta)];