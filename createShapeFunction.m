function [shape]=createShapeFunction(xi,eta)
shape=[(1/4)*(1-xi)*(1-eta);(1/4)*(1+xi)*(1-eta);(1/4)*(1+xi)*(1+eta);(1/4)*(1-xi)*(1+eta)];