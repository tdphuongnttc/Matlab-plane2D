function [rc,rf]=arrangeForce(gdof,bcDof,force)
%arrange rc and rf, the same 2D truss
activeDof=setdiff([1:gdof]',[bcDof]);
rf=force(activeDof);
rc=force(bcDof);