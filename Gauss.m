%...Generate Gauss points x(:) and weights w(:) of n-point formula
function [x,w]=Gauss(n)
eps = 3.d-14;
pi = acos(-1.0d0);
  
x1 = -1.0d0;
x2 = 1.0d0;
m = floor((n+1)/2);
xm = 0.5d0*(x2 + x1);
xl = 0.5d0*(x2 - x1);
for i = 1:m
    z = cos(pi*(i-0.25d0)/(n+0.5d0));
    p1 = 1.0d0;
    p2 = 0.0d0;
    for j=1:n
        p3 = p2;
        p2 = p1;
        p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j;
    end
    pp = n*(z*p1-p2)/(z*z-1.0d0);
    z1 = z;
    z = z1-p1/pp;
    while abs(z-z1) >= eps
      p1 = 1.0d0;
      p2 = 0.0d0;
      for j=1:n
        p3 = p2;
        p2 = p1;
        p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j;
      end
      pp = n*(z*p1-p2)/(z*z-1.0d0);
      z1 = z;
      z = z1-p1/pp;
    end
    x(i) = xm-xl*z;
    x(n+1-i) = xm+xl*z;
    w(i) = 2.0d0*xl/((1.0d0-z*z)*pp*pp);
    w(n+1-i) = w(i);
end