function [A,b,x] = phillips_alt(n)
%PHILLIPS test problem.
%
% [A,b,x] = phillips_alt(n)
%
% Discretization of the first-kind Fredholm integral
% equation discussed by D. L. Phillips by a Nystrom method
% based on the trapezoidal rule. Define the function
%
%    phi(t) = 1 + cos(t*pi/3),  |t| <  3, 
%    phi(t) = 0,                |t| >= 3.
%
% Then the kernel K, the solution f, and the right-hand side
% g are given by:
%
%    K(s,t) = phi(s-t),
%    f(t)   = phi(t),
%    g(s)   = (6-|s|)*(1+.5*cos(s*pi/3)) + 9/(2*pi)*sin(|s|*pi/3).
%
% Both s and t live in the interval [-6,6].
 
% Reference: D. L. Phillips, "A technique for the numerical solution
% of certain integral equations of the first kind", J. ACM 9
% (1962), 84-97.
 
% Check input.
if (n < 2), error('The order n must be at least 2'), end
%
% Compute the matrix A.
h = 12/(n-1); r1 = zeros(1,n);
if (rem((n-1),4) == 0)
    for i = 1:((n-1)/4)
        r1(i) = h*(1 + cos((i-1)*h*pi/3));
    end
else
    for i = 1:(fix((n-1)/4)+1)
        r1(i) = h*(1 + cos((i-1)*h*pi/3));
    end
end
A = toeplitz(r1);
for i = 1:(fix((n-1)/4)+1)
    A(i,1) = A(i,1)/2;
    A(n-i+1,n) = A(n-i+1,n)/2;
end

% Compute the right-hand side b.
if (nargout>1),
  b = zeros(n,1);
  for i = 1:n
      s = -6 + (i - 1)*h;
      b(i) = (6 - abs(s))*(1+0.5*cos(s*pi/3)) + 9/(2*pi)*sin(abs(s)*pi/3);
  end
end

% Compute the solution x.
if (nargout==3),
  x = zeros(n,1);
  i1 = fix((n-1)/4)+2;
  i2 = fix(3*(n-1)/4)+1;
  for i = i1:i2
      s = -6 + (i - 1)*h;
      x(i)= 1 + cos(s*pi/3);
  end
end
%cond(A)
