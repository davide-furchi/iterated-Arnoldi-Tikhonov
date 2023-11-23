% mainstIAT.m - driver for stiatik = stationary iterated arnoldi-tikhonov with identity as regularization 
% operator 
% initialize environmnet
clear all; close all
format short e
iptsetpref('ImshowBorder','tigh');

% generate linear discrete ill-posed problem 
n = input('Enter the order of the matrix: ');
type = input('Select matrix: 1/baart, 2/phillips, 3/blur ');
if (type == 1)
   [A,b,xexact]=baart(n);
elseif (type == 2)
   [A,b,xexact]=phillips_alt(n); 
elseif (type == 3)
   [A,b,xexact]=blur(n);
   n=n^2; %dimension of the problem is n^2
end

%compute image b and norm of solution xexact
b=A*xexact;
normxexact=norm(xexact);

% generate perturbed rhs
seed=11;
randn('state',seed); 
relerr = input ('Enter relative norm of noise: '); %choose noise
err = randn(length(b),1);
err = relerr*norm(b)*err/norm(err);
b = b + err; 
delta=norm(err)

%compute the Arnoldi decomposition
m=input('Number of Arnoldi steps: ');
[h,Lambda,U,Imq,V,H,nrmRb2]=Arnoldi(A,m,b); 

% check the condition of equation (15)
E=norm(xexact); %set constant E in equation (14)
%E=3*E;  %E used in "Error estimates for Arnoldi-Tikhonov" for the AT
%method
Eh=E*h;
rhs=(Eh+delta)^2;
if (rhs<=nrmRb2)
 display('rhs bd o.k.')
else
 display('rhs bd violated')
end
pause

% choose number of iterated Tikohonov steps
i = input('Enter the number of steps of Tikhonov iteration: ');


% determine alpha
[alpha]=parameter(Lambda,U,Imq,V,m,b,rhs,i)

% compute Tikhonov solution and ERR
xcomp=0;
for k=1:i
  xcomp=(H'*H+alpha*eye(m))\((H'*V'*b)+alpha*xcomp); 
end
xcomp=V(:,1:m)*xcomp;
xcompmxexact=norm(xcomp-xexact)/norm(xexact)

% plot
ans = input('plot: 1/yes, 0/no');
if (ans == 1)
    t=-6+12*(0:n-1)/(n-1);
    plot(t,xexact,'r','LineWidth',1.5)
    hold
    plot(t,xcomp,'k','LineWidth',1.5)
end