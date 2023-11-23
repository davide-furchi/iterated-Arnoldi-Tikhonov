function [alpha] = parameter(Lambda,U,Imq,V,m,b,rhs,i)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
D=zeros(m+1);
b2=Imq*U'*V(:,1:m+1)'*b;
alpha=1;
% loop for determining regularization parameter alpha 
for k=1:m
  D(k,k)=1/(alpha+Lambda(k,k))^(2*i+1);
end
lhs=alpha^(2*i+1)*b2'*D*b2;
if lhs<rhs
  while lhs<rhs
    alphaold=alpha;
    alpha=2*alpha;
    for k=1:m
      D(k,k)=1/(alpha+Lambda(k,k))^(2*i+1);
    end
    lhs=alpha^(2*i+1)*b2'*D*b2;
  end
else
  while lhs>rhs
    alphaold=alpha;
    alpha=alpha/2;
    for k=1:m
      D(k,k)=1/(alpha+Lambda(k,k))^(2*i+1);
    end
    lhs=alpha^(2*i+1)*b2'*D*b2;
  end
end
% desired alpha-value between alpha and alphaold. refine by bisection
alphamax=max(alpha,alphaold);
alphamin=min(alpha,alphaold);
while abs(alpha-alphaold)>1e-5*(alpha+alphaold)
  alphamid=(alphamax+alphamin)/2;
  for k=1:m
     D(k,k)=1/(alphamid+Lambda(k,k))^(2*i+1);
  end
  lhs=alphamid^(2*i+1)*b2'*D*b2;
  if rhs>lhs
    alphamin=alphamid;
  else
    alphamax=alphamid;
  end
  alphaold=alpha;
  alpha=alphamid;
end
end