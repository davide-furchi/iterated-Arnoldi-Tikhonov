function [h,Lambda,U,Imq,V,H,nrmRb2] = Arnoldi(A,m,b)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
bnorm=norm(b);
V=b/bnorm;
j = 0;
H=zeros(m+1,m);
% arnoldi with initial vector V
while (j<m)
  j=j+1;
  w=A*V(:,j); 
  for k=1:j
    H(k,j)=dot(V(:,k),w);
    w=w-H(k,j)*V(:,k);
  end
  H(j+1,j)=norm(w);
  if H(j+1,j)<1e-12
     warning('*** H is numerically singular ')
     condH=cond(H)
     pause
  end
  V(:,j+1)=w/H(j+1,j);
% reorthogonalization
  v=V(:,j+1);
  for k=j:-1:1
    v=v-dot(V(:,k),v)*V(:,k);
  end
  V(:,j+1)=v/norm(v);
end
% define h
Am=V(:,1:m+1)*H*(V(:,1:m))';
h=norm(A-Am);
% determine matrices in Proposition 5.
[U,S,W]=svd(H);
sigmamaxH=S(1,1);
sigmaminH=S(m,m);
rankH=rank(S);
% determine q = rank(H)
q=rank(S);
for k=q+1:m
  S(k,k)=0;
end
Lambda=zeros(m+1);
for k=1:m
  Lambda(k,k)=S(k,k)^2;
end
Imq=zeros(m+1);
Imq(1:q,1:q)=eye(q);
% define Rm
Rm=V(:,1:m+1)*U*Imq*U'*V(:,1:m+1)';
nrmRb2=norm(Rm*b)^2;
end