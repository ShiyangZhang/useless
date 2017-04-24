%invL=inv(L);%clear%load L5.mat
%Pcurl=sparse(Pcurl(:,1)+1,Pcurl(:,2)+1,Pcurl(:,3),n,2*m);  %save L5-pcurl.mat
BTLB=B'*(L\(B));
kvec=[0,1,1.55,1.6,2,4];
results=zeros(5,6);

for xunhuan=1:6
    k=kvec(1,xunhuan);
rho=k^2+1;
FinA(1:n,1:n)=A-k^2*MM;
%FinA=sparse(FinA);

%invL(3,8)=invL(3,8)+1e-1;
Ark2MM=A+(rho-k*k)*MM;



  x=sparse(n+m,1);
  maxit=70;
  tol=1e-6;
  
 % rhs=sparse(FinF);
 rhs=ones(n+m,1);

  normrhs=norm(rhs);
r=rhs-(FinA*x);
z=invP(r,L,B,C,A,MM,rho,k,m,n ,Pcurl,ML );
p=z;
H=blkdiag(Ark2MM ,diag(ones(m,1)));
tt1=z'*H*z;
for kk=1:maxit
alpha=tt1/(p'*H*(invP(FinA*p,L,B,C,A,MM,rho,k,m,n ,Pcurl,ML )));
x=x+alpha*p;


r=r-alpha*(FinA*p);
z=invP(r,L,B,C,A,MM,rho,k,m,n,Pcurl,ML );
nor=norm(r)/normrhs;
% kk
% nor
% norm(invP*rhs-invP*FinA*x(:,kk+1))/norm(invP*rhs)

if nor<=tol
%%norm(r(:,kk+1))<=1e-10
break

%plotnorm(1,kk)=norm(rhs-FinA*x(:,kk+1))/norm(rhs);

end
tt0=tt1;
tt1=z'*H*z;
beta=tt1/tt0;
p=z+beta*p;
end
% 
%norm(invP*rhs-invP*FinA*x(:,kk+1))/norm(invP*rhs)
%G4 k=4,迭代次数：30





%rhs=[rand(1,n) zeros(1,m)]';
%以下跑 W-minres
%FinA=sparse(FinA);

v0=zeros(n+m,1);


 om=zeros(n+m,1);om0 =zeros(n+m,1);
     minresx=zeros(n+m,1);
    r=rhs-FinA*minresx; 
v=invP(r,L,B,C,A,MM,rho,k,m,n ,Pcurl,ML );
 
H=blkdiag(Ark2MM ,diag(ones(m,1)));
  ga=sqrt( v'*H* v);
  et=ga;

  s0=0;s=0;
  c=1;c0=1;
  
al=zeros(4,1);

for j=2:maxit
v=v./ga;
de=v'*H*invP(FinA*v,L,B,C,A,MM,rho,k,m,n ,Pcurl,ML );
vs=v;
v=invP(FinA*v,L,B,C,A,MM,rho,k,m,n,Pcurl,ML )-de*v-ga*v0;
v0=vs;
ga0=ga;
 ga=sqrt( v'*H* v);

 al(1,1)=c*de-c0*s*ga0;
 al(2,1)=sqrt(al(1,1)*al(1,1)+ga*ga);
 al(3,1)=s*de+c0*c*ga0;
 al(4,1)=s0*ga0;
 
 c0=c;
   c=al(1,1)/al(2,1);
   s0=s;
      s=ga/al(2,1);
      oms=om;
	  om=(v0-al(4,1)* om0-al(3,1)* om )./al(2,1);
      om0=oms;
minresx=minresx+c*et*om;
et=-s*et;
   
 nor=norm(rhs-FinA*minresx)/normrhs;
if nor<=tol

break

end
end
%j-1
%nor



%block
%FinA(1:n,1:n)=A-k^2*MM;rho=k^2+1;


P1=blkdiag(Ark2MM ,1/rho*L);
[x1,flag,relres,iter,resvec]=minres(FinA,rhs,tol,maxit,@(x)[innnnvm(x(1:n),A,MM,n,rho,k );rho*innnnvL(x(n+1:n+m),L,m )]);

%eigz=real(eig(A+rho*B'*invL*B-k^2*MM));
%eigy=real(eig(invP*FinA));
%z(1,xunhuan)=min(eigz);
%y(1,xunhuan)=min(eigy);

results(1,xunhuan)=kk;
results(2,xunhuan)=j-1;
results(3,xunhuan)=iter;



zz=rho*BTLB+FinA(1:n,1:n);
invPA=blkdiag((A+rho*MM-k^2*MM)\(zz) ,eye(m));


 results(4,xunhuan)=min(real(eig(full(zz))));
 results(5,xunhuan)= min(real(eig(full(invPA))));

end


results(1:5,:)
