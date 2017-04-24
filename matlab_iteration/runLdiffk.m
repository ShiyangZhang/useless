
%out>matlab-out2.txt
clear

load L4.mat;
%load L1.mat;

%Pcurl=sparse(Pcurl(:,1)+1,Pcurl(:,2)+1,Pcurl(:,3),n,2*m);  
resultsp=zeros(5,5);
resultsm=zeros(5,5);

timesp=zeros(5,5);
timesm=zeros(5,5);

% qsub  -q qlong matlab-diffk.sh
% qsub  -q qshort matlab-diffk.sh
% qsub  -q qmedium matlab-diffk.sh
Ark2MM=A+(rho-k*k)*MM;
%%opts.type = 'nofill';opts.michol = 'on';
tic
L1 = ichol(L);
setupl=toc;





kvec=[0,1,1.2,1.25,2,4];



%invL=inv(L);

 results=zeros(2,6);
 times=zeros(2,6);
%cputimes=zeros(5,6);
tolams=[1e-5,1e-4,1e-3,1e-2,1e-1];
tolls=[1e-5,1e-4,1e-3,1e-2,1e-1];

for it_am=4:4
tolam=tolams(it_am);
for it_l=1:1
toll=tolls(it_l);


for xunhuan=1:6

    k=kvec(1,xunhuan);

if k^2<1
rho=1;
else rho=5.0/4.0*k^2;
end

LM=L+(rho-k^2)*ML;
tic
LM1=ichol(LM);
setuplm=toc;

FinA(1:n,1:n)=A-k^2*MM;

%FinA=sparse(FinA);



%invL(3,8)=invL(3,8)+1e-1;

Ark2MM=A+(rho-k*k)*MM;
H=blkdiag(Ark2MM ,diag(ones(m,1)));



 rhs=ones(n+m,1);
  normrhs=norm(rhs);

 
%times(1,xunhuan)=tic;
tic
%cputimes(1,xunhuan)=cputime;


  x=sparse(n+m,1);

  maxit=150;

  tol=1e-6;

  






r=rhs-(FinA*x);

z=invP(r,L,B,C,A,MM,rho,k,m,n,Pcurl ,L1  ,LM, LM1, Ark2MM,tolam,toll);

p=z;


tt1=z'*H*z;

for kk=1:maxit

%alpha=tt1/(p'*H*(invP(FinA*p,L,B,C,A,MM,rho,k,m,n ,Pcurl,ML )));

tempBp=B*p(1:n);
alpha=tt1/( (norm(p(n+1:n+m)))^2 + p(1:n)'*FinA(1:n,1:n)*p(1:n) +rho* tempBp'* innnnvL(tempBp,L,L1,m,toll)    );


x=x+alpha*p;





r=r-alpha*(FinA*p);

z=invP(r,L,B,C,A,MM,rho,k,m,n,Pcurl ,L1  ,LM, LM1, Ark2MM,tolam,toll);

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



times(1,xunhuan)=toc;
%timesp(it_am,it_l)=toc;
%cputimes(1,xunhuan)=cputime-cputimes(1,xunhuan);




%times(2,xunhuan)=tic;
% tic
% cputimes(2,xunhuan)=cputime;



%rhs=[rand(1,n) zeros(1,m)]';

%以下跑 W-minres

%FinA=sparse(FinA);



% v0=zeros(n+m,1);





 % om=zeros(n+m,1);om0 =zeros(n+m,1);

     % minresx=zeros(n+m,1);

    % r=rhs-FinA*minresx; 

% v=invP(r,L,B,C,A,MM,rho,k,m,n ,Pcurl,ML );

 

% %H=blkdiag(Ark2MM ,diag(ones(m,1)));

  % ga=sqrt( v'*H* v);

  % et=ga;



  % s0=0;s=0;

  % c=1;c0=1;

  

% al=zeros(4,1);



% for j=2:maxit

% v=v./ga;

% %de=v'*H*invP(FinA*v,L,B,C,A,MM,rho,k,m,n ,Pcurl,ML );

% tempBv=B*v(1:n);
% de=norm(v(n+1:n+m))^2 + v(1:n)'*FinA(1:n,1:n)*v(1:n) +rho* tempBv'* innnnvL(tempBv,L,m )    ;


% vs=v;

% v=invP(FinA*v,L,B,C,A,MM,rho,k,m,n,Pcurl,ML )-de*v-ga*v0;

% v0=vs;

% ga0=ga;

 % ga=sqrt( v'*H* v);



 % al(1,1)=c*de-c0*s*ga0;

 % al(2,1)=sqrt(al(1,1)*al(1,1)+ga*ga);

 % al(3,1)=s*de+c0*c*ga0;

 % al(4,1)=s0*ga0;

 

 % c0=c;

   % c=al(1,1)/al(2,1);

   % s0=s;

      % s=ga/al(2,1);

      % oms=om;

	  % om=(v0-al(4,1)* om0-al(3,1)* om )./al(2,1);

      % om0=oms;

% minresx=minresx+c*et*om;

% et=-s*et;

   

 % nor=norm(rhs-FinA*minresx)/normrhs;

% if nor<=tol



% break



% end

% end

% %j-1

% %nor

% times(2,xunhuan)=toc;
% cputimes(2,xunhuan)=cputime-cputimes(2,xunhuan);


%times(3,xunhuan)=tic;



%results(2,xunhuan)=j-1;



%block

%FinA(1:n,1:n)=A-k^2*MM;rho=k^2+1;




tic
[x1,flag,relres,iter,resvec]=minres(FinA,rhs,tol,maxit,@(x)[invAM(x(1:n), A,  MM,L , C,Pcurl, m ,rho,k  ,L1  ,LM, LM1, Ark2MM,tolam) ;rho*innnnvL(x(n+1:n+m),L,L1,m,toll )]);
times(2,xunhuan)=toc;
%timesm(it_am,it_l)=toc;
%cputimes(3,xunhuan)=cputime-cputimes(3,xunhuan);




% cputimes(4,xunhuan)=cputime;

% [x1,flag,relres,iter4,resvec]=gmres(FinA,rhs,20,tol,maxit,@(x)[innnnvm(x(1:n),A,MM,n,rho,k );rho*innnnvL(x(n+1:n+m),L,m )]);
% cputimes(4,xunhuan)=cputime-cputimes(4,xunhuan);


% cputimes(5,xunhuan)=cputime;

% [x1,flag,relres,iter5,resvec]=gmres(FinA,rhs,20,tol,maxit,@(x)(invP(x,L,B,C,A,MM,rho,k,m,n,Pcurl,ML )));
% cputimes(5,xunhuan)=cputime-cputimes(5,xunhuan);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

tic
%cputimes(2,xunhuan)=cputime;


v0=zeros(n+m,1);
z0=zeros(n+m,1);





 om=zeros(n+m,1);om0 =zeros(n+m,1);

     minresx=zeros(n+m,1);

    v=rhs-FinA*minresx; 

%v=[innnnvm(r(1:n),A,MM,n,rho,k );rho*innnnvL(r(n+1:n+m),L,m )];
z=invP(v,L,B,C,A,MM,rho,k,m,n,Pcurl ,L1  ,LM, LM1, Ark2MM);

 

%H=blkdiag(Ark2MM ,diag(ones(m,1)));

  ga=sqrt( z'*H* z);

  et=ga;



  s0=0;s=0;

  c=1;c0=1;

  

al=zeros(4,1);



for j=2:maxit

v=v./ga;
z=z./ga;

%de=v'*H*invP(FinA*v,L,B,C,A,MM,rho,k,m,n ,Pcurl,ML );

tempBz=B*z(1:n);
de=(norm(z(n+1:n+m)))^2 + z(1:n)'*FinA(1:n,1:n)*z(1:n) +rho* tempBz'* innnnvL(tempBz,L,L1,m )    ;


vs=v;

v=FinA*z-de*v-ga*v0;

v0=vs;


z0=z;
z=invP(v,L,B,C,A,MM,rho,k,m,n,Pcurl ,L1  ,LM, LM1, Ark2MM);


ga0=ga;

 ga=sqrt( z'*H* z);



 al(1,1)=c*de-c0*s*ga0;

 al(2,1)=sqrt(al(1,1)*al(1,1)+ga*ga);

 al(3,1)=s*de+c0*c*ga0;

 al(4,1)=s0*ga0;

 

 c0=c;

   c=al(1,1)/al(2,1);

   s0=s;

      s=ga/al(2,1);

      oms=om;

	  om=(z0-al(4,1)* om0-al(3,1)* om )./al(2,1);

      om0=oms;

minresx=minresx+c*et*om;

et=-s*et;


 nor=norm(rhs-FinA*minresx)/normrhs;

if nor<=tol



break



end






end  


   








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




times(2,xunhuan)=toc;
%cputimes(2,xunhuan)=cputime-cputimes(2,xunhuan);
results(2,xunhuan)=j-1;
%}

%eigz=real(eig(A+rho*B'*invL*B-k^2*MM));

%eigy=real(eig(invP*FinA));

%z(1,xunhuan)=min(eigz);

%y(1,xunhuan)=min(eigy);



% results(1,xunhuan)=kk;


% results(3,xunhuan)=iter;

results(1,xunhuan)=kk;
results(2,xunhuan)=iter;

% results(4,xunhuan)=iter4(1);
% results(5,xunhuan)=iter5(1);


end % endfor xunhuan 

end%  endfor it_am
end%  endfor it_l


results
times