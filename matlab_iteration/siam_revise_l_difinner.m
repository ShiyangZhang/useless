
%%to get a more reliable comparision, we run this script fot two times, and
%% only the last one is used
%%on the grid  L4 , we get the time for P-cg and W-Minres

clear
pack
load L4.mat;

%Pcurl=sparse(Pcurl(:,1)+1,Pcurl(:,2)+1,Pcurl(:,3),n,2*m);  

timesp=zeros(5,5);
timesm=zeros(5,5);

% qsub  -q qlong siam_revise.sh
% qsub  -q qshort matlab.sh
% qsub  -q qmedium matlab.sh
Ark2MM=A+MM;
%opts.type = 'nofill';opts.michol = 'on';
%L1 = ichol(L, opts);
tic
L1 = ichol(L);
setupl=toc;






%temph=zeros(n+m,1);

%invL=inv(L);

% results=zeros(5,6);
% times=zeros(5,6);
%cputimes=zeros(5,6);

itC=[1e-5,1e-4,1e-3,1e-2,1e-1];
itR=[1e-5,1e-4,1e-3,1e-2,1e-1];


%setuplm=zeros(1,6);
% LM=L+(rho-k^2)*ML;
% 
% tic
% LM1=ichol(LM);
% setuplm(xunhuan)=toc;


LM=L+ML;

tic
LM1=ichol(LM);
setuplm=toc;

    k=0;

% if 2*k^2<1
% rho=1;
% else rho=2*k^2;
% end
rho = k^2+1;

FinA(1:n,1:n)=A-k^2*MM;


for it_i=1:5
for it_j=1:5
tolam=itR(it_i);
toll=itC(it_j);
    






%Ark2MM=A+(rho-k*k)*MM;



 rhs=ones(n+m,1);
 
 tic

  normrhs=norm(rhs);

 


  x=sparse(n+m,1);

  maxit=350;

  tol=1e-6;

  






r=rhs-(FinA*x);

z=invP(r,L,B,C,A,MM,rho,k,m,n,Pcurl ,L1  ,LM, LM1, Ark2MM,tolam,toll);

p=z;

tt1=z(1:n)'*Ark2MM*z(1:n) + z(n+1:n+m)'*z(n+1:n+m);
%tt1=z'*H*z;

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

break;



end 

tt0=tt1;
tt1=z(1:n)'*Ark2MM*z(1:n) + z(n+1:n+m)'*z(n+1:n+m);
%tt1=z'*H*z;

beta=tt1/tt0;

p=z+beta*p;

end 

% 

%norm(invP*rhs-invP*FinA*x(:,kk+1))/norm(invP*rhs)

%G4 k=4,迭代次数�?0



%times(1,xunhuan)=toc;
timesp(it_i,it_j)=toc+setupl+setuplm;
%cputimes(1,xunhuan)=cputime-cputimes(1,xunhuan);




%times(2,xunhuan)=tic;
% tic
% cputimes(2,xunhuan)=cputime;



%rhs=[rand(1,n) zeros(1,m)]';

%以下�?W-minres

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


timesm(it_i,it_j)=toc+setupl+setuplm;







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

resultsp(it_i,it_j) =kk;
resultsm(it_i,it_j)=iter;


end 
end


times=timesm./timesp;

resultsp
resultsm

timesp
timesm


timesm./timesp


 
% LM=L+(rho-k^2)*ML;
% tic
% tt=ichol(LM);
% toc
% clear tt;
