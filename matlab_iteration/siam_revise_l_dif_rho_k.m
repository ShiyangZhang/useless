
%%to get a more reliable comparision, we run this script fot two times, and
%% only the last one is used
%%on the grid  L4 , we get the time for P-cg and W-Minres

clear
pack
load L4.mat;

%Pcurl=sparse(Pcurl(:,1)+1,Pcurl(:,2)+1,Pcurl(:,3),n,2*m);  

timespm=zeros(2,4);

% qsub  -q qlong siam_revise.sh
% qsub  -q qshort matlab.sh
% qsub  -q qmedium matlab.sh
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






% if 2*k^2<1
% rho=1;
% else rho=2*k^2;
% end
kvec=[0,1,2,4];

for xunhuan=1:4
for it_i=4:4
for it_j=2:2
tolam=itR(it_i);
toll=itC(it_j);
    
k=kvec(xunhuan);
rho=k^2+1;

FinA(1:n,1:n)=A-k^2*MM;

LM=L+(rho-k^2)*ML;
Ark2MM=A+(rho-k^2)*MM;

tic
LM1=ichol(LM);
setuplm=toc;


%Ark2MM=A+(rho-k*k)*MM;



 rhs=ones(n+m,1);

 
tic
  normrhs=norm(rhs);


  x=sparse(n+m,1);

  maxit=150;

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
timespm(1,xunhuan)=toc+setupl+setuplm;
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


timespm(2,xunhuan)=toc+setupl+setuplm;







resultspm(1,xunhuan) =kk;
resultspm(2,xunhuan)=iter;


end 
end
end


resultspm

timespm


timespm(2,:)./timespm(1,:)


 
% LM=L+(rho-k^2)*ML;
% tic
% tt=ichol(LM);
% toc
% clear tt;
