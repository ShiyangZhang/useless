function [ output_args ] = Copy_of_invP( ff,L,B,C,A,MM,rho,k,m,n,Pcurl ,L1  ,LM, LM1, Ark2MM,tolam,toll)

%INVP 此处显示有关此函数的摘要

%   此处显示详细说明







if k~=0

temp1=ff(n+1:n+m)-1/(rho-k^2)*C'*ff(1:n);
temp2=C'*ff(1:n)+k^2*ff(n+1:n+m);


else
    temp1=rho/(rho-k^2) *ff(n+1:n+m);
 temp2=zeros(m,1);
end
%Lcff=innnnvL(temp2,L,L1,m ,toll );

output_args=zeros(n+m,1);

output_args(1:n,1)=invAM(ff(1:n), A,  MM,L , C,Pcurl    , m ,rho,k  ,L1  ,LM, LM1, Ark2MM,tolam) +C*innnnvL(temp1,L,L1,m,toll   )  ;

output_args(n+1:n+m,1)=innnnvL(temp2,L,L1,m,toll   );









% function [inL]=innvL(args)

% if sum(abs(nonzeros(args)))==0

        % inL=sparse(m,1);

    % else

% tol = 1e-2;

% maxit = 200;

% [x11,fl1,~,~,~] = pcg(L,args,tol,maxit,L1,L1');

% if fl1~=0

% error('程序遇到错误，返回！');

% else

% inL=x11; 

% end

% end

% end



% function [inam]=inAM(args)

% if sum(abs(nonzeros(args)))==0

        % inam=sparse(n,1);

    % else

% tol = 1e-2;

% maxit = 200;

% [x11,fl1,~,~,~] = pcg(Ark2MM,args,tol,maxit,M1,M1');

     % if fl1~=0

% error('程序遇到错误，返回！');

        % else

% inam=x11; 

     % end

% end

% end



end

