function [ output_args ] = invP( ff,L,B,C,A,MM,rho,k,m,n,Pcurl ,L1  ,LM, LM1, Ark2MM,tolam,toll)
%INVP 此处显示有关此函数的摘要
%   此处显示详细说明





Lcff=innnnvL(C'*ff(1:n),L,L1,m ,toll );
output_args=zeros(n+m,1);temp=innnnvL(ff(n+1:n+m),L,L1,m,toll   );
output_args(1:n,1)=invAM(ff(1:n)-B'*Lcff, A,  MM,L , C,Pcurl    , m ,rho,k  ,L1  ,LM, LM1, Ark2MM,tolam) +C*temp  ;
output_args(n+1:n+m,1)=Lcff+k*k*temp;



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
