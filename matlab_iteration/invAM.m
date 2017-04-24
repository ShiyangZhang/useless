function [ output_args ] = invAM( tttt, A,  MM,L , C,Pcurl   , m ,rho,k ,L1  ,LM, LM1, Ark2MM,tolam)

%INVAM 此处显示有关此函数的摘要

%   此处显示详细说明



%Ark2MM=A+(rho-k*k)*MM;

%L1 = ichol(L);

%LM=L+(rho-k^2)*ML;

%LM1=ichol(LM);



tol = tolam;

tollin=tolam;
tollin=1e-1;


maxit = 250;


[x11,fl1,~,iter,~] = pcg(Ark2MM,tttt,tol,maxit,@(x)hx1j(x)+hx2(x)+hx3(x));

if fl1~=0

error('程序遇到错误，返回！');

else

output_args=x11; 



end





%output_args=hx1j(tttt)+hx3(tttt);







function [ o1utput_args ] = hx1j( ttt )

%HX1 J-光滑

  

t=diag(diag(Ark2MM));

o1utput_args=t\ttt;

end

% 

% function [ output_args ] = hx1g( x )

% %HX1 G-光滑

% global  A  MM

% invtr=inv(tril(A+MM));

% tt=(tril(A+MM))'\x;

% output_args=(tril(A))\x+tt-(tril(A+MM))\(A*tt);

% 

% 

% end



function [ output_args ] = hx1g( x )

%HX1 J-光滑

tt=(tril(Ark2MM))'\x;

output_args=(tril(Ark2MM))\x+tt-(tril(Ark2MM))\(Ark2MM*tt);





end



function [ output_args ] = hx2( t )

%HX2 此处显示有关此函数的摘要

%   此处显示详细说明



output_args=Pcurl(:,1:m)*(invML(Pcurl(:,1:m)'*t))+Pcurl(:,m+1:2*m)*(invML(Pcurl(:,1+m:2*m)'*t));



end



function [ output_args ] = hx3( input_args )

%HX3 此处显示有关此函数的摘要

%   此处显示详细说明



output_args=C*(innvL(C'*input_args));



end



function [inL]=invML(args)

if sum(abs(nonzeros(args)))==0

        inL=sparse(m,1);

    else


maxit = 400;

[x11,fl1,~,~,~] = pcg(LM,args,tollin,maxit,LM1,LM1');

if fl1~=0

error('程序遇到错误，返回！');

else

inL=x11; 

end

end

end





function [inL]=innvL(args)

if sum(abs(nonzeros(args)))==0

        inL=sparse(m,1);

    else

maxit = 400;

[x11,fl1,~,~,~] = pcg(L,args,tollin,maxit,L1,L1');

if fl1~=0

error('程序遇到错误，返回！');

else

inL=x11; 

end

end

end



end



