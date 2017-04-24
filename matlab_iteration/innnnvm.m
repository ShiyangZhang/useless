
function [inMM]=innnnvm(args,A,MM,n,rho,k )

if sum(abs(nonzeros(args)))==0
        inMM=sparse(n,1);
else
     Ark2MM=A+(rho-k*k)*MM;
M1 = ichol(Ark2MM);

tol = 1e-2;
maxit = 800;
[x11,fl1,~,~,~] = pcg(Ark2MM,args,tol,maxit,M1,M1');
if fl1~=0
error('³ÌÐòÓöµ½´íÎó£¬·µ»Ø£¡');
else
inMM=x11; 
end
end
end