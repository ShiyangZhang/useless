function [inL]=innnnvL(args,L,L1,m,toll )
if sum(abs(nonzeros(args)))==0
        inL=sparse(m,1);
else

tol = toll;
maxit = 200;
[x11,fl1,~,~,~] = pcg(L,args,tol,maxit,L1,L1');
if fl1~=0
error('³ÌÐòÓöµ½´íÎó£¬·µ»Ø£¡');
else
inL=x11; 
end
end
end