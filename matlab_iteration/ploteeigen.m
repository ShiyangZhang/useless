
%%pcgfail
clear
pack
load G3.mat;
k=4;
rho=k^2+1;
FinA(1:n,1:n)=A-k^2*MM;
Arho=blkdiag(A+rho*B'*(L\B)-k^2*MM,eye(m));
tempblue=(real(eig((Arho+Arho')/2.0)));
tempred=(eig(FinA));



vblue=sort(tempblue(tempblue<0.3),'descend');
vred=sort(tempred(tempred<0.3),'descend');
[size_blue,~]=size(vblue);
[size_red,~]=size(vred);

%Distributions of eigenvalues smaller than 0.3  of the coefficient matrix $\mathcal{K}$ (red part)  and the matrix $A_{\eta}$ 
h1=plot(1:size_blue,vblue,'b+');
hold on;
h2=plot(1:size_red,vred,'r.');

h3= legend([h2 h1], 'For $$\mathcal{K}$$', 'For $$A_{\eta}$$');
set(h3,'Interpreter','latex');


%%comparasion -- equal
clear
pack
load G3.mat;
k=1.3;
rho=k^2+1;
ep=-1/(rho-k^2);

FinA(1:n,1:n)=A-k^2*MM;
preM=blkdiag(A+(rho-k^2)*MM,ep*L);
preM(1:n,n+1:n+m)=(1-rho*ep)*B';
prepK=blkdiag(preM(1:n,1:n)\(A+rho*B'*(L\B)-k^2*MM),eye(m));

tempblue=(real(eig(full(prepK))));
tempred=(real(eig(full(preM\FinA))));



vblue=sort(tempblue,'descend');
vred=sort(tempred,'descend');
[size_blue,~]=size(vblue);
[size_red,~]=size(vred);


plot(1:size_blue,vblue,'b+');
hold on;
plot(1:size_red,vred,'r.');


%%puG3diffk

clear
pack
load G3.mat;
kvec=[0,1,2,4];
for xunhuan=1:4
k=kvec(xunhuan);
rho=k^2+1;
ep=-1/(rho-k^2);

%%FinA(1:n,1:n)=A-k^2*MM;
prepK=blkdiag((A+(rho-k^2)*MM)\(A+rho*B'*(L\B)-k^2*MM),eye(m));

tempblue=(real(eig(full(prepK))));



vblue=sort(tempblue,'descend');
[size_blue,~]=size(vblue);

axis([1750 1800 -4.1 1.1]);
subplot(1,4,xunhuan)
plot(1750:size_blue,vblue(1750:size_blue),'b.');
end


%%pcgfail meshsize
clear
pack
load L4.mat;
mineig=zeros(20,1);
for xunhuan=1:1:20
k=2;
rho=xunhuan/10+k^2;
clear FinA C Pcurl
tempblue=eigs(A+rho*B'*(L\B)-k^2*MM, 1,'sm');
mineig(xunhuan)=tempblue
end
plot(1:20,mineig,'*')
xlabel('\eta-k^2 (k=4)')  
ylabel('smallest magnitude')  


