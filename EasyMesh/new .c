#include<stdio.h>
#include<string.h>
#include<math.h>
double MT[2][2]={{1.0,2.0},{1.0,1.0}};int p=5,k,i,j;double C[3][2]={{5.0,5.0},{6.0,5.0},{5.0,8.0}};double V[3][2];

int main(){
//对MT进行lu分解
for(p=0;p<2-1;p++)//主元所在行
{for(k=p+1;k<2;k++) //消去行循环
     { MT[k][p]= MT[k][p]/MT[p][p];


         for(i=p+1;i<2;i++)//后半行,就是除le第一个元素;
          MT[k][i]=MT[k][i]-MT[p][i]*MT[k][p];
	
	




	
	 }
}
double UVT[p];
//LUVT=CT
for(i=0;i<3;i++)//对CT列循环，也就是对C行循环;
{
UVT[0]=C[i][0];
for(j=0;j<2;j++)//行循环 C[i][0]
     {UVT [j] =    C[i][j];for(k=0;k<j;k++) UVT [j] -=   UVT [k]* MT[j][k]    ;}
	 

	 
 V[i][1]=UVT[1]/MT[1][1];
 for(p=0;p>-1;p--)
	       {V[i][p]=UVT[p];for(k=1;k>p;k--) V[i][p]-=MT[p][k]*V[i][k];V[i][p]=V[i][p]/MT[p][p];}

}



for(i=0;i<3;i++)  {printf("\n");for(j=0;j<2;j++)
printf("%f  ",V[i][j]);}
return 0;
}