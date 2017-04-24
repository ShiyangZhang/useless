//编译环境：windows 10 Chinese;notepad++;MinGW.
//gcc -std=c99
//********************************************************************************************************************************************
/*本程序为论文 Preconditioners for the discretized time-harmonic Maxwell
equations in mixed form  的混合有限元实现

********************************************************************************************************************************************
本程序只能在 在 u（x,y）=1-y^2,1-x^2), p = (1 ? x^2 )(1 ? y^2 )和the unit square(-1=<x=<1,-1=<y=<1,)下获得正确的右端项.
请忽略err1;请注意更改程序中的k值

********************************************************************************************************************************************
首先由EasyMesh(https://zsy.wodemo.com/file/387961 ; https://zsy.wodemo.com/cat/4506此压缩包包含所用网格)获得网格信息（EasyMesh中 ，对更细的网格，需修改EasyMesh的大小限制，如  #define MAX_NODES 8000）。生成以下三个文件
EasyMesh使用方法：编译EasyMesh后 在EasyMesh.exe所在文件夹按Shift+右键，打开  在此处打开命令窗口（W） 选项  输入 EasyMesh L3 +dxf  (此处L3代表L3.d,为EasyMesh输入文件)
C:\\Users\\zhang\\Desktop\\easymesh\\L3.n  
C:\\Users\\zhang\\Desktop\\easymesh\\L3.s  
C:\\Users\\zhang\\Desktop\\easymesh\\L3.e  
请在下面程序中搜索L3，更改您的相应文件名和文件位置
并且更改参数kk的值。（论文中k值），
生成相应的A，C FinF（右端项），M（MM） ,L， FinA(系数矩阵),写入程序中。
可改进程序：本程序对数据结构完全没有优化，另外建议改用matlab生成网格，以稀疏形式存储矩阵。
********************************************************************************************************************************************
以下是提取数据的matlab程序；double U[8000]用作计算误差， 可以不用管它。如果想计算误差(err1)，需要先运行一遍程序，运行matlab计算

%MatLab， R2014a
clear;clc;
FinA=load('D:\FinA.txt');FinF=load('D:\FinF.txt');MM=load('D:\MM.txt');L=load('D:\L.txt');C=load('D:\C.txt');A=load('D:\A.txt');
[n,m]=size(C);
B=FinA(n+1:m+n,1:n);


%以下需要在正确生成FinF的前提下使用，参阅上端说明。
x=FinA\FinF;
u=x(1:n);
%p的验证:导入输出的文件Xn.txt，注意文件夹位置
[xn1,xn2]=textread('D:\Xn.txt');
perr=(1.-xn1.*xn1).*(1.-xn2.*xn2)-x((n+1):(m+n));

%u的验证；
%将u值写入u.txt;
dlmwrite('C:\Users\zhang\Desktop\easymesh\u.txt',u,'precision', '%.16f', ...
'newline', 'pc');
%dlmwrite('C:\Users\zhang\Desktop\easymesh\u.txt',u,'precision', '%.16f', 'newline', 'pc');
%或者  ：  save 'C:\Users\zhang\Desktop\easymesh\u.txt' u -ascii;
然后再次运行该程序，求得正确err1值。注意每次算误差都需要运行以上matlab程序（注意程序中以下语句：FileU=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\u.txt" , "r")  ;）


注意可在matlab中修改k值，然后依靠一下代码即可完全生成有效的除FinF以外有效的数据。这样仅仅对kk=0 运行程序即可。
右端项可采用随机值。
k=1.5;FinA(1:n,1:n)=A-k^2*MM;


MM=sparse(MM);L=sparse(L);A=sparse(A);B=sparse(B);C=sparse(C);invL=inv(L);invL=sparse(invL);

********************************************************************************************************************************************


********************************************************************************************************************************************
若希望详细理解本文算法：本文是根据 Understanding and Implementing the Finite Element Method[Mark_S._Gockenbach] 此书写出possion方程 的标准有限元离散，
然后在此基础上，运用  《电磁场有限元方法》-金建铭  书中的最低阶edge元和刚度矩阵生成技巧计算相关结果。建议您也依照此顺序即可理解程序。
********************************************************************************************************************************************
Writed by 张仕洋 ，（请在此处增加名字说明），
任何意见，改进和询问请发邮件至  hydzhang@hotmail.com ，
2015.10.

修改：2015.12  增加err1的计算
修改：2016.3   

*/
#include<math.h>
#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h> 
#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef PI
#define PI    3.14159265359
#endif
#define MAX_NODES 88000
#define FREE 0
int sign(int a){if (a>0)return 1;else return -1;}

/*=========================================================================*/


 
double U[8000];  

//以下三个函数提取文件中的字符，整数，和string, writed by  Bojan NICENO.
int load_i(FILE *in, int *numb)
{
 char dum, dummy[128];

 for(;;)
  {fscanf(in,"%s", dummy);
   if(dummy[0]=='#' && strlen(dummy)>1 && dummy[strlen(dummy)-1]=='#') {}
   else if(dummy[0]=='#') {do{fscanf(in,"%c", &dum);} while(dum!='#');}
   else                   {*numb=atoi(dummy); break;} }
 
   return 0;
}

int load_d(FILE *in, double *numb)
{
 char dum, dummy[128];

 for(;;)
  {fscanf(in,"%s", dummy);
   if(dummy[0]=='#' && strlen(dummy)>1 && dummy[strlen(dummy)-1]=='#') {}
   else if(dummy[0]=='#') {do{fscanf(in,"%c", &dum);} while(dum!='#');}
   else                   {*numb=atof(dummy); break;} } return 0;
}

int load_s(FILE *in, char *string)
{
 char dum, dummy[128];

 for(;;)
  {fscanf(in,"%s", dummy);
   if(dummy[0]=='#' && strlen(dummy)>1 && dummy[strlen(dummy)-1]=='#') {}
   else if(dummy[0]=='#') {do{fscanf(in,"%c", &dum);} while(dum!='#');}
   else                   {strcpy(string, dummy); break;} } return 0;
}
/*-------------------------------------------------------------------------*/



int main(/*int argc, char *argv[]*/)
{
FILE *FileU;
FileU=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\u.txt" , "r")  ;

if(FileU==NULL)                                                                       
  {                                                                             
   printf("can't open nodes file\n");                                                 
   return -1;                                                                   
  } 
  
  
  
for(int i=0;fscanf(FileU,"%lf",&U[i])!=EOF;i++);
fclose(FileU);
  
  
  
  

   
//double pp=0.0;//原方程p项的系数，0或1
   
   


char dummy1[80];
FILE *nodes_file=NULL;FILE *sides_file=NULL,*elements_file=NULL;
int num_free=0,num_con=0;
  
  
  
  
//以下为.n文件的提取.*******************************************************************************************************************************************************************************  

nodes_file=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\L3.n","r");

 if(NULL==nodes_file)                                                                       
  {                                                                             
   printf("can't open nodes file\n");                                                 
   return -1;                                                                   
  }        
  
  
int dummyea=0;int dummyeb=0 ,dummymark,dummyfbe=0;int num_nodes,i,j;
load_i(nodes_file,&num_nodes);
 // printf("%d",num_nodes);system("pause");
 
struct node
{double x,y;
	int ptr;
}nodes[num_nodes];


	  
	  
	  int num_freedummy=0;
	  for(i=0;i<num_nodes;i++)
	  {
	  load_s(nodes_file, dummy1);     
       load_d(nodes_file, &nodes[i].x);
       load_d(nodes_file, &nodes[i].y);
       load_i(nodes_file, &nodes[i].ptr);//记录节点坐标;
	
	   if(nodes[i].ptr==FREE) {num_free++;}
	  }
	  
	  num_con=num_nodes-num_free;
	   
	   
	   
 struct nptrs
{int fnodeptrs[num_free];
int cnodeptrs[num_con];
}nodeptrs;
double g[num_con];double F[num_free];//在获得num_free和num_con后申明相关量.
memset(F,0,sizeof(double)*num_free); 


FILE *outXn;


outXn = fopen("D:\\zsy\\Xn.txt", "w");
 if (!outXn)
    {
        perror("cannot open file  outXn");
        
    } //输出自由变量得坐标值
	 





	   
for(i=0;i<num_nodes;i++)
{
	if(nodes[i].ptr==FREE)
	      {nodes[i].ptr=num_freedummy+1;
          nodeptrs.fnodeptrs[num_freedummy]=i;
		  fprintf(outXn, "%16.12f  %16.12f\n ", nodes[i].x,nodes[i].y);   //输出自由变量的坐标值
          num_freedummy++;}
	  
	   
	      else  {nodes[i].ptr=-(i-num_freedummy+1);
	            nodeptrs.cnodeptrs[i-num_freedummy]=i;//记录节点指针和自由节点,强制节点;负值表示排序为多少位的强制节点
		        /* g[i-num_freedummy]=sin(PI*nodes[i].x)*sin(PI*nodes[i].y)+nodes[i].x; */}
}    //此处为Dirichlet条件g表达式的输入点.nodes[i].x node[i].y为坐标值.
	   

fclose(outXn);
	
		

			
		
		
fclose(nodes_file)	;
		
		
	







//以下为.s文件的提取.*******************************************************************************************************************************************************************************  

int num_sides;
sides_file=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\L3.s","r");
if(NULL==sides_file)                                                                       
  {                                                                             
   printf("can't open sides file\n");                                                 
   return -1;                                                                   
  } 





load_i(sides_file,&num_sides);
	  
  
struct sid
{ int x;int y;int p,q;int ptr;
}sides[num_sides];//p,q指示边的两边的elements.

int interiorsides[num_sides];


double length[num_sides];memset(length,0,sizeof(double)*num_sides); 

	 


FILE *outXs;


outXs= fopen("D:\\zsy\\Xs.txt", "w");
 if (!outXs)
    {
        perror("cannot open file  Xs");
        
    }


int fbndyedges[num_sides],dummy_fbe=0;
int interior=0;
for(i=0;i<num_sides;i++)
{
	  load_s(sides_file, dummy1);     
       load_i(sides_file, &sides[i].x);
       load_i(sides_file, &sides[i].y);//记录节点坐标;
	   
	   load_i(sides_file, &dummyea);
	   load_i(sides_file, &dummyeb);	
        load_i(sides_file, &dummymark);		   
	   if(dummyea!=-1&&dummyeb!=-1){fprintf(outXs, "%f  %f %f %f  ", nodes[sides[i].x].x,nodes[sides[i].x].y,nodes[sides[i].y].x,nodes[sides[i].y].y);  fputc('\n',outXs); // printf("dsfdsadfasdfasdf %d",dummymark);system("pause");
		   
		   sides[i].ptr=interior+1;interiorsides[interior]=i;sides[i].p=dummyea+1;sides[i].q=dummyeb+1;length[interior]=sqrt((nodes[sides[i].x].x-nodes[sides[i].y].x)*(nodes[sides[i].x].x-nodes[sides[i].y].x)+(nodes[sides[i].x].y-nodes[sides[i].y].y)*(nodes[sides[i].x].y-nodes[sides[i].y].y));/* printf("%f",length[interior]);system("pause"); */interior++;}
	   else 
	   {  // printf("22222222222222222222 %d",dummymark);
		   
		   sides[i].ptr=0;
		   if(dummyea=-1)
	       {sides[i].p=dummyeb+1;
             if(dummymark==FREE)
			   {
			    sides[i].q=-(dummyfbe+1);
			    fbndyedges[dummyfbe]=i;
			    dummyfbe++;
			   }
		   else 
             {sides[i].q=0;}
	      }
	     else 
		    {sides[i].p=dummyea+1;
             if(dummymark==FREE)
			   {
			    sides[i].q=-(dummyfbe+1);
		 	   fbndyedges[dummyfbe]=i;
		 	   dummyfbe++;
		 	   }
		    else 
            {sides[i].q=0;}
	       }
	    
        }
}
 int num_fbe=dummy_fbe;  double ss;

fclose(sides_file);

fclose(outXs);
 
 
 
 
 
 
 
 
 
 
//以下为.e文件的提取.*******************************************************************************************************************************************************************************  	  

int dummya ,dummyb,dummyc,dummyd,num_elements;
elements_file=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\L3.e","r");

       
 if(NULL==elements_file)                                                                       
  {                                                                             
   printf("can't open elements file\n");                                                 
   return -1;                                                                   
  }      


load_i(elements_file,&num_elements);
struct ele
 {
  int x,y,z;
 }
elem[num_elements];



	  for(i=0;i<num_elements;i++)
	  {
	  load_s(elements_file, dummy1); 
      load_i(elements_file, &dummyea);   	
      load_i(elements_file, &dummyea); 
      load_i(elements_file, &dummyea);  
      load_i(elements_file, &dummyea);   
      load_i(elements_file, &dummyea);   
      load_i(elements_file, &dummyea);     
	  
      load_i(elements_file,  &dummya);   
      load_i(elements_file,  &dummyb);   
      load_i(elements_file,  &dummyc); 
	  
	  
//以下用按逆时针相互衔接的顺序将三边记录于elem中,负数代表方向与边的本来的方向相反	  
	  
if(sides[dummya].y!=sides[dummyb].x&&sides[dummya].x!=sides[dummyb].x)
    {ss=(nodes[sides[dummya].x].x-nodes[sides[dummyb].x].x)*(nodes[sides[dummya].y].y-nodes[sides[dummyb].x].y)-(nodes[sides[dummya].x].y-nodes[sides[dummyb].x].y)*(nodes[sides[dummya].y].x-nodes[sides[dummyb].x].x);
	  
    if(ss>0)//判断点在向量的左侧
	        {if(sides[dummya].y==sides[dummyb].y)
				 { elem[i].x=dummya+1; elem[i].y=-dummyb-1;   if(sides[dummyb].x==sides[dummyc].x) elem[i].z=dummyc+1;else elem[i].z=-dummyc-1;}
			else { elem[i].x=dummya+1; elem[i].z=dummyb+1;   if(sides[dummyb].x==sides[dummyc].y) elem[i].y=dummyc+1;else elem[i].y=-dummyc-1; }
   			}
	else
	        {if(sides[dummya].x==sides[dummyb].y)
				 { elem[i].x=-dummya-1; elem[i].y=-dummyb-1;   if(sides[dummyb].x==sides[dummyc].x) elem[i].z=dummyc+1;else elem[i].z=-dummyc-1;}
			else { elem[i].x=-dummya-1; elem[i].z=dummyb+1;   if(sides[dummyb].x==sides[dummyc].y) elem[i].y=dummyc+1;else elem[i].y=-dummyc-1; }
   			}
			  
    }
else 
   {ss=(nodes[sides[dummya].x].x-nodes[sides[dummyb].y].x)*(nodes[sides[dummya].y].y-nodes[sides[dummyb].y].y)-(nodes[sides[dummya].x].y-nodes[sides[dummyb].y].y)*(nodes[sides[dummya].y].x-nodes[sides[dummyb].y].x);
	  
    if(ss>0)//判断点在向量的左侧
	        {if(sides[dummya].y==sides[dummyb].x)
				 { elem[i].x=dummya+1; elem[i].y=dummyb+1;   if(sides[dummyb].y==sides[dummyc].x) elem[i].z=dummyc+1;else elem[i].z=-dummyc-1;}
			else { elem[i].x=dummya+1; elem[i].z=-dummyb-1;   if(sides[dummyb].y==sides[dummyc].y) elem[i].y=dummyc+1;else elem[i].y=-dummyc-1; }
   			}
	else
	        {if(sides[dummya].x==sides[dummyb].x)
				 { elem[i].x=-dummya-1; elem[i].y=dummyb+1;   if(sides[dummyb].y==sides[dummyc].x) elem[i].z=dummyc+1;else elem[i].z=-dummyc-1;}
			else { elem[i].x=-dummya-1; elem[i].z=-dummyb-1;   if(sides[dummyb].y==sides[dummyc].y) elem[i].y=dummyc+1;else elem[i].y=-dummyc-1; }
   			}
			  
    }


	  
	 //单元上的边,正负标明方向.

	  
	  load_s(elements_file, dummy1); 
	  load_s(elements_file, dummy1); 
	  load_i(elements_file, &dummyd); 
 
  }
  fclose(elements_file);
  printf("共有点:%d\n边:%d\n元素:%d\n自由边界边:%d\n自由节点m=%d\n内部边n=%d\nn+m=%d\n",num_nodes,num_sides,num_elements,num_fbe,num_free,interior,num_free+interior);   system("pause");

  
  


//数据提取完毕,以下为算法第一部分:刚度矩阵的计算
//double A[interior][interior]; 



double* * A=NULL;
   A = (double **)malloc(sizeof(double *)*interior);
	for (i=0; i<interior; i++)
	{A[i] = (double *)malloc(sizeof(double)*interior); memset(A[i],0,interior*sizeof(double)); } 





















//double MM[interior][interior],B[num_free][interior];


double* * MM=NULL;
   MM = (double **)malloc(sizeof(double *)*interior);
	for (i=0; i<interior; i++)
	{MM[i] = (double *)malloc(sizeof(double)*interior); memset(MM[i],0,interior*sizeof(double)); } 



double* * B=NULL;
   B = (double **)malloc(sizeof(double *)*num_free);
	for (i=0; i<num_free; i++)
	{B[i] = (double *)malloc(sizeof(double)*interior); memset(B[i],0,interior*sizeof(double)); } 







double* * C=NULL;
   C = (double **)malloc(sizeof(double *)*interior);
	for (i=0; i<interior; i++)
	{C[i] = (double *)malloc(sizeof(double)*num_free); memset(C[i],0,num_free*sizeof(double)); } 









//memset(A,0,sizeof(double)*interior*interior);memset(MM,0,sizeof(double)*interior*interior);memset(B,0,sizeof(double)*num_free*interior);
  
int eptr[3]={0},indices[3]={0},ptrs[3]={0};
double delta,coords[3][2],I,x0,y0,invM[2][3];double det,G[3][3];
double kk=0;int s,r,ceshi=-1;double err1=0;

//double FinA[interior+num_free][interior+num_free],FinF[interior+num_free];memset(FinA,0,sizeof(double)*(interior+num_free)*(interior+num_free));memset(FinF,0,sizeof(double)*(interior+num_free));//输出最终矩阵,只是用于放入matlab验证算法.实际需要给出稀疏存储.





double* * FinA=NULL;
   FinA = (double **)malloc(sizeof(double *)*(interior+num_free));
	for (i=0; i<(interior+num_free); i++)
	{FinA[i] = (double *)malloc(sizeof(double)*(interior+num_free)); memset(FinA[i],0,(interior+num_free)*sizeof(double)); } 

double* * Pcurl=NULL;
   Pcurl = (double **)malloc(sizeof(double *)*(interior));
	for (i=0; i<interior; i++)
	{Pcurl[i] = (double *)malloc(sizeof(double)*2*(num_free)); memset(Pcurl[i],0,2*(num_free)*sizeof(double)); } 




double* * L=NULL;
   L = (double **)malloc(sizeof(double *)*num_free);
	for (i=0; i<num_free; i++)
	{L[i] = (double *)malloc(sizeof(double)*num_free); memset(L[i],0,num_free*sizeof(double)); } 






double* FinF=NULL;

	FinF = (double *)malloc(sizeof(double)*(interior+num_free)); memset(FinF,0,(interior+num_free)*sizeof(double)); 




/* 


double* truleul=NULL;

	FinF = (double *)malloc(sizeof(double)*(interior)); memset(trueul,0,(interior)*sizeof(double));  */





















int k;











double localength[3];
double aaa[3];//记录m逆矩阵的首行,也就是有限元基的常数项
double fff[3][3];double emm[3][3];
memset(emm,0,sizeof(emm));



double gg[interior];memset(gg,0,sizeof(double)*interior);	


double umid[2]={0.0,0.0};




for(k=0;k<num_elements;k++)  //对三角单元循环
//以下为algorithm 7.1 *************************
//此处代码只有上帝能看懂= - =
{	umid[0]=0.0;umid[1]=0.0;


eptr[0]=elem[k].x;//三条边,eptr指示真实的三角形逆时针三边}
eptr[1]=elem[k].y;
eptr[2]=elem[k].z;

if(eptr[0]>0)
	     {indices[0]=sides[eptr[0]-1].x;indices[1]=sides[eptr[0]-1].y;  if(sides[eptr[0]-1].ptr>0) localength[0]=length[sides[eptr[0]-1].ptr-1]; else localength[0]=-1;}
 else {indices[0]=sides[-eptr[0]-1].y;indices[1]=sides[-eptr[0]-1].x;if(sides[-eptr[0]-1].ptr>0) localength[0]=length[sides[-eptr[0]-1].ptr-1]; else localength[0]=-1;}
		   
if(eptr[1]>0)
		 {indices[2]=sides[eptr[1]-1].y;if(sides[eptr[1]-1].ptr>0) localength[1]=length[sides[eptr[1]-1].ptr-1]; else localength[1]=-1;}
else {indices[2]=sides[-eptr[1]-1].x;if(sides[-eptr[1]-1].ptr>0) localength[1]=length[sides[-eptr[1]-1].ptr-1]; else localength[1]=-1;}
	 
if(eptr[2]>0) {if(sides[eptr[2]-1].ptr>0) localength[2]=length[sides[eptr[2]-1].ptr-1]; else localength[2]=-1;}
else {if(sides[-eptr[2]-1].ptr>0) localength[2]=length[sides[-eptr[2]-1].ptr-1]; else localength[2]=-1;}

for(i=0;i<3;i++)
             {ptrs[i]=nodes[indices[i]].ptr;}//此处ptrs继承ptr得特性 没有更改.

for(i=0;i<3;i++)
             {coords[i][0]=nodes[indices[i]].x;
		 }//三条边}



for(i=0;i<3;i++)
{coords[i][1]=nodes[indices[i]].y;
//printf("\n  %d \n ",indices[i]);system("pause");
}	
	


//*****************************矩阵逆的计算

det=coords[1][0]*coords[2][1] -    coords[2][0]*coords[1][1]  + coords[2][0]*coords[0][1]-coords[0][0]*coords[2][1]  + coords[0][0]*coords[1][1] - coords[1][0]*coords[0][1] ;
//printf("%.12f aaaaa",det);system("pause");  

  invM[0][0]=(coords[1][1]-coords[2][1])/det;
  invM[1][0]=(coords[2][0]-coords[1][0])/det;
  
    invM[0][1]=(coords[2][1]-coords[0][1])/det;
	  invM[1][1]=(coords[0][0]-coords[2][0])/det;
	  
	    invM[0][2]=(coords[0][1]-coords[1][1])/det;
		  invM[1][2]=(coords[1][0]-coords[0][0])/det;
  

  
 
aaa[0]=(coords[1][0]*coords[2][1]-coords[1][1]*coords[2][0])/det;
aaa[1]=(coords[2][0]*coords[0][1]-coords[2][1]*coords[0][0])/det;
aaa[2]=(coords[0][0]*coords[1][1]-coords[0][1]*coords[1][0])/det;
  
  
  
  //if(abs(aaa[0]+aaa[2]+aaa[1]-1)>0.00000000000001) {printf("%f  s  s ",aaa[0]+aaa[2]+aaa[1]);system("pause");}
/* 

FILE *outn;


    outn = fopen("D:\\zsy\\length.txt", "w");
 if (!outn)
    {
        perror("cannot open file");
        
    } 
	  for (i = 0; i<3; i++)
   { for (j = 0; j<2; j++)
       
	fprintf(outn, "%f ", coords[i][j]); fputc('\n',outn);
	
	
	
	}
        
    
    fclose(outn);system("pause");  */



//ceshi :delta=0.5*(aaa[0]+aaa[1]+aaa[1]);
/*     system("pause");
	if(abs(1-(aaa[0]+aaa[1]+aaa[2]))<0.001) printf("right ");else printf("cuowu ");
	system("pause"); */

//******************************algorithm 7.2

  
  for(r=0;r<3;r++)
  {for(s=r;s<3;s++)	 
	  G[r][s]=invM[0][r]*invM[0][s]+invM[1][r]*invM[1][s];
  }

   
  for(r=1;r<3;r++)
  {for(s=0;s<r;s++)	 
	  G[r][s]=G[s][r];
  }
  
  
  
    

  
  
  
  
  
/*   
  printf("\n  %d Run read with the command: \n read namen.n names.s namee.e\n ",k);system("pause");

 
FILE *outn;


    outn = fopen("D:\\zsy\\Fn.txt", "w");
 if (!outn)
    {
        perror("cannot open file");
        
    } 
	  for (i = 0; i<3; i++)
    {for (j = 0; j<3; j++){
       
	fprintf(outn, "%f ",G[i][j]);}
        
        fputc('\n',outn);
    }
    fclose(outn);
system("pause");
   */
  
delta=0.5 * det;  
/* printf("%f ",delta);system("pause"); */

  x0=(coords[0][0]+coords[1][0]+ coords[2][0])/3.0       ;y0=(coords[0][1]+coords[1][1]+ coords[2][1])/3.0 ;
//系数k的输入点----------------------------------------------------------------------------------------------------------------------------------------------------------
I=delta;
   
  
  for(r=0;r<3;r++)
  {for(s=r;s<3;s++)	 
	  
	  {if(ptrs[r]>0&&ptrs[s]>0)
		  {		  i=min(ptrs[r],ptrs[s])-1;j=max(ptrs[r],ptrs[s])-1;
		  L[i][j]+=G[r][s]*I;}
	    
	  }
}  
  
  
for(r=0;r<3;r++)
{

	for(s=r;s<3;s++)	 
	  
	  {
		  if(localength[r]>=0&&localength[s]>=0)
		  {/* printf("r=%d;s= %d;",r,s); */i=min(sides[abs(eptr[r])-1].ptr-1,sides[abs(eptr[s])-1].ptr-1);j=max(sides[abs(eptr[r])-1].ptr-1,sides[abs(eptr[s])-1].ptr-1);/* printf("%d %d",i,j);system("pause"); */
		 A[i][j]+= sign(eptr[r])*sign(eptr[s])/(delta);/* printf("%f",A[i][j]); */}
	   
	  }
}  

 for(r=0;r<3;r++)
  {for(s=0;s<3;s++)	 
  fff[r][s]=invM[0][r]*invM[0][s]+invM[1][r]*invM[1][s];}



  
  
  
 emm[0][0]=(fff[1][1]-fff[0][1]+fff[0][0])*(delta)/6;

 emm[0][1]=sign(eptr[0])*sign(eptr[1])*(fff[1][2]-fff[1][1]-2*fff[0][2]+fff[0][1])*(delta)/12;
 
 emm[0][2]=sign(eptr[0])*sign(eptr[2])*(fff[1][0]-2*fff[1][2]-fff[0][0]+fff[0][2])*(delta)/12; 

 
 emm[1][1]=(fff[2][2]-fff[1][2]+fff[1][1])*(delta)/6;
   
 emm[1][2]=sign(eptr[2])*sign(eptr[1])*(fff[2][0]-fff[2][2]-2*fff[1][0]+fff[1][2])*(delta)/12;

 emm[2][2]=(fff[0][0]-fff[0][2]+fff[2][2])*(delta)/6;
  
  
    for(r=0;r<3;r++)
  {for(s=r;s<3;s++)	 
	  
      {if(localength[r]>=0&&localength[s]>=0)
		 
	 
		  {		  i=min(sides[abs(eptr[r])-1].ptr-1,sides[abs(eptr[s])-1].ptr-1);j=max(sides[abs(eptr[r])-1].ptr-1,sides[abs(eptr[s])-1].ptr-1);
  MM[i][j]+=emm[r][s];}
      }
  }
  
  
  
  


  
  
  

  

  
  
  
  
for(r=0;r<3;r++)
  {for(s=0;s<3;s++)	// (s=0 s=1 )
	  
      {if(localength[s]>=0&&ptrs[r]>0)
		  {j=sides[abs(eptr[s])-1].ptr-1;//printf("(%d)===%d,%d===",r,s,j);
	       i=ptrs[r]-1;
	       B[i][j]+= delta*sign(eptr[s])*(  invM[0][r]*(aaa[s]*invM[0][(1+s)%3]-aaa[(1+s)%3]*invM[0][s]+y0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s]))	+
	  invM[1][r]* (aaa[s]*invM[1][(1+s)%3]-aaa[(1+s)%3]*invM[1][s]+x0*(invM[1][(1+s)%3]*invM[0][s]-invM[1][s]*invM[0][(1+s)%3])	)   );
	      }
	  }
  }
  
 //  printf("%f,%f;%f,%f\n",umid[0],umid[1],1-y0*y0,1-x0*x0);system("pause");
  
for(s=0;s<3;s++)
{ if(localength[s]>=0)
	{j=sides[abs(eptr[s])-1].ptr-1;
	   umid[0] += sign(eptr[s])*localength[s]*  U[j]*(aaa[s]*invM[0][(1+s)%3]-aaa[(1+s)%3]*invM[0][s]+y0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s]));
	  
	  
	  umid[1] += sign(eptr[s])*localength[s]* U[j]* (aaa[s]*invM[1][(1+s)%3]-aaa[(1+s)%3]*invM[1][s]+x0*(invM[1][(1+s)%3]*invM[0][s]-invM[1][s]*invM[0][(1+s)%3]))  ;


/* printf("%d,%f,%f\n",j,umid[0],umid[1]);system("pause"); */
}
  

  
} 

 // err1+=(sqrt(pow(1-y0*y0-umid[0],2)+pow(1-x0*x0-umid[1],2)))/sqrt(pow(1-y0*y0,2)+pow(1-x0*x0,2))/(num_elements); // printf(";%f,%f;%f,%f;%f\n",umid[0],umid[1],1-y0*y0,1-x0*x0,(fabs(1-y0*y0-umid[0])+fabs(1-x0*x0-umid[1]))/(1-y0*y0+1-x0*x0));system("pause");


//三角基元对刚度矩阵循环完毕,但是同时要对load矩阵循环,所以不加大括号 
 

  
//以下为算法第二部分:load矩阵的初始化和计算.*******************************************************************************************************************************************************************************  	  



//以下为algorithm 7.5的第一部分************************************



 for(s=0;s<3;s++)

      {if(localength[s]>=0)
		  {i=sides[abs(eptr[s])-1].ptr-1;
	       
	          gg[i]+= delta*sign(eptr[s])*(( 2-(kk*kk+2*x0)*(1-y0*y0) )  *(aaa[s]*invM[0][(1+s)%3]-aaa[(1+s)%3]*invM[0][s]+y0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s]))	+
               ( 2-(kk*kk+2*y0)*(1-x0*x0) )             * (aaa[s]*invM[1][(1+s)%3]-aaa[(1+s)%3]*invM[1][s]-x0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s])	)   );
           }
     }
  
  
  
/* 	printf("k=%d   %f   gg=%f",k,invM[0][(1+2)%3],gg[k]);system("pause");
 */


}    //对三角基元循环完毕.

//以下计算C




 for(i=0;i<interior;i++)//边循环
 {
	//ptr>0指示自由节点
if(nodes[sides[interiorsides[i]].x].ptr>0)
 C[i][nodes[sides[interiorsides[i]].x].ptr-1]=-1.0;
if(nodes[sides[interiorsides[i]].y].ptr>0)
 C[i][nodes[sides[interiorsides[i]].y].ptr-1]=1.0;


}

//边循环
 for(i=0;i<interior;i++)
 {
if(nodes[sides[interiorsides[i]].x].ptr>0)
{Pcurl[i][nodes[sides[interiorsides[i]].x].ptr-1]=0.5*(nodes[sides[interiorsides[i]].y].x-nodes[sides[interiorsides[i]].x].x);
Pcurl[i][nodes[sides[interiorsides[i]].x].ptr-1+num_free]=0.5*(nodes[sides[interiorsides[i]].y].y-nodes[sides[interiorsides[i]].x].y);}

if(nodes[sides[interiorsides[i]].y].ptr>0)
{
	Pcurl[i][nodes[sides[interiorsides[i]].y].ptr-1]=0.5*(nodes[sides[interiorsides[i]].y].x-nodes[sides[interiorsides[i]].x].x);
Pcurl[i][nodes[sides[interiorsides[i]].x].ptr-1+num_free]=0.5*(nodes[sides[interiorsides[i]].y].y-nodes[sides[interiorsides[i]].x].y);}

 }
















/* for(i=0;i<interior;i++)
{
	
	if(nodes[sides[interiorsides[i]].x].ptr>0)
C[i][nodes[sides[interiorsides[i]].x].ptr-1]=-1.0/length[i];
if(nodes[sides[interiorsides[i]].y].ptr>0)
C[i][nodes[sides[interiorsides[i]].y].ptr-1]=1.0/length[i];


}
 */

/* 


for(i=1;i<num_free;i++){printf("\n");
for(j=1;j<interior;j++)
{	printf("%f",B[i][j]);
}
} */




 for (i = 0;i<interior;i++)
    {
		
 //	printf("i=%d",i);system("pause");

		
        for (j=i;j<interior;j++)
        { 	
FinA[i][j]=A[i][j]-kk*kk*MM[i][j];}
 	


        for (j = interior; j <num_free+interior; j++)
          {FinA[i][j]=B[j-interior][i];}

	}


	
	
	
	 for (i = 0;i<num_free;i++)
    {
		

		
        for (j=0;j<i;j++)
        { 	
         L[i][j]=L[j][i];}
 	



	}


	
	
	 for (i = 0;i<interior;i++)
    {
		

		
        for (j=0;j<i;j++)
        { 	
         MM[i][j]=MM[j][i];}
 	



	}
	
	
	
	
		

	
	
 for (i = interior; i <interior+num_free; i++)
    {for (j =i; j <num_free+interior; j++)
          {FinA[i][j-interior]=0;}
	
	}

	
	//以下输出矩阵到文件
	
	
	

FILE *outPcurl;


    outPcurl = fopen("D:\\zsy\\Pcurl.txt", "w");
 if (!outPcurl)
    {
        perror("cannot open file");
        
    }
	  for (i = 0; i <interior; i++)
		  	  {for (j = 0; j <2*num_free; j++)
    fprintf(outPcurl, "%16.12f ", Pcurl[i][j]);fputc('\n',outPcurl);}
fclose(outPcurl);
//printf("\n  %d Run read with the command: \n read namen.n names.s namee.e\n ",k);system("pause");


FILE *out;


    out = fopen("D:\\zsy\\FinA.txt", "w");
 if (!out)
    {
        perror("cannot open file");
       
    }
	  for (i = 0; i <interior+num_free; i++)
    {
		
        for (j = 0; j <interior+num_free; j++)
        {
           if(j<i) fprintf(out, " %20.12f", FinA[j][i]) ;
		   else fprintf(out, " %20.12f", FinA[i][j]) ;
    } fputc('\n',out);
	
}
    fclose(out);
    
	
	FILE *outA;
 outA = fopen("D:\\zsy\\A.txt", "w");
 if (!outA)
    {
        perror("cannot open file");
       
    }
	  for (i = 0; i <interior; i++)
    {
		
        for (j = 0; j <interior; j++)
        {
           if(j<i) fprintf(outA, " %20.12f", A[j][i]) ;
		   else fprintf(outA, " %20.12f", A[i][j]) ;
    } fputc('\n',outA);
	
}
    fclose(outA); 

FILE *out1;


    out1 = fopen("D:\\zsy\\FinF.txt", "w");
 if (!out1)
    {
        perror("cannot open file");
        
    }
	  for (i = 0; i <interior; i++)
    {fprintf(out1, "%16.12f ", gg[i]);fputc('\n',out1);}

      for (i = interior; i <interior+num_free; i++)
     {fputc('0',out1);fputc('\n',out1);} 
 
fclose(out1);


 
FILE *outL;


	
	
    outL = fopen("D:\\zsy\\L.txt", "w");
 if (!outL)
    {
        perror("cannot open file");
        
    }
	
	
	  for (i = 0; i <num_free; i++)
	  {for (j = 0; j <num_free; j++)
	  fprintf(outL, "%20.12f  ",L[i][j]);fputc('\n',outL);}

  
  
  
fclose(outL);






FILE *outC;
    outC = fopen("D:\\zsy\\C.txt", "w");
 if (!outC)
    {
        perror("cannot open file");
        
    }
	
	
	  for (i = 0; i <interior; i++)
	  {for (j = 0; j <num_free; j++)
	  fprintf(outC, "%20.12f ",C[i][j]);fputc('\n',outC);}

  
fclose(outC);





 
FILE *outnn;


    outnn = fopen("D:\\zsy\\MM.txt", "w");
 if (!outnn)
    {
        perror("cannot open file");
        
    } 
	  for(i = 0; i <interior; i++)
	  {	  for(j = 0; j <interior; j++)
       
	fprintf(outnn, "%12.12f ",MM[i][j]);  fputc('\n',outnn);}
        
      
    
    fclose(outnn);










for(i=0;i<interior;i++)

    free(A[i]);

free(A);
A=NULL;








for(i=0;i<interior;i++)

    free(MM[i]);

free(MM);
MM=NULL;




for(i=0;i<num_free;i++)

    free(B[i]);

free(B);
B=NULL;









for(i=0;i<interior+num_free;i++)

    free(FinA[i]);

free(FinA);
FinA=NULL;


for(i=0;i<interior;i++)

    free(C[i]);

free(C);
C=NULL;






for(i=0;i<num_free;i++)

    free(L[i]);

free(L);
L=NULL;





free(FinF);
FinF=NULL;







return 0;
}


