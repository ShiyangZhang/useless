//编译环境：windows 10 Chinese;notepad++;MinGW.
//gcc -std=c99
//********************************************************************************************************************************************
/*本程序为论文 Preconditioners for the discretized time-harmonic Maxwell
equations in mixed form  的混合有限元实现

********************************************************************************************************************************************
本程序只能在 在 u（x,y）=1-y^2,1-x^2), p = (1 ? x^2 )(1 ? y^2 )和the unit square(-1=<x=<1,-1=<y=<1,)下获得正确的右端项和err1（u值在三角形单元中点的平均相对误差，我算得约为2.5% <n+m=1777>）;请注意更改程序中的k值
其他情况不能生成正确的右端项和err1，但是可以生成其他值，

********************************************************************************************************************************************
首先由EasyMesh(https://zsy.wodemo.com/file/387961 ; https://zsy.wodemo.com/cat/4506此压缩包包含所用网格)获得网格信息（EasyMesh中 ，对更细的网格，需修改EasyMesh的大小限制，如  #define MAX_NODES 8000）。生成以下三个文件
EasyMesh使用方法：编译EasyMesh后 在EasyMesh.exe所在文件夹按Shift+右键，打开  在此处打开命令窗口（W） 选项  输入 EasyMesh G3 +dxf  (此处G3代表G3.d,为EasyMesh输入文件)
C:\\Users\\zhang\\Desktop\\easymesh\\G3.n  
C:\\Users\\zhang\\Desktop\\easymesh\\G3.s  
C:\\Users\\zhang\\Desktop\\easymesh\\G3.e  
请在下面程序中搜索G3，更改您的相应文件名和文件位置
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

nodes_file=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\L2.n","r");

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


/*   printf("\n  %d Run read with the command: \n read namen.n names.s namee.e\n ",i);system("pause");
//用于错误测试
 
FILE *outn;


    outn = fopen("D:\\zsy\\Fn.txt", "w");
 if (!outn)
    {
        perror("cannot open file");
        
    } 
	  for (i = 0; i<num_free; i++)
    {for (num_freedummy = 0; num_freedummy<num_free; num_freedummy++){
       】
	fprintf(outn, "%f ",L[i][num_freedummy]);}
        
        fputc('\n',outn);
    }
    fclose(outn);
system("pause");
    */
/*  
FILE *outn;


    outn = fopen("D:\\zsy\\Fn.txt", "w");
 if (!outn)
    {
        perror("cannot open file");
        
    } 
	  
    {for (i = 0; i<num_free; i++){
       
	fprintf(outn, "%f ",F[i]); fputc('\n',outn);}
        
       
    }
    fclose(outn);
system("pause");
    
 */
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
sides_file=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\L2.s","r");
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
elements_file=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\L2.e","r");

       
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
		 A[i][j]+= sign(eptr[r])*sign(eptr[s])*localength[r]*localength[s]/(delta*delta);/* printf("%f",A[i][j]); */}
	    
	  }
}  

 for(r=0;r<3;r++)
  {for(s=0;s<3;s++)	 
  fff[r][s]=invM[0][r]*invM[0][s]+invM[1][r]*invM[1][s];}



  
  
  
 emm[0][0]=localength[0]*localength[0]*(fff[1][1]-fff[0][1]+fff[0][0])*(delta)/6;

 emm[0][1]=sign(eptr[0])*sign(eptr[1])*localength[0]*localength[1]*(fff[1][2]-fff[1][1]-2*fff[0][2]+fff[0][1])*(delta)/12;
 
 emm[0][2]=sign(eptr[0])*sign(eptr[2])*localength[0]*localength[2]*(fff[1][0]-2*fff[1][2]-fff[0][0]+fff[0][2])*(delta)/12; 

 
 emm[1][1]=localength[1]*localength[1]*(fff[2][2]-fff[1][2]+fff[1][1])*(delta)/6;
   
 emm[1][2]=sign(eptr[2])*sign(eptr[1])*localength[2]*localength[1]*(fff[2][0]-fff[2][2]-2*fff[1][0]+fff[1][2])*(delta)/12;

 emm[2][2]=localength[2]*localength[2]*(fff[0][0]-fff[0][2]+fff[2][2])*(delta)/6;
  
  
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
	       B[i][j]+= delta*sign(eptr[s])*localength[s]*localength[s]*(  invM[0][r]*(aaa[s]*invM[0][(1+s)%3]-aaa[(1+s)%3]*invM[0][s]+y0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s]))	+
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

  err1+=(sqrt(pow(1-y0*y0-umid[0],2)+pow(1-x0*x0-umid[1],2)))/sqrt(pow(1-y0*y0,2)+pow(1-x0*x0,2))/(num_elements); // printf(";%f,%f;%f,%f;%f\n",umid[0],umid[1],1-y0*y0,1-x0*x0,(fabs(1-y0*y0-umid[0])+fabs(1-x0*x0-umid[1]))/(1-y0*y0+1-x0*x0));system("pause");

  //err1+=(fabs(1-y0*y0-umid[0])+fabs(1-x0*x0-umid[1]))/(1-y0*y0+1-x0*x0)/(num_elements); // printf(";%f,%f;%f,%f;%f\n",umid[0],umid[1],1-y0*y0,1-x0*x0,(fabs(1-y0*y0-umid[0])+fabs(1-x0*x0-umid[1]))/(1-y0*y0+1-x0*x0));system("pause");

//三角基元对刚度矩阵循环完毕,但是同时要对load矩阵循环,所以不加大括号 
 

  
//以下为算法第二部分:load矩阵的初始化和计算.*******************************************************************************************************************************************************************************  	  



//以下为algorithm 7.5的第一部分************************************



 for(s=0;s<3;s++)

      {if(localength[s]>=0)
		  {i=sides[abs(eptr[s])-1].ptr-1;
	       
	          gg[i]+= delta*sign(eptr[s])*localength[s]*(( 2-(kk*kk+2*x0)*(1-y0*y0) )  *(aaa[s]*invM[0][(1+s)%3]-aaa[(1+s)%3]*invM[0][s]+y0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s]))	+
               ( 2-(kk*kk+2*y0)*(1-x0*x0) )             * (aaa[s]*invM[1][(1+s)%3]-aaa[(1+s)%3]*invM[1][s]-x0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s])	)   );
           }
     }
  
  
  
/* 	printf("k=%d   %f   gg=%f",k,invM[0][(1+2)%3],gg[k]);system("pause");
 */



}    //对三角基元循环完毕.


//以下计算C




/* for(i=0;i<interior;i++)
{
	
	if(nodes[sides[interiorsides[i]].x].ptr>0)
C[i][nodes[sides[interiorsides[i]].x].ptr-1]=-1.0/length[i];
if(nodes[sides[interiorsides[i]].y].ptr>0)
C[i][nodes[sides[interiorsides[i]].y].ptr-1]=1.0/length[i];


}
 */


for(i=0;i<interior;i++)
{
	
	if(nodes[sides[interiorsides[i]].x].ptr>0)
C[i][nodes[sides[interiorsides[i]].x].ptr-1]=1.0;
if(nodes[sides[interiorsides[i]].y].ptr>0)
C[i][nodes[sides[interiorsides[i]].y].ptr-1]=-1.0;


}


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
//printf("\n  Run read with the command: \n read namen.n names.s namee.e\n ");

 
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
 if (!outL)
    {
        perror("cannot open file");
        
    }
	
	
	  for (i = 0; i <interior; i++)
	  {for (j = 0; j <num_free; j++)
	  fprintf(outC, "%20.12f ",C[i][j]);fputc('\n',outC);}

  
  
  
fclose(outC);

/* 
FILE *outAA;


    outAA = fopen("D:\\zsy\\AA.txt", "w");
 if (!outAA)
    {
        perror("cannot open file");
       
    }
	  for (i = 0; i <interior; i++)
    {
		
        for (j = 0; j <interior; j++)
        {
           if(j<i) fprintf(outAA, " %20.12f", FinA[j][i]) ;
		   else fprintf(outAA, " %20.12f", FinA[i][j]) ;
    } fputc('\n',outAA);
}
    fclose(outAA);
    
 
FILE *outFF;


    outFF = fopen("D:\\zsy\\FF.txt", "w");
 if (!outFF)
    {
        perror("cannot open file");
        
    }
	  for (i = 0; i <interior; i++)
    {fprintf(outFF, "%.12f ", gg[i]);

     fputc('\n',outFF);} 
 
fclose(outFF);





 */


  


//  printf("\n  %d Run read with the command: \n read namen.n names.s namee.e\n ",k);system("pause");

 
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




for(i=0;i<interior;i++)

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







//此误差仅在对数据处理并写入u.txt后第二次运行有效
printf("误差为err1=%12.12f（处理数据后再使用才有效）",err1);system("pause");










return 0;
}













//以下命令涉及到以下文章的实现
//On Nonsingular Saddle-Point Systems with a Maximally Rank Deficient Leading Block
//Modified block preconditioners for the discretized time-harmonic
//Maxwell equations in mixed form
//以及正在写还未发表的文章


/*其他matlab命令part1
 
 Z=null(B);

 %part 2 验证[2006]文章中的关系，需要都接近于0
 max(max(abs(FinA*x-FinF)))
max(abs(B*u))
max(abs(B*Z(:,1)))

 max(max(abs(MM*C-B')))

max(max(abs(B*C-L)))

max(max(abs(A*C)))
max(abs(B*u))



R=B;for i=(m+1):n R(i,1:n)=R(2,:); end 
R(n-(m-1):n)=B;for i=1:n-(m-2) R(i,1:n)=R(n,:); end 
R=B'*inv(L)*B;rank(A+R)
W=rand(m);rank(W)
R=B'*inv(W)*B;rank(A+R)

inv(A+R)*(eye(n)-B'*inv(L)*C');
inv(A+R)*A*inv(A+R)*A;%是对称的 
max(max(abs(ans-ans')))

rank(inv(A+R)*A+C*inv(L)*B)

eig(inv(A+R)*A+C*inv(L)*B)
hist(eig(inv(A+R)*A+C*inv(L)*B),1000)

%part4 [2014]文章的验证。
R=1e-6*eye(n);


 for i=1:1000
W=rand(m);
i
R=B'*inv(L)*B;
if rank(A+R)<n
break;end
end





R=1e-6*eye(n);
R=MM;


invP1=[inv(A+R)*(eye(n)-B'*inv(L)*C') C*inv(L);inv(L)*C' zeros(m)];

[x1,flag,relres,iter,resvec] =pcg(invP1*FinA,invP1*FinF,1e-10,[])

Z=null(B);V=Z*inv(Z'*A*Z)*Z';VA=V*A;max(max(abs(VA-eye(n)+C*inv(L)*B)))

invP2=[V C*inv(L);inv(L)*C' zeros(m)];
 hist(eig(invP2*FinA),100)

max(max(abs(invP2-invP1)));


max(max(abs(x2-x1)));

V=Z*Z';invP2=[V C*inv(L);inv(L)*C' zeros(m)];

[x2,flag,relres,iter,resvec] =pcg(invP2*FinA,invP2*FinF,1e-10,[])


f=FinF(1:n);
max(abs(A*V*f+B'*x((n+1):(m+n))-f));
(u-V*f)




FINDB0=zeros(m,2);
 for i=1:m
FINDB0(i,:)=[find(B(i,:),1) find(B(i,:),1,'last')];
end

sum(B~=0,2)  %统计每行非零的个数


dd=Z(:,2);
for i=1:100
ddp=inv(B*B')*B*(FinF(1:n)-A*dd);
dd=VA*dd+C*inv(L)*inv(L)*C'*ddp;
end









f=FinF(1:n);
u0=(A)\(f-B'*inv(L)*C'*f);max(max((A)*u0-(f-B'*inv(L)*C'*f)))
unew=(eye(n)-C*inv(L)*B)*u0;max(abs(x1-unew)) %+C*inv(L)*g


%part 5

%从此处开始
R=MM;

invP1=[inv(A+R)*(eye(n)-B'*inv(L)*C') C*inv(L);inv(L)*C' zeros(m)];

[x1,flag,relres,iter,resvec] =pcg(invP1*FinA,invP1*FinF,[],10)



f=FinF(1:n);

ff=f-B'*inv(L)*C'*f  ;
                %右端项，计算特解用
				
				
	
	
for i=1:n
lengthc(i)=length(nonzeros(C(i,:)));
end
hist(lengthc,1000)





spec=zeros(size(C));
for i=1:n
lengthc=length(nonzeros(C(i,:))); 
if lengthc==1
spec(i,find((C(i,:))))=1;
else spec(i,find(C(i,:),1))=1;
end
end
rank(spec)


% for j=1:m
% cc=find(spec(:,j),1);
% mark(j)=cc;
% end


for j=1:m
cc=find(spec(:,j));
ccc=randperm(length(cc),1);
mark(j)=cc(ccc);
end
%随机产生需要删除量

mark=sort(mark);
 Ar=A;
Ar(:,mark)=[];
Ar(mark,:)=[];
rank(Ar)%验证剩余线性无关，重要
cond(Ar)










for t=1:1000000000000000000000000000000000000000000000000000000

for j=1:m
cc=find(spec(:,j));
ccc=randperm(length(cc),1);
mark(j)=cc(ccc);
end
%随机产生需要删除量

mark=sort(mark);
 Ar=A;
Ar(:,mark)=[];
Ar(mark,:)=[];
o=rank(Ar)
if o<421
cccccc=o
break;
end

end









%选取特定sort


spec=zeros(size(C));
for i=1:n
lengthc=length(nonzeros(C(i,:))); 
if lengthc==1
spec(i,find((C(i,:))))=1;
else spec(i,find(C(i,:),1))=-1;
end
end
rank(spec)

for j=1:m
cc=find(spec(:,j)==-1);
if isempty(cc)
mark(j)=find(spec(:,j),1,'last');
else mark(j)=find(spec(:,j)==-1,1,'last');
end
end
%随机产生需要删除量

mark=sort(mark);
 Ar=A;
Ar(:,mark)=[];
Ar(mark,:)=[];
rank(Ar)%验证剩余线性无关，重要
cond(Ar)


%假设spec每列都有-1



%选取特定sort


spec=zeros(size(C));
for i=1:n
lengthc=length(nonzeros(C(i,:))); 
if lengthc==1
spec(i,find((C(i,:))))=1;
else spec(i,find(C(i,:),1))=-1;
end
end
rank(spec)

for j=1:m
cc=find(spec(:,j)==1);
if isempty(cc)
mark(j)=find(spec(:,j),1,'last');
else mark(j)=find(spec(:,j)==1,1,'last');
end
end
%随机产生需要删除量

mark=sort(mark);
 Ar=A;
Ar(:,mark)=[];
Ar(mark,:)=[];
rank(Ar)%验证剩余线性无关，重要
cond(Ar)











%全部统计，找最小的cond
count=sum(spec~=0);prod(count)
v=zeros(1,m);
for j=1:m
for i=1:count(j)
mark(j)=find(spec(:,j),i);















ff=f-B'*inv(L)*C'*f  ;
ff0=ff;
ff0(mark,:)=[];
u00=Ar\ff0;max(max((Ar)*u00-ff0));

for i=1:m
u00(mark(i)+1:numel(u00)+1) = u00(mark(i):end);
u00(mark(i)) = 0;
end%补齐特解，补为 0；
u0new=(eye(n)-C*inv(L)*B)*u00;max(abs(x1(1:n)-u0new)) %+C*inv(L)*g  
%比较此算法最终结果


rank(spec)







%
a=0.0010:0.0001:0.0045;length(a)
 b=a;
 
t=0;
for i=1:length(a);
b(i)=cond(A+a(i)*C*C');
end
20+t


plot(a,b)



a=22e-5:0.00001:27e-5;length(a)
 b=a;
 
t=0;
for i=1:length(a);

b(i)=cond(A+a(i)*CC);
end
20+t


plot(a,b)


%G4 只有无穷范数和1范数较好，他两相等








 hist(eig(A+0.01*C*C'),1000)

 cond(A+0.0001*C*C')
 ich=ichol(sparse(A+0.0001*C*C')，struct('type','ict'));
[x1,flag,relres,iter,resvec] =pcg(A+0.01*C*C',ff,1e-10,100,ich,ich')


 ich=ichol(sparse(L));
[x1,flag,relres,iter,resvec] =pcg(L,B(:,1),1e-8,100,ich,ich')


 ich=ichol(sparse(Ar));
[x1,flag,relres,iter,resvec] =pcg(Ar,ff0,1e-10,100,ich,ich')



 ich=ichol(sparse(A+MM));
[x1,flag,relres,iter,resvec] =pcg(A+MM,A(:,1),1e-10,100,ich,ich')

CC=C*C';
a=0.5*sum(sum(A))/sum(sum(CC))
cond(A+a*C*C')




A(:,1)'*A*A(:,1)/(A(:,1)'*(C*C')*A(:,1))

norm(A)/norm(C*C')






sada=A(:,1); sada*A*sada'/(sada*C'*C*sada);










ggu=B'*inv(L)*C'*C*C'*C*inv(L)*B;
 max(max(abs(C'*(ggu-C*C')*C)))
 
 
 
 
 
 %G3实验数据全记录
 
 sada=eig(A);
 
 atei=sada(m+1)
 atea=sada(n)
 
 atea/atei
 
 
sada=eig(C'*C);
 
 ctei=sada(1)
 ctea=sada(m)
  ctea/ctei
 
 alphaa=atea/ctea
 alphai=atei/ctei
 midalp= (alphaa+ alphai)/2
 
 norm(A)/norm(C'*C)
 
norm(A,1)/norm(C'*C,1)
norm(A,1)/norm(C*C',1)%偏小
sum(sum(A))/sum(sum(C'*C))%too bad

sum(sum(abs(A)))/sum(sum(abs(C'*C)))%too bad,either

sum(sum(A))/sum(sum(C*C'))%偏大
0.5*sum(sum(A))/sum(sum(C*C'))- midalp % 有时很近 good but not good sometiomes
(norm(A,1)/norm(C'*C,1)-alphaa )% 有时很近（G2，G4）


(norm(A,1)/norm(C'*C,1)+0.5*sum(sum(A))/sum(sum(C*C')))/2

sum(sum(abs(A)))/sum(sum(abs(C*C')))% bad,but ok.

norm(A,'fro')/norm(C'*C,'fro')
norm(A,'fro')/norm(C*C','fro') %相同的根据f范数性质；

 
 hist(eig(inv(A+1*MM)*(A-2500*MM+(2500)*B'*inv(L)*B)),1000)
 prem(1:n,1:n)=(A+1*MM)
 
 
 
 k=0;
 FinAm=FinA;
 FinAm(1:n,1:n)=A-k^2*MM;
 invL=inv(L);
 p=1;
Pm=[inv(A+p*MM)*(eye(n)-B'*invL*C') C*invL;invL*C' k^2*inv(L)];
[x1,flag,relres,iter,resvec] =pcg(FinAm,FinF+sort(FinF),1e-10,100,@(input_args)Pm*input_args)
 
 
 
 
 
 
 
  hist(eig(inv(A+1*MM)*(A-2500*MM+(2500)*B'*inv(L)*B)+C*inv(L)*B),-1:0.01:1.1)
  hist(eig(A+B'*inv(L)*B-MM),1000)
  
  
  
  sada=A-25*MM+(25)*B'*inv(L)*B;
  sadaa=C*inv(L)*B;
  cond(inv(A+20*MM)*(sada)+sadaa,1)



*/




/*其他matlab命令part2
rhs=ones(m+n,1);


k=1.5;FinA(1:n,1:n)=A-k^2*MM;
invP=[inv(A+rho*MM-k^2*MM)*(eye(n)-B'*inv(L)*C') C*inv(L)     ;
inv(L)*C'  k^2*inv(L)];
rho=k^2+1;


rhs=ones(m+n,1);

chaP=blkdiag(A+rho*MM-k^2*MM, eye(m));
chaA=blkdiag(A+rho*B'*inv(L)*B-k^2*MM ,eye(m));



%测试
k=0;rho=k^2+1;
invP=[inv(A+rho*MM-k^2*MM)*(eye(n)-B'*inv(L)*C') C*inv(L)     ;
inv(L)*C'  k^2*inv(L)];



%cg jianjieyong
[x1,flag,relres,iter,resvec] =pcg(chaA,chaP*invP*rhs,1e-6,100,chaP)
norm(FinA*x1-rhs)/norm(rhs)


norm(chaP*invP*rhs)\norm(chaA*x1-chaP*invP*rhs)


max(max(abs(chaA*x1-chaP*invP*rhs)))



invP*FinA-inv(chaP)*invP
max(max(abs(invP*FinA-inv(chaP)*chaA)))


以下是k=0.5 的结果
norm(FinA*x1-rhs)/norm(rhs)

7下-- 1.1558e-09    //前面迭代到很小很快
10下-- 4.5818e-10
11下 --   4.5818e-10



\\ bolck
R=MM;
P1=blkdiag(A+(rho-k^2)*R ,1/rho*L);

[x1,flag,relres,iter,resvec] =minres(FinA,rhs,1e-13,11,P1)

7--   1.0937e-07
10-- 8.2819e-11
11--   2.5932e-11

\\cg  zhijieyong

[x1,flag,relres,iter,resvec] =pcg(FinA,rhs,1e-100,10,@(input_args)invP*input_args)
7--  1.0606e-09
10--- 2.1482e-13
11--   2.8313e-14



\\minres  zhijieyong
[x1,flag,relres,iter,resvec] =minres(FinA,rhs,1e-10,100,@(input_args)invP*input_args)
7--   1.0349e-09    
10--   2.1170e-13
11---  2.7299e-14



以下是k=0
norm(FinA*x1-rhs)/norm(rhs)
6---          4.9222e-09
7----        4.2658e-10
10---      3.9499e-10
11---       3.9499e-10





\\ bolck

R=MM;
P1=blkdiag(A+(rho-k^2)*R ,1/rho*L);
[x1,flag,relres,iter,resvec] =minres(FinA,rhs,1e-10,100,P1)
%7--     5.6018e-08
%10--   6.6176e-11
%11--      3.7469e-12


\\cg  zhijieyong

[x1,flag,relres,iter,resvec] =pcg(FinA,rhs,1e-10,100,@(input_args)invP*input_args)
%7--           1.6071e-10 
%10--        3.0086e-14
%11 ---      2.5163e-14       



\\tri
P=[ A+(rho-k^2)*MM  (1+rho)*B';zeros(m,n)  -L];
[x1,flag,relres,iter,resvec] =gmres(FinA,rhs,20,1e-10,100,P)

tri=eig(inv(P)*FinA);
dia=eig(invP*FinA);
max(max(abs(tri-dia)))


%invP  gmres
[x1,flag,relres,iter,resvec] =gmres(FinA,rhs,20,1e-10,19,@(input_args)invP*input_args)


k=1.6;FinA(1:n,1:n)=A-k^2*MM;rho=k^2+1;
invP=[inv(A+rho*MM-k^2*MM)*(eye(n)-B'*inv(L)*C') C*inv(L)     ;
inv(L)*C'  k^2*inv(L)];
z=min(eig(A+rho*B'*inv(L)*B-k^2*MM))
y=min(eig(invP*FinA))

%zhijieyong pcg
[x1,flag,relres,iter,resvec] =pcg(FinA,rhs,1e-10,100,@(input_args)invP*input_args)


[x1,flag,relres,iter,resvec] =minres(FinA,rhs,1e-10,100,@(input_args)invP*input_args)

%bolck
R=MM;
P1=blkdiag(A+(rho-k^2)*R ,1/rho*L);
[x1,flag,relres,iter,resvec] =minres(FinA,rhs,1e-10,100,P1)

G1   (1e-10)

k                 0                   1.0     1.2       1.3      1.5            1.55             1.6        1.7           2            4                      
直接用pcg         8                   10       11        12       14             15             error       error        error        error                                  
间接用pcg                                                                                           error       error                           
z               0.0437                0.0437   0.0437   0.0372       0.0104        0.0030      -0.0047      -0.0209    -0.0754      -0.6670                                                                          
y              0.7113                0.4227     0.2956    0.2235     0.0618        0.0178      -0.0276    -0.1229      0.1569        0.1735                                                       
直接用minres      8                     10      11       12      14             15          15          14          15           31                                   
    
bolck-minres     10                    13        13      14       16             18          18         18            19          36                                          
 

tri gmres(20)    1,7                                                                                           2,14                 
invp gmres(20)   1,7                                                                                           2,14             




\begin{center}
\begin{table}[h]\footnotesize
\begin{tabular}{lccccccccc}  % {lccc} 表示各列元素对齐方式，left-l,right-r,center-c
\hline
k&     0      &    1.0        &    1.2      &     1.3  &       1.5  & 1.55           &        1.6     &      1.7       &      2                          \\ \hline  % \hline 在此行下面画一横线
pcg  &  8    &    10       &    11   &      12&  14&         15   &       error          & error     &    error           \\
minres  &     8        &              10      & 11   &     12  &     14     &         15    &       15     &      14      &     15              \\
b-minres &   10         &           13     &   13     & 14     &  16      &       18       &   18       &  18       &     19               \\
z &      0.0437   &             0.0437  & 0.0437&   0.0372  &     0.0104    &    0.0030   &   -0.0047   &   -0.0209    &-0.0754                       \\
y &      0.7113      &          0.4227   &  0.2956  &  0.2235  &   0.0618   &     0.0178    &  -0.0276 &   -0.1229  &    0.1569        \\
 \hline         

\end{tabular}
\caption{G1}
\end{table}
\end{center}










G2   (1e-10)        ones
k                 0                    1/8            1.0         1.2        1.5          1.55       1.6              1.7           2             4                     
直接用pcg         8                    8              10          11        14              15       error             error          error        error                                                                                                   
间接用pcg                                                                                                                                                  
z               0.0173               0.0173         0.0173       0.0173     0.0048        0.0014       -0.0021       -0.0096      -0.0348       -0.3090                                                                                                                                                       
y               0.7116               0.7071         0.4232      0.2963      0.0627          0.0187     -0.0267       -0.1218       0.1578      0.1796                                                                                                                                        
直接用minres     8                    8              10        11            14              15           15           14          15           31                                                                                   
                                                        
bolck-minres     11                  11              13         14           17            18             18           18         19           36                                                                                         
 



L2   (1e-10)  右端项  ones
k                 0              1.0         1.2       1.3          1.5           1.7        2             4                      
直接用pcg        8               10        12           error     error          error       error         error                                                                                                                            
间接用pcg                                                                                                                        
z              0.0045           0.0042     1.3048e-04     -0.0052    -0.0197       -0.0385    -0.0741      -0.6404                                                                                                                                                                                                        
y               0.5913           0.1827     0.0029      -0.0993    0.2839         0.1429    -0.1016      0.1875                                                                                                                                                                     
直接用minres   8                   10         15            13       14              15        16           29                                                                                                            
   
bolck-minres   10                  13           15           14      14              17       18            34                                                                                                                

tri gmres(20)                                                                                                                         
invp gmres(20)  

\begin{center}
\begin{table}[h]\footnotesize
\begin{tabular}{lcccccccccc}  % {lccc} 表示各列元素对齐方式，left-l,right-r,center-c
\hline
k&     0   &   1.0       &  1.2   &    1.3    &      1.5       &    1.7     &   2       &      4                               \\ \hline  % \hline 在此行下面画一横线
pcg  &   8& 10&         12         &   error  &    error    &       error  &      error  &        error              \\
minres  &  8  &                 10   &      15       &     13   &    14     &         15     &   16    &       29       \\
b-minres &   10   &    13       &     15    &        14  &     14     &          17   &     18      &       34                       \\
z &   0.0045    &       0.0042  &   1.3048e-04  &   -0.0052 &   -0.0197  &     -0.0385   & -0.0741     & -0.6404      \\
y &   0.5913      &     0.1827 &    0.0029    &  -0.0993  &  0.2839   &      0.1429&    -0.1016     & 0.1875        \\
 \hline         

\end{tabular}
\caption{L2}
\end{table}
\end{center}


G3   (1e-10)  右端项  ones
k                 0        1.0          1.2         1.3       1.55      1.5            1.6         1.7           2                
直接用pcg        8         10           11          12           15       14          error         error       error                                                                                                                                                              
间接用pcg                                                                                                                                                   
z               0.0046    0.0046     0.0046       0.0044     3.7141e-04    0.0012     -5.3503e-04   -0.0024      -0.0089                                                                                                                                                                                                                                                                                                  
y               0.7116    0.4232     0.2963       0.2242        0.0187    0.0627      -0.0267       -0.1219      0.1575                                                                                                                                                                                                                                        
直接用minres    8         10          11          12           15          14           15         15            15                                                                                                                                                 
   
bolck-minres   10        13            13        14             18          17          18          18           19                                                                                                                                                                            

tri gmres(20)                                                                                                                         
invp gmres(20)  








\begin{center}
\begin{table}[h]\footnotesize
\begin{tabular}{lcccccccccc}  % {lccc} 表示各列元素对齐方式，left-l,right-r,center-c
\hline
k&     0      &    1.0        &    1.2      &     1.3  &       1.5  & 1.55           &        1.6     &      1.7       &      2                          \\ \hline  % \hline 在此行下面画一横线
pcg  &  8    &    10       &    11   &      12&15&      &         14    &       error          & error     &    error           \\
minres  &   8    &     10       &   11        &  12     &      15        &  14     &      15     &    15   &         15          \\
b-minres &   10     &      13     &          13     &      14      &          18    &         17      &       18     &        18     &         19          \\
z &     0.0046   &  0.0046    &  0.0046    &    0.0044     & 3.7141e-04  &   0.0012   &   -5.3503e-04  &  -0.0024   &    -0.0089                       \\
y &     0.7116    &   0.4232    &    0.2963   &       0.2242     &      0.0187  &     0.0627       &  -0.0267   &       -0.1219    &     0.1575       \\
 \hline         

\end{tabular}
\caption{G3}
\end{table}
\end{center}









G4   (1e-10)  右端项  ones
k                0              1.0           1.2            1.3           1.5            1.55             1.6               1.7              2                
直接用pcg        8              10           11              12            13             14              error           error           error                                                                                                                                                                     
间接用pcg                                                                                                                                                              
z              0.0012          0.0012      0.0012          0.0011      3.1390e-04     9.3715e-05          -1.3383e-04     -6.1075e-04     -0.0022                                                                                                                                                                                                                                                                                                         
y              0.7116          0.4232       0.2963         0.2242      0.0627          0.0187              -0.0267         -0.1219       0.1575                                                                                                                                                                                                                                                  
直接用minres    8               9            11             12           13            14                   14              14              15                                                                                                           
   
bolck-minres    10             12             13             15          16           18                    18              17              18                                                                                                                                                       

tri gmres(20)                                                                                                                                            
invp gmres(20)  












L4  (1e-10)  右端项  ones
k                             0               1.0              1.2               1.3              1.5           1.7            2               4       
直接用pcg                   8                10                12             error              error           error        error           error                                                                                                                           
间接用pcg                                                                                                                                                                                                
z                        4.0389e-05        4.0389e-05          2.2624e-05    -1.7483e-04       -7.3466e-04        -0.0014      -0.0028       -0.0224                                                                                                                                                                                                                                                                                                                                                                                                                          
y                         0.5959         0.1917                 0.0139       -0.0871             0.2833           0.1421        -0.1027      0.1804             
直接用minres               8               10                   14              13               13                14             15          28                                                                                                                                                                                                                                                 
   
bolck-minres              10              14                    15            15                 15                 16            18          34                                                                                                                                                                                         

tri gmres(20)                                                                                                                         
invp gmres(20)  




\begin{center}
\begin{table}[h]\footnotesize
\begin{tabular}{lcccccccccc}  % {lccc} 表示各列元素对齐方式，left-l,right-r,center-c
\hline
k&     0   &   1.0       &  1.2   &    1.3    &      1.5       &    1.7     &   2       &      4                               \\ \hline  % \hline 在此行下面画一横线
pcg  &               8   &               10        &          12     &          error  &              error   &          error      &    error      &       error               \\
minres  &         8          &     10        &           14          &    13      &         13       &         14      &       15    &      28      \\
b-minres &        10     &          14      &               15        &     15      &            15        &          16     &        18       &    34          \\
z &   4.03e-05   &     4.03e-05    &      2.26e-05 &   -1.74e-04     &  -7.34e-04      &  -0.0014     & -0.0028 &      -0.0224              \\
y &   0.5959       &   0.1917            &      0.0139      &  -0.0871    &          0.2833       &     0.1421    &     -0.1027      & 0.1804         \\
 \hline         

\end{tabular}
\caption{L4}
\end{table}
\end{center}










\begin{center}
\begin{table}[h]\footnotesize
\begin{tabular}{lcccccccccc}  % {lccc} 表示各列元素对齐方式，left-l,right-r,center-c
\hline
k&       0        &      1.0      &     1.2     &       1.3     &      1.5     &       1.55         &    1.6          &     1.7        &      2                             \\ \hline  % \hline 在此行下面画一横线
pcg  &  8    &  10     &       11        &       12        &     13     &         14       &        error    &        error       &     error             \\
minres  &   8      &          9     &        11        &      12      &      13         &    14        &            14          &     14    &           15   \\
b-minres &    10     &         12           &   13      &        15       &    16       &     18       &              18       &        17        &       18        \\
z &   0.0012     &     0.0012   &   0.0012      &    0.0011   &   3.1390e-04  &   9.3715e-05     &     -1.3383e-04   &  -6.1075e-04   &  -0.0022     \\
y &    0.7116     &     0.4232  &     0.2963        & 0.2242    &  0.0627   &       0.0187     &         -0.0267       &  -0.1219     &  0.1575       \\
 \hline         

\end{tabular}
\caption{G4}
\end{table}
\end{center}










































L3  (1e-10)  右端项  ones
k                         0              1.0                   1.2            1.3           1.5              1.7           2            4           
直接用pcg                8              10                    12            error          error             error         error        error                                                                                                                                                                                                                                                                        
间接用pcg                                                                                                                                                                                                                                      
z                     2.7366e-04       2.7344e-04         9.1260e-05     -8.0887e-04     -0.0033           -0.0065        -0.0123      -0.0977                                                                                                                                                                                                                                                                                                                                                            
y                     0.5952             0.1904            0.0123         -0.0889         0.2834            0.1422         -0.1025       0.1809                                                                                                                                                                                                                                                                                                                      
直接用minres            8               10                  14             13              13               14              16           28                                                                                                                                                                                                                                                                                
                 
bolck-minres           10               12                  15            15               15               16              17          34                                                                                                                                                                                                                                                         
                        
tri gmres(20)                                                                                                                                       
invp gmres(20)  




\begin{center}
\begin{table}[h]\footnotesize
\begin{tabular}{lcccccccccc}  % {lccc} 表示各列元素对齐方式，left-l,right-r,center-c
\hline
k&     0   &   1.0       &  1.2   &    1.3    &      1.5       &    1.7     &   2       &      4                               \\ \hline  % \hline 在此行下面画一横线
pcg  &          8  &            10      &              12           & error    &      error     &        error  &       error  &      error            \\
minres  &       8     &          10     &             14       &      13    &          13       &        14        &      16        &   28   \\
b-minres &    10     &        13     &          16         &       14    &     16        &          17        &     18  &         34                 \\
z &    2.73e-04  &      2.73e-04    &      9.12e-05   &   -8.0887e-04   &   -0.0033   &         -0.0065    &     -0.0123&       -0.0977           \\
y &   0.5952    &          0.1904    &         0.0123   &       -0.0889   &       0.2834 &            0.1422    &      -0.1025 &        0.1809          \\
 \hline         

\end{tabular}
\caption{L3}
\end{table}
\end{center}


















L2  (1e-10)  右端项  ones
k                         0              1.0          1.2            1.3           1.5                 1.7         2           4          
直接用pcg                8               10          12              error         error             error        error       error                                                                                                                                                                                                                                                                                                                 
间接用pcg                                                                                                                                                                                                                                                                                
z                     8.9042e-04     8.8097e-04     1.7344e-04      -0.0022     -0.0091            -0.0181        -0.0349   -0.3001                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
y                       0.5940        0.1879            0.0093      -0.0922      0.2841             0.1431        -0.1014    0.1823                                                                                                                                                                                                                                                                                                                                                                                           
直接用minres             8            10                15             14        14                  15           16          28                                                                                                                                                                                                                                                                                                                                     
                                                                         
bolck-minres            10            13              16               14        16                 17            18          34                                                                                                                                                                                                                                                                           

tri gmres(20)                                                                                                                         
invp gmres(20)  


\begin{center}
\begin{table}[h]\footnotesize
\begin{tabular}{lcccccccccc}  % {lccc} 表示各列元素对齐方式，left-l,right-r,center-c
\hline
k&     0   &   1.0       &  1.2   &    1.3    &      1.5       &    1.7     &   2       &      4                               \\ \hline  % \hline 在此行下面画一横线
pcg  &    8          &     10  &        12      &        error  &       error       &      error      &  error     &  error           \\
minres  &  8   &          10            &     15           &   14     &    14    &               15   &         16     &      28         \\
b-minres &    10     &        13     &          16         &       14    &     16        &          17        &     18  &         34                 \\
z &   8.90e-04  &   8.80e-04    & 1.73e-04&      -0.0022    & -0.0091 &           -0.0181  &      -0.0349 &  -0.3001        \\
y &    0.5940 &       0.1879  &          0.0093 &     -0.0922   &   0.2841  &           0.1431      &  -0.1014&    0.1823          \\
 \hline         

\end{tabular}
\caption{L2}
\end{table}
\end{center}











MM=sparse(MM);L=sparse(L);A=sparse(A);B=sparse(B);C=sparse(C);invL=inv(L);invL=sparse(invL);


k=1.2;rho=k^2+1;

invP=[inv(A+rho*MM-k^2*MM)*(eye(n)-B'*invL*C') C*invL     ;
invL*C'  k^2*invL];
invP=sparse(invP);



lingshi=full(A+rho*B'*invL*B-k^2*MM);
z=min(eig(lingshi))
clear lingshi;


FinA(1:n,1:n)=A-k^2*MM;
FinA=sparse(FinA);

lingshi2=full(invP*FinA);
y=min(eig(lingshi2))
clear lingshi2
%%y=min(eig(lingshi2))



%zhijieyong pcg
[x1,flag,relres,iter1,resvec] =pcg(FinA,rhs,1e-10,100,@(input_args)invP*input_args)


[x1,flag,relres,iter2,resvec] =minres(FinA,rhs,1e-10,100,@(input_args)invP*input_args)

%bolck
P1=blkdiag(A+(rho-k^2)*MM ,1/rho*L);
P1=sparse(P1);
[x1,flag,relres,iter3,resvec] =minres(FinA,rhs,1e-10,100,P1)
z
y



FinA=A+rho*B'*invL*B-k^2*MM;
P=A+rho*MM-k^2*MM;   rhs=(eye(n)-B'*invL*C')*ones(n,1)+B'*invL*ones(m,1);

[x1,flag,relres,iter1,resvec] =pcg(FinA,rhs,1e-10,100,P)

*/













