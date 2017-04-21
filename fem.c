//编译环境：windows 10 Chinese;notepad++;MinGW.
//gcc -std=c99
//********************************************************************************************************************************************
/*本程序为论文 Preconditioners for the discretized time-harmonic Maxwell
equations in mixed form  的混合有限元实现

********************************************************************************************************************************************
本程序只能在 在 u（x,y）=1-y^2,1-x^2), p = (1 ? x^2 )(1 ? y^2 )和the unit square(-1=<x=<1,-1=<y=<1,)下获得正确的右端项.
请忽略err1;请注意更改程序中的k值
编译完成后，请重新在所有文件所在的同一个目录下打开命令行 程序。
********************************************************************************************************************************************
首先由EasyMesh(https://zsy.wodemo.com/file/387961 ; https://zsy.wodemo.com/cat/4506此压缩包包含所用网格)获得网格信息（EasyMesh中 ，对更细的网格，需修改EasyMesh的大小限制，如  #define MAX_NODES 8000）。生成以下三个文件
EasyMesh使用方法：编译EasyMesh后 在EasyMesh.exe所在文件夹按Shift+右键，打开  在此处打开命令窗口（W） 选项  输入 EasyMesh G6 +dxf  (此处G6代表G6.d,为EasyMesh输入文件)
C:\\Users\\zhang\\Desktop\\easymesh\\G6.n  
C:\\Users\\zhang\\Desktop\\easymesh\\G6.s  
C:\\Users\\zhang\\Desktop\\easymesh\\G6.e  
请在下面程序中搜索G6，更改您的相应文件名和文件位置
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

任何意见，改进和询问请发邮件至  hydzhang@hotmail.com ，
2015.10.

修改：2015.12  增加err1的计算
修改：2016.3   
修改：2017.2 --reconsidered after major revision 

*/
#include<math.h>
#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h> 
#include<errno.h>

#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef PI
#define PI    3.14159265359
#endif
#define MAX_NODES 8800
#define FREE 0
#define OK 1  
#define ERROR 0  
#define TRUE 1  
#define FALSE 0  

/* 
#define interior 15526
#define num_free 5087
 */
/*=========================================================================*/


int sign(int a){if (a>0)return 1;else return -1;}

double doubleelem=0.0;
#define MAXSIZE 42500*4  
  //MAXSIZE需要等于或稍微大于边的个数
typedef int Status;  
typedef double ElemType;  

typedef struct{//三元组结构  
    int i,j;//非零元素行下标和列下标  
    ElemType e;//非零元素值  
}Triple;  

typedef struct{  
    Triple *data;//三元组表，data[0]可用  
	
    int tu;//矩阵的行数、列数和非零元素个数  
}TSMatrix;  
  
  TSMatrix NewMatrix(int m,int n);    TSMatrix NewMatrixhalf(int m,int n);  
  TSMatrix NewMatrixhalfhalf(int m,int n);  

    //新建一个三元组表示的稀疏矩阵  
Status InsertElem(TSMatrix *M,int row,int col,double e);  
    //在三元组表示的稀疏矩阵M，第 row 行，第 col 列位置插入元素e  
    //插入成功，返回OK，否则返回ERROR  



 char name[40];
 char nname[40];
 char sname[40];
 char ename[40];



 

 
 
 
 
 
 
 


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



int main(int argc, char *argv[])
{
int arg;
 if(argc==1)
 {  	  
	 printf(  "\n#*******************************************************");
     printf(  "\n");
     printf(  "\n  本程序为论文 Preconditioners for the discretized time-harmonic .");   
     printf(  "\n   Maxwell  equations in mixed form  的混合有限元实现");

     printf(  "\n");
     printf(  "\n     +--------------------------|");
     printf(  "\n     |                          |");
     printf(  "\n     |   hydzhang@hotmail.com   |");
     printf(  "\n     |                          |");
     printf(  "\n     |                          |");
     printf(  "\n     |                          |");
     printf(  "\n     |                          |");
     printf(  "\n     +--------------------------|");
     printf(  "\n");
 	 printf(  "\n*******************************************************#\n\n\n\n\n\n");

	  
	  


 printf("编译完成后，请重新在命令行 下打开程序，格式如下：\n1) 输入程序名称之后，输入以 .n .s .e为后缀的文件的文件名（确保文件在当前命令行所在目录下）， 并用空格隔开;\n\n2) 如果三个文件同名，请不带后缀输入一次该名称即可。\n\n您没有提供EasyMesh生成的网格文件，程序将默认寻找名称为G1.n G1.s G1.e的文件 \n\n\n\n");
 getchar();
 getchar();

	 name[0]='G';
         name[1]='1';  name[2]='\0';


	  strcpy(nname,    name);    strcat(nname, ".n");
	  strcpy(sname,    name);    strcat(sname, ".s");
	  strcpy(ename,    name);    strcat(ename, ".e");


	  
	  
	  
	  
}else

 if(argc==2)
 {    strcpy(name,     argv[1]);
int len=strlen(name);
   if(name[len-2]=='.')
       name[len-2]='\0';

	  strcpy(nname,    name);    strcat(nname, ".n");
	  strcpy(sname,    name);    strcat(sname, ".s");
	  strcpy(ename,    name);    strcat(ename, ".e");
}
else
 if(argc!=4)
  {printf("编译完成后，请重新在命令行 下打开程序，格式如下：\n1) 输入程序名称之后，输入以 .n .s .e为后缀的文件的文件名（确保文件在当前命令行所在目录下）， 并用空格隔开;\n\n2) 如果三个文件同名，请不带后缀输入一次该名称即可。\n\n您没有提供EasyMesh生成的网格文件.\n");

	printf("\n缺少.n文件名参数，请输入:");
		scanf("%s",nname);
			printf("缺少.s文件名参数，请输入:");
		scanf("%s",sname);
			printf("缺少.e文件名参数，请输入:");
		scanf("%s",ename);

 }
 
 else
  {
    strcpy(nname,     argv[1]);
	  strcpy(sname,     argv[2]);
	    strcpy(ename,     argv[3]);
  }
   
//double pp=0.0;//原方程p项的系数，0或1
   



char dummy1[80];
FILE *nodes_file=NULL;FILE *sides_file=NULL,*elements_file=NULL;
int num_free=0,num_con=0;
  
  
 //以下转换网格格式给MFEM使用 
  FILE *Meshtrans=NULL;
  Meshtrans=fopen("D:\\zsy\\Meshtrans.mesh", "w");
 if (!Meshtrans)
    {
        perror("cannot open file Meshtrans");
    }
fprintf(Meshtrans, " MFEM mesh v1.0");fputc('\n',Meshtrans );fputc('\n',Meshtrans );fputc('\n',Meshtrans );fputc('\n',Meshtrans );

fprintf(Meshtrans, "dimension");fputc('\n',Meshtrans );
fprintf(Meshtrans, "%d",2);fputc('\n',Meshtrans );fputc('\n',Meshtrans );fputc('\n',Meshtrans );


fputc('\n',Meshtrans );fputc('\n',Meshtrans );fputc('\n',Meshtrans );fputc('\n',Meshtrans );


	fprintf(Meshtrans, "vertices");fputc('\n',Meshtrans );
	
//以下为.n文件的提取.*******************************************************************************************************************************************************************************  

nodes_file=fopen(nname,"r");

 if(!nodes_file)                                                                       
  {   fprintf(stderr, "%s \n", strerror(errno));
   printf("请确保是在网格文件所在目录下打开命令行；\nCannot open  nodes file  %s\n",nname);                                                 
   return -1;                                                                   
  }        
  
  
int dummyea=0;int dummyeb=0 ,dummymark,dummyfbe=0;int num_nodes,i,j;
load_i(nodes_file,&num_nodes);
 // printf("%d",num_nodes);system("pause");
 
 fprintf(Meshtrans, "%d",num_nodes);fputc('\n',Meshtrans );
  fprintf(Meshtrans, "%d",2);fputc('\n',Meshtrans );


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
	   	  fprintf(Meshtrans, "%f	%f",nodes[i].x,nodes[i].y);fputc('\n',Meshtrans );
       load_i(nodes_file, &nodes[i].ptr);//记录节点坐标;

	
	   if(nodes[i].ptr==FREE) {num_free++;}
	  }
	  
fputc('\n',Meshtrans );fputc('\n',Meshtrans );fputc('\n',Meshtrans );fputc('\n',Meshtrans );

	  num_con=num_nodes-num_free;
	   
	   
	   
 struct nptrs
{int fnodeptrs[num_free];
int cnodeptrs[num_con];
}nodeptrs;
double g[num_con];

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
		
		
	





	fprintf(Meshtrans, "boundary");fputc('\n',Meshtrans );
//以下为.s文件的提取.*******************************************************************************************************************************************************************************  

int num_sides;
sides_file=fopen(sname,"r");
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
       load_i(sides_file, &sides[i].y);//记录节点指标;
	   
	   load_i(sides_file, &dummyea);
	   load_i(sides_file, &dummyeb);	
        load_i(sides_file, &dummymark);		   
	   if(dummyea!=-1&&dummyeb!=-1){fprintf(outXs, "%18.12f	  %18.12f 	%18.12f 	%18.12f  ", nodes[sides[i].x].x,nodes[sides[i].x].y,nodes[sides[i].y].x,nodes[sides[i].y].y);  fputc('\n',outXs); // printf("dsfdsadfasdfasdf %d",dummymark);system("pause");
		   
		   sides[i].ptr=interior+1;interiorsides[interior]=i;sides[i].p=dummyea+1;sides[i].q=dummyeb+1;length[interior]=sqrt((nodes[sides[i].x].x-nodes[sides[i].y].x)*(nodes[sides[i].x].x-nodes[sides[i].y].x)+(nodes[sides[i].x].y-nodes[sides[i].y].y)*(nodes[sides[i].x].y-nodes[sides[i].y].y));/* printf("%f",length[interior]);system("pause"); */interior++;}
	   else 
	   {  // printf("22222222222222222222 %d",dummymark);
		   
		   sides[i].ptr=0;
		   if(dummyea==-1)
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
 
  fprintf(Meshtrans, "%d",num_sides-interior);fputc('\n',Meshtrans );
  
  for(i=0;i<num_sides;++i){
	  if(sides[i].ptr==0){
  fprintf(Meshtrans, "%d	%d	%d	%d",1,1,sides[i].x,sides[i].y);fputc('\n',Meshtrans );
	  }
  }
 fputc('\n',Meshtrans );fputc('\n',Meshtrans );fputc('\n',Meshtrans );fputc('\n',Meshtrans );

 
 

 
 
 
 
 
 	fprintf(Meshtrans, "elements");fputc('\n',Meshtrans );

//以下为.e文件的提取.*******************************************************************************************************************************************************************************  	  

int dummya ,dummyb,dummyc,dummyd,num_elements;
elements_file=fopen(ename,"r");

       
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
  fprintf(Meshtrans, "%d",num_elements);fputc('\n',Meshtrans );



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
	  
if(sides[dummya].y!=sides[dummyb].x&&sides[dummya].x!=sides[dummyb].x)//if  dummya  intersects with dummyb at dummyb.y
    {ss=(nodes[sides[dummya].x].x-nodes[sides[dummyb].x].x)*(nodes[sides[dummya].y].y-nodes[sides[dummyb].x].y)-(nodes[sides[dummya].x].y-nodes[sides[dummyb].x].y)*(nodes[sides[dummya].y].x-nodes[sides[dummyb].x].x);
    if(ss>0)//判断点在向量的左侧
	        {if(sides[dummya].y==sides[dummyb].y)
				 { elem[i].x=dummya+1; elem[i].y=-dummyb-1;   if(sides[dummyb].x==sides[dummyc].x) elem[i].z=dummyc+1;else elem[i].z=-dummyc-1;}
			else { elem[i].x=dummya+1; elem[i].z=dummyb+1;   if(sides[dummyb].x==sides[dummyc].y) elem[i].y=dummyc+1;else elem[i].y=-dummyc-1; }
   			}
	else
	        {if(sides[dummya].x==sides[dummyb].y)
				 { 	elem[i].x=-dummya-1; elem[i].y=-dummyb-1;   if(sides[dummyb].x==sides[dummyc].x) elem[i].z=dummyc+1;else elem[i].z=-dummyc-1;}
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
	
	
		if(elem[i].x>0)
			if(elem[i].y>0)
		{fprintf(Meshtrans, "%d	%d	%d	%d	%d	",2,2,sides[elem[i].x-1].x,sides[elem[i].x-1].y,sides[elem[i].y-1].y);fputc('\n',Meshtrans);}
	else		{fprintf(Meshtrans, "%d	%d	%d	%d	%d	",2,2,sides[elem[i].x-1].x,sides[elem[i].x-1].y,sides[-elem[i].y-1].x);fputc('\n',Meshtrans);	 }
	
	else    
		if(elem[i].y>0)
		{fprintf(Meshtrans, "%d	%d	%d	%d	%d	",2,2,sides[-elem[i].x-1].y,sides[-elem[i].x-1].x,sides[elem[i].y-1].y);fputc('\n',Meshtrans);	}

	else		{fprintf(Meshtrans, "%d	%d	%d	%d	%d	",2,2,sides[-elem[i].x-1].y,sides[-elem[i].x-1].x,sides[-elem[i].y-1].x);fputc('\n',Meshtrans); }



	
	 //单元上的边,正负标明方向.

	  
	  load_s(elements_file, dummy1); 
	  load_s(elements_file, dummy1); 
	  load_i(elements_file, &dummyd); 
 
  }  
  fclose(Meshtrans);

  fclose(elements_file);
  printf("共有点:%d\n边:%d\n元素:%d\n自由边界边:%d\n自由节点m=%d\n内部边n=%d\nn+m=%d\n",num_nodes,num_sides,num_elements,num_fbe,num_free,interior,num_free+interior);   system("pause");

   
//数据提取完毕,以下为算法第一部分:刚度矩阵的计算

//定义curl-curl稀疏矩阵，第一行为行标 第二行为列标，第三行为值,不排序，不允许重复。
//注意：A	L ML MM 计算的都是只有上三角部分
//而Pcurl, B计算的是完整部分
	



TSMatrix A=NewMatrix(interior,interior);


TSMatrix MM=NewMatrix(interior,interior);

TSMatrix B=NewMatrixhalf(num_free,interior);

TSMatrix  L=NewMatrixhalfhalf(num_free,num_free);
TSMatrix  ML=NewMatrixhalfhalf(num_free,num_free);





double* FinF=NULL;
	FinF = (double *)malloc(sizeof(double)*(interior+num_free)); memset(FinF,0,(interior+num_free)*sizeof(double)); 


int eptr[3]={0},indices[3]={0},ptrs[3]={0};
double delta,coords[3][2],I,x0,y0,invM[2][3];double det,G[3][3];
double kk=0;int s,r,ceshi=-1;double err1=0;
















int k;





double localength[3];
double aaa[3];//记录m逆矩阵的首行,也就是有限元基的常数项
double fff[3][3];double emm[3][3];
memset(emm,0,sizeof(emm));



double gg[interior]; memset(gg,0,sizeof(double)*interior);	






for(k=0;k<num_elements;k++)  //对三角单元循环
//以下为algorithm 7.1 *************************
{	


eptr[0]=elem[k].x;//三条边,eptr指示真实的三角形逆时针三边}
eptr[1]=elem[k].y;
eptr[2]=elem[k].z;

if(eptr[0]>0)
	     {indices[0]=sides[eptr[0]-1].x;indices[1]=sides[eptr[0]-1].y;  if(sides[eptr[0]-1].ptr>0) localength[0]=length[sides[eptr[0]-1].ptr-1]; else localength[0]=-1.0;}
 else {indices[0]=sides[-eptr[0]-1].y;indices[1]=sides[-eptr[0]-1].x;if(sides[-eptr[0]-1].ptr>0) localength[0]=length[sides[-eptr[0]-1].ptr-1]; else localength[0]=-1.0;}
		   
if(eptr[1]>0)
		 {indices[2]=sides[eptr[1]-1].y;if(sides[eptr[1]-1].ptr>0) localength[1]=length[sides[eptr[1]-1].ptr-1]; else localength[1]=-1.0;}
else {indices[2]=sides[-eptr[1]-1].x;if(sides[-eptr[1]-1].ptr>0) localength[1]=length[sides[-eptr[1]-1].ptr-1]; else localength[1]=-1.0;}
	 
if(eptr[2]>0) {if(sides[eptr[2]-1].ptr>0) localength[2]=length[sides[eptr[2]-1].ptr-1]; else localength[2]=-1.0;}
else {if(sides[-eptr[2]-1].ptr>0) localength[2]=length[sides[-eptr[2]-1].ptr-1]; else localength[2]=-1.0;}

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
  
  
  
 



//******************************algorithm 7.2

  
  for(r=0;r<3;r++)
  {for(s=r;s<3;s++)	 
	  G[r][s]=invM[0][r]*invM[0][s]+invM[1][r]*invM[1][s];
  }

   
  for(r=1;r<3;r++)
  {for(s=0;s<r;s++)	 
	  G[r][s]=G[s][r];
  }
  
  
		    
  
  
  
delta=0.5 * det;  
/* printf("%f ",delta);system("pause"); */

  x0=(coords[0][0]+coords[1][0]+ coords[2][0])/3.0       ;y0=(coords[0][1]+coords[1][1]+ coords[2][1])/3.0 ;
//系数k的输入点----------------------------------------------------------------------------------------------------------------------------------------------------------
I=delta;
   
  
  for(r=0;r<3;r++)
  {for(s=r;s<3;s++)	 
	  
	  {if(ptrs[r]>0&&ptrs[s]>0)
		  {		  i=min(ptrs[r],ptrs[s])-1;j=max(ptrs[r],ptrs[s])-1;
		        doubleelem=G[r][s]*I; if ( !InsertElem( &L,i,j, doubleelem) )
  { printf("\nError:L插入失败\n");  
        return ERROR;  
    }  }
	    
	  }
}  
  


  
  
  
    for(r=0;r<3;r++)
  {for(s=r;s<3;s++)	   
    	  {if(ptrs[r]>0&&ptrs[s]>0)
		  {		  i=min(ptrs[r],ptrs[s])-1;j=max(ptrs[r],ptrs[s])-1;  
  
         doubleelem= I*(
  1.0/3*(         (aaa[r]+invM[0][r]* (2.0/3*coords[0][0]+1.0/6*coords[1][0]+1.0/6*coords[2][0]) +invM[1][r]*(2.0/3*coords[0][1]+1.0/6*coords[1][1]+1.0/6*coords[2][1]) )*  
                  (aaa[s]+invM[0][s]* (2.0/3*coords[0][0]+1.0/6*coords[1][0]+1.0/6*coords[2][0]) +invM[1][s]*(2.0/3*coords[0][1]+1.0/6*coords[1][1]+1.0/6*coords[2][1]) )      ) +
 1.0/3*(         (aaa[r]+invM[0][r]* (1.0/6*coords[0][0]+2.0/3*coords[1][0]+1.0/6*coords[2][0]) +invM[1][r]*(1.0/6*coords[0][1]+2.0/3*coords[1][1]+1.0/6*coords[2][1]) )*  
                  (aaa[s]+invM[0][s]* (1.0/6*coords[0][0]+2.0/3*coords[1][0]+1.0/6*coords[2][0]) +invM[1][s]*(1.0/6*coords[0][1]+2.0/3*coords[1][1]+1.0/6*coords[2][1]) )	   )+
 1.0/3*(         (aaa[r]+invM[0][r]* (1.0/6*coords[0][0]+1.0/6*coords[1][0]+2.0/3*coords[2][0]) +invM[1][r]*(1.0/6*coords[0][1]+1.0/6*coords[1][1]+2.0/3*coords[2][1]) )*  
                  (aaa[s]+invM[0][s]* (1.0/6*coords[0][0]+1.0/6*coords[1][0]+2.0/3*coords[2][0]) +invM[1][s]*(1.0/6*coords[0][1]+1.0/6*coords[1][1]+2.0/3*coords[2][1]) )	   )	
                )	; if ( !InsertElem( &ML,i,j, doubleelem) )
  { printf("\nError:ML插入失败\n");  
        return ERROR;  
    }  
				  
	//	printf("%f", ML[i][j]);system("pause");
          }   
		  }
  }
  
  

  
  //load矩阵的初始化和计算.*******************************************************************************************************************************************************************************  	  



//以下为algorithm 7.5的第一部分************************************
  
  for(s=0;s<3;s++)

      {if(localength[s]>=0)
		  {
			  i=sides[abs(eptr[s])-1].ptr-1;
	          gg[i]+= delta*sign(eptr[s])*(( 2.0 )  *(aaa[s]*invM[0][(1+s)%3]-aaa[(1+s)%3]*invM[0][s]+y0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s]))	+
              2.0     * (aaa[s]*invM[1][(1+s)%3]-aaa[(1+s)%3]*invM[1][s]-x0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s])	)   );
           }
      }
  
  
  
  
for(r=0;r<3;r++)
{

	for(s=r;s<3;s++)	 
	  
	  {
		  if(localength[r]>=0&&localength[s]>=0)
		  {/* printf("r=%d;s= %d;",r,s); */i=min(sides[abs(eptr[r])-1].ptr-1,sides[abs(eptr[s])-1].ptr-1);j=max(sides[abs(eptr[r])-1].ptr-1,sides[abs(eptr[s])-1].ptr-1);/* printf("%d %d",i,j);system("pause"); */
		doubleelem= sign(eptr[r])*sign(eptr[s])/(delta);/* printf("%f",A[i][j]); */ if ( !InsertElem( &A,i,j, doubleelem) )
  { printf("\nError:A插入失败\n");  
        return ERROR;  
    }   }
	   
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
  doubleelem=emm[r][s]; if ( !InsertElem( &MM,i,j, doubleelem) )
  { printf("\nError:插入失败\n");  
        return ERROR;  
    }   }
      }
  }
  
  
  
  


  
  
  

  

  
  
  
  
for(r=0;r<3;r++)
  {for(s=0;s<3;s++)	// (s=0 s=1 )
	  
      {if(localength[s]>=0&&ptrs[r]>0)
		  {j=sides[abs(eptr[s])-1].ptr-1;//printf("(%d)===%d,%d===",r,s,j);
	       i=ptrs[r]-1;
	       doubleelem= delta*sign(eptr[s])*(  invM[0][r]*(aaa[s]*invM[0][(1+s)%3]-aaa[(1+s)%3]*invM[0][s]+y0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s]))	+
	  invM[1][r]* (aaa[s]*invM[1][(1+s)%3]-aaa[(1+s)%3]*invM[1][s]+x0*(invM[1][(1+s)%3]*invM[0][s]-invM[1][s]*invM[0][(1+s)%3])	)   ); if ( !InsertElem( &B,i,j, doubleelem) )
  { printf("\nError:插入失败\n");  
        return ERROR;  
    }   
	      }
	  }
  }
  


//三角基元对刚度矩阵循环完毕,但是同时要对load矩阵循环,所以不加大括号 
 

  


}    //对三角基元循环完毕.

//以下计算C




	

 
FILE *outA;

    outA = fopen("D:\\zsy\\A.txt", "w");
 if (!outA)
    {
        perror("cannot open file");
        
    }
	
	
	
	  for (i = 0; i <A.tu; i++)
	  {
	  fprintf(outA, "%d	 ",A.data[i].i);
	  
	  fprintf(outA, "%d	 ",A.data[i].j);
	  fprintf(outA, "%20.12f ",A.data[i].e);fputc('\n',outA );

	  
	  }
fclose(outA);

FILE *out1;


    out1 = fopen("D:\\zsy\\FinF.txt", "w");
 if (!out1)
    {
        perror("cannot open file");
        
    }
	  for (i = 0; i <interior; i++)
    {fprintf(out1, "%20.12f ", gg[i]);fputc('\n',out1);}

      for (i = interior; i <interior+num_free; i++)
     {fputc('0',out1);fputc('\n',out1);} 
 
fclose(out1);


 
FILE *outB;

    outB = fopen("D:\\zsy\\B.txt", "w");
 if (!outB)
    {
        perror("cannot open file");
        
    }
	
	
	
	  for (i = 0; i <B.tu; i++)
	  {
	  fprintf(outB, "%d	 ",B.data[i].i);
	  
	  fprintf(outB, "%d	 ",B.data[i].j);
	  fprintf(outB, "%20.12f ",B.data[i].e);fputc('\n',outB );

	  
	  }
fclose(outB);
	
	FILE *outL;

    outL = fopen("D:\\zsy\\L.txt", "w");
 if (!outL)
    {
        perror("cannot open file");
        
    }
	
	
	
	  for (i = 0; i <L.tu; i++)
	  {
	  fprintf(outL, "%d	 ",L.data[i].i);
	  
	  fprintf(outL, "%d	 ",L.data[i].j);
	  fprintf(outL, "%20.12f ",L.data[i].e);fputc('\n',outL);

	  
	  }
  
  
  
fclose(outL);


FILE *outML;

    outML = fopen("D:\\zsy\\ML.txt", "w");
 if (!outML)
    {
        perror("cannot open file");
        
    }
	
	
	  for (i = 0; i <ML.tu; i++)
	  {
	  fprintf(outML, "%d	 ",ML.data[i].i);
	  
	  fprintf(outML, "%d	 ",ML.data[i].j);
	  fprintf(outML, "%20.12f ",ML.data[i].e);fputc('\n',outML);

	  
	  }
  
fclose(outML);





 
FILE *outnn;


    outnn = fopen("D:\\zsy\\MM.txt", "w");
 if (!outnn)
    {
        perror("cannot open file");
        
    } 
	
	  for (i = 0; i <MM.tu; i++)
	  {
	  fprintf(outnn, "%d	 ",MM.data[i].i);
	  
	  fprintf(outnn, "%d	 ",MM.data[i].j);
	  fprintf(outnn, "%20.12f ",MM.data[i].e);fputc('\n',outnn);

	  
	  }
    
    fclose(outnn);




free(A.data);
free(MM.data);
free(L.data);
free(B.data);
free(ML.data);







TSMatrix C=NewMatrixhalfhalf(interior,num_free);
TSMatrix Pcurl=NewMatrix(interior,2*num_free);

 for(i=0;i<interior;i++)//边循环
 {
	//ptr>0指示自由节点
if(nodes[sides[interiorsides[i]].x].ptr>0)
{
 if ( !InsertElem( &C,i,nodes[sides[interiorsides[i]].x].ptr-1, -1.0) )
  { printf("\nError:插入失败\n");  
        return ERROR;  
    }   
}

if(nodes[sides[interiorsides[i]].y].ptr>0)
{
 if ( !InsertElem( &C,i,nodes[sides[interiorsides[i]].y].ptr-1, 1.0) )
  { printf("\nError:插入失败\n");  
        return ERROR;  
    }   
 }
}

//边循环
 for(i=0;i<interior;i++)
 {
if(nodes[sides[interiorsides[i]].x].ptr>0)
{doubleelem=0.5*(nodes[sides[interiorsides[i]].y].x-nodes[sides[interiorsides[i]].x].x);
 if ( !InsertElem( &Pcurl,i,nodes[sides[interiorsides[i]].x].ptr-1, doubleelem) )
  { printf("\nError:插入失败\n");  
        return ERROR;  
    }   
doubleelem=0.5*(nodes[sides[interiorsides[i]].y].y-nodes[sides[interiorsides[i]].x].y);

 if ( !InsertElem( &Pcurl,i,nodes[sides[interiorsides[i]].x].ptr-1+num_free, doubleelem) )
  { printf("\nError:插入失败\n");  
        return ERROR;  
    }   
}

if(nodes[sides[interiorsides[i]].y].ptr>0)
{
	

doubleelem=0.5*(nodes[sides[interiorsides[i]].y].x-nodes[sides[interiorsides[i]].x].x);
 if ( !InsertElem( &Pcurl,i,nodes[sides[interiorsides[i]].y].ptr-1, doubleelem) )
  { printf("\nError:插入失败\n");  
        return ERROR;  
    }   
doubleelem=0.5*(nodes[sides[interiorsides[i]].y].y-nodes[sides[interiorsides[i]].x].y);

 if ( !InsertElem( &Pcurl,i,nodes[sides[interiorsides[i]].x].ptr-1+num_free, doubleelem) )
  { printf("\nError:插入失败\n");  
        return ERROR;  
    }  
}



 }






FILE *outPcurl;


    outPcurl = fopen("D:\\zsy\\Pcurl.txt", "w");
 if (!outPcurl)
    {
        perror("cannot open file");
        
    }

	  for (i = 0; i <Pcurl.tu; i++)
	  {
	  fprintf(outPcurl, "%d	 ",Pcurl .data[i].i);
	  
	  fprintf(outPcurl, "%d	 ",Pcurl.data[i].j);
	  fprintf(outPcurl, "%20.12f ",Pcurl.data[i].e);fputc('\n',outPcurl );

	  }
fclose(outPcurl);
//printf("\n  %d Run read with the command: \n read namen.n names.s namee.e\n ",k);system("pause");





FILE *outC;
    outC = fopen("D:\\zsy\\C.txt", "w");
 if (!outC)
    {
        perror("cannot open file");
        
    }
	
	
	  for (i = 0; i <C.tu; i++)
	  {
	  fprintf(outC, "%d	 ",C.data[i].i);
	  
	  fprintf(outC, "%d	 ",C.data[i].j);
	  fprintf(outC, "%20.12f ",C.data[i].e);fputc('\n',outC);

	  
	  }
 
fclose(outC);





printf("Congratulations! end of program!");

return 0;
}

  
TSMatrix NewMatrix(int m,int n){  
    //新建一个三元组表示的稀疏矩阵  
Triple* adata=NULL;

	adata = (Triple *)malloc(sizeof(Triple)*MAXSIZE); memset(adata,0,MAXSIZE*sizeof(Triple)); 
	
    TSMatrix M; 
M.data=	adata;
 // M.mu=m;  
 //  M.nu=n;  
    M.tu=0;  
    return M;  
}  


  TSMatrix NewMatrixhalf(int m,int n){  
    //新建一个三元组表示的稀疏矩阵  
Triple* adata=NULL;

	adata = (Triple *)malloc((sizeof(Triple)*MAXSIZE/2)); memset(adata,0,(MAXSIZE*sizeof(Triple)/2)); 
	
    TSMatrix M; 
M.data=	adata;
  //  M.mu=m;  
 //  M.nu=n;  
    M.tu=0;  
    return M;  
}  
  
  TSMatrix NewMatrixhalfhalf(int m,int n){  
    //新建一个三元组表示的稀疏矩阵  
Triple* adata=NULL;

	adata = (Triple *)malloc((sizeof(Triple)*MAXSIZE/4)); memset(adata,0,(MAXSIZE*sizeof(Triple)/4)); 
	
    TSMatrix M; 
M.data=	adata;
   // M.mu=m;  
  // M.nu=n;  
    M.tu=0;  
    return M;  
}  
  
  
  
  
Status InsertElem(TSMatrix *M,int row,int col,ElemType e){  
    //在三元组表示的稀疏矩阵M，第 row+1 行，第 col+1 列位置插入元素e  
    //插入成功，返回OK，否则返回ERROR  
    int z,t;  
    if(M->tu>=MAXSIZE){//当前三元组表已满  
        printf("\nError:There is no space in the matrix;\n");  
        return ERROR;  
    }  
 // if(M->mu<=row||M->nu<=col){//当前三元组表已满  
        // printf("\nError:too many rows or cols;\n");  
        // return ERROR;  
    // }  

	
    for(z=0;z<M->tu;++z)  
        if(M->data[z].i==row&&M->data[z].j==col){  
         M->data[z].e+=e;  
            return TRUE;  
        }  
	   M->data[M->tu].i=row;  
        M->data[M->tu].j=col;  
        M->data[M->tu].e=e;  
        M->tu++;  
		 return OK;  
}  

  

 /*
  
%第一行为行标 第二行为列标，第三行为值,不排序，不允许重复。
%注意：A	L ML MM 计算的都是只有上三角部分
%而Pcurl, B计算的是完整部分
	
  clear
  load D:\zsy\A.txt; load D:\zsy\B.txt;load D:\zsy\C.txt;
  %load D:\zsy\FinA.txt;
  load D:\zsy\FinF.txt;
  load D:\zsy\L.txt;load D:\zsy\ML.txt;
  load D:\zsy\MM.txt;load D:\zsy\Pcurl.txt;
n=max(A(:,1))+1;
m=max(L(:,1))+1;
%m+n

  A=sparse(A(:,1)+1,A(:,2)+1,A(:,3),n,n); A=A+A'-diag(diag(A));
  B=sparse(B(:,1)+1,B(:,2)+1,B(:,3),m,n);  
  C=sparse(C(:,1)+1,C(:,2)+1,C(:,3),n,m); 
  L=sparse(L(:,1)+1,L(:,2)+1,L(:,3),m,m); L=L+L'-diag(diag(L));
  ML=sparse(ML(:,1)+1,ML(:,2)+1,ML(:,3),m,m); ML=ML+ML'-diag(diag(ML));
  MM=sparse(MM(:,1)+1,MM(:,2)+1,MM(:,3),n,n); MM=MM+MM'-diag(diag(MM));
  Pcurl=sparse(Pcurl(:,1)+1,Pcurl(:,2)+1,Pcurl(:,3),n,2*m); 
  

%invL=inv(L);
  
  
%max(max(abs(MM*C-B')))

%max(max(abs(B*C-L)))

%max(max(abs(A*C)))
 
 k=2; rho=6;
 FinA(1:n,1:n)=A-k*k*MM;
 FinA(1:n,n+1:n+m)=B';
 FinA(n+1:n+m,1:n)=B;
save
  
*/
