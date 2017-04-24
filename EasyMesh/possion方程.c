#include<math.h>
#include<stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef PI
#define PI    3.14159265359
#endif
#define MAX_NODES 3000

#define FREE 0
/*#define  FileAddressn  C:\\Users\\zhang\\Desktop\\easymesh\\example.n
#define  FileAddresss  C:\\Users\\zhang\\Desktop\\easymesh\\example.s
#define  FileAddresse  C:\\Users\\zhang\\Desktop\\\easymesh\\example.e*/


//由于是局部定义数组（全局变量数组中长度不能使用变量）网格数不能太密，否则可能溢出，需要调整程序，使用malloc对数组动态分配内存。

/*=========================================================================*/



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
	
	
/*int arg, ans, Nn0;
 if(argc<4)
  {printf("\n  Run read with the command: \n read namen.n names.s namee.e\n ");
 return 0;}
 else
  */
   
   
   
   
   


char dummy1[80];
FILE *nodes_file=NULL;FILE *sides_file=NULL,*elements_file=NULL;
int num_free=0,num_con=0;
  
  
  
  
//以下为.n文件的提取。*******************************************************************************************************************************************************************************  

nodes_file=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\example.n","r");

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
double g[num_con];double F[num_free],K[num_free][num_free];//在获得num_free和num_con后申明相关量。
memset(F,0,sizeof(double)*num_free); memset(K,0,sizeof(double)*num_free*num_free);


/*   printf("\n  %d Run read with the command: \n read namen.n names.s namee.e\n ",i);system("pause");
//用于错误测试
 
FILE *outn;


    outn = fopen("D:\\Fn.txt", "w");
 if (!outn)
    {
        perror("cannot open file");
        
    } 
	  for (i = 0; i<num_free; i++)
    {for (num_freedummy = 0; num_freedummy<num_free; num_freedummy++){
       
	fprintf(outn, "%f ",K[i][num_freedummy]);}
        
        fputc('\n',outn);
    }
    fclose(outn);
system("pause");
    */
/*  
FILE *outn;


    outn = fopen("D:\\Fn.txt", "w");
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
FILE *outX;


    outX = fopen("D:\\X.txt", "w");
 if (!outX)
    {
        perror("cannot open file  outX");
        
    } //输出自由变量得坐标值
	 






	   
for(i=0;i<num_nodes;i++)
{
	if(nodes[i].ptr==FREE)
	      {nodes[i].ptr=num_freedummy+1;
          nodeptrs.fnodeptrs[num_freedummy]=i;
		  fprintf(outX, "%f  %f\n ", nodes[i].x,nodes[i].y);   //输出自由变量的坐标值
          num_freedummy++;}
	  
	   
	      else  {nodes[i].ptr=-(i-num_freedummy+1);
	            nodeptrs.cnodeptrs[i-num_freedummy]=i;//记录节点指针和自由节点，强制节点;负值表示排序为多少位的强制节点
		        g[i-num_freedummy]=sin(PI*nodes[i].x)*sin(PI*nodes[i].y)+nodes[i].x;}
}    //此处为Dirichlet条件g表达式的输入点。nodes[i].x node[i].y为坐标值。
	   

fclose(outX);
	
		
	
		
		
fclose(nodes_file)	;
		
		
	







//以下为.s文件的提取。*******************************************************************************************************************************************************************************  

int num_sides;
sides_file=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\example.s","r");
if(NULL==sides_file)                                                                       
  {                                                                             
   printf("can't open sides file\n");                                                 
   return -1;                                                                   
  } 





load_i(sides_file,&num_sides);
	  
	  
	  
	  
	  
	  
struct sid
{ int x;int y;int p,q;
}sides[num_sides];//p,q指示边的两边的elements。








int fbndyedges[num_sides],dummy_fbe=0;

for(i=0;i<num_sides;i++)
{
	  load_s(sides_file, dummy1);     
       load_i(sides_file, &sides[i].x);
       load_i(sides_file, &sides[i].y);//记录节点坐标;
	   
	   load_i(sides_file, &dummyea);
	   load_i(sides_file, &dummyeb);	
        load_i(sides_file, &dummymark);		   
	   if(dummyea!=-1&&dummyeb!=-1){sides[i].p=dummyea+1;sides[i].q=dummyeb+1;}
	   else 
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
 int num_fbe=dummy_fbe;  double ss;
fclose(sides_file);







//以下为.e文件的提取。*******************************************************************************************************************************************************************************  	  

int dummya ,dummyb,dummyc,dummyd,num_elements;
elements_file=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\example.e","r");

       
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
	  
	  
//以下用按逆时针相互衔接的顺序将三边记录于elem中，负数代表方向与边的本来的方向相反	  
	  
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


	  
	  
	  
	  
	  
	  






  
	  
	 //单元上的边，正负标明方向。

	  
	  load_s(elements_file, dummy1); 
	  load_s(elements_file, dummy1); 
	  load_i(elements_file, &dummyd); 
 
  }
  fclose(elements_file);
  printf("共有点：%d\n边：%d\n元素：%d\n自由边：%d\n自由节点：%d\n",num_nodes,num_sides,num_elements,num_fbe,num_free);

  
  



  //数据提取完毕，以下为算法第一部分：刚度矩阵的计算
  
  
int eptr[3],indices[3],ptrs[3];double A,coords[3][2],I,x0,y0,invM[2][3];double det,G[3][3];
int k,s,r,ceshi=-1;


double f=0.0;//f的输入-------------------------------------------------------------------------------------------------------------------------------------------------------------

double w[3];//记录G的边界值


for(k=0;k<num_elements;k++)
//以下为algorithm 7.1 *************************
{ceshi=-1;
	
	
	
	
	
	eptr[0]=elem[k].x;//三条边}
eptr[1]=elem[k].y;
eptr[2]=elem[k].z;

if(eptr[0]>0)
	     {indices[0]=sides[eptr[0]-1].x;indices[1]=sides[eptr[0]-1].y;}
           else {indices[0]=sides[-eptr[0]-1].y;indices[1]=sides[-eptr[0]-1].x;}
 	 if(eptr[1]>0)
		 {indices[2]=sides[eptr[1]-1].y;}
	 else {indices[2]=sides[-eptr[1]-1].x;}


for(i=0;i<3;i++)
             {ptrs[i]=nodes[indices[i]].ptr;}//此处ptrs继承ptr得特性 没有更改。

	         for(i=0;i<3;i++)
             {coords[i][0]=nodes[indices[i]].x;
		 }//三条边}

/* for(i=0;i<3;i++){
 printf(" %d  %d  %d \n ",indices[0],indices[1],indices[2]);system("pause");



} */

for(i=0;i<3;i++)
{coords[i][1]=nodes[indices[i]].y;
/*printf("\n  %d \n ",indices[i]);system("pause");*/
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
  

  
  
 
FILE *outn;


    outn = fopen("D:\\length.txt", "w");
 if (!outn)
    {
        perror("cannot open file");
        
    } 
	  for (i = 0; i<2; i++)
   { for (j = 0; j<3; j++)
       
	fprintf(outn, "%f ", invM[i][j]); fputc('\n',outn);}
        
       
    
    fclose(outn);system("pause");

  
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


    outn = fopen("D:\\Fn.txt", "w");
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
  
A=0.5 * (      coords[0][0]*coords[1][1]+ coords[1][0]*coords[2][1] +coords[2][0]*coords[0][1]-coords[0][0]*coords[2][1]-coords[1][0]*coords[0][1]-coords[2][0]*coords[1][1] );
/* printf("%f ",A);system("pause"); */

  x0=(coords[0][0]+coords[1][0]+ coords[2][0])/3.0       ;y0=(coords[0][1]+coords[1][1]+ coords[2][1])/3.0 ;
//系数k的输入点----------------------------------------------------------------------------------------------------------------------------------------------------------
I=A*1;
//I=A*K(x0,y0);
   
  
  for(r=0;r<3;r++)
  {for(s=r;s<3;s++)	 
	  
	  {if(ptrs[r]>0&&ptrs[s]>0)
		  {		  i=min(ptrs[r],ptrs[s])-1;j=max(ptrs[r],ptrs[s])-1;
		  K[i][j]+=G[r][s]*I;}
	    
	  }
}  
  
/*   
  
  printf("\n  %d Run read with the command: \n read namen.n names.s namee.e\n ",k);system("pause");

 
FILE *outn;


    outn = fopen("D:\\.txt", "w");
 if (!outn)
    {
        perror("cannot open file");
        
    } 
	  for (i = 0; i<num_free; i++)
    {for (num_freedummy = 0; num_freedummy<num_free; num_freedummy++){
       
	fprintf(outn, "%f ",K[i][num_freedummy]);}
        
        fputc('\n',outn);
    }
    fclose(outn);
system("pause"); */
  
//｝三角基元对刚度矩阵循环完毕，但是同时要对load矩阵循环，所以不加大括号 
 

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
//以下为算法第二部分：load矩阵的初始化和计算。*******************************************************************************************************************************************************************************  	  



//以下为algorithm 7.5的第一部分************************************



//******************************
if(f!=1.1)//if(f!={0})
{
//计算面积
A=0.5 * (      coords[0][0]*coords[1][1]+ coords[1][0]*coords[2][1] +coords[2][0]*coords[0][1]-coords[0][0]*coords[2][1]-coords[1][0]*coords[0][1]-coords[2][0]*coords[1][1] );
x0=(coords[0][0]+coords[1][0]+ coords[2][0])/3.0       ;y0=(coords[0][1]+coords[1][1]+ coords[2][1])/3.0 ;
f=2*PI*PI*sin(PI*x0)*sin(PI*y0);//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
I=A*f/3.0;


ceshi=k;
for(r=0;r<3;r++)
{if(ptrs[r]>0){F[ptrs[r]-1]+=I;}
}
}

//以下为algorithm 7.5的第二部分************************************
if(/*g!={0}&&*/ptrs[0]<0||ptrs[1]<0||ptrs[2]<0)
{
if(ceshi!=k)
{A=0.5 * (      coords[0][0]*coords[1][1]+ coords[1][0]*coords[2][1] +coords[2][0]*coords[0][1]-coords[0][0]*coords[2][1]-coords[1][0]*coords[0][1]-coords[2][0]*coords[1][1] );
x0=(coords[0][0]+coords[1][0]+ coords[2][0])/3.0       ;y0=(coords[0][1]+coords[1][1]+ coords[2][1])/3.0 ;
}
I=A*1;
//I=A*K(x0,y0);------------------------------------------------------------------------------------------



for(r=0;r<3;r++)
{if(ptrs[r]<0)
{w[r]=g[-ptrs[r]-1];
}
else w[r]=0;
}



for(r=0;r<3;r++)
//计算梯度
{if(ptrs[r]>0)
{F[ptrs[r]-1]+=-(w[0]*G[0][r]+w[1]*G[1][r]+w[2]*G[2][r])*I;}   

/* if(!(F[ptrs[r]-1]<33.9&&F[ptrs[r]-1]>-66))
{printf("%d %f %f %f %f",k,G[0][r],G[1][r],G[2][r],F[ptrs[r]-1]);system("pause");} 
 */
}


}


}//对三角基元循环完毕。





FILE *out;


    out = fopen("D:\\A.txt", "w");
 if (!out)
    {
        perror("cannot open file");
       
    }
	  for (i = 0; i <num_free; i++)
    {
		
        for (j = 0; j <num_free; j++)
        {
           if(j>=i) fprintf(out, "%.12f ", K[i][j]);
		   else  fprintf(out, "%.12f ", K[j][i]);
        }
        fputc('\n',out);
    }
    fclose(out);
    
//printf("\n  Run read with the command: \n read namen.n names.s namee.e\n ");


FILE *out1;


    out1 = fopen("D:\\F.txt", "w");
 if (!out1)
    {
        perror("cannot open file");
        
    }
	  for (i = 0; i <num_free; i++)
    {
       
            fprintf(out1, "%.12f ", F[i]);
        
        fputc('\n',out1);
    }
    fclose(out1);
  return 0;
}









































