//���뻷����windows 10 Chinese;notepad++;MinGW.
//gcc -std=c99
//********************************************************************************************************************************************
/*������Ϊ���� Preconditioners for the discretized time-harmonic Maxwell
equations in mixed form  �Ļ������Ԫʵ��

********************************************************************************************************************************************
������ֻ���� �� u��x,y��=1-y^2,1-x^2), p = (1 ? x^2 )(1 ? y^2 )��the unit square(-1=<x=<1,-1=<y=<1,)�»����ȷ���Ҷ���.
�����err1;��ע����ĳ����е�kֵ
������ɺ��������������ļ����ڵ�ͬһ��Ŀ¼�´������� ����
********************************************************************************************************************************************
������EasyMesh(https://zsy.wodemo.com/file/387961 ; https://zsy.wodemo.com/cat/4506��ѹ����������������)���������Ϣ��EasyMesh�� ���Ը�ϸ���������޸�EasyMesh�Ĵ�С���ƣ���  #define MAX_NODES 8000�����������������ļ�
EasyMeshʹ�÷���������EasyMesh�� ��EasyMesh.exe�����ļ��а�Shift+�Ҽ�����  �ڴ˴�������ڣ�W�� ѡ��  ���� EasyMesh G6 +dxf  (�˴�G6����G6.d,ΪEasyMesh�����ļ�)
C:\\Users\\zhang\\Desktop\\easymesh\\G6.n  
C:\\Users\\zhang\\Desktop\\easymesh\\G6.s  
C:\\Users\\zhang\\Desktop\\easymesh\\G6.e  
�����������������G6������������Ӧ�ļ������ļ�λ��
���Ҹ��Ĳ���kk��ֵ����������kֵ����
������Ӧ��A��C FinF���Ҷ����M��MM�� ,L�� FinA(ϵ������),д������С�
�ɸĽ����򣺱���������ݽṹ��ȫû���Ż������⽨�����matlab����������ϡ����ʽ�洢����
********************************************************************************************************************************************
��������ȡ���ݵ�matlab����double U[8000]���������� ���Բ��ù����������������(err1)����Ҫ������һ���������matlab����

%MatLab�� R2014a
clear;clc;
FinA=load('D:\FinA.txt');FinF=load('D:\FinF.txt');MM=load('D:\MM.txt');L=load('D:\L.txt');C=load('D:\C.txt');A=load('D:\A.txt');
[n,m]=size(C);
B=FinA(n+1:m+n,1:n);


%������Ҫ����ȷ����FinF��ǰ����ʹ�ã������϶�˵����
x=FinA\FinF;
u=x(1:n);
%p����֤:����������ļ�Xn.txt��ע���ļ���λ��
[xn1,xn2]=textread('D:\Xn.txt');
perr=(1.-xn1.*xn1).*(1.-xn2.*xn2)-x((n+1):(m+n));

%u����֤��
%��uֵд��u.txt;
dlmwrite('C:\Users\zhang\Desktop\easymesh\u.txt',u,'precision', '%.16f', ...
'newline', 'pc');
%dlmwrite('C:\Users\zhang\Desktop\easymesh\u.txt',u,'precision', '%.16f', 'newline', 'pc');
%����  ��  save 'C:\Users\zhang\Desktop\easymesh\u.txt' u -ascii;
Ȼ���ٴ����иó��������ȷerr1ֵ��ע��ÿ��������Ҫ��������matlab����ע�������������䣺FileU=fopen("C:\\Users\\zhang\\Desktop\\easymesh\\u.txt" , "r")  ;��


ע�����matlab���޸�kֵ��Ȼ������һ�´��뼴����ȫ������Ч�ĳ�FinF������Ч�����ݡ�����������kk=0 ���г��򼴿ɡ�
�Ҷ���ɲ������ֵ��
k=1.5;FinA(1:n,1:n)=A-k^2*MM;


MM=sparse(MM);L=sparse(L);A=sparse(A);B=sparse(B);C=sparse(C);invL=inv(L);invL=sparse(invL);

********************************************************************************************************************************************


********************************************************************************************************************************************
��ϣ����ϸ��Ȿ���㷨�������Ǹ��� Understanding and Implementing the Finite Element Method[Mark_S._Gockenbach] ����д��possion���� �ı�׼����Ԫ��ɢ��
Ȼ���ڴ˻����ϣ�����  ����ų�����Ԫ������-����  ���е���ͽ�edgeԪ�͸նȾ������ɼ��ɼ�����ؽ����������Ҳ���մ�˳�򼴿�������
********************************************************************************************************************************************
Writed by ������ �������ڴ˴���������˵������
�κ�������Ľ���ѯ���뷢�ʼ���  hydzhang@hotmail.com ��
2015.10.

�޸ģ�2015.12  ����err1�ļ���
�޸ģ�2016.3   
�޸ģ�2017.2 --reconsidered after major revision 

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
#define MAX_NODES 8800
#define FREE 0
#define OK 1  
#define ERROR 0  
#define TRUE 1  
#define FALSE 0  
int sign(int a){if (a>0)return 1;else return -1;}

/* 
#define interior 15526
#define num_free 5087
 */
/*=========================================================================*/
#include<errno.h>
double doubleelem=0.0;
#define MAXSIZE 42500*4  
  //MAXSIZE��Ҫ���ڻ���΢���ڱߵĸ���
typedef int Status;  
typedef double ElemType;  
typedef struct{//��Ԫ��ṹ  
    int i,j;//����Ԫ�����±�����±�  
    ElemType e;//����Ԫ��ֵ  
}Triple;  
typedef struct{  
    Triple *data;//��Ԫ���data[0]����  
	
    int tu;//����������������ͷ���Ԫ�ظ���  
}TSMatrix;  
  
  TSMatrix NewMatrix(int m,int n);    TSMatrix NewMatrixhalf(int m,int n);  
  TSMatrix NewMatrixhalfhalf(int m,int n);  

    //�½�һ����Ԫ���ʾ��ϡ�����  
Status InsertElem(TSMatrix *M,int row,int col,double e);  
    //����Ԫ���ʾ��ϡ�����M���� row �У��� col ��λ�ò���Ԫ��e  
    //����ɹ�������OK�����򷵻�ERROR  



 char name[40];

 char nname[40];
 char sname[40];
 char ename[40];



 

 
 
 
 
 
 
 


//��������������ȡ�ļ��е��ַ�����������string, writed by  Bojan NICENO.
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
     printf(  "\n  ������Ϊ���� Preconditioners for the discretized time-harmonic .");   
     printf(  "\n   Maxwell  equations in mixed form  �Ļ������Ԫʵ��");

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

	  
	  


  printf("������ɺ��������������� �´򿪳��򣬸�ʽ���£�\n1) �����������֮�������� .n .s .eΪ��׺���ļ����ļ�����ȷ���ļ��ڵ�ǰ����������Ŀ¼�£��� ���ÿո����;\n\n2) ��������ļ�ͬ�����벻����׺����һ�θ����Ƽ��ɡ�\n\n��û���ṩEasyMesh���ɵ������ļ�������Ĭ��Ѱ������ΪG1.n G1.s G1.e���ļ� \n\n\n\n");
  system("pause");
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
  {printf("������ɺ��������������� �´򿪳��򣬸�ʽ���£�\n1) �����������֮�������� .n .s .eΪ��׺���ļ����ļ�����ȷ���ļ��ڵ�ǰ����������Ŀ¼�£��� ���ÿո����;\n\n2) ��������ļ�ͬ�����벻����׺����һ�θ����Ƽ��ɡ�\n\n��û���ṩEasyMesh���ɵ������ļ�.\n");

	printf("\nȱ��.n�ļ���������������:");
		scanf("%s",nname);
			printf("ȱ��.s�ļ���������������:");
		scanf("%s",sname);
			printf("ȱ��.e�ļ���������������:");
		scanf("%s",ename);

 }
 
 else
  {
    strcpy(nname,     argv[1]);
	  strcpy(sname,     argv[2]);
	    strcpy(ename,     argv[3]);
  }
   
//double pp=0.0;//ԭ����p���ϵ����0��1
   



char dummy1[80];
FILE *nodes_file=NULL;FILE *sides_file=NULL,*elements_file=NULL;
int num_free=0,num_con=0;
  
  
 //����ת�������ʽ��MFEMʹ�� 
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
	
//����Ϊ.n�ļ�����ȡ.*******************************************************************************************************************************************************************************  

nodes_file=fopen(nname,"r");

 if(!nodes_file)                                                                       
  {   fprintf(stderr, "%s \n", strerror(errno));
   printf("��ȷ�����������ļ�����Ŀ¼�´������У�\nCannot open  nodes file  %s\n",nname);                                                 
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
       load_i(nodes_file, &nodes[i].ptr);//��¼�ڵ�����;

	
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
        
    } //������ɱ���������ֵ
	 





	   
for(i=0;i<num_nodes;i++)
{
	if(nodes[i].ptr==FREE)
	      {nodes[i].ptr=num_freedummy+1;
          nodeptrs.fnodeptrs[num_freedummy]=i;
		  fprintf(outXn, "%16.12f  %16.12f\n ", nodes[i].x,nodes[i].y);   //������ɱ���������ֵ
          num_freedummy++;}
	  
	   
	      else  {nodes[i].ptr=-(i-num_freedummy+1);
	            nodeptrs.cnodeptrs[i-num_freedummy]=i;//��¼�ڵ�ָ������ɽڵ�,ǿ�ƽڵ�;��ֵ��ʾ����Ϊ����λ��ǿ�ƽڵ�
		        /* g[i-num_freedummy]=sin(PI*nodes[i].x)*sin(PI*nodes[i].y)+nodes[i].x; */}
}    //�˴�ΪDirichlet����g���ʽ�������.nodes[i].x node[i].yΪ����ֵ.
	   

fclose(outXn);
	
		

			
		
		
fclose(nodes_file)	;
		
		
	





	fprintf(Meshtrans, "boundary");fputc('\n',Meshtrans );
//����Ϊ.s�ļ�����ȡ.*******************************************************************************************************************************************************************************  

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
}sides[num_sides];//p,qָʾ�ߵ����ߵ�elements.

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
       load_i(sides_file, &sides[i].y);//��¼�ڵ�ָ��;
	   
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

//����Ϊ.e�ļ�����ȡ.*******************************************************************************************************************************************************************************  	  

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
	  
//�����ð���ʱ���໥�νӵ�˳�����߼�¼��elem��,������������ߵı����ķ����෴	  
	  
if(sides[dummya].y!=sides[dummyb].x&&sides[dummya].x!=sides[dummyb].x)//if  dummya  intersects with dummyb at dummyb.y
    {ss=(nodes[sides[dummya].x].x-nodes[sides[dummyb].x].x)*(nodes[sides[dummya].y].y-nodes[sides[dummyb].x].y)-(nodes[sides[dummya].x].y-nodes[sides[dummyb].x].y)*(nodes[sides[dummya].y].x-nodes[sides[dummyb].x].x);
    if(ss>0)//�жϵ������������
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
    if(ss>0)//�жϵ������������
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



	
	 //��Ԫ�ϵı�,������������.

	  
	  load_s(elements_file, dummy1); 
	  load_s(elements_file, dummy1); 
	  load_i(elements_file, &dummyd); 
 
  }  
  fclose(Meshtrans);

  fclose(elements_file);
  printf("���е�:%d\n��:%d\nԪ��:%d\n���ɱ߽��:%d\n���ɽڵ�m=%d\n�ڲ���n=%d\nn+m=%d\n",num_nodes,num_sides,num_elements,num_fbe,num_free,interior,num_free+interior);   system("pause");

   
//������ȡ���,����Ϊ�㷨��һ����:�նȾ���ļ���

//����curl-curlϡ����󣬵�һ��Ϊ�б� �ڶ���Ϊ�б꣬������Ϊֵ,�����򣬲������ظ���
//ע�⣺A	L ML MM ����Ķ���ֻ�������ǲ���
//��Pcurl, B���������������
	



TSMatrix A=NewMatrix(interior,interior);


TSMatrix MM=NewMatrix(interior,interior);

TSMatrix B=NewMatrixhalf(num_free,interior);

TSMatrix  L=NewMatrixhalfhalf(num_free,num_free);
TSMatrix  ML=NewMatrixhalfhalf(num_free,num_free);





double* FinF=NULL;
	FinF = (double *)malloc(sizeof(double)*(interior+num_free)); memset(FinF,0,(interior+num_free)*sizeof(double)); 
//double FinA[interior+num_free][interior+num_free],FinF[interior+num_free];memset(FinA,0,sizeof(double)*(interior+num_free)*(interior+num_free));memset(FinF,0,sizeof(double)*(interior+num_free));//������վ���,ֻ�����ڷ���matlab��֤�㷨.ʵ����Ҫ����ϡ��洢.




//memset(A,0,sizeof(double)*interior*interior);memset(MM,0,sizeof(double)*interior*interior);memset(B,0,sizeof(double)*num_free*interior);
  
int eptr[3]={0},indices[3]={0},ptrs[3]={0};
double delta,coords[3][2],I,x0,y0,invM[2][3];double det,G[3][3];
double kk=0;int s,r,ceshi=-1;double err1=0;
















int k;











double localength[3];
double aaa[3];//��¼m����������,Ҳ��������Ԫ���ĳ�����
double fff[3][3];double emm[3][3];
memset(emm,0,sizeof(emm));



double gg[interior]; memset(gg,0,sizeof(double)*interior);	






for(k=0;k<num_elements;k++)  //�����ǵ�Ԫѭ��
//����Ϊalgorithm 7.1 *************************
{	


eptr[0]=elem[k].x;//������,eptrָʾ��ʵ����������ʱ������}
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
             {ptrs[i]=nodes[indices[i]].ptr;}//�˴�ptrs�̳�ptr������ û�и���.

for(i=0;i<3;i++)
             {coords[i][0]=nodes[indices[i]].x;
		 }//������}



for(i=0;i<3;i++)
{coords[i][1]=nodes[indices[i]].y;
//printf("\n  %d \n ",indices[i]);system("pause");
}	
	


//*****************************������ļ���

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
//ϵ��k�������----------------------------------------------------------------------------------------------------------------------------------------------------------
I=delta;
   
  
  for(r=0;r<3;r++)
  {for(s=r;s<3;s++)	 
	  
	  {if(ptrs[r]>0&&ptrs[s]>0)
		  {		  i=min(ptrs[r],ptrs[s])-1;j=max(ptrs[r],ptrs[s])-1;
		        doubleelem=G[r][s]*I; if ( !InsertElem( &L,i,j, doubleelem) )
  { printf("\nError:L����ʧ��\n");  
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
  { printf("\nError:ML����ʧ��\n");  
        return ERROR;  
    }  
				  
	//	printf("%f", ML[i][j]);system("pause");
          }   
		  }
  }
  
  

  
  //load����ĳ�ʼ���ͼ���.*******************************************************************************************************************************************************************************  	  



//����Ϊalgorithm 7.5�ĵ�һ����************************************
/* for(s=0;s<3;s++)

      {if(localength[s]>=0)
		  {i=sides[abs(eptr[s])-1].ptr-1;
	       
	          gg[i]+= delta*sign(eptr[s])*(( 2-(kk*kk+2*x0)*(1-y0*y0) )  *(aaa[s]*invM[0][(1+s)%3]-aaa[(1+s)%3]*invM[0][s]+y0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s]))	+
               ( 2-(kk*kk+2*y0)*(1-x0*x0) )             * (aaa[s]*invM[1][(1+s)%3]-aaa[(1+s)%3]*invM[1][s]-x0*(invM[1][s]*invM[0][(1+s)%3]-invM[1][(1+s)%3]*invM[0][s])	)   );
           }
     } */
  
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
  { printf("\nError:A����ʧ��\n");  
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
  { printf("\nError:����ʧ��\n");  
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
  { printf("\nError:����ʧ��\n");  
        return ERROR;  
    }   
	      }
	  }
  }
  


 // err1+=(sqrt(pow(1-y0*y0-umid[0],2)+pow(1-x0*x0-umid[1],2)))/sqrt(pow(1-y0*y0,2)+pow(1-x0*x0,2))/(num_elements); // printf(";%f,%f;%f,%f;%f\n",umid[0],umid[1],1-y0*y0,1-x0*x0,(fabs(1-y0*y0-umid[0])+fabs(1-x0*x0-umid[1]))/(1-y0*y0+1-x0*x0));system("pause");


//���ǻ�Ԫ�ԸնȾ���ѭ�����,����ͬʱҪ��load����ѭ��,���Բ��Ӵ����� 
 

  

  
  
/* 	printf("k=%d   %f   gg=%f",k,invM[0][(1+2)%3],gg[k]);system("pause");
 */


}    //�����ǻ�Ԫѭ�����.

//���¼���C




	

 
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

 for(i=0;i<interior;i++)//��ѭ��
 {
	//ptr>0ָʾ���ɽڵ�
if(nodes[sides[interiorsides[i]].x].ptr>0)
{
 if ( !InsertElem( &C,i,nodes[sides[interiorsides[i]].x].ptr-1, -1.0) )
  { printf("\nError:����ʧ��\n");  
        return ERROR;  
    }   
}

if(nodes[sides[interiorsides[i]].y].ptr>0)
{
 if ( !InsertElem( &C,i,nodes[sides[interiorsides[i]].y].ptr-1, 1.0) )
  { printf("\nError:����ʧ��\n");  
        return ERROR;  
    }   
 }
}

//��ѭ��
 for(i=0;i<interior;i++)
 {
if(nodes[sides[interiorsides[i]].x].ptr>0)
{doubleelem=0.5*(nodes[sides[interiorsides[i]].y].x-nodes[sides[interiorsides[i]].x].x);
 if ( !InsertElem( &Pcurl,i,nodes[sides[interiorsides[i]].x].ptr-1, doubleelem) )
  { printf("\nError:����ʧ��\n");  
        return ERROR;  
    }   
doubleelem=0.5*(nodes[sides[interiorsides[i]].y].y-nodes[sides[interiorsides[i]].x].y);

 if ( !InsertElem( &Pcurl,i,nodes[sides[interiorsides[i]].x].ptr-1+num_free, doubleelem) )
  { printf("\nError:����ʧ��\n");  
        return ERROR;  
    }   
}

if(nodes[sides[interiorsides[i]].y].ptr>0)
{
	

doubleelem=0.5*(nodes[sides[interiorsides[i]].y].x-nodes[sides[interiorsides[i]].x].x);
 if ( !InsertElem( &Pcurl,i,nodes[sides[interiorsides[i]].y].ptr-1, doubleelem) )
  { printf("\nError:����ʧ��\n");  
        return ERROR;  
    }   
doubleelem=0.5*(nodes[sides[interiorsides[i]].y].y-nodes[sides[interiorsides[i]].x].y);

 if ( !InsertElem( &Pcurl,i,nodes[sides[interiorsides[i]].x].ptr-1+num_free, doubleelem) )
  { printf("\nError:����ʧ��\n");  
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
    //�½�һ����Ԫ���ʾ��ϡ�����  
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
    //�½�һ����Ԫ���ʾ��ϡ�����  
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
    //�½�һ����Ԫ���ʾ��ϡ�����  
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
    //����Ԫ���ʾ��ϡ�����M���� row+1 �У��� col+1 ��λ�ò���Ԫ��e  
    //����ɹ�������OK�����򷵻�ERROR  
    int z,t;  
    if(M->tu>=MAXSIZE){//��ǰ��Ԫ�������  
        printf("\nError:There is no space in the matrix;\n");  
        return ERROR;  
    }  
 // if(M->mu<=row||M->nu<=col){//��ǰ��Ԫ�������  
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
  
%��һ��Ϊ�б� �ڶ���Ϊ�б꣬������Ϊֵ,�����򣬲������ظ���
%ע�⣺A	L ML MM ����Ķ���ֻ�������ǲ���
%��Pcurl, B���������������
	
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
 
save L1
  
*/