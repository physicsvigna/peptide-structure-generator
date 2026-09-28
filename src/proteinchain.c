/*
 * Peptide Structure Generator from Ramachandran Angles
 *
 * Generates peptide-backbone Cartesian coordinates from prescribed
 * phi, psi, and omega dihedral angles using recursive homogeneous
 * transformation matrices.
 *
 * Geometric / algorithmic basis:
 *   H. Stamati, A. Shehu, and L. Kavraki,
 *   "Computing Forward Kinematics for Protein-like linear systems
 *   using Denavit-Hartenberg Local Frames," Department of Computer
 *   Science, 2007.
 *
 * See references.bib in the repository for the BibTeX entry.
 */

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

int main(int argc, char *argv[])
{
    FILE *fp;
    FILE *fq;

    if(argc != 2)
    {
        fprintf(stderr,
                "Usage: %s <input_angle_file>\n",
                argv[0]);
        return 1;
    }

    fp=fopen("output_structure.pdb","w");

    if(fp == NULL)
    {
        perror("Error opening output file");
        return 1;
    }

//fh=fopen("newangles","r");
//initializing all the variables from rotational matrix to other integer variables to use as counter in constructing the peptide chain 
//we also inialize the oxygen coordinate array matrix.
int i,j,p,k,q,pe,l;

const char s[]=", ";
char *token;

float sum,sm,dis,theta,alpha,coord[5][2];
float mul[5][5],oxcoord[5][5],rot[4][4];
float rotcoord[50][5][5][5],mat1,mat2,tempang[50][4];

sum=0;
sm=0;

//We declare the values of geometry of the peptide such as bond length (bndlen) and bond angles (bndang)
//We also intialize the ramachandran angles matrix and origin colomn vector for the local frames attached.

float bndlen[4]={1.46,1.51,1.24,1.33};
float bndang[4]={1.012,1.204,1.038,1.117};
float ramang[50][4];

      coord[0][0]=0.0;
      coord[1][0]=0.0;
      coord[2][0]=0.0;
      coord[3][0]=1.0;                 // {{0.0},{0.0},{0.0},{1.0}};

float rotrec[4][4]={{1.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0},{0.0,0.0,1.0,0.0},{0.0,0.0,0.0,1.0}};


printf("\n");

//for(p=0;p<4;p++)
//{
// for(q=0;q<4;q++)
//  {
//   printf("\t%.2f",rotrec[p][q]);
//  }
//   printf("\n");
//}

//we also inialize the oxygen coordinate in the local frame at the atom C of the given residue.

oxcoord[0][0]=0; //{{0},{1.24*(sin(1.038))},{1.24*(cos(1.038))},{1}};
oxcoord[1][0]=1.24*(sin(1.038));
oxcoord[2][0]=1.24*(cos(1.038));
oxcoord[3][0]=1;

//printf("\nEnter the number of peptides\n");
//scanf("%d",&pe);

//printf("\n Enter the ramachandran angle \n");
//fq=fopen("input_backbone_angles.txt","r");


fq=fopen(argv[1],"r");


k=0;
l=0;
 if(fq != NULL)
    {
	char line[200];
	while(fgets(line, sizeof line, fq) != NULL)
	{
		if(line[0] == '\n' || line[0] == '\r')
			{
				continue;
			}
	    token = strtok(line, s);
	    for(i=0;i<2;i++)
	    {
		if(i==0)
		{
		    //fprintf(fp,"\nfrom file: %s \t, K= %d",token,k);
			printf("\nfrom file: %s \t, K= %d",token,k);
		   // if(atoi(token)!=0){
			 ramang[k][l]=atof(token);
			 l++;
		    if(l==3)
			{
			  k++;
			  l=0;
			 }
		     //}
		    token = strtok(NULL,s);
		} else {
		  // printf("from else: %d\n",atoi(token));
		  }
	    }

	}
	fclose(fq);
    }
    else
      {
         perror("Error opening input angle file");
         fclose(fp);
         return 1;
      }
       pe=k;
       //fprintf(fp,"\nThe number of peptides is = %d",pe);
	   printf("\nThe number of residues is = %d",pe);
    //printf("\nK Value: %d",k);
    for(i=0;i<pe;i++){
      for(int j=0;j<3;j++){
       //fprintf(fp,"\nramang[%d][%d] = %f",i,j,ramang[i][j]);
	   printf("\nramang[%d][%d] = %f",i,j,ramang[i][j]);
      }
    }
//for(i=0;i<pe;i++)
//  {
//    for(j=0;j<3;j++)
//     {
//       scanf("%f",&ramang[i][j]);

//     }
//  }
for(i=0;i<pe;i++)
 {
  for(j=0;j<3;j++)
   {
    ramang[i][j]=180-ramang[i][j];
    tempang[i][j]= (3.142/180)*ramang[i][j];
    ramang[i][j]=tempang[i][j];
    }
 }
//for(i=0;i<pe;i++)
//  {
//    for(j=0;j<3;j++)
//     {
//       fprintf(fp,"\n%f",ramang[i][j]);
//     }
//  }



for(i=0;i<pe;i++)
  {
   for(j=0;j<4;j++)
     {
     if(j==1)
       {
	mat1=ramang[i][j];
       }

      else if(j==2)
       {
	mat2=ramang[i][j];
	ramang[i][j]=mat1;
       }

      else if(j==3)
       {
	ramang[i][j]=mat2;
       }
  }
}


//for(i=0;i<pe;i++)
//{
//  for(j=0;j<4;j++)
//  {
//     fprintf(fp,"\n Ramachandran angle[%d][%d]= %.2f",i,j,ramang[i][j]);
//    }
//    printf("\n\n");
//}


//for(i=0;i<4;i++)
//  {
//  printf("\n The bond length[%d] = %.2f",i,bndlen[i]);
//  }

// printf("\n");

//for(i=0;i<4;i++)
//{
//  printf("\n The bond angle[%d] =%.2f",i,bndang[i]);
//}

for(i=0;i<pe;i++)
{
 for(j=0;j<4;j++)
 {

  dis=bndlen[j];
  alpha=bndang[j];
  theta=ramang[i][j];

//  float ct=cos(theta);
//  float st=sin(theta);
//  float ca=cos(alpha);
//  float sa=sin(alpha);
//  float rot[4][4]={{ct,st,0,0},{(st*ca),(ct*ca),(-sa),(-d*sa)},{(st*sa),(ct*sa),ca,(d*ca)},{0,0,0,1}};


if(j!=2)
{

 if(i%2==0)
  {
    if(i==0)
     {
       if(j==0)
	 {
	  alpha=0;
//	  printf("\nThe alpha of [%d][%d] = %f",i,j,alpha);
	 }
       if(j==3)
	{
	 alpha=(-1)*bndang[j];
//	 printf("\nThe alpha of [%d][%d] = %f",i,j,alpha);
	}

      if(j==1)
	{
	  alpha=bndang[j];
//	  printf("\nThe alpha of [%d][%d] = %f",i,j,alpha);
	}
      }
    else
     {
       if(j==0||j==3)
	   {
	    alpha=(-1)*bndang[j];
//	    printf("\nThe alpha of [%d][%d] = %f",i,j,alpha);
	   }
	  else
	   {
	    alpha=bndang[j];
//	    printf("\nThe alpha of [%d][%d] = %f",i,j,alpha);
	   }
     }
 }


 else
  {
    if(j==1)
     {
       alpha=(-1)*bndang[j];
//       printf("\nThe alpha of [%d][%d] = %f",i,j,alpha);
     }
    else
    {
       alpha=bndang[j];
//       printf("\nThe alpha of [%d][%d] = %f",i,j,alpha);
    }
  }

rot[0][0]= (cos(theta));
rot[0][1]= (sin(theta));
rot[0][2]= 0;
rot[0][3]= 0;

rot[1][0]= (-sin(theta))*(cos(alpha));
rot[1][1]= (cos(theta))*(cos(alpha));
rot[1][2]= (sin(alpha));
rot[1][3]=  ((sin(alpha)))*dis;

rot[2][0]= (sin(theta))*(sin(alpha));
rot[2][1]= (cos(theta))*(-sin(alpha));
rot[2][2]= (cos(alpha));
rot[2][3]= (cos(alpha))*dis;

rot[3][0]= 0;
rot[3][1]= 0;
rot[3][2]= 0;
rot[3][3]= 1;

//printf("\n");
//for(p=0;p<4;p++)
// {
//  for(q=0;q<4;q++)
//   {
//    printf("\t%.2f",rot[p][q]);
//   }
//   printf("\n");
// }

for(p=0;p<4;p++)
 {
  for(q=0;q<4;q++)
   {
    for(k=0;k<4;k++)
      {
	sum=sum+(rotrec[p][k]*rot[k][q]);
      }
      mul[p][q]=sum;
      sum=0;
   }
 }




//printf("\n The matrix is multiplied\n");
//for(p=0;p<4;p++)
// {
//  for(q=0;q<4;q++)
//   {
//    printf("\t%.2f",mul[p][q]);
//   }
//   printf("\n");
// }





for(p=0;p<4;p++)
 {
  for(q=0;q<4;q++)
   {
    rotrec[p][q]=mul[p][q];
    }
 }




for(p=0;p<4;p++)
 {
  for(q=0;q<1;q++)
   {
    for(k=0;k<4;k++)
      {
	sm=sm+(rotrec[p][k]*coord[k][q]);
      }
      rotcoord[i][j][p][q]=sm;
//      printf("\t %.2f",rotcoord[i][j][p][q]);
      sm=0.0;
   }

 }

// printf("\n\n");
//   for(p=0;p<4;p++)
//     {
//      for(q=0;q<4;q++)
//       {
//	 printf("\t %.2f",rot[p][q]);
//       }
//       printf("\n");
//     }
//   }

 }

else
  {
 //   printf("\n");
//    for(p=0;p<4;p++)
//    {
//     for(q=0;q<4;q++)
//       {
//	 printf("\t %.2f",rotrec[p][q]);
//      }
//      printf("\n");
//    }

   if(i%2==0)
      {
       alpha=bndang[j];
//       printf("\nThe alpha of [%d][%d] = %f",i,j,alpha);
       oxcoord[0][0]=0; //{{0},{1.24*(sin(1.038))},{1.24*(cos(1.038))},{1}};
       oxcoord[1][0]=1.24*(sin(1.038));
       oxcoord[2][0]=1.24*(cos(1.038));
       oxcoord[3][0]=1;
      }
     else
      {
	     alpha=(-1)*bndang[j];
//       printf("\nThe alpha of [%d][%d] = %f",i,j,alpha);
       oxcoord[0][0]=0; //{{0},{1.24*(sin(1.038))},{1.24*(cos(1.038))},{1}};
       oxcoord[1][0]=-1.24*(sin(1.038));
       oxcoord[2][0]=1.24*(cos(1.038));
       oxcoord[3][0]=1;
      }


for(p=0;p<4;p++)
 {
  for(q=0;q<1;q++)
   {
    for(k=0;k<4;k++)
      {
	sum=sum+(rotrec[p][k]*oxcoord[k][q]);
      }
      rotcoord[i][j][p][q]=sum;
//      printf("\t%.2f",rotcoord[i][j][p][q]);
      sum=0.0;
   }
 }

}
// printf("\n");

}

// else
//    {
//    for(p=0;p<4;p++)
//     {
//      for(q=0;q<4;q++)
//       {
//	 printf("\t %f",rot[p][q]);
//       }
//       printf("\n");
//     }
//  }
}
//printf("\n The set of coordinates in the order of x,y,z,t respectively in a row is\n");
//for(i=0;i<4;i++)
// {
//  for(j=0;j<1;j++)
//   {
//    printf("\n %f",coord[i][j]);
//    }
//  }



//fprintf(fp,"\n");
//fprintf(fp,"         1         2         3         4         5         6         7         8");
//fprintf(fp,"\n12345678901234567890123456789012345678901234567890123456789012345678901234567890\n");

int u=1;

char ele1[3]="CA";
char ele2[3]="C";
char ele3[3]="O";
char ele4[3]="N";
char ami[6]="ALA A";
//fprintf(fp,"ATOM%6d%5s%7s%4d",u,ele4,ami,u);

fprintf(fp,
        "ATOM%6d%5s%7s%4d%12.3f%8.3f%8.3f\n",
        u, ele4, ami, 1,
        0.0, 0.0, 0.0);

u++;

for(i=0;i<pe;i++)
 {
  for(j=0;j<4;j++)
   {
	   
	
   if(i == pe-1 && j == 3)
    {
       continue;
    }   
	   
   if(j==0)
   {
   fprintf(fp,"ATOM%6d%5s%7s%4d",u,ele1,ami,i+1);
   }
   if(j==1)
   {
   fprintf(fp,"ATOM%6d%5s%7s%4d",u,ele2,ami,i+1);
   }
   if(j==2)
   {
   fprintf(fp,"ATOM%6d%5s%7s%4d",u,ele3,ami,i+1);
   }
   if(j==3)
	{
    fprintf(fp,"ATOM%6d%5s%7s%4d",u,ele4,ami,i+2);
	}
    for(p=0;p<4;p++)
     {

      for(q=0;q<1;q++)
      {
      //fprintf(fp,"[%d][%d]",i,j);
       if(p==0)
       {
	fprintf(fp,"%12.3f",rotcoord[i][j][p][q]);
       }
       if(p==1)
       {
	fprintf(fp,"%8.3f",rotcoord[i][j][p][q]);
       }
       if(p==2)
	{
	fprintf(fp,"%8.3f",rotcoord[i][j][p][q]);
	}
      //   fprintf(fp,"\t");
     }

   }
   fprintf(fp,"\n") ;
   u++;
   }
 }
fprintf(fp,"TER\n");
fprintf(fp,"END\n"); 
fclose(fp);// close output file
return 0; // successful completion
}

