
#include<iostream>
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include "DataStructure.h"

//as for the test data:graph-16
#define N 234
#define M 442
#define T 10

using namespace std;

void CalculateForg_function(NeighbourMatrix *g,double s[N][T],List ListIn[N],List ListOut[N]);
void swap(double &a,double &b);
double max(double a,double b);
double min(double a,double b);
void MatrixOutput(int **a,int h,int w);
void MatrixOutput(double **a,int h,int w);


int main()
{
	FILE *EDGE=fopen("edge.txt","r");
	FILE *DIS=fopen("distribution.txt","r");
	List ListIn[N],ListOut[N];//sparse storage
	//diffcount represents the number of nodes that have different scores
	int temp[3],i,j,k,times=0,notcount=0,diffcount=0;
	double s[N][T]={0},damper=0.5;
	NeighbourMatrix  g[N],b[N],r[N],a[N];
	int yold[N][T]={0};
	for(i=0;i<M;i++)
	{
		fscanf(EDGE,"%d%d%d",&temp[0],&temp[1],&temp[2]);//temp[0]->temp[1]
		ListIn[temp[1]].AddNode(temp[0],temp[2]);
		ListOut[temp[0]].AddNode(temp[1],temp[2]);
	}
	for(i=0;i<N;i++)
		for(j=0;j<T;j++) 
			fscanf(DIS,"%lf",&s[i][j]);

	fclose(EDGE);
	fclose(DIS);

	//1.1 calculate g(vi,yi,z), tentatively no correlation between topics is considered here
	CalculateForg_function(g,s,ListIn,ListOut);
	/*
	MatrixOutput(g[2].matrix,ListOut[2].len+1,T);
    */
	//1.2 Eq8, calculate bz,ij
	for(i=0;i<N;i++)
	{
		b[i].BuildMatrix(ListOut[i].len+1,T);
		r[i].BuildMatrix(ListOut[i].len+1,T);
		a[i].BuildMatrix(ListOut[i].len+1,T);//1.3 initialize a and r
		double sum[T]={0};
		for(j=0;j<=ListOut[i].len;j++)
			for(k=0;k<T;k++) 
				sum[k]=sum[k]+g[i].matrix[j][k];
		for(j=0;j<=ListOut[i].len;j++)
			for(k=0;k<T;k++) 
				b[i].matrix[j][k]=log(g[i].matrix[j][k]/sum[k]);
	}
	//MatrixOutput(b[2].matrix,ListOut[2].len+1,T);

	//the main loop in AP
	while(times<=500)
	{
		times++;
		//接下来的工作和上一次的比较像，区别是有多个topic和并不是每个点之间都有关联
		
		//step1: update r according to Eq 5
			
		for(i=0;i<N;i++)
		{
			//handle special cases, e.g., those nodes with no neighborhoods
			double firstmax[T],secondmax[T],temp=0,maxk[T];
			if(ListOut[i].len<1)
				for(k=0;k<T;k++)
					r[i].matrix[0][k]=b[i].matrix[0][k];
			else
			{
				for(k=0;k<T;k++)
				{
					firstmax[k]=b[i].matrix[0][k]+a[i].matrix[0][k];
					secondmax[k]=b[i].matrix[1][k]+a[i].matrix[1][k];
					maxk[k]=0;
					if(secondmax[k]>firstmax[k]) {swap(secondmax[k],firstmax[k]);maxk[k]=1;}
				}	
				for(j=2;j<=ListOut[i].len;j++)
				{
					for(k=0;k<T;k++)
					{
						if((temp=a[i].matrix[j][k]+b[i].matrix[j][k])>secondmax[k]) swap(secondmax[k],temp);
						if(secondmax[k]>firstmax[k]) swap(firstmax[k],secondmax[k]),maxk[k]=j;
					}
				}
				for(j=0;j<=ListOut[i].len;j++)
				{
					for(k=0;k<T;k++)
					{
						if(j==maxk[k]) r[i].matrix[j][k]=(b[i].matrix[j][k]-secondmax[k])*(1-damper)+r[i].matrix[j][k]*damper;
						else r[i].matrix[j][k]=(b[i].matrix[j][k]-firstmax[k])*(1-damper)+r[i].matrix[j][k]*damper;
					}
				}
			}
		}
		
	/*	MatrixOutput(r[0].matrix,ListOut[0].len+1,T);
		MatrixOutput(r[1].matrix,ListOut[1].len+1,T);
		MatrixOutput(r[2].matrix,ListOut[2].len+1,T);*/

		//second step: update matrix a. First update a[j].matrix[j][1...k], here we do not consider the summation.
		
		
		double firstmax[N][T]={0},secondmax[N][T]={0},temp;
		int maxk[N][T]={0};//calculate max min{r z,kj,0}
		memset(maxk,-1,N*T*sizeof(int));//initialization
		for(j=0;j<N;j++)//enumarte all inlinks
		{
			//where maxk[N] records the maximum value of min{r z,kj,0}
			if(ListIn[j].len==0)
			{
				for(k=0;k<T;k++) firstmax[j][k]=0;
				continue;
			}
		


			int neighbour=ListIn[j].GetNeighbour(0);
			int pos=ListOut[neighbour].Searchpos(j); //according to outlinks
//			ListOut[neighbour].Output();
			for(k=0;k<T;k++)
			{
				firstmax[j][k]=min(r[neighbour].matrix[pos][k],0);
				maxk[j][k]=neighbour;
			}

			if(ListIn[j].len>=2)
			{
				neighbour=ListIn[j].GetNeighbour(1);
				pos=ListOut[neighbour].Searchpos(j);
				for(k=0;k<T;k++)
				{
					secondmax[j][k]=min(r[neighbour].matrix[pos][k],0);
					if(secondmax[j][k]>firstmax[j][k])
					{
						swap(firstmax[j][k],secondmax[j][k]);
						maxk[j][k]=neighbour;
					}
				}

				for(i=2;i<ListIn[j].len;i++)
				{
					neighbour=ListIn[j].GetNeighbour(i);
					pos=ListOut[neighbour].Searchpos(j);
					for(k=0;k<T;k++)
					{
						temp=min(r[neighbour].matrix[pos][k],0);
						if(temp>secondmax[j][k]) swap(secondmax[j][k],temp);
						if(secondmax[j][k]>firstmax[j][k])
						{
							swap(firstmax[j][k],secondmax[j][k]);
							maxk[j][k]=neighbour;
						}
					}
				}
			}	
		}

		for(i=0;i<N;i++)
		{
			//first update aii
			for(k=0;k<T;k++) a[i].matrix[ListOut[i].len][k]=firstmax[i][k];
			//then update aij
			for(j=0;j<ListOut[i].len;j++)  
			{
				int neighbour=ListOut[i].GetNeighbour(j);//get all neighbourhood nodes' IDs
				for(k=0;k<T;k++)
				{
					if(i==maxk[neighbour][k]) 
						a[i].matrix[j][k]=(min(max(r[neighbour].matrix[ListOut[neighbour].len][k],0),-min(r[neighbour].matrix[ListOut[neighbour].len][k],0)-secondmax[i][k]))*(1-damper)+a[i].matrix[j][k]*damper;
					else 
						a[i].matrix[j][k]=(min(max(r[neighbour].matrix[ListOut[neighbour].len][k],0),-min(r[neighbour].matrix[ListOut[neighbour].len][k],0)-firstmax[i][k]))*(1-damper)+a[i].matrix[j][k]*damper;
				}
			}
		}


		/*MatrixOutput(a[0].matrix,ListOut[0].len+1,T);
		MatrixOutput(r[0].matrix,ListOut[0].len+1,T);
		ListOut[0].Output();
		MatrixOutput(r[64].matrix,ListOut[64].len+1,T);
		ListOut[64].Output();
		MatrixOutput(r[73].matrix,ListOut[73].len+1,T);
		ListOut[73].Output();
		MatrixOutput(r[232].matrix,ListOut[232].len+1,T);
		ListOut[232].Output();*/


		//third step:check for convergence
		if(times==20)//
		{
			double temp;
			for(i=0;i<N;i++)
			{
				for(k=0;k<T;k++)
				{
					double firstmax=r[i].matrix[ListOut[i].len][k]+a[i].matrix[ListOut[i].len][k];
					int rep=-1;//represents that the nodes itself is a representative node
					for(j=0;j<ListOut[i].len;j++)
						if((temp=r[i].matrix[j][k]+a[i].matrix[j][k])>firstmax) 
						{
							swap(temp,firstmax);
							rep=j;
						}
					if(rep==-1) rep=i;
					else rep=ListOut[i].GetNeighbour(rep);
					yold[i][k]=rep;
				}
			}
		}
		diffcount=0;
		if(times>=21)
		{
			double temp;
			for(i=0;i<N;i++)
			{
				for(k=0;k<T;k++)
				{
					double firstmax=r[i].matrix[ListOut[i].len][k]+a[i].matrix[ListOut[i].len][k];
					int rep=-1;//represents that the nodes itself is a representative node
					for(j=0;j<ListOut[i].len;j++)
					{
						temp=r[i].matrix[j][k]+a[i].matrix[j][k];
						if(temp>firstmax) 
						{
							rep=j;
							swap(temp,firstmax);
						}
					}
					if(rep==-1) rep=i;
					else rep=ListOut[i].GetNeighbour(rep);
					if(yold[i][k]!=rep) diffcount++;
					yold[i][k]=rep;
				}
			}
			if(diffcount==0) notcount++;
			else notcount=0;
			cout<<"times:"<<times<<" difference:"<<diffcount<<endl;
			if(notcount==100) break;

		}
	}
	
	//output the representative node (the node having the highest influence on the current) of each node on each topic
	FILE *out=fopen("Result.txt","w");
	fprintf(out,"Times:%d\n",times);
	for(i=0;i<N;i++)
	{
		fprintf(out,"%d:\t",i+1);
		for(k=0;k<T;k++)
			fprintf(out,"%d\t",yold[i][k]);
		fprintf(out,"\n");
	}
	fclose(out);
	system("pause");
	return 0;
}

void CalculateForg_function(NeighbourMatrix *g,double s[N][T],List ListIn[N],List ListOut[N])
{
	/*
		//debug
		cout<<"In:";
		ListIn[2].Output();
		cout<<endl<<"Out:";
		ListOut[2].Output();
	*/
	int i,j,k;
	LinkNode *current;
	//according to (1), considering directed graph, only when vi's target nodes or itself can be representative 
	for(i=0;i<N;i++)
	{
		g[i].BuildMatrix(ListOut[i].len+1,T);
		//calcualte denominator in eq. (1)
		
		//the sum of T topics's indegree W*i, 
		//the sum of T topics's outdegree W*i, 
		double Sumin[T]={0},Sumout[T]={0};
		//the sum of weighted outdegree
		for(current=ListOut[i].first;current!=NULL;current=current->link)
			for(k=0;k<T;k++)
				Sumout[k]=Sumout[k]+current->edgeweight*s[current->neighbour][k];
		//the sum of weighted indegree
		for(current=ListIn[i].first;current!=NULL;current=current->link)
			for(k=0;k<T;k++)
			{
				Sumin[k]=Sumin[k]+current->edgeweight*s[i][k];
				//calculate y z,i=i
				g[i].matrix[ListOut[i].len][k]=Sumin[k]/(Sumin[k]+Sumout[k]);
			}
		//calculate function g
		for(current=ListOut[i].first,j=0;current!=NULL;current=current->link,j++)
			for(k=0;k<T;k++)
				g[i].matrix[j][k]=current->edgeweight*s[current->neighbour][k]/(Sumin[k]+Sumout[k]);//邻域的计算
		
	}
}

void swap(double &a,double &b)
{
	double t;
	t=a;
	a=b;
	b=t;  
}

double max(double a,double b)
{
	return a>b?a:b;
}

double min(double a,double b)
{
	return a>b?b:a;
}

void MatrixOutput(int **a,int h,int w)
{
	cout<<"Matrix Output:"<<endl;
	for(int i=0;i<h;i++)
	{
		for(int j=0;j<w;j++)
			cout<<a[i][j]<<" ";
		cout<<endl;
	}
}

void MatrixOutput(double **a,int h,int w)
{
	cout<<"Matrix Output:"<<endl;
	for(int i=0;i<h;i++)
	{
		for(int j=0;j<w;j++)
			cout<<a[i][j]<<" ";
		cout<<endl;
	}
}