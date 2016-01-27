#ifndef DATASTRUCTURE_H
#define DATASTRUCTURE_H
#include <cstdio>
#include<cstring>
#include <iostream>

using namespace std;


struct LinkNode												//node
{
	//data
	int neighbour;
	double edgeweight; 
	LinkNode *link;											
	
	LinkNode(int neighbour=0,double edgeweight=0,LinkNode *ptr=NULL) 
	{
		this->neighbour=neighbour;
		this->edgeweight=edgeweight;
		link=ptr;
	}				
};

class List
{
public:
	List(){first=NULL;len=0;}																							
	List(int neighbour,double edgeweight){first=new LinkNode(neighbour,edgeweight);}	
	~List(){MakeEmpty();}													
	void MakeEmpty();														
	bool IsEmpty() const {return first==NULL?true:false;}					
	void Output();															//output
	int Searchpos(int neighbour);											//search the position of neighbour
	int GetNeighbour(int pos);												
	void AddNode(int neighbour,double edgeweight);
public:																		
	LinkNode *first;														
	int len;																	

};

void List::MakeEmpty()										
{
	LinkNode *q,*p;
	p=first;
	for(int i=0;i<len;i++)
	{
		q=p->link;
		delete p;
		p=q;
	}
}

void List::Output()
{
	LinkNode *current=first;
	while (current!=NULL)
	{
		cout << current->neighbour << " " <<current->edgeweight << endl;
		current=current->link;
	}
}

void List::AddNode(int neighbour,double edgeweight)
{
	if(len==0) first=new LinkNode(neighbour,edgeweight);
	else
	{
		LinkNode *current=first;
		while(current->link!=NULL) current=current->link;
		current->link=new LinkNode(neighbour,edgeweight);
	}
	len++;
}

int List::Searchpos(int neighbour)
{
	LinkNode *current=first;
	int pos=-1;
	while(current!=NULL)
	{
		pos++;
		if(current->neighbour==neighbour) return pos;
		current=current->link;
	}
	return -1;
}

int List::GetNeighbour(int pos)
{
	if(pos>=len) return -1;
	int count=0;
	LinkNode *current=first;
	while(true)
	{
		if(count==pos) return current->neighbour;
		count++;
		current=current->link;
	}
}


//g¾ØÕó
struct NeighbourMatrix
{
	int h,w;
	double **matrix;
	NeighbourMatrix(int h=0,int w=0)
	{
		this->h=h;
		this->w=w;
	}
	~NeighbourMatrix()
	{
		if(matrix==NULL);
		else
		{
			for(int i=0;i<h;i++) delete [] matrix[i];
			delete [] matrix;

		}

	}
	void BuildMatrix(int h,int w)
	{
		this->h=h;
		this->w=w;
		matrix=new double *[h];
		for(int i=0;i<h;i++)
		{
			matrix[i]=new double[w];
			memset(matrix[i],0,w*sizeof(double));//intitialization
		}
	}
};


#endif
