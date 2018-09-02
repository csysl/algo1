#pragma once
#include<stdio.h>
#include<time.h>       //time
#include<math.h>       //power
#include<algorithm>    //max
#include<string>       //string
#include<set>          //set
#include<vector>       //vector
#include<random>       //random_device;uniform_int_distribution<>
#include<fstream>      //ftream
#include<iostream>     //cin;cout
#include<iomanip>      //control cin cout
#include<glpk.h>               //gMIP;gLP;gstFlow
#include<ilcplex\ilocplex.h>   //cMIP;cLP;CstFlow

using namespace std;

double *rand_weight(int n, int m);    //set weight,range[1,10]
set<set<int>> randSS(int n, int k);   //set SS  |SS|=m |SS|<3*n;|SS|>q
int randq(set<set<int>> SS, int k);   //get q
set<set<int>> *randGG(set<set<int>> SS, int qq);  //get a partition of SS which named GG  |GG|=q; 0.15*m<=|GG|<=0.25*m;|GG|>k


void gMIP(int n, int k, double *weight, set<set<int>> SS, set<set<int>> *GG, int q);       //solve 0-1_LP problem by GLPK
void cMIP(int n, int k, double *weight, set<set<int>> SS, set<set<int>> *GG, int q);       //solve 0-1_LP problem by CPLEX
double *gLP(int n, int k, double *weight, set<set<int>> SS, set<set<int>> *GG, int q);     //solve LP problem by GLPK
double *cLP(int n, int k, double *weight, set<set<int>> SS, set<set<int>> *GG, int q);     //solve LP problem by CPLEX


struct Edge {
	int no;
	int u, v;
	double weight;
};
//������¼��ǰ�ڵ���ߣ�������ߣ�����һ���� �ߵ���� �ߵ�ֵ
struct V_Weight {
	int v;
	int no;
	double weight;
};
class DirectedGraph {
public:
	vector<V_Weight> *ulist;    //��¼ÿ���ڵ�ĳ��ߵ����ԣ��ڵ����Ϊ���±�
	vector<V_Weight> *vlist;    //��¼ÿ���ڵ����ߵ����ԣ��ڵ����Ϊ���±�
	vector<Edge> edge;          //��¼�߼�

	//constructor 
	DirectedGraph(vector<Edge> const&edges, int n, bool judge) {
		/*
		edges:   the set of edges
		n:       the num of points
		judge:   ��ȡ0,1��0:������¼��1:weightΪ0��ʱ�򲻼�¼
		ulist:   ��¼��u�����ıߣ������±�u����¼ÿһ��v���Ѿ���Ӧ�ߵ���ź�Ȩֵ
		vlist:   ��¼�����v�ıߣ������±�v����¼ÿһ��u���Ѿ���Ӧ�ߵ���ź�Ȩֵ
		*/

		//allocate memory
		ulist = new vector<V_Weight>[n];
		vlist = new vector<V_Weight>[n];
		V_Weight vw;

		//initialization
		edge = edges;
		if (judge == false) {
			
			//add edges to the Graph
			for (int i = 0; i < edges.size(); i++) {
				int u = edges[i].u;
				vw.v = edges[i].v;
				vw.no = edges[i].no;
				vw.weight = edges[i].weight;
				//for every u,add all its v
				ulist[u].push_back(vw);
			}
			cout << endl;
			for (int i = 0; i < edges.size(); i++) {
				int u = edges[i].v;
				vw.v = edges[i].u;
				vw.no = edges[i].no;
				vw.weight = edges[i].weight;
				//for every v,add all its u
				vlist[u].push_back(vw);
			}
		}
		else if(judge==true){

			//add edges to the Graph
			for (int i = 0; i < edges.size(); i++) {
				if (edges[i].weight == 0.0) {
					continue;
				}
				int u = edges[i].u;
				vw.v = edges[i].v;
				vw.no = edges[i].no;
				vw.weight = edges[i].weight;
				//for every u,add all its v
				ulist[u].push_back(vw);
			}
			for (int i = 0; i < edges.size(); i++) {
				if (edges[i].weight == 0.0) {
					continue;
				}
				int u = edges[i].v;
				vw.v = edges[i].u;
				vw.no = edges[i].no;
				vw.weight = edges[i].weight;
				//for every v,add all its u
				vlist[u].push_back(vw);
			}
		}
		
	}

	//constructor
	~DirectedGraph() {
		int u_size = ulist->size();
		int v_size = vlist->size();
		for (int i = 0; i <= u_size; i++) {
			ulist[i].clear();
		}
		for (int i = 0; i <= v_size; i++) {
			vlist[i].clear();
		}
	}
};
vector<Edge> Edgess(int q, int k);//get the set of edges
DirectedGraph Gcstflow(DirectedGraph G, vector<Edge> edges, double *Z, int q, int k, int n); //get an stFlow with value 1 by cplex
set<set<double>> GFf(DirectedGraph G, int q, int k);//get a collections of path-flows


//double *gstFlow(double *Z, int q, int k);//Compute an st-flow of value 1 in G by GLPK
//double *cstFlow(double *Z, int q, int k);//Compute an st-flow of value 1 in G by CPLEX
//set<set<double>> Ff(double *Zx, int q, int k);//get a collections of path-flows


double *ff(set<set<double>> FF, int k);//select a path from path-flows

set<set<double>> GGGFf(DirectedGraph G, int q, int k);
set<set<double>> GGFf(DirectedGraph G, int q, int k);