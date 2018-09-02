/*
|U|=n  |SS|=m  |GG|=q
|S|<=0.4*n
n<=m<=3*n
k<=q<=0.4*m
*/
#include"func.h"

int main() {
	system("chcp 437");
	int n, k;
	while (1) {
		cin >> n >> k;  //nΪ����U�ķ�Χ��kΪ���󼯺ϵĸ���,qΪ���ָ���

#pragma region SS_GG
	SS:
		//�����ݼ����Ӽ�SS
		set<set<int>> SS;
		SS = randSS(n, k);
		cout << "********************SS set!" << endl;

		//��SS�����setSS.txt
		ofstream coutSS("setSS.txt");
		set<set<int>>::iterator it1 = SS.begin();
		while (it1 != SS.end()) {
			set<int>::iterator it2 = (*it1).begin();
			while (it2 != (*it1).end()) {
				coutSS << *it2 << " ";
				it2++;
			}
			coutSS << endl;
			it1++;
		}
	GG:
		//����SS�Ļ���GG
		int q = randq(SS, k);
		//int q = k;
		set<set<int>> *GG = new set<set<int>>[q + 1];
		GG = randGG(SS, q);
		cout << "********************GG set!" << endl;
		//��GG�����setGG.txt
		ofstream coutGG("setGG.txt");
		for (int i = 1; i <= q; i++) {
			set<set<int>>::iterator it1 = GG[i].begin();
			while (it1 != GG[i].end()) {
				set<int>::iterator it2 = (*it1).begin();
				while (it2 != (*it1).end()) {
					coutGG << *it2 << " ";
					it2++;
				}
				coutGG << endl;
				it1++;
			}
			coutGG << endl;
		}
#pragma endregion
	weight:
		//����Ȩ��wi
		double *weight = new double[n + SS.size() + 1];
		weight = rand_weight(n, SS.size());
		cout << "********************weight set!\n" << endl;
	cLP:
		double *X = new double[SS.size() + 1];
		//ʹ��GLPK����LP
		//X = gLP(n, k, weight, SS, GG, q);
		//cout << "***************" << endl;
		//ʹ��cplex����LP
		X = cLP(n, k, weight, SS, GG, q);
		cout << "********************" << endl;
	Z:
		//get Z*
		double *Z = new double[q + 1];
		int tag_x = 1;
		for (int i = 1; i <= q; i++) {
			int s_g = GG[i].size();
			Z[i] = 0.0;
			for (int j = 1; j <= s_g; j++) {
				Z[i] += X[tag_x];
				//cout << "X" << tag_x << "  " << X[tag_x] << endl;  //�õ���
				++tag_x;
			}
			//cout << "        Z" << i << ": " << Z[i] << endl;   //�õ�Zl*
		}
	Graph:
		vector<Edge> edges = Edgess(q, k);
		int nn = 2 * k * q + 3;   //��ʱ����n�����������
		DirectedGraph G1(edges, nn, 0);
	GstFlow:
		//�������ͼ
		DirectedGraph G = Gcstflow(G1, edges, Z, q, k, nn);
		//���������txt
		ofstream gout("GstFlow.txt");
		for (int i = 0; i < G.edge.size(); i++) {
			if (G.edge[i].weight != 0.0) {
				gout << G.edge[i].u << " " << G.edge[i].v << " " << G.edge[i].weight << endl;
			}
		}
		gout << k << " " << q << " " << 0.0 << endl;
		//����������ն�
		for (int t = 0; t < k*q; t++) {
			if (G.edge[t].weight != 0.0) {
				cout << "Z" << G.edge[t].no << ":  " << G.edge[t].weight << endl;
			}
			if (t%q == q - 1) {
				cout << endl;
			}
		}
		cout << "********************" << endl;
	GFF:
		//�õ��ֽ�·��
		set<set<double>> FF;
		FF = GGFf(G, q, k);

		//��FF�����FF.txt
		ofstream coutFF("FF.txt");
		double tmpp = 0.0;
		set<set<double>>::iterator itt1 = FF.begin();
		while (itt1 != FF.end()) {
			set<double>::iterator itt2 = (*itt1).begin();
			tmpp += (*itt2);
			while (itt2 != (*itt1).end()) {
				coutFF << *itt2 << " ";
				itt2++;
			}
			coutFF << endl;
			itt1++;
		}
		cout << "*******************FF get!" << endl;
		cout << "the sum of FF: " << tmpp << endl;

	f:
		//select a path from path-flows
		double *f = new double[k + 2];
		f = ff(FF, k);
		//cout << "The ordinal number of the selected edge: ";
		//for (int i = 1; i <= k; i++) {
		//	cout << f[i] << "\t";
		//}
		//cout<<endl;

		cout << "********************f get!" << endl;
		ofstream coutf("f.txt");
		for (int t = 0; t < k; t++) {
			coutf << f[t + 1] << endl;
		}

	CC:
		//get result CC
		set<set<int>> *CC = new set<set<int>>[k + 1];
		for (int i = 1; i <= k; i++) {
			int tmp_cc = int(f[i]) % q;
			CC[i] = GG[tmp_cc];
		}
		//���CC��outCC.txt
		ofstream coutCC("outCC.txt");
		for (int i = 1; i <= k; i++) {
			auto it_c = CC[i].begin();
			while (it_c != CC[i].end()) {
				auto it_c2 = (*it_c).begin();
				while (it_c2 != (*it_c).end()) {
					coutCC << *it_c2 << " ";
					++it_c2;
				}
				coutCC << endl;
				++it_c;
			}
			coutCC << endl;
		}
		cout << "********************CC get!" << endl;
		cout << endl;

		//erase the memory
		
		for (int i = 1; i <= k; i++) {
			CC[i].clear();
		}
		SS.clear();
		for (int i = 1; i <= q; i++) {
			GG[i].clear();
		}
		delete[]weight, X, Z;
			//Ze;
		cout << "----------**************************************----------" << endl << endl;
	}

	system("pause");
	return 0;
}