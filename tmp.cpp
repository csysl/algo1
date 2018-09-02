#include"func.h"
//get a collections of path-flows
//��δ�����ʼ�Ĵ���
set<set<double>> GGGFf(DirectedGraph G, int q, int k) {
	set<set<double>> FF;
	random_device r;
	int tag_no;
	double tag_w;//��ŵ�ǰ�ҵ�����С�ߵ���ź�Ȩֵ
	while (1) {
		set<double> f;   //��һ��ֵ�洢·��f�ĸ��ʣ������k��ֵ�洢k����
		double minv = INFINITY;  //��¼��СȨֵ
		tag_no = -1;
		//�õ�Ȩֵ��С�ı�
		for (int i = 0; i < k*q; i++) {
			if (G.edge[i].weight < minv&&G.edge[i].weight != 0.0) {
				minv = G.edge[i].weight;
				tag_no = G.edge[i].no;
			}
		}
		cout << "no: " << tag_no << "  min: " << minv << endl;
		if (tag_no == -1)
			break;
		f.insert(minv);
		f.insert(tag_no + 1);

		int t = tag_no / q + 1;  //��ǰ��ǰ·�����ڵ���
		int u, v, tmp_tag_no = tag_no;
		tag_w = minv;
		//��ǰ��·��
		for (int i = t - 1; i >= 1; i--) {
			u = 2 * (tmp_tag_no + 1) - 1;  //��ǰ�ߵĳ�����
			int ss = G.vlist[u].size();
			if (ss == 0)//�����vlist��ɾ�����������Լ����жϣ�����õ�û����ߣ���ѳ�ʼ����Ϊ0 
				goto goodbye;
			uniform_int_distribution<> ui(0, ss - 1);
			int tmp = ui(r);
			v = G.vlist[u][tmp].v;  //����õ�����������һ��������

			G.edge[G.vlist[u][tmp].no].weight -= tag_w; //������
			//�������
			if (G.edge[G.vlist[u][tmp].no].weight <= 1.0e-7) {
				//if (G.edge[G.vlist[u][tmp].no].weight <=0.0) {
				G.edge[G.vlist[u][tmp].no].weight = 0.0;//�����ߵ�ȨֵС��1.0e-7ʱ����0
				//��������Ȩֵ����0ʱ����vlist��ɾ��
				auto iter1 = G.vlist[u].begin();
				while (iter1 != G.vlist[u].end()) {
					if ((*iter1).v == v) {
						G.vlist[u].erase(iter1);
						break;
					}
					++iter1;
				}
				//��������Ȩֵ����0ʱ����ulist��ɾ��
				auto iter2 = G.ulist[v].begin();
				while (iter2 != G.ulist[v].end()) {
					if ((*iter2).v == u) {
						G.ulist[v].erase(iter2);
						break;
					}
					++iter2;
				}
			}

			tmp_tag_no = v / 2 - 1; //�ҵ���ǰ�ߵ����
			G.edge[tmp_tag_no].weight -= tag_w;
			if (G.edge[tmp_tag_no].weight <= 1.0e-7)
				G.edge[tmp_tag_no].weight = 0.0;
			f.insert(tmp_tag_no + 1);
		}
		//������·��
		tmp_tag_no = tag_no;
		for (int i = t + 1; i <= k; i++) {
			u = 2 * (tmp_tag_no + 1);
			int ss = G.ulist[u].size();
			if (ss == 0)
				goto goodbye;
			uniform_int_distribution<> ui(0, ss - 1);
			int tmp = ui(r);
			v = G.ulist[u][tmp].v;

			G.edge[G.ulist[u][tmp].no].weight -= tag_w;
			//�������
			if (G.edge[G.ulist[u][tmp].no].weight <= 1.0e-7) {
				//if (G.edge[G.ulist[u][tmp].no].weight <=0.0) {
				G.edge[G.ulist[u][tmp].no].weight = 0.0;
				//��������Ȩֵ����0ʱ����ulist��ɾ��
				auto iter1 = G.ulist[u].begin();
				while (iter1 != G.ulist[u].end()) {
					if ((*iter1).v == v) {
						G.ulist[u].erase(iter1);
						break;
					}
					++iter1;
				}
				//��������Ȩֵ����0ʱ����vlist��ɾ��
				auto iter2 = G.vlist[v].begin();
				while (iter2 != G.vlist[v].end()) {
					if ((*iter2).v == u) {
						G.vlist[v].erase(iter2);
						break;
					}
					++iter2;
				}
			}

			tmp_tag_no = (v + 1) / 2 - 1;
			G.edge[tmp_tag_no].weight -= tag_w;
			if (G.edge[tmp_tag_no].weight <= 1.0e-7)
				G.edge[tmp_tag_no].weight = 0.0;
			f.insert(tmp_tag_no + 1);
		}
		FF.insert(f);
		f.clear();
	goodbye:
		G.edge[tag_no].weight = 0.0;
	}

	return FF;
}

//get a collections of path-flows
//��δ�������ʦ���ķ����㷨����������   ����������Сʱ����һ���������������ĸ�����
set<set<double>> GFf(DirectedGraph G, int q, int k) {
	set<set<double>> FF;
	random_device r;
	int tag_no;
	double tag_w;//��ŵ�ǰ�ҵ�����С�ߵ���ź�Ȩֵ
	vector<int> no_set;
	while (1) {

	goodbye:
		set<double> f;   //��һ��ֵ�洢·��f�ĸ��ʣ������k��ֵ�洢k����
		double minv = INFINITY;  //��¼��СȨֵ
		tag_no = -1;
		//�õ�Ȩֵ��С�ı�
		for (int i = 0; i < k*q; i++) {
			if (G.edge[i].weight < minv&&G.edge[i].weight != 0.0) {
				minv = G.edge[i].weight;
				tag_no = G.edge[i].no;
			}
		}
		cout << endl << "no: " << tag_no << "  min: " << minv;
		if (tag_no == -1)
			break;
		f.insert(minv);
		f.insert(tag_no + 1);

		int t = tag_no / q + 1;  //��ǰ��ǰ·�����ڵ���
		int u, v, tmp_tag_no = tag_no;
		tag_w = minv;
		int debug = 0;
		//��ǰ��·��
		for (int i = t - 1; i >= 1; i--) {
			/*
			1���õ���ǰ�ߣ��Ǹ����ߣ�
			2����鵱ǰ���Ƿ�����ߣ�û�У��Ѹñ���0��ͬʱ�ص����   ���⣺��ǰ��֮ǰ�ı��Ѿ�������ѡ�еĳ�ʼ�ߵ�Ȩֵ
			3������ߣ��������ѡһ��������ߣ���ȥ��ʼ�ߵ�Ȩ�أ���������С��1.0e-7����0�������ø����ߴ�vlist��ulist��ɾ��
			4��
			*/
			u = 2 * (tmp_tag_no + 1) - 1;  //��ǰ�ߵĳ�����
			int ss = G.vlist[u].size();

			//�����vlist��ɾ�����������Լ����жϣ�����õ�û����ߣ���ѳ�ʼ����Ϊ0
			if (ss == 0) {
				G.edge[tmp_tag_no].weight = 0.0;
				cout << "+++++++++++" << endl;
				goto goodbye;
			}
			uniform_int_distribution<> ui(0, ss - 1);
		go_tag1:
			int tmp = ui(r);
			auto vlis = G.vlist[u][tmp];
			v = vlis.v;  //����õ�����������һ��������
			double fuck;
			if (ss == 1)
				fuck = tag_w - (1.0e-10);
			if (G.edge[vlis.no].weight < fuck) {
				cout << "weight:" << G.edge[vlis.no].weight << endl;
				goto go_tag1;
			}
			G.edge[vlis.no].weight -= tag_w; //������

			//�������
			if (G.edge[vlis.no].weight <= 1.0e-10) {
				//if (G.edge[vlis.no].weight <=0.0) {
				G.edge[vlis.no].weight = 0.0;            //�����ߵ�ȨֵС��1.0e-10ʱ����0

				//��������Ȩֵ����0ʱ����vlist��ɾ��
				auto iter1 = G.vlist[u].begin();
				while (iter1 != G.vlist[u].end()) {
					if ((*iter1).v == v) {
						G.vlist[u].erase(iter1);
						break;
					}
					++iter1;
				}
				//��������Ȩֵ����0ʱ����ulist��ɾ��
				auto iter2 = G.ulist[v].begin();
				while (iter2 != G.ulist[v].end()) {
					if ((*iter2).v == u) {
						G.ulist[v].erase(iter2);
						break;
					}
					++iter2;
				}
			}


			/*if (G.edge[v/2-1].weight <= 0.0) {
				goto go_tag1;
			}*/
			tmp_tag_no = v / 2 - 1; //�ҵ���ǰ�ߵ����
			G.edge[tmp_tag_no].weight -= tag_w;
			if (G.edge[tmp_tag_no].weight <= 1.0e-10)
				G.edge[tmp_tag_no].weight = 0.0;
			f.insert(tmp_tag_no + 1);
			cout << " " << tmp_tag_no;
		}
		//������·��
		tmp_tag_no = tag_no;
		for (int i = t + 1; i <= k; i++) {
			u = 2 * (tmp_tag_no + 1);
			int ss = G.ulist[u].size();
			//�����ulist��ɾ�����������Լ����жϣ�����õ�û�г��ߣ���ѳ�ʼ����Ϊ0

			if (ss == 0) {
				G.edge[tmp_tag_no].weight = 0.0;
				cout << "xxxxxxxxxxxxxxxx" << endl;
				goto goodbye;
			}
			uniform_int_distribution<> ui(0, ss - 1);
		go_tag2:
			int tmp = ui(r);
			auto ulis = G.ulist[u][tmp];
			v = ulis.v;

			double fuck;
			if (ss == 1)
				fuck = tag_w - (1.0e-10);
			if (G.edge[ulis.no].weight < fuck) {
				cout << "weight:" << G.edge[ulis.no].weight << endl;
				goto go_tag2;
			}
			G.edge[ulis.no].weight -= tag_w;

			//�������
			if (G.edge[ulis.no].weight <= 1.0e-10) {
				//if (G.edge[ulis.no].weight <=0.0) {
				G.edge[ulis.no].weight = 0.0;
				//��������Ȩֵ����0ʱ����ulist��ɾ��
				auto iter1 = G.ulist[u].begin();
				while (iter1 != G.ulist[u].end()) {
					if ((*iter1).v == v) {
						G.ulist[u].erase(iter1);
						break;
					}
					++iter1;
				}
				//��������Ȩֵ����0ʱ����vlist��ɾ��
				auto iter2 = G.vlist[v].begin();
				while (iter2 != G.vlist[v].end()) {
					if ((*iter2).v == u) {
						G.vlist[v].erase(iter2);
						break;
					}
					++iter2;
				}
			}

			/*if (G.edge[(v + 1) / 2 - 1].weight <= 0.0) {
				goto go_tag2;
			}*/
			tmp_tag_no = (v + 1) / 2 - 1;
			G.edge[tmp_tag_no].weight -= tag_w;
			if (G.edge[tmp_tag_no].weight <= 1.0e-10)
				G.edge[tmp_tag_no].weight = 0.0;
			f.insert(tmp_tag_no + 1);
			cout << " " << tmp_tag_no;
		}
		FF.insert(f);
		f.clear();

		G.edge[tag_no].weight = 0.0;

	}

	return FF;
}

//get a collections of path-flows   
//����Ǹĺ�Ĵ��룬��������ʦ���ķ����㷨��̫һ��,����õ������
set<set<double>> GGFf(DirectedGraph G, int q, int k) {
	set<set<double>> FF;
	random_device r;
	int tag_no, tag_u, tag_v;
	double tag_w;//��ŵ�ǰ�ҵ�����С�ߵ���ź�Ȩֵ
	vector<int> no_set;
	while (1) {
		set<double> f;   //��һ��ֵ�洢·��f�ĸ��ʣ������k��ֵ�洢k����
		double minv = INFINITY;  //��¼��СȨֵ
		tag_no = -1;
		//�õ�Ȩֵ��С�ı�
		for (int i = 0; i < G.edge.size()-2*q; i++) {
			if (G.edge[i].weight < minv&&G.edge[i].weight != 0.0) {
				minv = G.edge[i].weight;
				tag_no = G.edge[i].no;
				tag_u = G.edge[i].u;
				tag_v = G.edge[i].v;
			}
		}
		cout << endl << "no: " << tag_no << "  min: " << minv;
		if (tag_no == -1)
			break;
		f.insert(minv);
		f.insert(tag_no + 1);

		//int t = tag_no / q + 1;  //��ǰ��ǰ·�����ڵ���
		int u = tag_u, v = tag_v;
		tag_w = minv;
		int debug = 0;
		//��ǰѰ��·��
		while (1) {
			//cout << "u:   " << u << endl;
			//����ǰ�ҵ��ĵ�һ������Դ��sʱ������ѭ��
			auto iter = G.vlist[u].begin();
			if ((*iter).v==2*k*q+1) {
				break;
			}
			///�㷨��ȷ���е�����£���������������
			int ss = G.vlist[u].size();
			if (ss == 0) {
				//G.edge[tmp_tag_no].weight = 0.0;
				cout << "+++++++++++" << endl;
				goto goodbye;
			}
			uniform_int_distribution<> ui(0, ss - 1);
		go_tag1:
			int tmp = ui(r);
			auto vlis = G.vlist[u][tmp];    //��ǰ�ҵ��ı�
			//cout << "no:" << vlis.no << endl;
			double fuck;

			if (abs(G.edge[vlis.no].weight - tag_w) < 1.0e-7|| G.edge[vlis.no].weight> tag_w) {
				G.edge[vlis.no].weight -= tag_w; //��ǰ�ߵ�Ȩ�ؼ�ȥ��ʼ�ߵ�Ȩ��
				goto go_tag11;
			}
			goto go_tag1;
		go_tag11:
			/*if (G.edge[vlis.no].weight < fuck) {
				cout << "weight:" << G.edge[vlis.no].weight << endl;
				goto go_tag1;
			}*/
			int tav = G.edge[vlis.no].u;
			
			//�������
			if (G.edge[vlis.no].weight <= 1.0e-10) {
				//if (G.edge[vlis.no].weight <=0.0) {
				G.edge[vlis.no].weight = 0.0;            //�����ߵ�ȨֵС��1.0e-10ʱ����0

				//��������Ȩֵ����0ʱ����vlist��ɾ��
				auto iter1 = G.vlist[u].begin();
				while (iter1 != G.vlist[u].end()) {
					if ((*iter1).v == tav) {
						G.vlist[u].erase(iter1);
						break;
					}
					++iter1;
				}
				//��������Ȩֵ����0ʱ����ulist��ɾ��
				auto iter2 = G.ulist[tav].begin();
				while (iter2 != G.ulist[tav].end()) {
					if ((*iter2).v == u) {
						G.ulist[tav].erase(iter2);
						break;
					}
					++iter2;
				}
			}
			u = G.edge[vlis.no].u;
			if (vlis.no < k*q) {
				f.insert(vlis.no + 1);
			}
			
		}
		u = tag_u, v = tag_v;
		while (1) {
			//cout << "v:   " << v << endl;
			auto iter1 = G.ulist[v].begin();
			if ((*iter1).v == 2 * k*q + 2) {
				break;
			}
			int ss = G.ulist[v].size();
			//�����ulist��ɾ�����������Լ����жϣ�����õ�û�г��ߣ���ѳ�ʼ����Ϊ0

			if (ss == 0) {
				//G.edge[tmp_tag_no].weight = 0.0;
				cout << "xxxxxxxxxxxxxxxx" << endl;
				goto goodbye;
			}
			uniform_int_distribution<> ui(0, ss - 1);
		go_tag2:
			int tmp = ui(r);
			auto ulis = G.ulist[v][tmp];
			//cout << "no:" << ulis.no << endl;
			double fuck;
			if (ss == 1&&ulis.no < k*q)
				fuck = tag_w - (1.0e-7);
			if (abs(G.edge[ulis.no].weight - tag_w) < 1.0e-7 || G.edge[ulis.no].weight > tag_w) {
				G.edge[ulis.no].weight -= tag_w; //��ǰ�ߵ�Ȩ�ؼ�ȥ��ʼ�ߵ�Ȩ��
				goto go_tag22;
			}
			goto go_tag2;
		go_tag22:
			//if(abs(G.edge[ulis.no].weight - tag_w)<1.0e-7){
			////if (G.edge[ulis.no].weight < fuck) {
			//	cout << "weight:" << G.edge[ulis.no].weight << endl;
			//	goto go_tag2;
			//}
			int tau = G.edge[ulis.no].v;

			//�������
			if (G.edge[ulis.no].weight <= 1.0e-10) {
				//if (G.edge[ulis.no].weight <=0.0) {
				G.edge[ulis.no].weight = 0.0;
				//��������Ȩֵ����0ʱ����ulist��ɾ��
				auto iter1 = G.ulist[v].begin();
				while (iter1 != G.ulist[v].end()) {
					if ((*iter1).v == v) {
						G.ulist[v].erase(iter1);
						break;
					}
					++iter1;
				}
				//��������Ȩֵ����0ʱ����vlist��ɾ��
				auto iter2 = G.vlist[tau].begin();
				while (iter2 != G.vlist[tau].end()) {
					if ((*iter2).v == v) {
						G.vlist[tau].erase(iter2);
						break;
					}
					++iter2;
				}
			}
			v = G.edge[ulis.no].v;
			if (ulis.no < k*q) {
				f.insert(ulis.no + 1);
			}
			
		}
		//����·��
		FF.insert(f);
		f.clear();
	goodbye:
		G.edge[tag_no].weight = 0.0;

	}

	return FF;
}
