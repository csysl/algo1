#include"func.h"
//get a collections of path-flows
//这段存放了最开始的代码
set<set<double>> GGGFf(DirectedGraph G, int q, int k) {
	set<set<double>> FF;
	random_device r;
	int tag_no;
	double tag_w;//存放当前找到的最小边的序号和权值
	while (1) {
		set<double> f;   //第一个值存储路径f的概率，后面的k个值存储k条边
		double minv = INFINITY;  //记录最小权值
		tag_no = -1;
		//得到权值最小的边
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

		int t = tag_no / q + 1;  //当前当前路径所在的列
		int u, v, tmp_tag_no = tag_no;
		tag_w = minv;
		//往前找路径
		for (int i = t - 1; i >= 1; i--) {
			u = 2 * (tmp_tag_no + 1) - 1;  //当前边的出发点
			int ss = G.vlist[u].size();
			if (ss == 0)//后面对vlist有删除操作，所以加入判断，如果该点没有入边，则把初始边置为0 
				goto goodbye;
			uniform_int_distribution<> ui(0, ss - 1);
			int tmp = ui(r);
			v = G.vlist[u][tmp].v;  //随机得到与其相连的一个出发点

			G.edge[G.vlist[u][tmp].no].weight -= tag_w; //辅助边
			//消除误差
			if (G.edge[G.vlist[u][tmp].no].weight <= 1.0e-7) {
				//if (G.edge[G.vlist[u][tmp].no].weight <=0.0) {
				G.edge[G.vlist[u][tmp].no].weight = 0.0;//辅助边的权值小于1.0e-7时，置0
				//当辅助边权值减到0时，从vlist中删除
				auto iter1 = G.vlist[u].begin();
				while (iter1 != G.vlist[u].end()) {
					if ((*iter1).v == v) {
						G.vlist[u].erase(iter1);
						break;
					}
					++iter1;
				}
				//当辅助边权值减到0时。从ulist中删除
				auto iter2 = G.ulist[v].begin();
				while (iter2 != G.ulist[v].end()) {
					if ((*iter2).v == u) {
						G.ulist[v].erase(iter2);
						break;
					}
					++iter2;
				}
			}

			tmp_tag_no = v / 2 - 1; //找到的前边的序号
			G.edge[tmp_tag_no].weight -= tag_w;
			if (G.edge[tmp_tag_no].weight <= 1.0e-7)
				G.edge[tmp_tag_no].weight = 0.0;
			f.insert(tmp_tag_no + 1);
		}
		//往后着路径
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
			//消除误差
			if (G.edge[G.ulist[u][tmp].no].weight <= 1.0e-7) {
				//if (G.edge[G.ulist[u][tmp].no].weight <=0.0) {
				G.edge[G.ulist[u][tmp].no].weight = 0.0;
				//当辅助边权值减到0时，从ulist中删除
				auto iter1 = G.ulist[u].begin();
				while (iter1 != G.ulist[u].end()) {
					if ((*iter1).v == v) {
						G.ulist[u].erase(iter1);
						break;
					}
					++iter1;
				}
				//当辅助边权值减到0时。从vlist中删除
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
//这段代码是老师给的分流算法，存在问题   即当主边最小时，不一定存在满足条件的辅助边
set<set<double>> GFf(DirectedGraph G, int q, int k) {
	set<set<double>> FF;
	random_device r;
	int tag_no;
	double tag_w;//存放当前找到的最小边的序号和权值
	vector<int> no_set;
	while (1) {

	goodbye:
		set<double> f;   //第一个值存储路径f的概率，后面的k个值存储k条边
		double minv = INFINITY;  //记录最小权值
		tag_no = -1;
		//得到权值最小的边
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

		int t = tag_no / q + 1;  //当前当前路径所在的列
		int u, v, tmp_tag_no = tag_no;
		tag_w = minv;
		int debug = 0;
		//往前找路径
		for (int i = t - 1; i >= 1; i--) {
			/*
			1、得到当前边（非辅助边）
			2、检查当前边是否有入边，没有，把该边置0，同时回到起点   问题：当前边之前的边已经减过被选中的初始边的权值
			3、有入边，则随机挑选一条辅助入边，减去初始边的权重，如果减完后小于1.0e-7，置0，并将该辅助边从vlist和ulist中删除
			4、
			*/
			u = 2 * (tmp_tag_no + 1) - 1;  //当前边的出发点
			int ss = G.vlist[u].size();

			//后面对vlist有删除操作，所以加入判断，如果该点没有入边，则把初始边置为0
			if (ss == 0) {
				G.edge[tmp_tag_no].weight = 0.0;
				cout << "+++++++++++" << endl;
				goto goodbye;
			}
			uniform_int_distribution<> ui(0, ss - 1);
		go_tag1:
			int tmp = ui(r);
			auto vlis = G.vlist[u][tmp];
			v = vlis.v;  //随机得到与其相连的一个出发点
			double fuck;
			if (ss == 1)
				fuck = tag_w - (1.0e-10);
			if (G.edge[vlis.no].weight < fuck) {
				cout << "weight:" << G.edge[vlis.no].weight << endl;
				goto go_tag1;
			}
			G.edge[vlis.no].weight -= tag_w; //辅助边

			//消除误差
			if (G.edge[vlis.no].weight <= 1.0e-10) {
				//if (G.edge[vlis.no].weight <=0.0) {
				G.edge[vlis.no].weight = 0.0;            //辅助边的权值小于1.0e-10时，置0

				//当辅助边权值减到0时，从vlist中删除
				auto iter1 = G.vlist[u].begin();
				while (iter1 != G.vlist[u].end()) {
					if ((*iter1).v == v) {
						G.vlist[u].erase(iter1);
						break;
					}
					++iter1;
				}
				//当辅助边权值减到0时。从ulist中删除
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
			tmp_tag_no = v / 2 - 1; //找到的前边的序号
			G.edge[tmp_tag_no].weight -= tag_w;
			if (G.edge[tmp_tag_no].weight <= 1.0e-10)
				G.edge[tmp_tag_no].weight = 0.0;
			f.insert(tmp_tag_no + 1);
			cout << " " << tmp_tag_no;
		}
		//往后着路径
		tmp_tag_no = tag_no;
		for (int i = t + 1; i <= k; i++) {
			u = 2 * (tmp_tag_no + 1);
			int ss = G.ulist[u].size();
			//后面对ulist有删除操作，所以加入判断，如果该点没有出边，则把初始边置为0

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

			//消除误差
			if (G.edge[ulis.no].weight <= 1.0e-10) {
				//if (G.edge[ulis.no].weight <=0.0) {
				G.edge[ulis.no].weight = 0.0;
				//当辅助边权值减到0时，从ulist中删除
				auto iter1 = G.ulist[u].begin();
				while (iter1 != G.ulist[u].end()) {
					if ((*iter1).v == v) {
						G.ulist[u].erase(iter1);
						break;
					}
					++iter1;
				}
				//当辅助边权值减到0时。从vlist中删除
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
//这段是改后的代码，方法与老师给的分流算法不太一样,最后用的是这段
set<set<double>> GGFf(DirectedGraph G, int q, int k) {
	set<set<double>> FF;
	random_device r;
	int tag_no, tag_u, tag_v;
	double tag_w;//存放当前找到的最小边的序号和权值
	vector<int> no_set;
	while (1) {
		set<double> f;   //第一个值存储路径f的概率，后面的k个值存储k条边
		double minv = INFINITY;  //记录最小权值
		tag_no = -1;
		//得到权值最小的边
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

		//int t = tag_no / q + 1;  //当前当前路径所在的列
		int u = tag_u, v = tag_v;
		tag_w = minv;
		int debug = 0;
		//向前寻找路径
		while (1) {
			//cout << "u:   " << u << endl;
			//当往前找到的第一个点是源点s时，跳出循环
			auto iter = G.vlist[u].begin();
			if ((*iter).v==2*k*q+1) {
				break;
			}
			///算法正确运行的情况下，不会有下面的输出
			int ss = G.vlist[u].size();
			if (ss == 0) {
				//G.edge[tmp_tag_no].weight = 0.0;
				cout << "+++++++++++" << endl;
				goto goodbye;
			}
			uniform_int_distribution<> ui(0, ss - 1);
		go_tag1:
			int tmp = ui(r);
			auto vlis = G.vlist[u][tmp];    //当前找到的边
			//cout << "no:" << vlis.no << endl;
			double fuck;

			if (abs(G.edge[vlis.no].weight - tag_w) < 1.0e-7|| G.edge[vlis.no].weight> tag_w) {
				G.edge[vlis.no].weight -= tag_w; //当前边的权重减去初始边的权重
				goto go_tag11;
			}
			goto go_tag1;
		go_tag11:
			/*if (G.edge[vlis.no].weight < fuck) {
				cout << "weight:" << G.edge[vlis.no].weight << endl;
				goto go_tag1;
			}*/
			int tav = G.edge[vlis.no].u;
			
			//消除误差
			if (G.edge[vlis.no].weight <= 1.0e-10) {
				//if (G.edge[vlis.no].weight <=0.0) {
				G.edge[vlis.no].weight = 0.0;            //辅助边的权值小于1.0e-10时，置0

				//当辅助边权值减到0时，从vlist中删除
				auto iter1 = G.vlist[u].begin();
				while (iter1 != G.vlist[u].end()) {
					if ((*iter1).v == tav) {
						G.vlist[u].erase(iter1);
						break;
					}
					++iter1;
				}
				//当辅助边权值减到0时。从ulist中删除
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
			//后面对ulist有删除操作，所以加入判断，如果该点没有出边，则把初始边置为0

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
				G.edge[ulis.no].weight -= tag_w; //当前边的权重减去初始边的权重
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

			//消除误差
			if (G.edge[ulis.no].weight <= 1.0e-10) {
				//if (G.edge[ulis.no].weight <=0.0) {
				G.edge[ulis.no].weight = 0.0;
				//当辅助边权值减到0时，从ulist中删除
				auto iter1 = G.ulist[v].begin();
				while (iter1 != G.ulist[v].end()) {
					if ((*iter1).v == v) {
						G.ulist[v].erase(iter1);
						break;
					}
					++iter1;
				}
				//当辅助边权值减到0时。从vlist中删除
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
		//插入路径
		FF.insert(f);
		f.clear();
	goodbye:
		G.edge[tag_no].weight = 0.0;

	}

	return FF;
}
