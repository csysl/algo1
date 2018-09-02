#include"func.h"

//get the auxiliary graph
vector<Edge> Edgess(int q, int k) {
	/*
	edges:   存放所有的边，每条边形式为{u,v,no,weight}  no:0,1,2,....
	u,v：    共2*k*q个
	s,t:     源点和汇点，序号分别为2*k*q+1和2*k*q+2
	*/
	vector<Edge> edges;
	Edge edge;
	int tag = 0;
	//添加主要的边
	for (int i = 1; i <= k * q; i++) {
		edge.u = 2 * i - 1;
		edge.v = 2 * i;
		edge.no = tag;
		edge.weight = 0.0;
		edges.push_back(edge);
		tag++;
	}
	//添加辅助的边
	for (int s = 1; s <= q - 1; s++) {
		for (int g = s + 1; g <= q; g++) {
			for (int t = 0; t < k - 1; t++) {
				edge.u = t * 2 * q + 2 * s;
				edge.v = (t + 1) * 2 * q + 2 * g - 1;
				edge.no = tag;
				edge.weight = 0.0;
				edges.push_back(edge);
				tag++;
			}
		}
	}
	//添加连接源点的边
	for (int i = 1; i <= q; i++) {
		edge.u = 2 * k*q + 1;
		edge.v = 2 * i - 1;
		edge.no = tag;
		edge.weight = 0.0;
		edges.push_back(edge);
		tag++;
	}
	//添加连接汇点的边
	for (int i = 1; i <= q; i++) {
		edge.u = 2 * (k - 1)*q + 2 * i;
		edge.v = 2 * k*q + 2;
		edge.no = tag;
		edge.weight = 0.0;
		edges.push_back(edge);
		tag++;
	}
	//cout << "the num of edges:" << tag << endl;

	return edges;
}

//get an stFlow with value 1 by cplex
DirectedGraph Gcstflow(DirectedGraph G, vector<Edge> edges, double *Z, int q, int k, int n) {
	IloEnv env;
	int num_col = k * q + (k - 1)*q*(q - 1) / 2 + 2 * q;
	cout << "the num of edges:" << num_col << endl;
	cout << "the num of edges:" << edges.size() << endl;
	try {
		IloModel st(env);
		IloNumVarArray vars(env);
		//添加列
		for (int i = 0; i < num_col; i++) {
			vars.add(IloNumVar(env, 0.0, 1.0));
		}
		st.add(vars);
		//添加目标
		IloExpr exp(env);
		for (int i = 0, t = 1; i < k*q; i++, t++) {
			exp += t * vars[i];
			if (t > q)
				t = 1;
		}
		IloObjective obj(IloMaximize(env, exp, "DcstFlow"));
		st.add(obj);
		//添加约束
		IloRangeArray cons(env);
		//前q个约束
		for (int i = 0; i < q; i++) {
			IloExpr texp(env);
			for (int j = 0; j < k; j++) {
				texp += vars[j*q + i];
			}
			cons.add(IloRange(env, Z[i + 1], texp, Z[i + 1]));
			texp.end();
		}
		//出边的约束
		for (int i = 1; i <= k * q; i++) {
			int u = 2 * i;
			IloExpr texp(env);
			for (auto iter : G.ulist[u]) {
				texp += vars[iter.no];
			}
			texp -= vars[i - 1];
			cons.add(IloRange(env, 0.0, texp, 0.0));
			texp.end();
		}
		//入边的约束
		for (int i = 1; i <= k * q; i++) {
			int v = 2 * i - 1;
			IloExpr texp(env);
			for (auto iter : G.vlist[v]) {
				texp += vars[iter.no];
			}
			texp -= vars[i - 1];
			cons.add(IloRange(env, 0.0, texp, 0.0));
			texp.end();
		}
		//value = 1
		IloExpr texp(env);
		for (auto iter : G.ulist[2 * k*q + 1]) {
			texp += vars[iter.no];
		}
		cons.add(IloRange(env, 1.0, texp, 1.0));
		texp.end();

		st.add(cons);
		//计算结果
		IloCplex cplex(st);
		cplex.exportModel("GcstFlow.lp");
		if (!cplex.solve()) {
			env.error() << "Failed to get GstFlow." << endl;
			throw(-1);
		}
		IloNumArray vals(env);
		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value = " << cplex.getObjValue() << endl;
		cplex.getValues(vals, vars);

		for (int i = 0; i < num_col; i++) {
			edges[i].weight = vals[i];
		}

	}
	catch (IloException& e) {
		cerr << "Concert Exception: " << e << endl;
	} //用于检测cplex的异常
	catch (...) {
		cerr << "Other Exception" << endl;
	} //用于检测非cplex的异常
	env.end();
	DirectedGraph GG(edges, n, 1);
	return GG;
}

//select a path from path-flows
double *ff(set<set<double>> FF, int k) {
	double *f = new double[k + 2];  //用来存放选出的路径的序号
	int size = FF.size();
	auto it1 = FF.begin();
	set<double>::iterator it2;

	random_device r;
	//理想情况下tmpp应该取1
	double tmpp = 0;
	for (auto iter : FF) {
		tmpp += (*iter.begin());
	}
	uniform_real_distribution<double> u(0.0, tmpp);
	double tag = u(r);      //随机产生一个随机数，落在哪个区间，就选择哪个f
	double tmp = 0.0;
	//cout << "debug1" << endl;
	while (it1 != FF.end()) {
		it2 = (*it1).begin();
		tmp += (*it2);
		if (tag <= tmp)
			break;
		++it1;
	}
	//cout << "debug2" << endl;

	if ((*it1).size() != k)//使用set存放路径，如果概率是1的话，且起始路径的序号是1，则set会只保留一个1；当两者不重合时，从第二个元素开始取
		++it2;
	int i = 1;
	while (it2 != (*it1).end()) {

		f[i] = (*it2);
		//cout << f[i] << "\t";
		++i;
		++it2;
	}
	cout << endl;
	//cout << tag << "\t" << tmp << endl;

	return f;
}
