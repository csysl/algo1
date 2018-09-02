#include"func.h"

//使用GLPK计算mip问题
void gMIP(int n, int k, double *weight, set<set<int>> SS, set<set<int>> *GG, int q) {


	int m = SS.size();//获取SS的长度

	char colName[10], rowName[10];
	int *ind = new int[12 * n + 1];  //m最大为10*n
	double *val = new double[12 * n + 1];

	glp_prob *lp;
	glp_iocp parm;
	lp = glp_create_prob();
	glp_set_prob_name(lp, "CMCG");
	glp_set_obj_name(lp, "max");
	glp_set_obj_dir(lp, GLP_MAX);
	//添加列
	glp_add_cols(lp, m + n);
	//前n列变量为y1-yn
	for (int i = 1; i <= n; i++) {
		sprintf_s(colName, "y%1d", i);
		glp_set_col_name(lp, i, colName);
		glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
		glp_set_col_kind(lp, i, GLP_BV);
		glp_set_obj_coef(lp, i, weight[i]);
	}
	//后m列变量为x1-xm
	for (int i = n + 1; i <= m + n; i++) {
		sprintf_s(colName, "x%1d", i - n);
		glp_set_col_name(lp, i, colName);
		glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
		glp_set_col_kind(lp, i, GLP_BV);
	}
	//添加行，行数n+q+1
	glp_add_rows(lp, n + q + 1);
	//前n行，s.t.1
	for (int i = 1; i <= n; i++) {
		sprintf_s(rowName, "st1_%ld", i);
		glp_set_row_name(lp, i, rowName);
		glp_set_row_bnds(lp, i, GLP_LO, 0.0, 0.0);
		//s.t.1:y[i]<=sum(x[j]),j: i∈Sj
		ind[1] = i;
		val[1] = -1;
		set<set<int>>::iterator itSS = SS.begin();
		for (int j = n + 1; j <= n + m; j++) {
			ind[j - n + 1] = j;
			//i属于S[j-n]
			if ((*itSS).count(i)) {
				val[j - n + 1] = 1;
			}
			else
				val[j - n + 1] = 0;
			itSS++;
		}
		glp_set_mat_row(lp, i, m + 1, ind, val);
	}
	//中间q行，s.t.2
	for (int i = 1; i <= q; i++) {
		sprintf_s(rowName, "st2_%ld", i);
		glp_set_row_name(lp, i + n, rowName);
		glp_set_row_bnds(lp, i + n, GLP_UP, 0.0, 1.0);
		set<set<int>>::iterator itSS = SS.begin();

		//s.t.2:sum(x[j])<=n[l],j: Sj∈Gl
		for (int j = n + 1; j <= n + m; j++) {
			ind[j - n] = j;

			if (GG[i].count(*itSS)) {
				val[j - n] = 1;
			}
			else
				val[j - n] = 0;
			itSS++;
			//cout << "va" << j - n << ":" << val[j - n] << "   ";
		}
		//cout << endl;
		glp_set_mat_row(lp, i + n, m, ind, val);
	}
	//最后一行，s.t.3
	glp_set_row_name(lp, n + q + 1, "st3");
	glp_set_row_bnds(lp, n + q + 1, GLP_UP, 0.0, k);
	for (int j = n + 1; j <= n + m; j++) {
		ind[j - n] = j;
		val[j - n] = 1;
	}
	glp_set_mat_row(lp, n + q + 1, m, ind, val);

	//初始化MIP
	glp_init_iocp(&parm);
	parm.presolve = GLP_ON;
	//计算时间
	clock_t startTime, endTime;
	startTime = clock();
	//计算MIP
	glp_intopt(lp, &parm);
	endTime = clock();
	cout << "Time usage of MIP by GLPK: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

	//结果输出到gMIP.txt
	glp_print_mip(lp, "gMIP.txt");

	//判断结果
	int ret = glp_mip_status(lp);
	switch (ret) {
	case GLP_OPT:
		printf("Success\n");
		break;
	default:
		printf("Failed\n");
		// delete problem
		glp_delete_prob(lp);
		return;
	}
	//输出x的取值
	for (int i = 1; i <= m; i++) {
		int xn = glp_mip_col_val(lp, i + n);
		if (xn == 1) {
			printf("S%ld:{ ", i);
			set<set<int>>::iterator item = SS.begin();
			for (int j = 1; j < i; j++) {
				item++;
			}
			set<int>::iterator item2 = (*item).begin();
			while (item2 != (*item).end()) {
				cout << *item2 << " ";
				item2++;
			}
			cout << "}" << endl;
		}
	}
	cout << "Elements are not selected:";
	for (int i = 1; i <= n; i++) {
		int yn = glp_mip_col_val(lp, i);
		//y[n]为0时输出
		if (yn == 0) {
			printf("%ld ", i);
		}
	}
	cout << endl;

	//释放内存
	glp_delete_prob(lp);
	delete[]ind;
	delete[]val;
}

//使用CPLEX计算mip问题
void cMIP(int n, int k, double *weight, set<set<int>> SS, set<set<int>> *GG, int q) {

}

//get Z* by GLPK
double *gLP(int n, int k, double *weight, set<set<int>> SS, set<set<int>> *GG, int q) {
	int m = SS.size();//获取SS的长度

	char colName[10], rowName[10];

	glp_prob *lp;
	glp_smcp parm;
	lp = glp_create_prob();
	glp_set_prob_name(lp, "gLP");
	glp_set_obj_name(lp, "max_gLP");
	glp_set_obj_dir(lp, GLP_MAX);
	//添加列
	glp_add_cols(lp, 2 * m + n);
	//前n列变量为y1-yn
	for (int i = 1; i <= n; i++) {
		sprintf_s(colName, "y%1d", i);
		glp_set_col_name(lp, i, colName);
		glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(lp, i, weight[i]);
	}
	//后m列变量为x1-xm
	for (int i = n + 1; i <= m + n; i++) {
		sprintf_s(colName, "x%1d", i - n);
		glp_set_col_name(lp, i, colName);
		glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
		//glp_set_obj_coef(lp, i, 0.0);
	}
	//最后m列为辅助参数yn+1到yn+m   小于0.001
	for (int i = m + n + 1; i <= 2 * m + n; i++) {
		sprintf_s(colName, "y%1d", i - m);
		glp_set_col_name(lp, i, colName);
		glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(lp, i, weight[i - m]);
	}
	//cout << "debug1" << endl;
	//添加行，行数n+q+1
	glp_add_rows(lp, n + q + 1 + m);
	//前n行，s.t.1
	set<set<int>>::iterator itSS = SS.begin();
	for (int i = 1; i <= n; i++) {
		int *ind = new int[2 * m + n + 2];
		double *val = new double[2 * m + n + 2];
		sprintf_s(rowName, "st1_%ld", i);
		glp_set_row_name(lp, i, rowName);
		glp_set_row_bnds(lp, i, GLP_LO, 0.0, 0.0);
		//s.t.1:y[i]<=sum(x[j]),j: i∈Sj
		ind[1] = i;
		val[1] = -1.0;
		itSS = SS.begin();
		for (int j = n + 1; j <= n + m; j++) {
			ind[j - n + 1] = j;
			//i属于S[j-n]
			if ((*itSS).count(i)) {
				val[j - n + 1] = 1.0;
			}
			else
				val[j - n + 1] = 0.0;
			itSS++;
		}
		glp_set_mat_row(lp, i, m + 1, ind, val);
		delete[]ind, val;
	}
	//cout << "debug2" << endl;
	//中间q行，s.t.2
	for (int i = 1; i <= q; i++) {
		int *ind = new int[2 * m + n + 2];
		double *val = new double[2 * m + n + 2];
		sprintf_s(rowName, "st2_%ld", i);
		glp_set_row_name(lp, i + n, rowName);
		glp_set_row_bnds(lp, i + n, GLP_UP, 0.0, 1.0);
		itSS = SS.begin();
		//s.t.2:sum(x[j])<=n[l],j: Sj∈Gl
		for (int j = n + 1; j <= n + m; j++) {
			ind[j - n] = j;

			if (GG[i].count(*itSS)) {
				val[j - n] = 1.0;
			}
			else
				val[j - n] = 0.0;
			itSS++;
			//cout << "va" << j - n << ":" << val[j - n] << "   ";
		}
		//cout << endl;
		glp_set_mat_row(lp, i + n, m, ind, val);
		delete[]ind, val;
	}
	//cout << "debug3" << endl;
	//最后一行，s.t.3
	int *indd = new int[2 * m + n + 1];
	double *vall = new double[2 * m + n + 1];
	glp_set_row_name(lp, n + q + 1, "st3");
	glp_set_row_bnds(lp, n + q + 1, GLP_UP, 0.0, k);
	for (int j = n + 1; j <= n + m; j++) {
		indd[j - n] = j;
		vall[j - n] = 1;
	}
	glp_set_mat_row(lp, n + q + 1, m, indd, vall);
	delete[]indd, vall;

	int tag = n + q + 2;
	//cout << "debug4" << endl;
	//辅助的m行
	for (int i = 1; i <= m; i++) {
		int *ind = new int[2 * m + n + 2];
		double *val = new double[2 * m + n + 2];
		sprintf_s(rowName, "st1_%ld", tag);
		glp_set_row_name(lp, tag, rowName);
		glp_set_row_bnds(lp, tag, GLP_LO, 0.0, 0.0);
		//s.t.1:y[i]<=sum(x[j]),j: i∈Sj
		ind[1] = m + n + i;
		val[1] = -1.0;
		itSS = SS.begin();
		for (int j = n + 1; j <= n + m; j++) {
			ind[j - n + 1] = j;
			//n+i属于S[j-n]
			if ((*itSS).count(n + i)) {
				val[j - n + 1] = 1.0;
				//cout << "+++++++++++++++" << (n+i)<< endl;
			}
			else
				val[j - n + 1] = 0.0;
			itSS++;
		}
		glp_set_mat_row(lp, tag, m + 1, ind, val);
		tag++;
		delete[]ind, val;
	}
	//cout << "debug5" << endl;
	//初始化LP
	glp_init_smcp(&parm);
	parm.meth = GLP_PRIMAL;
	//计算时间
	clock_t startTime, endTime;
	startTime = clock();
	//计算LP
	glp_simplex(lp, &parm);
	endTime = clock();
	cout << "Time usage to get Z* by GLPK: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

	//结果输出到gLP.txt
	glp_print_sol(lp, "gLP.txt");

	//判断结果
	int ret = glp_get_status(lp);
	switch (ret) {
	case GLP_OPT:
		cout << "success to get Z*!" << endl;
		break;
	default:
		cout << "failed to get Z*" << endl;
		break;
	}
	//输出x的取值
	/*double *Y = new double[n + 1];
	for (int i = 1; i <= n; i++) {
		double xn = glp_get_col_prim(lp, i);
		Y[i] = xn;
		//if (xn != 0.0)
		cout << "y" << i << ":" << xn << endl;
	}*/
	double *X = new double[m + 1];
	for (int i = 1; i <= m; i++) {
		double xn = glp_get_col_prim(lp, i + n);
		X[i] = xn;

	}
	//cout << "debug6" << endl;

	cout << endl;

	//释放内存
	glp_delete_prob(lp);

	return X;
}
//使用CPLEX计算lp问题
double *cLP(int n, int k, double *weight, set<set<int>> SS, set<set<int>> *GG, int q) {
	int m = SS.size();
	IloEnv env;
	double *X = new double[m + 1];
	try {
		IloModel lp(env);
		//添加参数,前n个为y，中间m个为辅助参数,最后n个为x
		IloNumVarArray var(env);
		for (int i = 1; i <= 2 * m + n; i++) {
			var.add(IloNumVar(env, 0.0, 1.0));
		}
		lp.add(var);
		//添加目标函数
		IloExpr expr(env);
		for (int i = 0; i < m + n; i++) {
			expr += weight[i + 1] * var[i];
		}
		IloObjective obj(IloMaximize(env, expr, "cLP"));
		lp.add(obj);
		//添加行约束
		IloRangeArray con(env);
		//st1添加前n+m行约束
		auto itSS = SS.begin();
		for (int i = 0; i < m + n; i++) {
			itSS = SS.begin();
			IloExpr exp(env);
			exp = -var[i];
			for (int j = 0; j < m; j++) {
				if ((*itSS).count(i + 1))
					exp += var[m + n + j];
				itSS++;
			}
			con.add(IloRange(env, 0.0, exp, INFINITY));
			exp.end();
		}
		//添加中间q行
		int tag = n + m;
		for (int i = 1; i <= q; i++) {
			itSS = SS.begin();
			IloExpr exp(env);
			int tag_g = GG[i].size();
			for (int j = 0; j < tag_g; j++) {
				exp += var[tag];
				tag++;
			}
			con.add(IloRange(env, 0.0, exp, 1.0));
			exp.end();
		}
		//添加最后一行
		IloExpr exp(env);
		for (int j = 0; j < m; j++) {
			exp += var[j + m + n];
		}
		con.add(IloRange(env, 0.0, exp, double(k)));
		lp.add(con);
		//求解
		IloCplex cplex(lp);
		cplex.exportModel("cLP.lp");
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << endl;
			throw(-1);
		}
		IloNumArray vals(env);
		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value = " << cplex.getObjValue() << endl;
		cplex.getValues(vals, var);
		for (int i = m + n; i < 2 * m + n; i++) {
			X[i - m - n + 1] = vals[i];
		}

	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}


	env.end();
	return X;
}