#include"func.h"

/*
���cpp�ĺ���û���õ�
*/

//Compute an st-flow of value 1 in G by GLPK
double *gstFlow(double *Z, int q, int k) {
	glp_prob *lp;
	lp = glp_create_prob();
	glp_smcp parm;
	glp_set_prob_name(lp, "stFlow");
	glp_set_obj_name(lp, "st_v1");
	glp_set_obj_dir(lp, GLP_MAX);
	//store ind and value of auxiliary 

	int tmp_col = k * q + ((k - 1)*(q - k + 2)*(q - k + 1)) / 2 + 2 * (q - k + 1);

	int *ind = new int[tmp_col + 2];
	double *val = new double[tmp_col + 2];
	glp_add_cols(lp, tmp_col);
	//add columns:edges    k*q
	int tag = 1;
	for (int i = 1; i <= k; i++) {
		for (int j = 1; j <= q; j++) {
			glp_set_col_bnds(lp, tag, GLP_DB, 0.0, 1.0);
			//find an st-flow with value 1 meet the objective
			//objective
			glp_set_obj_coef(lp, tag, j);
			tag++;
		}
	}
	//cout << "debug1" << endl;
	//add columns:auxiliary edgs   0.5*(k-1)*(q-k+2)*(q-k+1)
	for (int i = 1; i <= k - 1; i++) {
		int tmpp = q - k + 1;
		for (int j = 1; j <= (q - k + 1); j++) {
			for (int l = 1; l <= tmpp; l++) {
				glp_set_col_bnds(lp, tag, GLP_DB, 0.0, 1.0);
				tag++;
			}
			tmpp--;
		}
	}
	//cout << "debug2" << endl;
	//add columns:s- and -t  2*(q-k+1)
	for (int i = 1; i <= q - k + 1; i++) {
		glp_set_col_bnds(lp, tag, GLP_DB, 0.0, 1.0);
		tag++;
	}
	for (int i = 1; i <= q - k + 1; i++) {
		glp_set_col_bnds(lp, tag, GLP_DB, 0.0, 1.0);
		tag++;
	}
	//cout << "debug3" << endl;

	tag = 1;
	//add rows
	//int tmp_row = q + 2 * (k - 1)*(q - k + 1) + 1 + (q - k + 1)*2;
	int tmp_row = q + 2 * k*(q - k + 1) + 1;
	glp_add_rows(lp, tmp_row);
	//add q rows:s.t.2
	for (int i = 1; i <= q; i++) {
		glp_set_row_bnds(lp, i, GLP_FX, Z[i], Z[i]);
		for (int j = 1; j <= k; j++) {
			ind[j] = (j - 1)*q + i;
			val[j] = 1;
		}
		glp_set_mat_row(lp, i, k, ind, val);
		tag++;
	}
	delete[]ind, val;
	//cout << "debug2" << endl;

	//add 2 * (k - 1)*(q - k + 1) rows:s.t.1
	//cout << (2 * (k - 1)*(q - k + 1)) << endl;
	for (int i = 1; i <= k - 1; i++) {
		int tmpp = q - k + 1;
		int tmppp = 1;
		for (int j = 1; j <= q - k + 1; j++) {
			int *indd = new int[tmp_col + 2];
			double *vall = new double[tmp_col + 2];
			glp_set_row_bnds(lp, tag, GLP_FX, 0.0, 0.0);
			for (int l = 1; l <= tmpp; l++) {
				indd[l] = k * q + (i - 1)*(q - k + 2)*(q - k + 1) / 2 + tmppp;
				vall[l] = -1.0;
				tmppp++;
			}
			//cout << "    debug" << endl;
			indd[tmpp + 1] = (i - 1)*q + i - 1 + j;//�����ıߵ����
			vall[tmpp + 1] = 1.0;
			//cout << "debug+++++++" << endl;
			glp_set_mat_row(lp, tag, tmpp + 1, indd, vall);
			//cout << "debug-------" << endl;
			tmpp--;
			tag++;
			delete[]indd, vall;
			//cout << "debug" << (i - 1)*(q - k + 1) + j << endl;
		}
	}
	//cout << "debug3:  " << tag << endl;
	int tmpp, tmpp2, tmppp;
	for (int i = 1; i <= k - 1; i++) {
		tmpp = q - k + 1;
		tmpp2 = k * q + i * (tmpp + 1)*tmpp / 2 + 1;
		for (int j = q - k + 1, tt = 1; j >= 1; j--, tt++) {
			int *indd = new int[tmp_col + 2];
			double *vall = new double[tmp_col + 2];
			glp_set_row_bnds(lp, tag, GLP_FX, 0.0, 0.0);

			tmpp2 = tmpp2 - tt;
			tmppp = tmpp2;
			//cout << "debug+++++++" << endl;
			for (int l = 1; l <= tmpp; l++) {
				indd[l] = tmppp;
				vall[l] = -1.0;
				tmppp = tmppp - (l + tt - 1);
			}
			//cout << "debug--------" << endl;
			indd[tmpp + 1] = i * q + i + j;//�յ��ıߵ����
			vall[tmpp + 1] = 1.0;
			glp_set_mat_row(lp, tag, tmpp + 1, indd, vall);
			tmpp--;
			tag++;
			delete[]indd, vall;
		}
	}
	//cout << "debug4:  " << tag << endl;
	//add 2*(q-k+1) rows 
	int enu = k * q + ((k - 1)*(q - k + 2)*(q - k + 1)) / 2;
	for (int i = 1; i <= q - k + 1; i++) {
		glp_set_row_bnds(lp, tag, GLP_FX, 0.0, 0.0);
		int *indd = new int[tmp_col + 2];
		double *vall = new double[tmp_col + 2];
		indd[1] = i;
		vall[1] = 1.0;
		indd[2] = enu + i;
		vall[2] = -1.0;
		glp_set_mat_row(lp, tag, 2, indd, vall);
		delete[]indd, vall;
		tag++;
	}
	//cout << "debug5:  " << tag << endl;
	for (int i = 1; i <= q - k + 1; i++) {
		glp_set_row_bnds(lp, tag, GLP_FX, 0.0, 0.0);
		int *indd = new int[tmp_col + 2];
		double *vall = new double[tmp_col + 2];
		indd[1] = k * q - (i - 1);
		vall[1] = 1.0;
		indd[2] = tmp_col - (i - 1);
		vall[2] = -1.0;
		glp_set_mat_row(lp, tag, 2, indd, vall);
		delete[]indd, vall;
		tag++;
	}
	//cout << "debug6:  " << tag << endl;
	//add 1 rows
	int *indd = new int[tmp_col + 2];
	double *vall = new double[tmp_col + 2];

	glp_set_row_bnds(lp, tag, GLP_FX, 1.0, 1.0);
	for (int i = 1; i <= q - k + 1; i++) {
		indd[i] = enu + i;
		vall[i] = 1.0;
	}
	glp_set_mat_row(lp, tag, q - k + 1, indd, vall);

	//cout << "debug7:  " << tag << endl;
	//��ʼ��LP
	glp_init_smcp(&parm);
	parm.meth = GLP_PRIMAL;
	//����ʱ��
	clock_t startTime, endTime;
	startTime = clock();
	//����LP
	glp_simplex(lp, &parm);
	endTime = clock();
	cout << "Time usage to calculate stFlow by GLPK: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

	//��������gLP.txt
	glp_print_sol(lp, "stflow.txt");

	//���x��ȡֵ
	double *ZX = new double[tmp_col + 1];
	for (int i = 1; i <= tmp_col; i++) {
		double xn = glp_get_col_prim(lp, i);
		ZX[i] = xn;
		//ZX[i] = RRound(xn, 7);
///////////
		if (xn < 1.0e-10&&xn != 0.0)       //����С���ɺ��Ե�����ֱ����Ϊ0
			ZX[i] = 0.0;
	}

	//�жϽ��
	int ret = glp_get_status(lp);
	switch (ret) {
	case GLP_OPT:
		cout << "Success to get the st-flow!" << endl;
		break;
	default:
		//cout << "ret  " << ret << endl;
		//cout << "q:  " << q << endl;
		cout << "Failed to get the st-flow!" << endl;
		break;
	}

	cout << endl;

	//erase the memory
	glp_delete_prob(lp);
	delete[]indd;
	delete[]vall;

	return ZX;
}

//Compute an st-flow of value 1 in G by CPLEX
double *cstFlow(double *Z, int q, int k) {
	IloEnv env;
	int tmp_col = k * q + ((k - 1)*(q - k + 2)*(q - k + 1)) / 2 + 2 * (q - k + 1);
	int tmp_row = q + 2 * (k - 1)*(q - k + 1) + 1 + (q - k + 1) * 2;
	double *Zx = new double[tmp_col + 1];
	try {
		IloModel st(env);
		IloNumVarArray vars(env);
		//�����
		for (int i = 0; i < tmp_col; i++) {
			vars.add(IloNumVar(env, 0.0, 1.0));
		}
		st.add(vars);
		//���Ŀ�꺯��
		IloExpr exp(env);
		for (int i = 0, t = 1; i < k*q; i++, t++) {
			exp += t*vars[i];
			if (t > q)
				t = 1;
		}
		IloObjective obj(IloMaximize(env, exp, "cstFlow"));
		st.add(obj);
		//cout << "debug1" << endl;

		//���Լ��
		IloRangeArray cons(env);
		//ǰq��Լ��
		for (int i = 0; i < q; i++) {
			IloExpr texp(env);
			for (int j = 0; j < k; j++) {
				texp += vars[j*q + i];
			}
			cons.add(IloRange(env, Z[i + 1], texp, Z[i + 1]));
			texp.end();
		}
		//cout << "debug2" << endl;
		//���ߵ�(k - 1)*(q - k + 1)��Լ��
		for (int i = 0; i < k - 1; i++) {
			int tag = k * q + i * (q - k + 2)*(q - k + 1) / 2;   //��ǰ�е�һ�������ߵ����
			int tnum = q - k + 1;
			for (int j = i * q + i; j <= i * q + i + q - k; j++) {  //jΪ��ǰ�ߵ���ţ�jС�ڵ�ֵΪ��ǰ�п�ȡ�����һ���ߵ����
				IloExpr texp(env);
				for (int t = 0; t < tnum; t++) {
					texp += vars[tag];
					++tag;
				}
				texp -= vars[j];
				cons.add(IloRange(env, 0.0, texp, 0.0));
				texp.end();
				--tnum;
			}
		}
		
		//��ߵ�(k - 1)*(q - k + 1)��Լ��
		for (int i = 1; i <= k - 1; i++) {
			int ttag = k * q + i * (q - k + 2)*(q - k + 1) / 2;
			int tnum = q - k + 1;
			for (int j = i * q + i + q - k, tt = 1; j >= i * q + i; j--, tt++) {  //��ǰ�е����һ���߿�ʼ���Լ��
				int tag = ttag - (tt*(tt + 1) / 2);   //ǰһ����ñ����������һ�������ߵ����
				IloExpr texp(env);
				for (int t = 0; t < tnum; t++) {
					texp += vars[tag];
					tag = tag - (tt + t);
				}
				texp -= vars[j];
				cons.add(IloRange(env, 0.0, texp, 0.0));
				texp.end();
				--tnum;
			}
		}

		//���Դ��s��(q-k+1)��Լ��
		int tmpp = k * q + ((k - 1)*(q - k + 2)*(q - k + 1)) / 2;
		for (int i = 0; i < q - k + 1; i++) {
			cons.add(IloRange(env, 0.0, vars[i] - vars[tmpp + i], 0.0));
		}
		//��ӻ��t��(q-k+1)��Լ��
		for (int i = 0; i < q - k + 1; i++) {
			cons.add(IloRange(env, 0.0, vars[tmpp + q - k + 1 + i] - vars[tmpp - ((q - k + 1) - i)], 0.0));
		}
		
		//���Լ����ʹvalue=1
		IloExpr texp(env);
		for (int i = 0; i < q - k + 1; i++) {
			texp += vars[tmpp + i];
		}
		cons.add(IloRange(env, 1.0, texp, 1.0));
		texp.end();

		//���Լ����ʹÿһ��ֵΪ1
		for (int i = 0; i < k; i++) {
			IloExpr texp(env);
			for (int j = i * q + i; j < i*q + i + q - k + 1; j++) {
				texp += vars[j];
			}
			cons.add(IloRange(env, 1.0, texp, 1.0));
			texp.end();
		}

		st.add(cons);
		//cout << "debug7" << endl;

		IloCplex cplex(st);
		cplex.exportModel("cstFlow.lp");
		if (!cplex.solve()) {
			env.error() << "Failed to get stFlow." << endl;
			throw(-1);
		}
		IloNumArray vals(env);
		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value = " << cplex.getObjValue() << endl;
		cplex.getValues(vals, vars);
		for (int i = 1; i <= tmp_col; i++) {
			Zx[i] = vals[i - 1];
			/*if (Zx[i] < 1.0e-10&&Zx[i] != 0.0)       
				Zx[i] = 0.0;*/
		}
		//cout << "debug8" << endl;

	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}

	env.end();
	return Zx;
}

//get a collections of path-flows
//Ff��ÿһ��Ԫ�ش���һ��·����ÿһ��·���ĵ�һ��Ԫ�ش�������·���ĸ��ʣ�����Ԫ�ش���·���ıߵ����
set<set<double>> Ff(double *Zx, int q, int k) {
	set<set<double>> FF;
	random_device r;
	int tag;
	while (1) {
		set<double> f;   //��һ��ֵ�洢·��f�ĸ��ʣ������k��ֵ�洢k����
		double tmp = INFINITY;
		tag = 0;
		//�õ���С��Zx
		for (int i = 1; i <= k * q; i++) {
			if (Zx[i] < tmp&&Zx[i] != 0.0) {
				tmp = Zx[i];
				tag = i;
			}
		}//Zx[tag]Ϊ��С�ı�
		if (tag == 0)
			break;

		f.insert(Zx[tag]);    //�õ�·��f�ĸ���
		f.insert(tag);

		//cout << "tag:   " << tag << endl;
		int t = (tag - 1) / q + 1;  //��ǰZx�����·�����ڵ���
		//��ǰ��·��
		int tmp_h, tmpp, indd, numm, tmp_hx, tmp_p = tag;//tmp_h  tmp_hx  tmp_p
		for (int i = t - 1; i >= 1; i--) {  //iΪ��i��
			uniform_int_distribution<> u((i - 1)*q + i, tmp_p - q - 1);//��ǰһ�������ѡһ����
			while (1) {
				tmp_h = u(r);    //tmp_hΪǰһ��·������ţ�tagΪ��ʼ·�����
				if (Zx[tmp_h] == 0.0)
					continue;
				tmpp = (i - 1)*q + i;   //tmpp��¼ǰһ�е�һ���ߵı��
				indd = tmp_h - tmpp;    //indd��¼ǰһ�б�ѡ�б���ͬ�е�һ���ߵ���Ų�
				numm = k * q + ((i - 1)*(q - k + 2)*(q - k + 1)) / 2; //numm��¼��ѡ�еıߵĵ�һ������ǰһ���ߵ����
				for (int j = 0; j < indd; j++) {
					numm = numm + (q - k + 1) - j;
				}
				tmp_hx = numm + ((tmp_p - q) - tmp_h);//��¼���ӵ�ǰ�ߺ�ǰһ���ߵıߵ����
				if (Zx[tmp_hx] == 0.0)
					continue;
				break;
			}
			Zx[tmp_h] = Zx[tmp_h] - Zx[tag];
			///////////
			if (Zx[tmp_h] <= 1.0e-7)
				Zx[tmp_h] = 0.0;
			tmp_p = tmp_h;
			f.insert(tmp_h);
			cout << "head:   " << Zx[tmp_h] << "   " << tmp_h << endl;
		}
		//������·��
		int tmp_s = tag;
		for (int i = t + 1; i <= k; i++) {
			uniform_int_distribution<> u(tmp_s + q + 1, (i - 1)*q + (q - k + 1) + (i - 1));
			while (1) {
				tmp_h = u(r);
				if (Zx[tmp_h] == 0.0)
					continue;
				tmpp = (i - 2)*q + i - 1;   //tmpp��¼ǰһ�е�һ���ߵı��
				indd = tmp_s - tmpp;    //indd��¼ǰһ�б�ѡ�б���ͬ�е�һ���ߵ���Ų�
				numm = k * q + ((i - 2)*(q - k + 2)*(q - k + 1)) / 2; //numm��¼��ѡ�еıߵĵ�һ������ǰһ���ߵ����
				for (int j = 0; j < indd; j++) {
					numm = numm + (q - k + 1) - j;
				}
				tmp_hx = numm + ((tmp_h - q) - tmp_s);//��¼���ӵ�ǰ�ߺ�ǰһ���ߵıߵ����
				if (Zx[tmp_hx] == 0.0)
					continue;
				break;
			}
			Zx[tmp_h] = Zx[tmp_h] - Zx[tag];
			///////////
			if (Zx[tmp_h] <= 1.0e-7)
				Zx[tmp_h] = 0.0;
			tmp_s = tmp_h;
			f.insert(tmp_h);
			cout << "end:   " << Zx[tmp_h] << "  " << tmp_h << endl;
		}
		FF.insert(f);
		f.clear();
		Zx[tag] = 0.0;
		//break;
	}
	return FF;
}
