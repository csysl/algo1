#include"func.h"

using namespace std;

//生成N个权重,范围为1-10
double *rand_weight(int n, int m) {
	/*
	weight:存放权重，从weight[1]开始存储，范围为1-10；
	*/
	double *weight = new double[m + n + 1];
	random_device r;
	uniform_real_distribution<double> u(1.0, 10.0);
	for (int i = 1; i <= n; i++) {
		weight[i] = u(r);
	}
	//后m个weight为辅助权重
	uniform_real_distribution<> v(0.000001, 0.000009);
	for (int i = n + 1; i <= m + n; i++) {
		weight[i] = v(r);
	}
	return weight;
}

//生成集合SS，含有m个S,每个S含有f个元素,
set<set<int>> randSS(int n, int k) {
	/*
	parm:n  集合U的元素；
	parm:k  生成的SS的元素不少于k；
	parm:m  生成的SS的元素数量，n<=m<=2*n;
	算法：生成随机数m；
	对于SS的每一个元素，先生成一个随机数f，f<0.3*n；
	从n中随机挑选f个随机数，生成S，S含有f个元素；
	若SS中不存在S，将S插入SS；
	循环，知道SS含有m个S。
	*/
	random_device r;
	uniform_int_distribution<unsigned> u(n, 2 * n);
	int m = u(r);  //m个集合，m<=n*3,m>=k
	int f;    //标记每个集合中元素的个数
	set<int> S;
	set<set<int>> SS;
	int tag = 0;

	for (int i = 1; i <= m; i++) {
		while (1) {
			f = r() % (int(0.3*n));//限制S的范围，S最多含有0.3*n个元素
			if (f != 0)break;
		}

		//从n个数中挑选f个不同的数
		for (int j = 1; j <= n; j++) {
			if (r() % (n + 1 - j) < unsigned(f)) {
				S.insert(j);
				f--;
			}
		}
		/*
		n+i为辅助参数，为了再用GLPK求解时能把每个S都计算
		*/
		S.insert(n + i);
		SS.insert(S);
		tag++;
		if (tag != SS.size()) {
			tag = SS.size();
			m++;
		}
		S.clear();
	}
	cout << "the size of SS:" << SS.size() << endl;
	return SS;
}

//产生随机数q,不小于k，不大于m
int randq(set<set<int>> SS, int k) {
	int num = SS.size();
	int q, t = int(0.4*num);
	random_device r;
	uniform_int_distribution<unsigned> u(k, t);
	while (1) {
		q = u(r);//用于产生k-num之间的随机数
		if (q >= k)
			break;
	}

	//cout << "partation num:" << q << endl;
	return q;
}

//生成SS的划分GG,含有q个划分，q最小为k
set<set<int>> *randGG(set<set<int>> SS, int qq) {
	int num = SS.size();
	int q = qq;
	cout << "the size of GG:" << q << endl;

	//产生q个划分集对应的S个数
	set<int>tmp;
	tmp.insert(0);
	int tag = 1, m = q;
	random_device r;
	//生成q-1个中间随机数,首尾在插入0和num
	for (int i = 0; i < m - 1; i++) {
		int tmpa = r() % num;
		tmp.insert(tmpa);
		tag++;
		if (tag != tmp.size()) {
			tag--;
			m++;
		}

	}
	tmp.insert(num);
	//cout << "size of tmp:" << tmp.size() << endl;

	//生成GG[],1~q+1
	set<set<int>> *GG = new set<set<int>>[q + 1];
	auto it = tmp.begin();
	auto it1 = SS.begin();
	int tmp1 = 0, tmp2, tmp3;
	for (int i = 1; i <= q; i++) {
		if (it == tmp.end())break;
		tmp1 = *it;
		it++;
		tmp2 = *it;
		tmp3 = tmp2 - tmp1;
		//cout << "tmp3:" << tmp3 << "   tmp2:" << tmp2 << endl;
		for (int j = 0; j < tmp3; j++) {
			GG[i].insert(*it1);
			it1++;
		}
		//cout << GG[i].size() << "   " << i << endl;
	}
	return GG;
}