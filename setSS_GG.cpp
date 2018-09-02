#include"func.h"

using namespace std;

//����N��Ȩ��,��ΧΪ1-10
double *rand_weight(int n, int m) {
	/*
	weight:���Ȩ�أ���weight[1]��ʼ�洢����ΧΪ1-10��
	*/
	double *weight = new double[m + n + 1];
	random_device r;
	uniform_real_distribution<double> u(1.0, 10.0);
	for (int i = 1; i <= n; i++) {
		weight[i] = u(r);
	}
	//��m��weightΪ����Ȩ��
	uniform_real_distribution<> v(0.000001, 0.000009);
	for (int i = n + 1; i <= m + n; i++) {
		weight[i] = v(r);
	}
	return weight;
}

//���ɼ���SS������m��S,ÿ��S����f��Ԫ��,
set<set<int>> randSS(int n, int k) {
	/*
	parm:n  ����U��Ԫ�أ�
	parm:k  ���ɵ�SS��Ԫ�ز�����k��
	parm:m  ���ɵ�SS��Ԫ��������n<=m<=2*n;
	�㷨�����������m��
	����SS��ÿһ��Ԫ�أ�������һ�������f��f<0.3*n��
	��n�������ѡf�������������S��S����f��Ԫ�أ�
	��SS�в�����S����S����SS��
	ѭ����֪��SS����m��S��
	*/
	random_device r;
	uniform_int_distribution<unsigned> u(n, 2 * n);
	int m = u(r);  //m�����ϣ�m<=n*3,m>=k
	int f;    //���ÿ��������Ԫ�صĸ���
	set<int> S;
	set<set<int>> SS;
	int tag = 0;

	for (int i = 1; i <= m; i++) {
		while (1) {
			f = r() % (int(0.3*n));//����S�ķ�Χ��S��ຬ��0.3*n��Ԫ��
			if (f != 0)break;
		}

		//��n��������ѡf����ͬ����
		for (int j = 1; j <= n; j++) {
			if (r() % (n + 1 - j) < unsigned(f)) {
				S.insert(j);
				f--;
			}
		}
		/*
		n+iΪ����������Ϊ������GLPK���ʱ�ܰ�ÿ��S������
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

//���������q,��С��k��������m
int randq(set<set<int>> SS, int k) {
	int num = SS.size();
	int q, t = int(0.4*num);
	random_device r;
	uniform_int_distribution<unsigned> u(k, t);
	while (1) {
		q = u(r);//���ڲ���k-num֮��������
		if (q >= k)
			break;
	}

	//cout << "partation num:" << q << endl;
	return q;
}

//����SS�Ļ���GG,����q�����֣�q��СΪk
set<set<int>> *randGG(set<set<int>> SS, int qq) {
	int num = SS.size();
	int q = qq;
	cout << "the size of GG:" << q << endl;

	//����q�����ּ���Ӧ��S����
	set<int>tmp;
	tmp.insert(0);
	int tag = 1, m = q;
	random_device r;
	//����q-1���м������,��β�ڲ���0��num
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

	//����GG[],1~q+1
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