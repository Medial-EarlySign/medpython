#include "Utils.h"
#include <math.h>
#include <iostream>
#include <algorithm>

using namespace std;

struct val_idx {
	int idx;
	float val;
};


inline int auccompare( const void *arg1, const void *arg2 )
{
	//return((int)(10000*(nfc_lp[*(int *)arg2]-nfc_lp[*(int *)arg1])+0.5));

	if (((struct val_idx *)arg1)->val < ((struct val_idx *)arg2)->val) return 1;
	if (((struct val_idx *)arg1)->val > ((struct val_idx *)arg2)->val) return -1;
 
	return(0);
}

double Get_AUC(vector<float> &scores, vector<char> &Y, vector<double> &spe, vector<double> &sen)
{
	int i;

	int total0 = 0;
	int total1 = 0;
	int N = (int)scores.size();

	vector<val_idx> indexes(N);
	for(int i=0; i<N; i++) {
		if (Y[i]) total1++; else total0++;
		indexes[i].idx = i;
		indexes[i].val = (float)scores[i];
	}


	sort(indexes.begin(), indexes.end(), [](const val_idx &v1, const val_idx &v2){return v1.val > v2.val;});

	int current0 = total0;
	int sum = 0;

	for(i=0; i<N ; i++) {
		//cout << "i= " << i << " idx= " << indexes[i].idx << " val = " << indexes[i].val << " Y= " << Y[indexes[i].idx] << "\n";
		if (Y[indexes[i].idx]==1) 
			sum += current0; 
		else 
			current0--;
	}

	double auc = (double)sum/(double)(total0*total1);

	if (sen.size() > 0) {
		int t1[2] = {0,0};
		int t0[2] = {0,0};
		t0[0] = total0;
		t0[1] = total1;
		double spec,sens;
		spe.resize(sen.size());
		int k = 0;
		for (i=0; i<N && k<sen.size(); i++) {
			//cout << "i= " << i << " idx= " << indexes[i].idx << " val = " << indexes[i].val << " Y= " << Y[indexes[i].idx] << "\n";
			
			if (Y[indexes[i].idx] == 1) {
				t1[1]++;
				t0[1]--;
			} else {
				t1[0]++;
				t0[0]--;
			}
			spec = (double)t0[0]/(double)total0;
			sens = (double)t1[1]/(double)total1;
			if (sens >= sen[k]) {
				spe[k] = spec;
				sen[k] = sens;
				k++;
			}
		}
		
	}


	return auc;
}


double Get_AUC(vector<float> &scores, vector<char> &Y)
{
	vector<double> a,b;
	return( Get_AUC(scores, Y, a, b) );
}
