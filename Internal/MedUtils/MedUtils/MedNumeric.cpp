//
// MeddNumeric
//

#include "MedNumeric.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/math/distributions/skew_normal.hpp>

#include "Logger/Logger/Logger.h"
#define LOCAL_SECTION LOG_MEDMAT
#define LOCAL_LEVEL	LOG_DEF_LEVEL
extern MedLogger global_logger;

//...............................................................................................................
double square_dist_for_fit(float (*f_to_fit)(vector<float> &p, float xp), vector<float> &p, vector<float> &x, vector<float> &y, vector<float> &w)
{
	double err = 0;
	float d;
	for (int i=0; i<x.size(); i++) {
		float yp = (*f_to_fit)(p,x[i]);
		if (w.size() > 0)
			d = w[i]*(yp - y[i]);
		else
			d = (yp - y[i]);
		err += (double)(d*d);
	}

	return err;
}

//..............................................................................................................
// finding the p[0] scale factor to get the best fit
int get_fit_scale_factor(float (*f_to_fit)(vector<float> &p, float xp), vector<float> &p, vector<float> &x, vector<float> &y, vector<float> &w)
{
	double sxy,sxx;

	sxy = 0;
	sxx = 0;
	p[0] = 1;
	for (int i=0; i<x.size(); i++) {
		float yp = (*f_to_fit)(p,x[i]);
		sxy += (double)w[i]*yp*y[i];
		sxx += (double)w[i]*yp*yp;
	}
	if (sxx > 0)
		p[0] = (float)(sxy/sxx);
	else {
	}

	return 0;
}
//...............................................................................................................
// p[0] is ALWAYS a scaling factor that is fitted for each enumarated option
int full_enumaration_fit(float (*f_to_fit)(vector<float> &p, float xp), MedFitParams &params, vector<float> &x, vector<float> &y, vector<float> &w, vector<float> &res, double &res_err)
{
	int i,j;

	//MLOG("Starting fit of %d parameters\n",params.n_params);
	//for (i=0; i<params.n_params; i++) 
	//	MLOG("Param %d :: range %6.3f - %6.3f :: jump %6.4f\n",i,params.range[i].first,params.range[i].second,params.jump[i]);

	if (w.size() == 0) {
		w.resize(x.size());
		fill(w.begin(),w.end(),(float)1);
	}

	// initialize with minimal values
	vector<float> p(params.n_params);
	for (i=0; i<params.n_params; i++)
		p[i] = params.range[i].first;

	double best_err = square_dist_for_fit(f_to_fit,p,x,y,w);
	//for(i=0; i<x.size(); i++) {
	//	float yp = (*f_to_fit)(p,x[i]);
	//	MLOG("Fitting i=%d x %f y %f w %f yp %f\n",i,x[i],y[i],w[i],yp);
	//}

	vector<float> bestp = p;

	int go_on = 1;

	float epsilon = (float)1e-5;

	while (go_on) {

		// advance p
		int i=1; //p[0] will be estimated later as a scaling factor
		while (i<p.size()) {
			if (p[i] < params.range[i].second) {
				p[i] += max(epsilon,params.jump[i]);
				for (j=0; j<i; j++)
					p[j] = params.range[j].first;
				break;
			} else
				i++;
		}

		if (i == p.size()) go_on = 0;

		// test p
		get_fit_scale_factor(f_to_fit,p,x,y,w);
		double err = square_dist_for_fit(f_to_fit,p,x,y,w);
		if (err < best_err) {
			best_err = err;
			bestp = p;
			//MLOG("Best so far: p: %f %f %f , err %f\n",p[0],p[1],p[2],best_err);
		}
	}


	res_err = best_err;
	res = bestp;

	//MLOG("Best so far: p: %f %f %f , err %f\n",res[0],res[1],res[2],best_err);
	return 0;
}

float normal_cdf(float mean, float std, float x)
{
	boost::math::normal_distribution<float> nd(mean,std);
	return (boost::math::cdf(nd,x));
}

float log_normal_cdf(float mean, float std, float x)
{
	boost::math::lognormal_distribution<float> nd(mean,std);
	return (boost::math::cdf(nd,x));
}

//...............................................................................................................
// p[0] - scale p[1] - mean , p[2] - std , p[3] width - we want the probability to fall in [xp-p[3],xp+p[3]]
float scaled_normal_dist(vector<float> &p, float xp)
{

	boost::math::normal_distribution<float> nd(p[1],p[2]);

//	return(p[2]*boost::math::pdf(nd,xp));
	return(p[0]*(boost::math::cdf(nd,xp+p[3]) - boost::math::cdf(nd,xp-p[3])));

}

//...............................................................................................................
// p[0] - scale , p[1] - mean , p[2] - std , p[3] - shape , p[4]-width
float scaled_skewed_normal_dist(vector<float> &p, float xp)
{

	boost::math::skew_normal_distribution<double> nd(p[1],p[2],p[3]);

	return((float)((double)p[0]*(boost::math::cdf(nd,xp+p[4]) - boost::math::cdf(nd,xp-p[4]))));

}

//...............................................................................................................
float scaled_log_normal_dist(vector<float> &p, float xp)
{

	boost::math::lognormal_distribution<double> nd(p[1],p[2]);

	return((float)((double)p[0]*(boost::math::cdf(nd,xp+p[3]) - boost::math::cdf(nd,xp-p[3]))));

}

//...............................................................................................................
float scaled_inv_gauss_dist(vector<float> &p, float xp)
{

	boost::math::inverse_gaussian_distribution<double> nd(p[1],p[2]);

	return((float)((double)p[0]*(boost::math::cdf(nd,xp+p[3]) - boost::math::cdf(nd,xp-p[3]))));

}

//...............................................................................................................
int get_normal_dist_quantiles(float mean, float sd, vector<float> &p, vector<float> &v)
{
	boost::math::normal_distribution<float> nd(mean,sd);

	v.resize(p.size());
	for (int i=0; i<p.size(); i++)
		v[i] = boost::math::quantile(nd,p[i]);

	return 0;
}

//...............................................................................................................
int get_skewed_normal_dist_quantiles(float mean, float sd, float shape, vector<float> &p, vector<float> &v)
{
	boost::math::skew_normal_distribution<float> nd(mean,sd,shape);

	v.resize(p.size());
	for (int i=0; i<p.size(); i++)
		v[i] = boost::math::quantile(nd,p[i]);

	return 0;
}

//...............................................................................................................
int get_log_normal_dist_quantiles(float mean, float sd, vector<float> &p, vector<float> &v)
{
	boost::math::lognormal_distribution<float> nd(mean,sd);

	v.resize(p.size());
	for (int i=0; i<p.size(); i++)
		v[i] = boost::math::quantile(nd,p[i]);

	return 0;
}

//...............................................................................................................
int get_inv_gauss_dist_quantiles(float mean, float sd, vector<float> &p, vector<float> &v)
{
	boost::math::inverse_gaussian_distribution<float> nd(mean,sd);

	v.resize(p.size());
	for (int i=0; i<p.size(); i++)
		v[i] = boost::math::quantile(nd,p[i]);

	return 0;
}




