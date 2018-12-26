//
// MedUtils: Collection of several needed utilities (Files IO, Matrix data type, Data structures, etc)
//

// This is a simple h file to include all library parts

#ifndef __MED_UTILS_H__
#define __MED_UTILS_H__

//#include "MedGenUtils.h"
//#include "MedMat.h"
//#include "MedMedical.h"
//#include "MedDataStructures.h"
#include <Logger/Logger/Logger.h>
#include "assert.h"
#include "MedPlot.h"
#include <boost/program_options.hpp>
#include <boost/spirit/home/support/detail/hold_any.hpp>
#include <boost/algorithm/string.hpp>
#include <unordered_map>
#include <MedMat/MedMat/MedMatConstants.h>

using namespace std;

enum MedCorrelationType {
	CORR_PEARSON,
	CORR_LAST
};

enum MedBinningType {
	BIN_EQUIDIST,
	BIN_EQUISIZE,
	BIN_LAST
};

// Correlation of vectors
template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, int& n);
template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, int& n, float missing_value);
template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, int& n, MedCorrelationType type);
template <class S> int get_corr(vector<S>& _x, vector<S>& _y, double& corr, int& n, float missing_value, MedCorrelationType type);


// Pearson correlation of two vectors
template <class S> double get_corr_pearson(vector<S>& x, vector<S>& y);
double get_corr_pearson(float *v1, float *v2, int len); // A COPY OF THE TEMPLATED VERSION FOR C INTERFACING

template <class S> double get_squared_dist(vector<S>& x, vector<S>& y);
template <class S> double get_abs_avg_dist(vector<S>& x, vector<S>& y);
template <class S> double get_abs_relative_avg_dist(vector<S>& x, vector<S>& y);
template <class S> double get_vecs_accuracy(vector<S>& x, vector<S>& y, double epsilon);

// Discretization
template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins);
template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins, float missing_value);
template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins, MedBinningType binning);
template <class S> int discretize(vector<S>& x, vector<int>& binned_x, int& nbins, int max_bins, float missing_value, MedBinningType binning);

// Counts from binned vectors
void get_counts(vector<int>& x, map<int, int>& counts, int& n);
int get_co_counts(vector<int>& x, vector<int>& y, map<pair<int, int>, int>& counts, int& n);

// Mutual Information
double get_mutual_information(map<int, int>& x_count, int nx, map<int, int>& y_count, int ny, map<pair<int, int>, int>& xy_count, int n);
int get_mutual_information(vector<int>& x, vector<int>& y, double& mi, int &n);
double get_mutual_information(map<int, int>& x_count, int nx, map<int, int>& y_count, int ny, map<pair<int, int>, int>& xy_count, int n);

// Moments
//template <class S> int get_moments(vector<S>& values, vector<float>& wgts, S missing_value, S& mean, S& sd);
int get_moments(vector<float>& values, const vector<float>& wgts, float missing_value, float& mean, float& sd, bool do_missing = true);
int get_moments(float *values, const float *wgts, int n, float missing_value, float& mean, float& sd, bool do_missing = true);



namespace po = boost::program_options;

/**
* \brief medial namespace for function
*/
namespace medial {
	/*!
	*  \brief print
	*/
	namespace print {
		/// \brief general print to string woth format
		template<class T> string print_obj(T obj, const string &format);
		/// \brief printing vector elements in list [] with title to MLOG
		template<class T> void print_vec(const vector<T> &vec, const string &title, const string &format, const string &delimeter = ", ");
		/// \brief printing vector elements in list [] with title to MLOG
		void print_vec(const vector<string> &vec, const string &title, const string &delimeter = ", ");
		/// \brief printing vector elements hist prctiles in list [] with title to MLOG
		template<class T> void print_hist_vec(const vector<T> &vec, const string &title, const string &format, const vector<double> *prctile_samples = NULL);
		/// \brief print boost program options object
		string print_any(po::variable_value &a);

		void log_with_file(ofstream &fw, const char *format_str, ...);
	}
	/*!
	*  \brief process
	*/
	namespace process {
		/// \brief calc prctile
		template<class T> void prctils(const vector<T> &x, const vector<double> &prc, vector<T> &res, const vector<float> *weights = NULL);
		/// \brief binary search for index. -1 if not found
		template<typename T> int binary_search_index(const T *begin, const T *end, T val);
		/// \brief binary search for position to add new element in sorted manner (first position if equal elements found).
		template<typename T> int binary_search_position(const T *begin, const T *end, T val, bool reversed = false);
		/// \brief binary search for position to add new element in sorted manner (last position if equal elements found).
		template<typename T> int binary_search_position_last(const T *begin, const T *end, T val, bool reversed = false);
	}
	/*!
	*  \brief io
	*/
	namespace io {
		/// \brief reads file with codes name to vector
		void read_codes_file(const string &file_path, vector<string> &tokens);
		template<class T> string get_list(const unordered_map<string, T> &ls, const string &delimeter = ",");
		template<class T> string get_list_op(const unordered_map<T, string> &ls, const string &delimeter = ",");
		template<class ContainerType> string get_list(const ContainerType &ls, const string &delimeter = ",");
		/**
		* A basic class wrapper to parse command args
		* has default "h", "help", "debug" and "base_config" for reading all arguments from file
		* and prinring help. You just need to implement the Ctor of inheritence class and call init function
		* to use this class. you may also override post_process hook for setting some variable after all
		* arguments were set by the program_options. to use in main call "parse_parameters" function
		*/
		class ProgramArgs_base {
		private:
			string base_config; ///< config file with all arguments - in CMD we override those settings
			bool init_called; ///< mark for calling init function

			po::options_description desc; ///< the program_options object
			/// converts string arguments to enums if the program has some. 
			/// keep the raw string params from the user input as private and keep
			/// the Enum result of the converted as public.
			virtual void post_process() {};
			/// finds module section help in full help message
			string get_section(const string &full_help, const string &search);
		protected:
			/// an init function
			void init(po::options_description &prg_options, const string &app_l = "");
			/// the ctor of base class
			ProgramArgs_base() { init_called = false; debug = false; }
		public:
			bool debug; ///< a debug flag for verbose printing. will be init from command args
			string app_logo = "\
##     ## ######## ########  ####    ###    ## \n\
###   ### ##       ##     ##  ##    ## ##   ##    \n\
#### #### ##       ##     ##  ##   ##   ##  ##       \n\
## ### ## ######   ##     ##  ##  ##     ## ##       \n\
##     ## ##       ##     ##  ##  ######### ##       \n\
##     ## ##       ##     ##  ##  ##     ## ##       \n\
##     ## ######## ########  #### ##     ## ######## "; ///< the application logo/name



			/// the main function to parse the command arguments
			virtual int parse_parameters(int argc, char *argv[]);

		};
	}
}

#include "MedUtils_imp.h"

#endif