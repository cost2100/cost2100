/*!
 * \file
 * \addtogroup COST2100
 * \brief COST2100 Channel Model Specifications
 * \author K.K.
 */

#ifndef COST2100_SPECIFICAION_H_
#define COST2100_SPECIFICAION_H_

#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <iostream>
#include <string>
#include <libxml/parser.h>
#include <vector>

#include "debug.h"
#include "mat3.h"

#define DEBUG
#undef DBG_LVL
#define DBG_LVL (DBG_ERR|DBG_INFO)//DBG_TRACE|DBG_WARN|DBG_INFO|DBG_ENTER|DBG_LEAVE)

namespace penux {

//! \addtogroup COST2100
//@{

using namespace itpp;
using std::complex;
using std::string;

//! Predefined COST2100 channel model profiles: including Macrocell, Microcell, Picocell, AALTO or Customized Test Scenarios
enum CHANNEL_PROFILE {
	COST2100_MACRO, COST2100_MICRO, COST2100_PICO, COST2100_AALTO, COST2100_TEST
};

//! LoS type: Line-of-Sight or Non-Line-of-Sight
enum LINE_OF_SIGHT {
	LOS, NLOS
};
//! Side type: BS Side, MS Side or Both BS/MS Sides
enum SIDE_TYPE {
	BS, MS, BS_MS
};
//! Cluster type: Single Cluster, Twin Cluster, Local BS Cluster or Local MS Cluster
enum CLUSTER_TYPE {
	SINGLE, TWIN, LOCAL_AT_BS, LOCAL_AT_MS
};
//! Spread type: Delay Spread, AoD Angular Spread or EoD Angular Spread to BS/MS side
enum SPREAD_TYPE {
	DELAY_AT_BS, DELAY_AT_MS, AOD_AT_BS, AOD_AT_MS, EOD_AT_BS, EOD_AT_MS
};
//! Angle type: Azimuth or Elevation Angle at BS/MS side
enum ANGLE_TYPE {
	AZIMUTH_AT_BS, AZIMUTH_AT_MS, ELEVATION_AT_BS, ELEVATION_AT_MS
};
//! Coordinate-Transformation type: Cartesian to Spherical,  Spherical to Cartesian, Polar to Cartesian or Cartesian to Polar
enum COOR_TRANS_TYPE {
	CART2SPH, SPH2CART, POL2CART, CART2POL
};
//! Distribution type: Normal Distribution, Uniform Distribution or Poisson Distribution
enum PDF_TYPE {
	NORMAL, UNIFORM, POISSON
};
//! Channel band type: Wide-band or Narrow-band
enum BAND_TYPE {
	WIDEBAND, NARROWBAND
};

/*!
 * \brief MPC(Multiple Paths Component) Class
 *
 * MPC Class stands for a set of MPCs occur in one cluster.
 *
 * Each MPC has positions with reference to both BS/MS sides, and its amplitude is calculated as Rayleigh distributed.
 */
class MPC {
public:
	//! Empty constructor
	MPC() {}
	//! Default constructor
	MPC(const mat &mpc_bs_pos, const mat &mpc_ms_pos, const cvec &mpc_amp) { MPC_BS_pos = mpc_bs_pos; MPC_MS_pos = mpc_ms_pos; MPC_amp = mpc_amp; }
	//! Destructor
	virtual ~MPC() {}

	//! Update single MPC position to BS/MS side
	void update_MPC_pos(const SIDE_TYPE sd_type, int mpc_idx, const vec pos);

	//! Get MPC positions to BS/MS side
	mat get_MPC_pos(const SIDE_TYPE sd_type) { switch (sd_type) { case BS: return MPC_BS_pos; case MS: return MPC_MS_pos; default: return mat(0); }	}
	//! Get MPC amplitudes
	cvec get_MPC_amplitude() { return MPC_amp; }

protected:
	mat MPC_BS_pos;							//!< MPCs positions to BS side
	mat MPC_MS_pos;							//!< MPCs Generatepositions to MS side
	cvec MPC_amp;							//!< MPCs amplitudes
};

/*!
 * \brief Cluster Class
 *
 * For clusters, there can be three types: local cluster, single cluster or twin clusters.
 *
 * Each cluster has a set of delay and angular spread parameters, together with its link delay,
 * shadow fading and all the MPCs occur in the cluster.
 */
class Cluster {
public:
	//! Empty constructor
	Cluster() {}
	//! Default constructor
	Cluster(const CLUSTER_TYPE c_type, const vec &c_bs_pos, const vec &c_ms_pos, const vec &spread, double c_sf, double c_tau_link, const vec &angle, const MPC mpc);
	//! Destructor
	virtual ~Cluster() {}

	//! Update cluster MPCs
	void update_cluster_MPC(const MPC mpc) { C_MPC = mpc; }

	//! Get cluster type
	CLUSTER_TYPE get_cluster_type() { return C_type; };
	//! Get cluster positions to BS/MS side
	vec get_cluster_pos(const SIDE_TYPE sd_type) { switch (sd_type) { case BS: return C_BS_pos;	case MS: return C_MS_pos; default: return zeros(3); } }
	//! Get cluster shadow fading
	double get_cluster_shadowing_fading() { return C_shadow_f; }
	//! Get cluster link delay
	double get_cluster_link_delay() { return C_tau_link; }
	//! Get cluster MPCs
	MPC get_cluster_MPC() { return C_MPC; }
	//! Get cluster spread
	double get_cluster_spread(const SPREAD_TYPE s_type);
	//! Get cluster spread vector to BS/MS side
	vec get_cluster_spread_vec(const SIDE_TYPE sd_type);
	//! Get cluster angle
	double get_cluster_angle(const ANGLE_TYPE a_type);
	//! Get cluster angle vector to BS/MS side
	vec get_cluster_angle_vec(const SIDE_TYPE sd_type);

protected:
	CLUSTER_TYPE C_type;					//!< cluster type
	vec C_BS_pos;							//!< cluster position to BS side
	vec C_MS_pos;							//!< cluster position to MS side
	double C_delay_BS;						//!< cluster spatial delay spread to BS side
	double C_delay_MS;						//!< cluster spatial delay spread to MS side
	double C_AoD_BS;						//!< cluster AoD to BS side
	double C_AoD_MS;						//!< cluster AoD to MS side
	double C_EoD_BS;						//!< cluster EoD to BS side
	double C_EoD_MS;						//!< cluster EoD to MS side
	double C_shadow_f;						//!< cluster shadowing fading.
	double C_tau_link;						//!< cluster link delay.
	double C_Phi_BS;						//!< cluster azimuth angle to BS side
	double C_Theta_BS;						//!< cluster elevation angle to BS side
	double C_Phi_MS;						//!< cluster azimuth angle to MS side
	double C_Theta_MS;						//!< cluster elevation angle to MS side
	MPC C_MPC;								//!< cluster MPCs
};

/*!
 * \brief BS Information Class
 *
 * BS_Info Class stands for the information of all BSs in the simulation case.
 *
 * BS Information contains position, antenna, common clusters,
 * visibility regions and local cluster information of every BSs.
 */
class BS_Info {
public:
	//! Empty constructor
	BS_Info() {}
	//! Default constructor
	BS_Info(const mat &pos, const vec &ratio, const vec &antenna) { set_BS_pos(pos); set_BS_common_ratio(ratio); set_BS_antenna_num(antenna); BS_num = pos.rows(); }
	//! Destructor
	virtual ~BS_Info() {}

	//! Set BS positions
	void set_BS_pos(const mat &pos) { if (pos.cols() != 3) ERROR("BS_Info::set_BS_pos(): invalid BS position dimension!"); else BS_pos = pos; }
	//! Set BS common cluster ratios
	void set_BS_common_ratio(const vec &ratio) { common_ratio = ratio; }
	//! Set BS antenna numbers
	void set_BS_antenna_num(const vec &antenna) { antenna_num = antenna; }
	//! Set BS visibility regions
	void set_BS_visibility_region(const mat &vr, const mat &vr_los) { VR = vr; VR_LOS = vr_los; }
	//! Set BS local clusters
	void set_BS_local_cluster(const Array<Cluster> &cluster) { local_cluster = cluster; }

	//! Get BS numbers
	int get_BS_num() { return BS_num; }
	//! Get BS positions
	mat get_BS_pos() { return BS_pos; }
	//! Get BS common cluster ratios
	vec get_BS_common_ratio() { return common_ratio; }
	//! Get BS antenna numbers
	vec get_BS_antenna_num() { return antenna_num; }
	//! Get BS visibility regions
	mat get_VR() { return VR; }
	//! Get BS LoS visibility regions
	mat get_VR_LOS() { return VR_LOS; }
	//! Get BS local clusters
	Array<Cluster> get_BS_local_cluster() { return local_cluster; }

protected:
	int BS_num;								//!< BS number
	mat BS_pos;								//!< BS positions (x, y, z)				size: BS_num x 3
	vec common_ratio;						//!< BS common ratios					size: BS_num
	vec antenna_num;						//!< BS antenna number					size: BS_num
	mat VR;									//!< (BS)visibility regions (x, y)		size: VR_num x 2
	mat VR_LOS;								//!< (BS)LOS visibility regions (x, y)	size: BS_num x 2
	Array<Cluster> local_cluster;			//!< BS local clusters
};

/*!
 * \brief MS Information Class
 *
 * MS_Info Class stands for the information of all MSs in the simulation case.
 *
 * MS Information contains position, velocity, antenna, common clusters,
 * and local cluster information of every MSs.
 */
class MS_Info {
public:
	//! Empty constructor
	MS_Info() {}
	//! Default constructor
	MS_Info(const mat &pos, const mat &velo, const double mscc, const vec &antenna) { set_MS_pos(pos); set_MS_velo(velo); set_MS_common(mscc); set_MS_antenna_num(antenna); MS_num = pos.rows();}
	//! Destructor
	virtual ~MS_Info() {}

	//! Set MS positions
	void set_MS_pos(const mat &pos) { if (pos.cols() != 3) ERROR("MS_Info::set_MS_pos(): invalid MS position dimension!"); else MS_pos = pos; }
	//! Set MS velocities
	void set_MS_velo(const mat &velo) { if (velo.cols() != 3) ERROR("MS_Info::set_MS_velo(): invalid MS velocity dimension!"); else MS_velo = velo; }
	//! Set MS common cluster groups
	void set_MS_common(const double mscc) { MS_common = mscc; }
	//! Set MS antenna numbers
	void set_MS_antenna_num(const vec &antenna) { antenna_num = antenna; }
	//! Set MS local clusters
	void set_MS_local_cluster(const Array<Cluster> &cluster) { local_cluster = cluster; }
	//! Update specified MS local cluster
	void update_MS_local_cluster(const int idx, const Cluster &cluster) { local_cluster(idx) = cluster; }
	//! Update specified MS position
	void update_MS_pos(const mat &pos) { set_MS_pos(pos); }

	//! Get MS numbers
	int get_MS_num() { return MS_num; }
	//! Get MS positions
	mat get_MS_pos() { return MS_pos; }
	//! Get MS velocities
	mat get_MS_velo() { return MS_velo; }
	//! Get MS common cluster groups
	double get_MS_common() { return MS_common; }
	//! Get MS antenna numbers
	vec get_MS_antenna_num() { return antenna_num; }
	//! Get MS local clusters
	Array<Cluster> get_MS_local_cluster() { return local_cluster; }

protected:
	int MS_num;								//!< MS number
	mat MS_pos;								//!< MS positions						size: MS_num x 3
	mat MS_velo;							//!< MS velocities						size: MS_num x 3
	double MS_common;						//!< MS group common clusters
	vec antenna_num;						//!< MS antenna number					size: MS_num
	Array<Cluster> local_cluster;			//!< MS local clusters
};

/*
class Impulse_Response {
public:
	Impulse_Response() {}
	Impulse_Response(int bs, int ms, cmat3 &h) { set_bs_idx(bs); set_ms_idx(ms); set_impulse_response(h); }
	virtual ~Impulse_Response() {}

	void set_bs_idx(int bs) { bs_idx = bs; }
	void set_ms_idx(int ms) { ms_idx = ms; }
	void set_impulse_response(cmat3 &h) { impulse = h; }

	int get_bs_idx() { return bs_idx; }
	int get_ms_idx() { return ms_idx; }
	cmat3 get_impulse_response() { return impulse; }

protected:
	int bs_idx,ms_idx;
	cmat3 impulse;
};
*/

/*!
 * \brief Transfer Function Class
 *
 * Transfer functions of a MIMO channel link between a pair of BS and MS at each frequency bin.
 */
class Transfer_Function {
public:
	//! Empty constructor
	Transfer_Function() {}
	//! Default constructor
	Transfer_Function(int bs, int ms, vec &f, Array<cmat> trans) { set_bs_idx(bs); set_ms_idx(ms); set_frequency(f); set_transfer_function(trans); }
	//! Destructor
	virtual ~Transfer_Function() {}

	//! Set BS index
	void set_bs_idx(int bs) { bs_idx = bs; }
	//! Set MS index
	void set_ms_idx(int ms) { ms_idx = ms; }
	//! Set frequency bins
	void set_frequency(vec &f) { freq = f; }
	//! Set channel matrices
	void set_transfer_function(Array<cmat> trans) { H = trans; }

	//! Get BS index
	int get_bs_idx() { return bs_idx; }
	//! Get MS index
	int get_ms_idx() { return ms_idx; }
	//! Get frequency bins
	vec get_frequency() { return freq; }
	//! Get channel matrices
	Array<cmat> get_transfer_function() { return H; }

protected:
	int bs_idx;								//!< BS index
	int ms_idx;								//!< MS index
	Array<cmat> H;							//!< array of channel matrices
	vec freq;								//!< frequency vector
};

/*!
 * \brief COST2100 Channel Specification Class
 *
 * Initialize and simulate the COST2100 channel model specifications with the parameters imported from configuration files.
 */
class Channel_Specification {
public:
	//! Empty constructor
	Channel_Specification() {}
	//! Default constructor
	Channel_Specification(string fileName);
	//! Destructor
	virtual ~Channel_Specification() {}

	//! Get parameters from the configuration file
	void getParamsFromXML(string fileName);

	//!Initialize channel model
	void init_channel();
	//! Update MS information
	void update_MS_info();

	//! Generate BS common cluster assignment
	void init_BSCC();
	//! Generate MS common cluster assignment
	void init_MSCC();
	//! Generate visibility regions
	void init_VR();
	//! Generate clusters
	void init_cluster();
	//! Initialize cluster random parameters
	void init_cluster_params(int num_in, mat &corr_randn_out, vec &shadow_factor_out, vec &tau_C_out, vec &theta_C_BS_out, vec &phi_C_BS_out,
							vec &theta_C_MS_out, vec &phi_C_MS_out, vec &d_tau_out);
	//! Generate far clusters
	Array<Cluster> init_far_clusters();
	//! Generate local clusters
	Array<Cluster> init_local_clusters(SIDE_TYPE sd_type);
	//! Get reference BS and visibility region positions for cluster
	void get_ref_BS_VR(int cluster_idx, vec &ref_bs, vec &ref_vr);
	//! Generate MPCs for cluster
	MPC get_MPC(const CLUSTER_TYPE c_type, const vec &c_bs_pos, const vec &c_ms_pos, const vec& spread, const vec &bs_angle, const vec &ms_angle);

	//! Get all channel transfer functions of the simulation
	vector<vector<Transfer_Function> > get_all_transfer_functions() { return transfer_function; }
	//! Get all channel clusters of the simulation
	Array<Cluster> get_all_channel_cluster() { return channel_cluster; }

	//! Get frequency bins
	int get_freq_bins() { return N_freq; }
	//! Get snapshot numbers
	int get_snap_num() { return snapNum; }
	//! Get TX/RX antenna number vector
	vec get_antenna_nums(int bs_idx, int ms_idx);
	//!	Get BS/MS number vector
	vec get_BS_MS_nums();

	//! Calculate channel path loss
	double calc_path_loss(const double dist);

	void get_DMC() {}

protected:
	CHANNEL_PROFILE profile;				//!< channel profile

	// Channel parameters
	// External parameters of all scenarios
	double freq_start;						//!< start frequency			(Hz)
	double freq_stop;						//!< stop frequency 			(Hz)
	double freq_div;						//!< frequency division			(Hz)
	double bandwidth;						//!< bandwidth					(Hz)
	double freq_c;							//!< center frequency			(Hz)
	int N_freq;								//!< number of frequency bins
	double cellRadius;						//!< cell radius				(m)
	double delay_max;						//!< maximum delay				(sec)
	double tau_idx_max;						//!< maximum delay bin index
	double snapRate;						//!< snapshot rate
	int snapNum;							//!< number of snapshots
	double sampleRate;						//!< sample rate
	double overSampleRate;					//!< oversample rate
	BAND_TYPE bandType;						//!< channel band type

	// External parameters of macro and micro cells
	double h_rooftop;						//!< rooftop height				(m)
	double w_road;							//!< road width					(m)
	double w_street;						//!< street width				(m) or building separation
	double phi_road;						//!< road orientation 			(degree)
	double h_BS;							//!< BS height					(m)
	double h_BS_min;						//!< minimum BS height			(m)
	double h_BS_max;						//!< maximum BS height			(m)
	double h_MS;							//!< MS height					(m)

	// External parameters of picocells
	double l_room;							//!< room length				(m)
	double w_room;							//!< room width					(m)
	double n_floor;							//!< floors between BS and MS

	// Stochastic parameters of all scenariosBS
	double R_C;								//!< visibility region radius	(m)
	double L_C;								//!< transition region width	(m)
	double k_tau;							//!< cluster power				(dB/us)
	double tau_B;							//!< excess delay				(s)
	// LoS occurrence
	double d_co;							//!< LoS cut-off distance		(m)
	double R_L;								//!< LoS VR radius				(m)
	double L_L;								//!< LoS TR width				(m)
	double EPL;								//!< Excess path loss
	double mu_K;							//!< mean of LoS power factor	(dB)
	double sigma_K;							//!< std. of LoS power factor	(dB)
	double power_factor;					//!< LoS power factor
	double K_sel; 							//!< rate of single interacting clusters
	double N_C_local;						//!< (average) number of local clusters, assume = 1
	double mu_N_C_add;						//!< mean of additional(far) clusters, poisson distributed
	int N_C_add;							//!< poisson distributed number of additional(far) clusters
	int N_C_far;							//!< number of additional(far) clusters for simulation
	int N_C_single;							//!< number of single clusters for simulation
	int N_C_twin;							//!< number of twin clusters for simulation
	bool BS_local;							//!< activeness of BS local clusters
	bool MS_local;							//!< activeness of MS local clusters
	double N_MPC;							//!< number of MPCs per cluster
	double K_MPC;							//!< Rice factor of additional clusters
	double mu_diff;							//!< mean of diffuse radiation
	double sigma_diff;						//!< std. of diffuse radiation
	double PDP_diff;						//!< PDP of diffuse radiation, uniform in azimuth exp(-t/tau), tau = 0.5us
	// dealy spread
	double mu_tau;							//!< mean of delay spread		(us)
	double sigma_tau;						//!< std. of delay spread		(dB)
	// angular spread
	double mu_phi_BS;						//!< mean of AoD spread			(degree)
	double sigma_phi_BS;					//!< std. of AoD spread			(degree)
	double mu_theta_BS;					Cluster	//!< mean of EoD spread			(degree)
	double sigma_theta_BS;					//!< std. of EoD spread			(degree)
	double mu_phi_MS;						//!< mean of AoA spread			(degree)
	double sigma_phi_MS;					//!< std. of AoA spread			(degree)
	double mu_theta_MS;						//!< mean of EoA spread			(degree)
	double sigma_theta_MS;					//!< std. of EoA spread			(degree)
	double theta_MS_min;					//!< minimum EoA spread(uniform)(degree)
	double theta_MS_max;					//!< maximum EoA spread(uniform)(degree)
	double sigma_sf;						//!< std. of shadow fading 		(dB)
	// autocorrelation distances (m)
	double L_tau;							//!< autocorrelation of delay spread
	double L_phi_BS;						//!< autocorrelation of AoD spread
	double L_theta_BS;						//!< autocorrelation of EoD spread
	double L_phi_MS;						//!< autocorrelation of AoA spread
	double L_theta_MS;						//!< autocorrelation of EoA spread
	double L_sf;							//!< autocorrelation of shadow fading
	// Cross-correlations
	mat x_corr_matrix;						//!< cross-correlation matrix (Cholesky factorized)
	double rho_BS_BS;						//!< cross-correlation coefficient between BSs
	// Polarization XPD
	double mu_xpdv;							//!< mean of XPDV
	double sigma_xpdv;						//!< std. of XPDV
	double mu_xpdh;							//!< mean of XPDH
	double sigma_xpdh;						//!< std. of XPDH
	double mu_cpr;							//!< mean of CPR
	double sigma_cpr;						//!< std. of CPR
	// DMC
	double mu_spread_dmc;					//!< mean of DMC spatial spread	(m)
	double sigma_spread_dmc;				//!< std. of DMC spatial spread	(m)
	double beta_dmc;						//!< DMC delay power decaying	(dB)

	// other parameters
	// cluster to visibility region parameters
	PDF_TYPE pdf_phi_C;						//!< PDF of azimuth angle of cluster to visibility region
	PDF_TYPE pdf_theta_C;					//!< PDF of elevation angle of cluster to visibility region
	PDF_TYPE pdf_r_C;						//!< PDF of distance between cluster to visibility region
	vec phi_C;								//!< distribution factor vector of cluster azimuth angle
	vec theta_C;							//!< distribution factor vector of cluster elevation angle
	vec r_C;								//!< distribution factor vector of cluster distance
	double mu_tau_C;						//!< mean of cluster delay
	double sigma_tau_C;						//!< std. of cluster delay

	BS_Info BS_info;						//!< BS information
	MS_Info MS_info;						//!< MS information

	imat BSCC;								//!< BS common cluster table
	ivec MSCC;								//!< MS common cluster table

	Array<Cluster> channel_cluster;			//!< set of all channel clusters

	vector<Transfer_Function> transfer;		//!< channel transfer functions
	//! all simulated channel transfer functions at different snapshots
	vector<vector<Transfer_Function> > transfer_function;

	// static parameters
	int c0;									//!< wave speed					(m/s)
	double lambda;							//!< wavelength c0 / freq		(m)

};
//! Generates a Poisson distributed random number
inline double get_Poisson_number(double lambda) { double L = exp(-1 * lambda); double k = 0; double p = 1.0; while (p > L) { k ++; p = p * randu(); } k --; return (k == 0) ? 1 : k; }
//! Generates a sized index vector
inline ivec get_idx_vec(const int size) { ivec a(size); for (int i = 0; i < size; i ++) { a.set(i, i); } return a; }
//! Find elements in the vector
inline ivec find_element(const ivec &v, const int ele) { ivec temp; int pos = 0; for (int i = 0; i < v.length(); i ++) { if (v.get(i) == ele) { temp.ins(pos, i); pos++; } } return temp; }
//! Delete a index in the index vector
inline bool delete_idx_element(ivec &v, const int idx) { for (int i = 0; i < v.length(); i ++) { if (v.get(i) == idx) { v.del(i); return true; } } return false; }
//! Delete elements in the vector
inline ivec delete_elements(ivec &v, const double ele) { ivec temp; int pos = 0; for (int i = 0; i < v.length(); i ++) { if (v.get(i) != ele) { temp.ins(pos, i); pos++; } } return temp; }
//! Delete repeated elements in the vector (for simplicity, ivec is chosen)
inline ivec delete_repeated_elements(const ivec &v) { ivec vec_sorted_idx = sort_index(v); ivec temp = v; for (int i = 1; i < v.length(); i ++) { if (v.get(vec_sorted_idx.get(i)) == v.get(vec_sorted_idx.get(i-1))) { temp.set(vec_sorted_idx.get(i), 0); } } return delete_elements(temp, 0); }
/*! \brief Transform the coordinates according to the coordinate transformation type
 *
 * 	CART2SPH: \f$[x,y,z]\f$ to \f$[\varphi, \theta, r]\f$
 *
 *  SPH2CART: \f$[\varphi, \theta, r]\f$ to \f$[x, y, z]\f$
 *
 *  POL2CART: \f$[\varphi, r, z]\f$ to \f$[x, y, z]\f$
 *
 *  CART2POL: \f$[x, y, z]\f$ to \f$[\varphi, r, z]\f$
 */
inline vec coordinate_transformation(double a, double b, double c, const COOR_TRANS_TYPE ct_type) { vec output = zeros(3); switch (ct_type) { case CART2SPH: output.set(0, atan(b / a)); output.set(1, atan(c / sqrt(pow(a, 2) + pow(b, 2)))); output.set(2, sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2))); break; case SPH2CART: output.set(0, c * cos(b) * cos(a)); output.set(1, c * cos(b) * sin(a)); output.set(2, c * sin(b)); break; case POL2CART: output.set(0, b * cos(a)); output.set(1, b * sin(a)); output.set(2, c); break; case CART2POL: output.set(0, atan(b / a)); output.set(1, sqrt(pow(a, 2) + pow(b, 2))); output.set(2, c); break; } return output; }
//! Transform the coordinates according to the coordinate transformation type
inline vec coordinate_transformation(const vec &input, const COOR_TRANS_TYPE ct_type) { if (input.length() != 3) return zeros(3); return coordinate_transformation(input.get(0), input.get(1), input.get(2), ct_type); }
//! Generate random number according to the distribution parameters
inline double get_random_number(const PDF_TYPE pdf_type, const vec &para) { double output; switch (pdf_type) { case NORMAL: output = para.get(0) + randn() * para.get(1); break; case UNIFORM: output = randu() * (para.get(1) - para.get(0)) + para.get(0); break; case POISSON: output = get_Poisson_number(para.get(0)); break;} return output; }
//! Calculate the distance between two positions
inline double calc_dist(const vec &pos1, const vec &pos2) { if (pos1.length() != pos2.length()) { ERROR("calc_dist(): Two positions cannot be of different dimensions!"); return 0; } vec pos = pos1 - pos2; double dist = 0; for (int i = 0; i < pos.length(); i ++)	dist += pow(pos.get(i), 2); return sqrt(dist); }
//! Calculate the distance between one position and a set of other positions
inline vec calc_dist(const vec &pos1, const Vec<vec> &pos2) { vec dists; for (int i = 0; i < pos2.length(); i ++) dists.ins(i, calc_dist(pos1, pos2.get(i))); return dists; }
//! Calculate the distance between one position and a set of other positions
inline vec calc_dist(const vec &pos1, const mat &pos2) { vec dists; for (int i = 0; i < pos2.rows(); i ++) dists.ins(i, calc_dist(pos1, pos2.get_row(i))); return dists; }
/*! Generate a 3x3 rotation matrix  by given rotation in azimuth \f$\varphi\f$ and elevation \f$\theta\f$, assume \f$\delta = 0\f$
 *
 * The rotation matrix is

 * \f{eqnarray*}{
 * 		\cos(\varphi)\cos(\theta)	&	-\sin(\varphi)	&	\cos(\varphi)\sin(\theta) \\
 *		\sin(\varphi)\cos(\theta)	&	\cos(\varphi)	&	\sin(\varphi)\cos(\theta) \\
 *		-\sin(\theta)				&	0				&	\cos(\theta)
 * \f}
 */
inline mat rotate_matrix(const vec &phi_theta) { double delta = 0; double phi = phi_theta.get(0); double theta = phi_theta.get(1); mat m = repmat(zeros(3), 1, 3); double sin_delta = sin(delta); double cos_delta = cos(delta); mat Tx = m; Tx.set(0, 0, 1); Tx.set(1, 1, cos_delta); Tx.set(1, 2, -1 * sin_delta); Tx.set(2, 1, sin_delta); Tx.set(2, 2, cos_delta); double sin_theta = sin(theta);double cos_theta = cos(theta); mat Ty = m; Ty.set(0, 0, cos_theta); Ty.set(0, 2, sin_theta); Ty.set(1, 1, 1); Ty.set(2, 0, -1 * sin_theta); Ty.set(2, 2, cos_theta); double sin_phi = sin(phi); double cos_phi = cos(phi); mat Tz = m;	Tz.set(0, 0, cos_phi); Tz.set(0, 1, sin_phi); Tz.set(1, 0, -1 * sin_phi); Tz.set(1, 1, cos_phi); Tz.set(2, 2, 1); return Tx * Ty * Tz; }

//@}

} // namespace penux

#endif /* COST2100_SPECIFICATION_H_ */
