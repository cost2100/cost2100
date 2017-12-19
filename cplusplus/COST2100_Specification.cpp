/*
 * COST2100_Specification.cpp
 *
 *  Created on: Mar 30, 2010
 *      Author: K.K.
 */

#include "COST2100_Specification.h"
#include "debug.h"

#define DEBUG
#undef DBG_LVL
#define DBG_LVL (DBG_ERR|DBG_INFO)//DBG_TRACE|DBG_WARN|DBG_INFO|DBG_ENTER|DBG_LEAVE)


#include <string>
#include <vector>

namespace penux {

using namespace itpp;

void MPC::update_MPC_pos(const SIDE_TYPE sd_type, int mpc_idx, const vec pos) {
	switch (sd_type) {
		case BS:
			MPC_BS_pos.set_row(mpc_idx, pos);
			break;
		case MS:
			MPC_MS_pos.set_row(mpc_idx, pos);
			break;
		case BS_MS:
			MPC_BS_pos.set_row(mpc_idx, pos);
			MPC_MS_pos.set_row(mpc_idx, pos);
			break;
	}
}

Cluster::Cluster(const CLUSTER_TYPE c_type, const vec &c_bs_pos, const vec &c_ms_pos, const vec &spread, double c_sf, double c_tau_link, const vec &angle, const MPC mpc) {
	C_type = c_type;
	C_BS_pos = c_bs_pos;
	C_MS_pos = c_ms_pos;
	C_delay_BS = spread.get(0);
	C_AoD_BS = spread.get(1);
	C_EoD_BS = spread.get(2);
	C_delay_MS = spread.get(3);
	C_AoD_MS = spread.get(4);
	C_EoD_MS = spread.get(5);
	C_shadow_f = c_sf;
	C_tau_link = c_tau_link;
	C_Phi_BS = angle.get(0);
	C_Theta_BS = angle.get(1);
	C_Phi_MS = angle.get(2);
	C_Theta_MS = angle.get(3);
	C_MPC = mpc;
}

double Cluster::get_cluster_spread(const SPREAD_TYPE s_type) {
	double spread;
	switch (s_type) {
		case DELAY_AT_BS:
			spread = C_delay_BS;
			break;

		case DELAY_AT_MS:
			spread = C_delay_MS;
			break;

		case AOD_AT_BS:
			spread = C_AoD_BS;
			break;

		case AOD_AT_MS:
			spread = C_AoD_MS;

		case EOD_AT_BS:
			spread = C_EoD_BS;
			break;

		case EOD_AT_MS:
			spread = C_EoD_MS;
			break;

	}
	return spread;
}

vec Cluster::get_cluster_spread_vec(const SIDE_TYPE sd_type) {
	vec s_vec;
	switch (sd_type) {
		case BS:
			s_vec.ins(0, C_delay_BS);
			s_vec.ins(1, C_AoD_BS);
			s_vec.ins(2, C_EoD_BS);
			break;
		case MS:
			s_vec.ins(0, C_delay_MS);
			s_vec.ins(1, C_AoD_MS);
			s_vec.ins(2, C_EoD_MS);
			break;
		case BS_MS:
			s_vec.ins(0, C_delay_BS);
			s_vec.ins(1, C_AoD_BS);
			s_vec.ins(2, C_EoD_BS);
			s_vec.ins(3, C_delay_MS);
			s_vec.ins(4, C_AoD_MS);
			s_vec.ins(5, C_EoD_MS);
			break;
	}
	return s_vec;
};

double Cluster::get_cluster_angle(const ANGLE_TYPE a_type) {
	double angle;
	switch (a_type) {
		case AZIMUTH_AT_BS:
			angle = C_Phi_BS;
			break;

		case AZIMUTH_AT_MS:
			angle = C_Phi_MS;
			break;

		case ELEVATION_AT_BS:
			angle = C_Theta_BS;
			break;

		case ELEVATION_AT_MS:
			angle = C_Theta_MS;
			break;

	}
	return angle;
}

vec Cluster::get_cluster_angle_vec(const SIDE_TYPE sd_type) {
	vec a_vec;
	switch (sd_type) {
		case BS:
			a_vec.ins(0, C_Phi_BS);
			a_vec.ins(1, C_Theta_BS);
			break;

		case MS:
			a_vec.ins(0, C_Phi_MS);
			a_vec.ins(1, C_Theta_MS);
			break;

		case BS_MS:
			a_vec.ins(0, C_Phi_BS);
			a_vec.ins(1, C_Theta_BS);
			a_vec.ins(2, C_Phi_MS);
			a_vec.ins(3, C_Theta_MS);
			break;

	}
	return a_vec;
}

Channel_Specification::Channel_Specification(string fileName) {
	// initialize the random seeds
	RNG_randomize ();
	getParamsFromXML(fileName);
	init_BSCC();
	init_MSCC();
	init_VR();
	init_cluster();
	init_channel();
}

void Channel_Specification::init_channel() {
	// for each snapshot
	for (int t = 0; t < snapNum; t ++) {
		int bs_num = BS_info.get_BS_num();
		int ms_num = MS_info.get_MS_num();
		// for each BS
		for (int i = 0; i < bs_num; i ++) {
			vec bs_pos = BS_info.get_BS_pos().get_row(i);
			ivec bs_vr_idx = find_element(BSCC.get_row(i), 1);
			mat bs_vr = BS_info.get_VR().get_rows(bs_vr_idx);
			vec bs_vr_los = BS_info.get_VR_LOS().get_row(i);
			// for each MS
			for (int j = 0; j < ms_num; j ++) {
				mat channel_matrix;
				vec ms_pos = MS_info.get_MS_pos().get_row(j);
				double d_ms_bs = calc_dist(bs_pos, ms_pos); // distance between BS and MS
				if (d_ms_bs == 0) {
					ERROR("init_channel(): BS/MS at the same position");
					return;
				} else {
					// active far clusters
					double pathloss = calc_path_loss(d_ms_bs);
					ivec active_vr_idx;
					vec vr_gain, vr_amp;
					int pos = 0;
					for (int k = 0; k < bs_vr_idx.length(); k ++) {
						double dist_tmp = calc_dist(ms_pos.get(0, 1), bs_vr.get_row(k).get(0, 1));
						if (dist_tmp < R_C) { // active VR
							double y = dist_tmp + L_C - R_C;
							vr_gain.ins(pos,  1 / 2 - atan(2 * sqrt(2) * y / sqrt(lambda * L_C)) / pi);
							active_vr_idx.ins(pos, bs_vr_idx.get(k));
							pos ++;
						}
					}
					if (active_vr_idx.length() > 0) {
						double tau_0 = d_ms_bs / c0;
						imat mscc = MSCC;
						ivec active_cluster_idx = mscc.get_rows(active_vr_idx).get_col(0);
						ivec active_cluster = delete_repeated_elements(active_cluster_idx);
						mat vrgain = vr_gain;
						for (int l = 0; l < active_cluster.length(); l ++) {
							int active_cluster_idx = active_cluster.get(l);
							vr_amp.ins(l, max(vrgain.get_rows(find_element(MSCC, active_cluster_idx)).get_col(0)));
							Cluster current_active_cluster = channel_cluster(active_cluster_idx);
							double d_bs_c_ms = calc_dist(bs_pos, current_active_cluster.get_cluster_pos(BS))
											+ calc_dist(ms_pos, current_active_cluster.get_cluster_pos(MS));
							double tau = d_bs_c_ms / c0 + current_active_cluster.get_cluster_link_delay();
							double cluster_att = std::max(exp(-k_tau * (tau - tau_0) * 1e6), exp(-k_tau * (tau_B - tau_0) * 1e6)); // attenuation of clusters
							MPC cluster_MPC = current_active_cluster.get_cluster_MPC();
							double cluster_sf = current_active_cluster.get_cluster_shadowing_fading();
							for (int m = 0; m < N_MPC; m ++) {
								vec mpc_bs_sph = coordinate_transformation(cluster_MPC.get_MPC_pos(BS).get_row(m) - bs_pos, CART2SPH);
								vec mpc_ms_sph = coordinate_transformation(cluster_MPC.get_MPC_pos(MS).get_row(m) - ms_pos, CART2SPH);
								double tau_mpc = (mpc_bs_sph.get(2) + mpc_ms_sph.get(2)) / c0 + current_active_cluster.get_cluster_link_delay();
								complex<double> phase(0, -2 * pi * freq_c * tau_mpc);
								complex<double> channel_amp = sqrt(cluster_sf * cluster_att)
														* cluster_MPC.get_MPC_amplitude().get(m) * pathloss * exp(phase);
								double channel_amp_real = channel_amp.real();
								double channel_amp_imag = channel_amp.imag();
								// DDIR
								vec channel = concat(mpc_bs_sph.get(0, 1), mpc_ms_sph.get(0, 1), to_vec(tau_mpc), to_vec(channel_amp_real), to_vec(channel_amp_imag));
								channel_matrix.append_row(channel);
							}
						}
					}

					std::cout << "Info: init_channel(): channel initializing...!" << std::endl;

					// active BS local cluster
					if (BS_local) {
						Cluster local_cluster = (BS_info.get_BS_local_cluster())(i);
						MPC local_MPC = local_cluster.get_cluster_MPC();
						double cluster_sf = local_cluster.get_cluster_shadowing_fading();
						for (int m = 0; m < N_MPC; m ++) {
							vec mpc_bs_sph = coordinate_transformation(local_MPC.get_MPC_pos(BS).get_row(m) - bs_pos, CART2SPH);
							vec mpc_ms_sph = coordinate_transformation(local_MPC.get_MPC_pos(MS).get_row(m) - ms_pos, CART2SPH);
							double tau_mpc = (mpc_bs_sph.get(2) + mpc_ms_sph.get(2)) / c0 + local_cluster.get_cluster_link_delay();
							complex<double> phase(0, -2 * pi * freq_c * tau_mpc);
							complex<double> channel_amp = sqrt(cluster_sf) * local_MPC.get_MPC_amplitude().get(m) * pathloss * exp(phase);
							double channel_amp_real = channel_amp.real();
							double channel_amp_imag = channel_amp.imag();
							// DDIR
							vec channel = concat(mpc_bs_sph.get(0, 1), mpc_ms_sph.get(0, 1), to_vec(tau_mpc), to_vec(channel_amp_real), to_vec(channel_amp_imag));
							channel_matrix.append_row(channel);
						}
					}
					// active MS local cluster
					if (MS_local) {
						Cluster local_cluster = (MS_info.get_MS_local_cluster())(j);
						MPC local_MPC = local_cluster.get_cluster_MPC();
						double cluster_sf = local_cluster.get_cluster_shadowing_fading();
						for (int m = 0; m < N_MPC; m ++) {
							vec mpc_bs_sph = coordinate_transformation(local_MPC.get_MPC_pos(BS).get_row(m) - bs_pos, CART2SPH);
							vec mpc_ms_sph = coordinate_transformation(local_MPC.get_MPC_pos(MS).get_row(m) - ms_pos, CART2SPH);
							double tau_mpc = (mpc_bs_sph.get(2) + mpc_ms_sph.get(2)) / c0 + local_cluster.get_cluster_link_delay();
							complex<double> phase(0, -2 * pi * freq_c * tau_mpc);
							complex<double> channel_amp = sqrt(cluster_sf) * local_MPC.get_MPC_amplitude().get(m) * pathloss * exp(phase);
							double channel_amp_real = channel_amp.real();
							double channel_amp_imag = channel_amp.imag();
							// DDIR
							vec channel = concat(mpc_bs_sph.get(0, 1), mpc_ms_sph.get(0, 1), to_vec(tau_mpc), to_vec(channel_amp_real), to_vec(channel_amp_imag));
							channel_matrix.append_row(channel);
						}
					}

					// LOS component
					double power_los, power_factor_los;
					if (d_ms_bs > d_co) {
						power_los = 0;
						power_factor_los = 0;
					} else {
						double d_ms_vr_los = calc_dist(ms_pos.get(0, 1), bs_vr_los.get(0, 1));
						if (d_ms_vr_los > R_L) {
							power_los = 0;
							power_factor_los = 0;
						} else {
							double y = d_ms_vr_los + L_L - R_L;
							double vr_los_gain = 1 / 2 - atan(2 * sqrt(2) * y / sqrt(lambda * R_L)) / pi;
							double power_other = 0;
							power_factor_los = power_factor;
							power_other = pow(sum(channel_matrix.get_col(5)), 2) + pow(sum(channel_matrix.get_col(6)), 2);
							power_los = abs(power_factor_los) * power_other;
						}
					}

					vec los_bs_sph = coordinate_transformation(ms_pos - bs_pos, CART2SPH);
					vec los_ms_sph = coordinate_transformation(bs_pos - ms_pos, CART2SPH);

					double d_bs_ms_xy = calc_dist(bs_pos.get(0, 1), ms_pos.get(0, 1));// distance BS to MS, in X-Y plane
					double tau_0_xy = d_bs_ms_xy / c0; // delay of the LOS
					complex<double> phase_los(0,  -2 * pi * freq_c * tau_0_xy);
					complex<double> channel_amp_los = sqrt(power_los) * exp(phase_los);
					double channel_amp_los_real = channel_amp_los.real();
					double channel_amp_los_imag = channel_amp_los.imag();
					// DDIR
					vec channel_los = concat(los_bs_sph.get(0, 1), los_ms_sph.get(0, 1), to_vec(tau_0_xy), to_vec(channel_amp_los_real), to_vec(channel_amp_los_imag));
					channel_matrix.append_row(channel_los);

					std::cout << "Info: init_channel(): " << channel_matrix.size() << std::endl;



					/********** Creating frequency response with antenna settings: start **********/
					// FIXME: change it to be compatible with external antenna settings

					// Reserved TX, RX numbers TODO: get from external parameters
					int tx_num = 2;
					int rx_num = 2;

					// Get transfer function with DDIRs
					Array<cmat> H_all(N_freq);
					vec freqs(N_freq );
					for (int f = 0; f < N_freq; f++) {
						cmat H(rx_num, tx_num);
						double freq = freq_start + f * freq_div;
						freqs.set(f, freq);
						for (int m = 0; m < channel_matrix.rows(); m++ ) {
							int tau_idx = ceil(channel_matrix.get(m, 4) / sampleRate);
							if (tau_idx > 0 && tau_idx < tau_idx_max ) {
								complex<double> delay_phase(0,  -2 * pi * freq * channel_matrix.get(m, 4));
								complex<double> delay_response = exp(delay_phase);

								double aod = channel_matrix.get(m, 0);
								double eod = channel_matrix.get(m, 1);
								mat bsv(1, 3);
								bsv.set(0, 0, sin(eod) * cos(aod));
								bsv.set(0, 1, sin(eod)*sin(aod));
								bsv.set(0, 2, cos(eod));
								bsv = bsv * rotate_matrix(vec("0, 0"));//txRot(1), txRot(2)
								aod = acos(bsv(0, 2));
								eod = atan2(bsv(0, 1), bsv(0, 0));

								double aoa = channel_matrix.get(m, 2);
								double eoa = channel_matrix.get(m, 3);
								mat msv(1, 3);
								msv.set(0, 0, sin(eoa)*cos(aoa));
								msv.set(0, 1, sin(eoa)*sin(aoa));
								msv.set(0, 2, cos(eoa));
								msv = msv * rotate_matrix(vec("0, 0")); //rxRot(1), rxRot(2)
								aod=acos(msv(0, 2));
								eod=atan2(msv(0, 1),msv(0, 0));

								cmat antenna_response_tx(tx_num, 1);
								cmat antenna_response_rx(rx_num, 1);

								/*
								// Antenna repsonse general version
								// load from data
								vec txAziRange, txEleRange, rxAziRange, rxEleRange;
								Array<cmat> txResponses, rxResponses;
								int txAziIdx = min_index(abs(txAziRange - aod * 180 / pi));
								int txEleIdx = min_index(abs(txEleRange - eod * 180 / pi));
								int rxAziIdx = min_index(abs(rxAziRange - aoa * 180 / pi));
								int rxEleIdx = min_index(abs(rxEleRange - eoa * 180 / pi));
								for (int tx = 0; tx < tx_num; tx++) {
									antenna_response_tx.set(tx, 0, txResponses(tx).get(txAziIdx, txEleIdx));
								}
								for (int rx = 0; rx < rx_num; rx++) {
									antenna_response_rx.set(rx, 0, rxResponses(rx).get(rxAziIdx, rxEleIdx));
								}
								*/

								// 2x2 antenna system test, assume half wavelength distance
								complex<double> antenna_shift_tx(0,  -pi * cos(aod));
								complex<double> antenna_shift_rx(0,  -pi * cos(aoa));
								antenna_response_tx.set(0, 0, 1);
								antenna_response_tx.set(1, 0, exp(antenna_shift_tx));
								antenna_response_rx.set(0, 0, 1);
								antenna_response_rx.set(1, 0, exp(antenna_shift_rx));

								complex<double> channel_amplitude(channel_matrix.get(m, 5), channel_matrix(m, 6));
								H += channel_amplitude * delay_response * antenna_response_rx * antenna_response_tx.transpose();
							} else {// out of the delay range
								continue;
							}
						}
						H_all(f) = H;
						//TODO: channel normalization with filters
					}
					
					/********** Creating frequency response with antenna settings: end **********/



					transfer.push_back(Transfer_Function(i, j, freqs, H_all));
				}
			}
		}
		update_MS_info();
		transfer_function.push_back(transfer);
	}

	std::cout << "Info: init_channel(): channel initialized!" << std::endl;
}



void Channel_Specification::getParamsFromXML(string fileName) {
	char *fn=new char[fileName.size()+1];
	fn[fileName.size()] = 0;
	memcpy(fn, fileName.c_str(), fileName.size());
	xmlDocPtr doc = xmlReadFile(fn,"UTF-8",XML_PARSE_RECOVER); // parse the file
	xmlNodePtr curNode, bsNodes, msNodes, channelNodes, exChannelNodes, stChannelNodes, spreadNodes, clusterNodes;
	xmlChar *xmlKey;

	if (NULL == doc) {
		ERROR("getParamsFromXML(): file cannot be read!");
		return;
	}
	curNode = xmlDocGetRootElement(doc);
	if (NULL == curNode) {
		xmlFreeDoc(doc);
		ERROR("getParamsFromXML(): file is blank!");
		return;
	}
	if (xmlStrcmp(curNode->name, BAD_CAST "ChannelSpecification")) { // check root element name
		xmlFreeDoc(doc);
		ERROR("getParamsFromXML(): file is invalid!");
		return;
	}
	// get children nodes
	curNode = curNode->xmlChildrenNode;
	while (curNode != NULL) {
		// set the parameters from the given xml

		// load BS information
		if ((!xmlStrcmp(curNode->name, BAD_CAST "BSInfo"))) {
			bsNodes = curNode->xmlChildrenNode;
			mat bs_pos;
			vec bscc_ratio;
			vec antenna_num;
			while (bsNodes != NULL) {
				if ((!xmlStrcmp(bsNodes->name, BAD_CAST "position"))) {
					xmlKey = xmlNodeGetContent(bsNodes);
					string pos((char*) xmlKey);
					bs_pos = mat(pos);
					xmlFree(xmlKey);
				} else if ((!xmlStrcmp(bsNodes->name, BAD_CAST "commonRatio"))) {
					xmlKey = xmlNodeGetContent(bsNodes);
					string ratio((char*) xmlKey);
					bscc_ratio = vec(ratio);
					xmlFree(xmlKey);
				} else if ((!xmlStrcmp(bsNodes->name, BAD_CAST "antennaNum"))) {
					xmlKey = xmlNodeGetContent(bsNodes);
					string antenna((char*) xmlKey);
					antenna_num = vec(antenna);
					xmlFree(xmlKey);
				}
				// TODO: add here for more parameter nodes if necessary
				bsNodes = bsNodes->next;
			}
			xmlFree(bsNodes);
			if (bs_pos.rows() != bscc_ratio.length()) {
				ERROR("getParamsFromXML(): invalid BS information!");
				return;
			}
			BS_info = BS_Info(bs_pos, bscc_ratio, antenna_num); // initialize the BS_Info object
			std::cout << "Info: getParamsFromXML(): BS_Info loaded!" << std::endl;
		}
		// load MS information
		else if ((!xmlStrcmp(curNode->name, BAD_CAST "MSInfo"))) {
			msNodes = curNode->xmlChildrenNode;
			mat ms_pos, ms_velo;
			double mscc;
			vec antenna_num;
			while (msNodes != NULL) {
				if ((!xmlStrcmp(msNodes->name, BAD_CAST "position"))) {
					xmlKey = xmlNodeGetContent(msNodes);
					string pos((char*) xmlKey);
					ms_pos = mat(pos);
					xmlFree(xmlKey);
				} else if ((!xmlStrcmp(msNodes->name, BAD_CAST "velocity"))) {
					xmlKey = xmlNodeGetContent(msNodes);
					string velo((char*) xmlKey);
					ms_velo = mat(velo);
					xmlFree(xmlKey);
				} else if ((!xmlStrcmp(msNodes->name, BAD_CAST "commonCluster"))) {
					xmlKey = xmlNodeGetContent(msNodes);
					mscc = atof((char*) xmlKey);
					xmlFree(xmlKey);
				} else if ((!xmlStrcmp(msNodes->name, BAD_CAST "antennaNUm"))) {
					xmlKey = xmlNodeGetContent(msNodes);
					string antenna((char*) xmlKey);
					antenna_num = vec(antenna);
					xmlFree(xmlKey);
				}
				// TODO: add here for more parameter nodes if necessary
				msNodes = msNodes->next;
			}
			xmlFree(msNodes);
			if (ms_pos.rows() != ms_velo.rows()) {
				ERROR("getParamsFromXML(): invalid MS information!");
				return;
			}
			MS_info = MS_Info(ms_pos, ms_velo, mscc, antenna_num);
			std::cout << "Info: getParamsFromXML(): MS_Info loaded!" << std::endl;
		}
		// load channel information
		else if ((!xmlStrcmp(curNode->name, BAD_CAST "Channel"))) {
			channelNodes = curNode->xmlChildrenNode;
			while (channelNodes != NULL) {
				if ((!xmlStrcmp(channelNodes->name, BAD_CAST "profile"))) {
					xmlKey = xmlNodeGetContent(channelNodes);
					string prof((char*) xmlKey);
					if (prof == "macro")
						profile = COST2100_MACRO;
					else if (prof == "micro")
						profile = COST2100_MICRO;
					else if (prof == "pico")
						profile = COST2100_PICO;
					else if (prof == "aalto")
						profile = COST2100_AALTO;
					else
						profile = COST2100_TEST;
					xmlFree(xmlKey);
					std::cout << "Info: getParamsFromXML(): profile loaded!" << std::endl;
				} else if ((!xmlStrcmp(channelNodes->name, BAD_CAST "frequency"))) {
					xmlKey = xmlNodeGetContent(channelNodes);
					string freq_str((char*) xmlKey);
					vec freq = vec(freq_str);
					if (freq.length() != 3 && freq.get(0) > freq.get(1) ) {
						ERROR("getParamsFromXML(): incorrect frequency range!");
						return;
					} else {
						freq_start = freq.get(0);
						freq_stop = freq.get(1);
						freq_div = freq.get(2);
						bandwidth = freq_stop - freq_start;
						N_freq = bandwidth / freq_div;
						freq_c = (freq_start + freq_stop) / 2;
						sampleRate = 1 / bandwidth;
						c0 = 3e8;
						lambda = c0 / freq_c;
					}
					xmlFree(xmlKey);
					std::cout << "Info: getParamsFromXML(): frequency loaded!" << std::endl;
				} else if ((!xmlStrcmp(channelNodes->name, BAD_CAST "snapRate"))) {
					xmlKey = xmlNodeGetContent(channelNodes);
					snapRate = atof((char*) xmlKey);
					xmlFree(xmlKey);
					std::cout << "Info: getParamsFromXML(): snapRate loaded!" << std::endl;
				} else if ((!xmlStrcmp(channelNodes->name, BAD_CAST "snapNum"))) {
					xmlKey = xmlNodeGetContent(channelNodes);
					snapNum = atof((char*) xmlKey);
					xmlFree(xmlKey);
					std::cout << "Info: getParamsFromXML(): snapNum loaded!" << std::endl;
				} else if ((!xmlStrcmp(channelNodes->name, BAD_CAST "overSampleRate"))) {
					xmlKey = xmlNodeGetContent(channelNodes);
					overSampleRate = atof((char*) xmlKey);
					xmlFree(xmlKey);
					std::cout << "Info: getParamsFromXML(): overSampleRate loaded!" << std::endl;
				} else if ((!xmlStrcmp(channelNodes->name, BAD_CAST "bandType"))) {
					xmlKey = xmlNodeGetContent(channelNodes);
					string band((char*) xmlKey);
					if (band == "wide") {
						bandType = WIDEBAND;
					} else if (band == "narrow")  {
						bandType = NARROWBAND;
					}
					xmlFree(xmlKey);
					std::cout << "Info: getParamsFromXML(): bandType loaded!" << std::endl;
				} else if ((!xmlStrcmp(channelNodes->name, BAD_CAST "external"))) { // external channel parameters
					exChannelNodes = channelNodes->xmlChildrenNode;
					while (exChannelNodes != NULL) {
						if ((!xmlStrcmp(exChannelNodes->name, BAD_CAST "bsHeight"))) { // TODO: check the range height expression
							xmlKey = xmlNodeGetContent(exChannelNodes);
							h_BS = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(exChannelNodes->name, BAD_CAST "msHeight"))) { // TODO: check the range height expression
							xmlKey = xmlNodeGetContent(exChannelNodes);
							h_MS = atof((char*) xmlKey);
							xmlFree(xmlKey);
						}  else if ((!xmlStrcmp(exChannelNodes->name, BAD_CAST "cellRadius"))) {
							xmlKey = xmlNodeGetContent(exChannelNodes);
							cellRadius = atof((char*) xmlKey);
							delay_max = cellRadius * 5; // maximum delay, 5 times of cell_radius [s]
							tau_idx_max = ceil(delay_max / sampleRate);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(exChannelNodes->name, BAD_CAST "rooftopHeight"))) {
							xmlKey = xmlNodeGetContent(exChannelNodes);
							h_rooftop = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(exChannelNodes->name, BAD_CAST "roadWidth"))) {
							xmlKey = xmlNodeGetContent(exChannelNodes);
							w_road = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(exChannelNodes->name, BAD_CAST "buildingSeparation"))) {
							xmlKey = xmlNodeGetContent(exChannelNodes);
							w_street = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(exChannelNodes->name, BAD_CAST "roadOrientation"))) {
							xmlKey = xmlNodeGetContent(exChannelNodes);
							phi_road = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(exChannelNodes->name, BAD_CAST "roomSize"))) {
							xmlKey = xmlNodeGetContent(exChannelNodes);
							string room_str((char*) xmlKey);
							vec room = vec(room_str);
							if (room.length() != 2) {
								ERROR("getParamsFromXML(): incorrect room size!");
								return;
							} else {
								l_room = room.get(0);
								w_room = room.get(1);
							}
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(exChannelNodes->name, BAD_CAST "floorNum"))) {
							xmlKey = xmlNodeGetContent(exChannelNodes);
							n_floor = atof((char*) xmlKey);
							xmlFree(xmlKey);
						}
						exChannelNodes = exChannelNodes->next;
					}
					xmlFree(exChannelNodes);
					std::cout << "Info: getParamsFromXML(): other external parameters loaded!" << std::endl;
				} else if ((!xmlStrcmp(channelNodes->name, BAD_CAST "stochastic"))) { // stochastic channel parameters
					stChannelNodes = channelNodes->xmlChildrenNode;
					while (stChannelNodes != NULL) {
						if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "vrRadius"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							R_C = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "trRadius"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							L_C = atof((char*) xmlKey);
							xmlFree(xmlKey);
						}  else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "clusterPower"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							k_tau = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "excessDelay"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							tau_B = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "cutoffDistLOS"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							d_co = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "vrRadiusLOS"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							R_L = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "trRadiusLOS"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							L_L = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "factorLOS"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							string factor_str((char*) xmlKey);
							vec factor_LOS = vec(factor_str);
							if (factor_LOS.length() != 2)
								ERROR("getParamsFromXML(): incorrect LOS factors!");
							else {
								mu_K = factor_LOS.get(0);
								sigma_K = factor_LOS.get(1);
								power_factor = mu_K * pow(10, (randn() * sigma_K / 10));
							}
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "singleClusterRatio"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							K_sel = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "averageLocalCluster"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							N_C_local = atof((char*) xmlKey);
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "averageFarCluster"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							mu_N_C_add = atof((char*) xmlKey);
							N_C_add = get_Poisson_number(mu_N_C_add);
							double rho_C = (N_C_add) / (pi * pow(R_C - L_C, 2));
							N_C_far = round_i(rho_C * (pi * pow(cellRadius, 2)));
							xmlFree(xmlKey);
							if (N_C_far == 0) {
								ERROR("getParamsFromXML(): incorrect average far cluster number! (unknown exception)");
								return;
							}
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "activeLocalCluster"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							string local_str((char*) xmlKey);
							vec local_active = vec(local_str);
							if (local_active.length() != 2) {
								ERROR("getParamsFromXML(): incorrect local cluster activeness!");
								return;
							} else {
								BS_local = (local_active.get(0) == 1) ? true : false;
								MS_local = (local_active.get(1) == 1) ? true : false;
							}
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "mpc"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							string mpc_str((char*) xmlKey);
							vec mpc_factor = vec(mpc_str);
							if (mpc_factor.length() != 2) {
								ERROR("getParamsFromXML(): incorrect MPC factors!");
								return;
							} else {
								N_MPC = mpc_factor.get(0);
								K_MPC = mpc_factor.get(1);
							}
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "diffuseRadiation"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							string diffuse_str((char*) xmlKey);
							vec diffuse_radiation = vec(diffuse_str);
							if (diffuse_radiation.length() != 3) {
								ERROR("getParamsFromXML(): incorrect diffuse radiation factor!");
								return;
							} else {
								mu_diff = diffuse_radiation.get(0);
								sigma_diff = diffuse_radiation.get(1);
								PDP_diff = diffuse_radiation.get(2);
							}
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "shadowing"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							string shadow_str((char*) xmlKey);
							vec shadowing = vec(shadow_str);
							if (shadowing.length() != 2) {
								ERROR("getParamsFromXML(): incorrect shadowing parameters!");
								return;
							} else {
								sigma_sf = shadowing.get(0);
								L_sf = shadowing.get(1);
							}
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "spread"))) {
							spreadNodes = stChannelNodes->xmlChildrenNode;
							while (spreadNodes != NULL) {
								if ((!xmlStrcmp(spreadNodes->name, BAD_CAST "delay"))) {
									xmlKey = xmlNodeGetContent(spreadNodes);
									string delay_str((char*) xmlKey);
									vec delay = vec(delay_str);
									if (delay.length() != 3) {
										ERROR("getParamsFromXML(): incorrect delay spread!");
										return;
									} else {
										mu_tau = delay.get(0);
										sigma_tau = delay.get(1);
										L_tau = delay.get(2);
									}
									xmlFree(xmlKey);
								} else if ((!xmlStrcmp(spreadNodes->name, BAD_CAST "AoD"))) {
									xmlKey = xmlNodeGetContent(spreadNodes);
									string aod_str((char*) xmlKey);
									vec aod = vec(aod_str);
									if (aod.length() != 3) {
										ERROR("getParamsFromXML(): incorrect AoD spread!");
										return;
									} else {
										mu_phi_BS = aod.get(0);
										sigma_phi_BS = aod.get(1);
										L_phi_BS = aod.get(2);
									}
									xmlFree(xmlKey);
								} else if ((!xmlStrcmp(spreadNodes->name, BAD_CAST "EoD"))) {
									xmlKey = xmlNodeGetContent(spreadNodes);
									string eod_str((char*) xmlKey);
									vec eod = vec(eod_str);
									if (eod.length() != 3) {
										ERROR("getParamsFromXML(): incorrect EoD spread!");
										return;
									} else {
										mu_theta_BS = eod.get(0);
										sigma_theta_BS = eod.get(1);
										L_theta_BS = eod.get(2);
									}
									xmlFree(xmlKey);
								} else if ((!xmlStrcmp(spreadNodes->name, BAD_CAST "AoA"))) {
									xmlKey = xmlNodeGetContent(spreadNodes);
									string aoa_str((char*) xmlKey);
									vec aoa = vec(aoa_str);
									if (aoa.length() != 3) {
										ERROR("getParamsFromXML(): incorrect AoA spread!");
										return;
									} else {
										mu_phi_MS = aoa.get(0);
										sigma_phi_MS = aoa.get(1);
										L_phi_MS = aoa.get(2);
									}
									xmlFree(xmlKey);
								} else if ((!xmlStrcmp(spreadNodes->name, BAD_CAST "EoA"))) {
									xmlKey = xmlNodeGetContent(spreadNodes);
									string eoa_str((char*) xmlKey);
									vec eoa = vec(eoa_str);
									if (eoa.length() != 3) {
										ERROR("getParamsFromXML(): incorrect EoA spread!");
										return;
									} else {
										mu_theta_MS = eoa.get(0);
										sigma_theta_MS = eoa.get(1);
										L_theta_MS = eoa.get(2);
									}
									xmlFree(xmlKey);
								}
								spreadNodes = spreadNodes->next;
							}
							xmlFree(spreadNodes);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "crossCorrelation"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							string x_corr_str((char*) xmlKey);
							x_corr_matrix = chol(mat(x_corr_str));
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "cluster"))) {
							clusterNodes = stChannelNodes->xmlChildrenNode;
							while (clusterNodes != NULL) {
								if ((!xmlStrcmp(clusterNodes->name, BAD_CAST "azimuth"))) {
									xmlKey = xmlNodeGetContent(clusterNodes);
									string azimuth_str((char*) xmlKey);
									vec azimuth = vec(azimuth_str);
									if (azimuth.length() != 3) {
										ERROR("getParamsFromXML(): incorrect azimuth parameters of cluster to visibility region!");
										return;
									} else {
										phi_C = azimuth.get(0, 1);
										pdf_phi_C = (azimuth.get(2)==0) ? NORMAL : UNIFORM;
									}
									xmlFree(xmlKey);
								} else if ((!xmlStrcmp(clusterNodes->name, BAD_CAST "elevation"))) {
									xmlKey = xmlNodeGetContent(clusterNodes);
									string elevation_str((char*) xmlKey);
									vec elevation = vec(elevation_str);
									if (elevation.length() != 3) {
										ERROR("getParamsFromXML(): incorrect elevation parameters of cluster to visibility region!");
										return;
									} else {
										theta_C = elevation.get(0, 1);
										pdf_theta_C = (elevation.get(2)==0) ? NORMAL : UNIFORM;
									}
									xmlFree(xmlKey);
								} else if ((!xmlStrcmp(clusterNodes->name, BAD_CAST "distance"))) {
									xmlKey = xmlNodeGetContent(clusterNodes);
									string distance_str((char*) xmlKey);
									vec distance = vec(distance_str);
									if (distance.length() != 3) {
										ERROR("getParamsFromXML(): incorrect distance parameters of cluster to visibility region!");
										return;
									} else {
										r_C = distance.get(0, 1);
										pdf_r_C = (distance.get(2)==0) ? NORMAL : UNIFORM;
									}
									xmlFree(xmlKey);
								} else if ((!xmlStrcmp(clusterNodes->name, BAD_CAST "linkDelay"))) {
									xmlKey = xmlNodeGetContent(clusterNodes);
									string link_str((char*) xmlKey);
									vec link = vec(link_str);
									if (link.length() != 2) {
										ERROR("getParamsFromXML(): incorrect link delay parameters of cluster to visibility region!");
										return;
									} else {
										mu_tau_C = link.get(0);
										sigma_tau_C = link.get(1);
									}
									xmlFree(xmlKey);
								}
								clusterNodes = clusterNodes->next;
							}
							xmlFree(clusterNodes);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "polarization"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							string xpd_str((char*) xmlKey);
							vec xpd = vec(xpd_str);
							if (xpd.length() != 6) {
								ERROR("getParamsFromXML(): incorrect XPD parameters!");
								return;
							} else {
								mu_xpdv = xpd.get(0);
								sigma_xpdv = xpd.get(1);
								mu_xpdh = xpd.get(2);
								sigma_xpdh = xpd.get(3);
								mu_cpr = xpd.get(4);
								sigma_cpr = xpd.get(5);
							}
							xmlFree(xmlKey);
						} else if ((!xmlStrcmp(stChannelNodes->name, BAD_CAST "dmc"))) {
							xmlKey = xmlNodeGetContent(stChannelNodes);
							string dmc_str((char*) xmlKey);
							vec dmc = vec(dmc_str);
							if (dmc.length() != 3) {
								ERROR("getParamsFromXML(): incorrect DMC parameters!");
								return;
							} else {
								mu_spread_dmc = dmc.get(0);
								sigma_spread_dmc = dmc.get(1);
								beta_dmc = dmc.get(2);
							}
							xmlFree(xmlKey);
						}
						stChannelNodes = stChannelNodes->next;
					}
					xmlFree(stChannelNodes);
					std::cout << "Info: getParamsFromXML(): stochastic parameters loaded!" << std::endl;
				}// TODO: add here for more parameter nodes
				channelNodes = channelNodes->next;
			}
			xmlFree(channelNodes);
			std::cout << "Info: getParamsFromXML(): Channel parameters loaded!" << std::endl;
		}
		curNode = curNode->next;
	}
	xmlFree(curNode);
	xmlFreeDoc(doc);
	std::cout << "Info: getParamsFromXML(): all file loaded!" << std::endl;

}


void Channel_Specification::update_MS_info() {
	// calculate new MS positions
	mat ms_pos_new = MS_info.get_MS_pos() + MS_info.get_MS_velo() * snapRate;
	for (int i = 0; i < MS_info.get_MS_num(); i ++) {
		Cluster local_cluster = (MS_info.get_MS_local_cluster())(i);
		MPC local_mpc = local_cluster.get_cluster_MPC();
		vec d_local_mpc_ms_new = calc_dist(ms_pos_new.get_row(i).get(0, 1), local_mpc.get_MPC_pos(MS).get_cols(0, 1));
		// check for MPC re-calculation
		for (int j = 0; j < d_local_mpc_ms_new.length(); j ++) {
			vec new_mpc_tmp(3);
			double d_mpc_tmp = d_local_mpc_ms_new.get(j);
			if (d_mpc_tmp > local_cluster.get_cluster_spread(DELAY_AT_MS)) {
				new_mpc_tmp = (randn(1, 3) * diag(local_cluster.get_cluster_spread_vec(MS) / 3)).get_row(0) + ms_pos_new.get_row(i);
				d_mpc_tmp= calc_dist(new_mpc_tmp, ms_pos_new.get_row(i));
				while (d_mpc_tmp > local_cluster.get_cluster_spread(DELAY_AT_MS)) {
					new_mpc_tmp = (randn(1, 3) * diag(local_cluster.get_cluster_spread_vec(MS) / 3)).get_row(0) + ms_pos_new.get_row(i);
					d_mpc_tmp= calc_dist(new_mpc_tmp, ms_pos_new.get_row(i));
				}
				local_mpc.update_MPC_pos(BS_MS, j, new_mpc_tmp);
			}
		}
		local_cluster.update_cluster_MPC(local_mpc);
		MS_info.update_MS_local_cluster(i, local_cluster);
	}
	MS_info.update_MS_pos(ms_pos_new);
	std::cout << "Info: update_MS_info(): MS updated!" << std::endl;
}

Array<Cluster> Channel_Specification::init_local_clusters(SIDE_TYPE sd_type) {
	int cluster_num;
	CLUSTER_TYPE c_type;
	mat pos;
	switch (sd_type) {
		case BS:
			c_type = LOCAL_AT_BS;
			cluster_num = BS_info.get_BS_num();
			pos = BS_info.get_BS_pos();
			break;

		case MS:
			c_type = LOCAL_AT_MS;
			cluster_num = MS_info.get_MS_num();
			pos = MS_info.get_MS_pos();
			break;

		default:
			ERROR ("Error local cluster side type!");
			break;

	}
	Array<Cluster> clusters(cluster_num);
	mat corr_randn;
	vec shadow_factor, tau_C, theta_C_BS, phi_C_BS, theta_C_MS, phi_C_MS, d_tau; // spatial delay spread
	init_cluster_params(cluster_num, corr_randn, shadow_factor, tau_C, theta_C_BS, phi_C_BS, theta_C_MS, phi_C_MS, d_tau);
	// calculate parameters of each cluster
	for (int i = 0; i < cluster_num; i ++) {
		vec current_pos = pos.get_row(i);
		double current_d_tau = d_tau.get(i);
		vec local_spread;
		local_spread.ins(0, current_d_tau);
		local_spread.ins(1, current_d_tau);
		local_spread.ins(2, current_d_tau / 20);
		MPC local_mpc = get_MPC(c_type, current_pos, current_pos, local_spread, zeros(2), zeros(2));
		clusters(i) = Cluster(c_type, current_pos, current_pos, repmat(local_spread, 2), shadow_factor.get(i), 0, zeros(4), local_mpc);
	}
	return clusters;
}

Array<Cluster> Channel_Specification::init_far_clusters() {
	int cluster_num = max(MSCC);
	Array<Cluster> clusters(cluster_num);
	mat corr_randn;
    vec shadow_factor, tau_C, theta_C_BS, phi_C_BS, theta_C_MS, phi_C_MS, d_tau; // spatial delay spread
    init_cluster_params(cluster_num, corr_randn, shadow_factor, tau_C, theta_C_BS, phi_C_BS, theta_C_MS, phi_C_MS, d_tau);
    // calculate parameters of each cluster
	for (int i = 0; i < cluster_num; i ++) {
		CLUSTER_TYPE c_type;
		vec pos_C_BS, pos_C_MS, spread, vr_pos, bs_pos;// activity
		get_ref_BS_VR(i, bs_pos, vr_pos);
		vec VR_direct =  vr_pos - bs_pos;
		VR_direct.set(2, 0);
		vec VR_sph1 = coordinate_transformation(VR_direct, CART2SPH);
		vec VR_sph2;
		VR_sph2.ins(0, VR_sph1.get(0) + get_random_number(pdf_phi_C , phi_C));
		VR_sph2.ins(1, VR_sph1.get(1) + get_random_number(pdf_theta_C , theta_C));
		VR_sph2.ins(2, get_random_number(pdf_r_C, r_C));

		pos_C_BS = bs_pos + coordinate_transformation(VR_sph2, SPH2CART); // cluster to VR at BS side
		double d_C_BS = calc_dist(pos_C_BS, bs_pos);
		// single or twin cluster
		int cluster_type = (randu() > K_sel) ? 1 : 0;
		double tau_C_link;
		switch (cluster_type) {
			case 0:
			c_type = SINGLE;
			pos_C_MS = pos_C_BS;
			tau_C_link = 0;
			break;

		case 1:
			c_type = TWIN;
			VR_sph2.set(0, VR_sph1.get(0) + get_random_number(pdf_phi_C , phi_C));
			VR_sph2.set(1, VR_sph1.get(1) + get_random_number(pdf_theta_C , theta_C));
			VR_sph2.set(2, d_C_BS * tan(phi_C_BS(i) / 180 * pi / 2) / tan(phi_C_MS(i) / 180 * pi / 2));
			pos_C_MS = vr_pos + coordinate_transformation(VR_sph2, SPH2CART); // cluster to VR at MS side
			tau_C_link = (calc_dist(pos_C_BS, pos_C_MS) / c0) + mu_tau_C * pow(10,(0.1 * sigma_tau_C * randn()));
			break;
		}

		spread.ins(0, d_tau(i));
		spread.ins(1, d_C_BS * tan(phi_C_BS(i) / 180 * pi / 2));
		spread.ins(2, d_C_BS * tan(theta_C_BS(i) / 180 * pi / 2));
		vec BS_angle = coordinate_transformation(bs_pos - pos_C_BS, CART2SPH).get(0, 1);
		vec MS_angle = coordinate_transformation(vr_pos - pos_C_MS, CART2SPH).get(0, 1);
		MPC mpc = get_MPC(c_type, pos_C_BS, pos_C_MS, spread, BS_angle, MS_angle);
		clusters(i) = Cluster(c_type, pos_C_BS, pos_C_MS, repmat(spread, 2),shadow_factor.get(i), tau_C_link, concat(BS_angle, MS_angle), mpc);
	}
	return clusters;
}

MPC Channel_Specification::get_MPC(const CLUSTER_TYPE c_type, const vec &c_bs_pos, const vec &c_ms_pos, const vec& spread, const vec &bs_angle, const vec &ms_angle) {
	mat MPC_BS_pos, MPC_MS_pos;
	mat tmp = randn(N_MPC, 3) * diag(spread / 3);
	switch (c_type) {
    	case LOCAL_AT_BS: // local cluster at BS
    		MPC_BS_pos = tmp + repmat(c_bs_pos.transpose(), N_MPC, 1);
    		MPC_MS_pos = MPC_BS_pos;
    		break;

    	case SINGLE: // single cluster
    		MPC_BS_pos = tmp * rotate_matrix(bs_angle) + repmat(c_bs_pos.transpose(), N_MPC, 1);
    		MPC_MS_pos = MPC_BS_pos;
    		break;

    	case TWIN: // twin cluster
	   		MPC_BS_pos = tmp * rotate_matrix(bs_angle) + repmat(c_bs_pos.transpose(), N_MPC, 1);
	   		MPC_MS_pos = tmp * rotate_matrix(ms_angle) + repmat(c_ms_pos.transpose(), N_MPC, 1);
	   		break;

    	case LOCAL_AT_MS: // local cluster at MS
	   		MPC_MS_pos = tmp + repmat(c_ms_pos.transpose(), N_MPC, 1);
	   		MPC_BS_pos = MPC_MS_pos;
	   		break;
	}
	cvec P_S = to_cvec(randn(N_MPC), randn(N_MPC)); // complex Rayleigh distribution
	return MPC(MPC_BS_pos, MPC_MS_pos, P_S / abs(sum(P_S))); //normalized attenuation of each MPC's complex Rayleigh amplitude

}

void Channel_Specification::init_cluster() {
	BS_info.set_BS_local_cluster(init_local_clusters(BS));
	MS_info.set_MS_local_cluster(init_local_clusters(MS));
	channel_cluster = init_far_clusters();

	std::cout << "Info: init_cluster(): cluster initialized!" << std::endl;
}

void Channel_Specification::init_BSCC() {
	int bs_num = BS_info.get_BS_num();
	vec bs_common = BS_info.get_BS_common_ratio();
	int vr_avg = N_C_far;
	imat bscc = to_imat(zeros(1,1));
	if (bs_common.length() != bs_num) {
		ERROR("Please specify correct common ratio!");
	} else {
		int vr_all = bs_num * vr_avg;
		ivec bs_idx = get_idx_vec(bs_num);
		ivec bs_nums = bs_idx + to_ivec(ones(bs_num));
		ivec bscc_vr_nums = to_ivec(round(vr_all * bs_common / (bs_nums * bs_common)));
		bscc = to_imat(zeros(bs_num, sum(bscc_vr_nums)));
		int pos = 0;
		int col_counter = 0;
		// assign BS common clusters
		for (int i = bs_num; i > 0 ; i --) {
			int bscc_vr_num = bscc_vr_nums.get(i - 1);
			for (int j = 0; j < bscc_vr_num; j ++) {
				ivec bscc_bs_idx = bs_idx;
				ivec bscc_tmp = to_ivec(zeros(bs_num));
				if (bscc_bs_idx.length() == 0) {
					ERROR("Empty vector for common cluster loaded!");
					break;
				} else {
					for(int k = 0; k < i; k ++) {
						pos = randi(0, bscc_bs_idx.length() - 1);
						int current_cc_bs_idx = bscc_bs_idx.get(pos);
						while (sum(bscc.get_row(current_cc_bs_idx)) >= vr_avg) {
							if(delete_idx_element(bs_idx, current_cc_bs_idx) && delete_idx_element(bscc_bs_idx, current_cc_bs_idx)) {
								if (bscc_bs_idx.length() == 0) {
									ERROR("Empty vector for common cluster loaded!");
									pos = -1;
									current_cc_bs_idx = -1;
									break;
								} else {
									pos = randi(0, bscc_bs_idx.length() - 1);
									current_cc_bs_idx = bscc_bs_idx.get(pos);
								}
							} else {
								ERROR("BS index not found!");
								break;
							}
						}
						if (pos != -1 && current_cc_bs_idx != -1) {
							bscc_tmp.set(current_cc_bs_idx, 1);
							bscc_bs_idx.del(pos);
						} else {
							break;
						}
					}
					if (sum(bscc_tmp) > 0) {
						bscc.set_col(col_counter, bscc_tmp);
						col_counter ++;
					}
				}
			}
		}
		bscc = bscc.get_cols(0, col_counter - 1);
	}
	BSCC = bscc;
	std::cout << "Info: init_BSCC(): BSCC initialized!" << std::endl;
}

void Channel_Specification::init_MSCC() {
	int vr_num = BSCC.cols();
	double ms_common = MS_info.get_MS_common();
	ivec vr_idx = get_idx_vec(vr_num);
	ivec mscc = to_ivec(zeros(vr_num));
	int cluster_counter = 0;
	int pos = 0;
	while (vr_idx.length() > 0) {
		// randomly assign MS common clusters
		int mscc_vr_group = get_Poisson_number(ms_common);
		while (vr_idx.length() < mscc_vr_group || mscc_vr_group == 0) {
			mscc_vr_group = get_Poisson_number(ms_common);
		}
		int current_cc_vr_idx;
		for (int i = 0; i < mscc_vr_group; i ++) {
			pos = randi(0, vr_idx.length() - 1);
			current_cc_vr_idx = vr_idx.get(pos);
			vr_idx.del(pos);
			mscc.set(current_cc_vr_idx, cluster_counter);
		}
		cluster_counter ++;
	}
	MSCC = mscc;
	std::cout << "Info: init_MSCC(): MSCC initialized!" << std::endl;
}

void Channel_Specification::init_VR() {
	mat bs_pos = BS_info.get_BS_pos();
	int vr_num = BSCC.cols();
	int bs_num = BS_info.get_BS_num();

	// calculate visibility regions
	double x_min = min(bs_pos.get_col(0)) - cellRadius;
	double x_max = max(bs_pos.get_col(0)) + cellRadius;
	double y_min = min(bs_pos.get_col(1)) - cellRadius;
	double y_max = max(bs_pos.get_col(1)) + cellRadius;

	mat vr(vr_num, 3);
	vr.set_col(0, (x_max - x_min) * randu(vr_num) + x_min * ones(vr_num));
	vr.set_col(1, (y_max - y_min) * randu(vr_num) + y_min * ones(vr_num));
	vr.set_col(2, zeros(vr_num));
	//TODO: check validity

	// calculate LoS visibility regions
	mat vr_los(bs_num, 3);
	double x, y;
	for (int i = 0; i < bs_num; i++) {
		x = randu() * (d_co - R_L);
		y = randu() * (d_co - R_L);
		while (sqrt(pow(x, 2) + pow(y, 2)) > d_co) {
			x = randu() * (d_co - R_L);
			y = randu() * (d_co - R_L);
		}
		vr_los.set(i, 0, x);
		vr_los.set(i, 1, y);
	}
	vr_los +=  bs_pos.get_cols(0, 1);
	vr_los.set_col(2, zeros(bs_num));
	BS_info.set_BS_visibility_region(vr, vr_los);
	std::cout << "Info: init_VR(): VR initialized!" << std::endl;
}

void Channel_Specification::init_cluster_params(int num_in, mat &corr_randn_out, vec &shadow_factor_out, vec &tau_C_out,
											vec &theta_C_BS_out, vec &phi_C_BS_out, vec &theta_C_MS_out, vec &phi_C_MS_out, vec &d_tau_out) {
	corr_randn_out = randn(num_in, 6) * x_corr_matrix; // get the random delay angular spread for each cluster
	shadow_factor_out = pow(10, (0.1 * sigma_sf * corr_randn_out.get_col(0))); // shadow fading
	tau_C_out = mu_tau * pow(10, (0.1 * sigma_tau * corr_randn_out.get_col(1))); // delay spread
	theta_C_BS_out = mu_theta_BS * pow(10, (0.1 * sigma_theta_BS * corr_randn_out.get_col(2))); // elevation spread to BS
	phi_C_BS_out = mu_phi_BS * pow(10, (0.1 * sigma_phi_BS * corr_randn_out.get_col(3))); // azimuth spread to BS
	theta_C_MS_out = mu_theta_MS * pow(10, (0.1 * sigma_theta_MS * corr_randn_out.get_col(4))); // elevation spread to MS
	phi_C_MS_out = mu_phi_MS * pow(10, (0.1 * sigma_phi_MS * corr_randn_out.get_col(5))); //azimuth spread to MS
	d_tau_out = tau_C_out * c0 / 2; // spatial delay spread
}

void Channel_Specification::get_ref_BS_VR(int cluster_idx, vec &ref_bs, vec &ref_vr) {
	// get reference BS position
	ivec vr_grp = find_element(MSCC, cluster_idx);
	imat bs_vr_grp = BSCC.get_cols(vr_grp);
	ivec bs_weight = sum(bs_vr_grp, 2);
	ivec ref_bs_ind = find_element(bs_weight, max(bs_weight));
	int ref_bs_idx = ref_bs_ind.get(randi(0, ref_bs_ind.length() - 1));
	ref_bs = BS_info.get_BS_pos().get_row(ref_bs_idx);

	// get reference visibility region position
	ivec ref_bs_vr_grp = bs_vr_grp.get_row(ref_bs_idx);
	ivec vr_weight = sum(BSCC.get_cols(ref_bs_vr_grp));
	ivec ref_vr_ind = find_element(vr_weight, max(vr_weight));
	int ref_vr_idx = ref_vr_ind.get(randi(0, ref_vr_ind.length() - 1));
	ref_vr = BS_info.get_VR().get_row(ref_vr_idx);
}

double Channel_Specification::calc_path_loss(const double dist) {
	double pathloss;
	switch(profile) {
		case COST2100_MACRO: { //COST 231 Walfisch- Ikegami Model
				double delta_base = h_BS - h_rooftop;
				double pl_0 = 32.45 + 20 * log10(dist / 1000) + 20 * log10(freq_c / 1e6); // free space loss
				double pl_ori;
				if (phi_road >= 0 && phi_road < 35) {
					pl_ori = -10+0.354*phi_road;
				} else if (phi_road >= 35 && phi_road < 55) {
					pl_ori = 2.5+0.075*(phi_road-35);
				} else if (phi_road >= 55 && phi_road < 90) {
					pl_ori = 4.0-0.114*(phi_road-55);
				}
				double pl_lts = -16.9 - 10*log10(w_road) + 10 * log10(freq_c / 1e6) + 20 * log10(h_rooftop - h_MS) + pl_ori;
				double pl_bsh, k_a, k_d;
				if (h_BS > h_rooftop) {
					pl_bsh = -18 * log10(1 + delta_base);
					k_a = 54;
					k_d = 18;
				} else {
					pl_bsh = 0;
					k_d = 18 - 15 * delta_base / h_rooftop;
					if (dist >= 500) {
						k_a = 54 - 0.8 * delta_base;
					} else {
						k_a = 54 - 0.8 * delta_base * dist / 500;
					}
				}
				double k_f = -4 + 0.7 * (freq_c / 1e6 / 925-1); // k_f = -4 + 1.5 * (freq_c / 1e6 / 925-1);
				double pl_msd = pl_bsh + k_a + k_d * log10(dist / 1000) + k_f * log10(freq_c / 1e6) - 9 * log10(w_street);
				if (pl_msd + pl_lts >= 0 ) {
					pathloss = pl_0 + pl_msd + pl_lts;
				} else {
					pathloss = pl_0;
				}
				pathloss = pow(10, (-pathloss/20));
			}
			break;
		case COST2100_MICRO: {
				double pl = 10 * 2.6 * log10(dist) + 20 * log10(4 * pi * 1/(c0 / freq_c));
				pathloss = pow(10, -pl / 20);
			}
			break;
		case COST2100_PICO:
		case COST2100_AALTO:
		default: {
				double pl_0 = 4 * pi * dist * freq_c / c0; //Free space pathloss
				double pl_d = n_floor * 30; //Average number of floors * 30 dB extra pathloss
				if (pl_0==0) {
					pathloss = pow(10, -pl_d/20);
				} else {
					pathloss = 1 / pl_0 * pow(10, -pl_d/20);
				}
			}
			break;

	}
	return pathloss;
}

vec Channel_Specification::get_antenna_nums(int bs_idx, int ms_idx) {
	vec a = zeros(2);
	a.set(0, BS_info.get_BS_antenna_num().get(bs_idx)); a.set(1, MS_info.get_MS_antenna_num().get(ms_idx));
	return a;
}

vec Channel_Specification::get_BS_MS_nums() {
	vec a(2);
	a.set(0, BS_info.get_BS_num());
	a.set(1, MS_info.get_MS_num());
	return a;
}

}
