/*
 * COST2100_Channel.cpp
 *
 *  Created on: May 24, 2010
 *      Author: K.K.
 */

#include "COST2100_Channel.h"
#include "port.h"
#include "debug.h"

#define DEBUG
#undef DBG_LVL
#define DBG_LVL (DBG_ERR|DBG_INFO)//DBG_TRACE|DBG_WARN|DBG_INFO|DBG_ENTER|DBG_LEAVE)


#include <itpp/base/vec.h>
#include <itpp/stat/misc_stat.h>
#include <sys/stat.h>
#include <assert.h>
#include <vector>
#include <string>

namespace penux {

using namespace itpp;

COST2100_Channel::COST2100_Channel(string blockName, string fileName, int idxBS, int idxMS, int numTx, int numRx, int numFreqBins, int numTime, double SNRdB)
: BasicProcessBlock(blockName), controlPort(this, duration), fileName(fileName), idxBS(idxBS), idxMS(idxMS),
  numTx(numTx), numRx(numRx), numFreqBins(numFreqBins), SNRdB(SNRdB)
{
	ENTER("COST2100_Channel::COST2100_Channel()");

	if (exist(fileName)) {
		//initialize channel
		channel = Channel_Specification (fileName);
	    INFO("Load channel parameters!");
	} else {
	    ERROR("Unable to open file named "<<fileName);
	    assert(0);
	}

	vec numTxRx = channel.get_antenna_nums(idxBS, idxMS);
	vec numBSMS = channel.get_BS_MS_nums();

	if(this->numTx != numTxRx.get(0) || this->numRx != numTxRx.get(1) || this->numFreqBins != channel.get_freq_bins()
			|| this->idxBS >= numBSMS.get(0) || this->idxMS >= numBSMS.get(1) || this->numTime > channel.get_snap_num()) {
	    ERROR("COST2100_Channel::COST2100_Channel() - Proposed parameter values does not match channel specification on file.");
	    assert(0);
	}
	/*
	// reset the proposed parameters
	else {
		numTx = numTxRx.get(0);
		numRx = numTxRx.get(1);
		numFreqBins = channel.get_freq_bins();
		numTime = channel.get_snap_num() - 1;
	}
	*/

	currentSnapshot = 0;

	vector<vector<Transfer_Function> > v = channel.get_all_transfer_functions();
	vector<vector<Transfer_Function> >::iterator it;
	for(it = v.begin(); it != v.end(); it++) {
		vector<Transfer_Function>::iterator itt;
		for(itt = (*it).begin(); itt != (*it).end(); itt++) {
			Transfer_Function t = *itt;
			if (t.get_bs_idx() == idxBS && t.get_ms_idx() == idxMS) {
				transfer.push_back(t.get_transfer_function());
			}
		}
	}

	//Create and add the data input ports
	for (int i = 0; i < numTx; i++)
	{
	  InPort<cvec>* pDataInPort =new InPort<cvec>(this, duration);
	  vpDataInPort.push_back(pDataInPort); //Add to objects variable vector
	  vInPort.push_back(pDataInPort); //Add to "hidden" interface vector
	}

	//Create and add the channel output ports
	for (int i=0; i < numRx; i++)
	{
	  OutPort<cvec>* pDataOutPort = new OutPort<cvec>(this, duration);
	  vpDataOutPort.push_back(pDataOutPort); //Add to objects variable vector
	  vOutPort.push_back(pDataOutPort); //Add to "hidden" interface vector
	}

	//Add control port, to be used for setting when new channel realization
	vInPort.push_back(&controlPort);

	LEAVE("COST2100_Channel::COST2100_Channel()");
}

COST2100_Channel::~COST2100_Channel()
{
	ENTER("COST2100_Channel::~COST2100_Channel()");
	for(int i = 0; i < numTx; i++) {
		delete vpDataInPort[i];
	}
	for(int i = 0; i < numRx; i++) {
		delete vpDataOutPort[i];
	}
	LEAVE("COST2100_Channel::~COST2100_Channel()");

}

void COST2100_Channel::run(int time)
{
  ENTER("COST2100_Channel::run()");

	//Set meta-data for the out-port stream

	//=========== Extract and build matrix from indata ================
	cmat inputDataMat(numTx, numFreqBins);

	//Go through the data inports and extract data
	int currentRow = 0;
	vector<InPort<cvec>*>::iterator ipInPort;
	for (ipInPort = vpDataInPort.begin(); ipInPort != vpDataInPort.end(); ++ipInPort, ++currentRow) {
		shared_ptr<const cvec> pInData = (*ipInPort)->getConstData();

		if (!pInData) {
			//If input vector empty, set all zeros as received vector
			shared_ptr<cvec> pTemp(new cvec(numFreqBins));
			pTemp->zeros();
			pInData = pTemp;
		}
		if (length(*pInData) != numFreqBins) {
			ERROR(name <<": Number of input samples does not match the number of channel taps!");
		} else {
			//Extract data
			inputDataMat.set_row(currentRow, *pInData);
		}
	}

	if (controlPort.isConnected()) {
		if (controlPort.getConstData()) {
			loadChannel();
			INFO("COST2100_Channel::run() - Channel created");
		}
	} else if (H.length() == 0) { //when no ctrl ch, and H empty
		loadChannel();
		INFO("COST2100_Channel::run() - Channel created");
	}

	//=========== CREATE OUTPUT USING THE INPUT AND CHANNEL ================
	INFO("COST2100_Channel::run() - Create output \n");

	int currentRx = 0; //index variable
	vector<OutPort<cvec>*>::iterator ipOutPort; //iterate over outports
	for (ipOutPort = vpDataOutPort.begin(); ipOutPort != vpDataOutPort.end(); ++ipOutPort, ++currentRx) {
		cvec RxDataVec(numFreqBins);
		vec chPow(numFreqBins);

		//Create output data
		for (int f = 0; f < numFreqBins; f ++) {
			RxDataVec(f) = H(f).get_row(currentRx) * inputDataMat.get_col(f);
		}
		//Add noise
		double var_noise = pow(10, (-1) * SNRdB / 10);
		cvec noise = randn_c(numFreqBins) * sqrt(2 * var_noise);
		RxDataVec = RxDataVec + noise;

		//Create data pointer, set metadata and write data to OutPort
		shared_ptr<cvec> pData(new cvec(RxDataVec));
		(*ipOutPort)->sendData(pData);
	}

	LEAVE("COST2100_Channel::run()");
}

void COST2100_Channel::loadChannel() {

	ENTER("COST2100_Channel::loadChannel()");
	// Load channel matrix from created channel simulation
	vector<Array<cmat> >::iterator it;
	int snapCount = 0;
	for(it = transfer.begin(); it != transfer.end(); it ++, snapCount ++) {
		if (currentSnapshot == snapCount) {
			H = *it;
			break;
		}
	}

	currentSnapshot++;

	LEAVE("COST2100_Channel::loadChannel()");
}

}
