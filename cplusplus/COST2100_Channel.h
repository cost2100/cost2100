/*
 * \file
 * \brief COST2100 Channel Model Block
 *
 * \author K.K.
 */

#ifndef COST2100_CHANNEL_H_
#define COST2100_CHANNEL_H_

#include <itpp/itbase.h>
#include <itpp/itcomm.h>

#include "COST2100_Specification.h"
#include "processblock.h"

#include <string>


namespace penux {

using namespace itpp;
using namespace std;

class COST2100_Channel : public BasicProcessBlock
{
public:
	//! Default constructor
	COST2100_Channel(string blockNname = "COST2100_Channel",
					string fileName = "channel_spec.xml",
					int idxBS = 0,
					int idxMS = 0,
					int numTx = 1,
					int numRx = 1,
					int numFreqBins = 0,
					int numTime = 0,
					double SNRdB = 10);
	//! Destructor
	virtual ~COST2100_Channel();

	//! run function of the block
	virtual void run(int time);

	//! Set SNR in dB scale
	inline void setSNRdB(double snrdb){ SNRdB=snrdb; };

	//! Load current channel realization
	void loadChannel();

	//! Get channel matrices
	inline Array<cmat> getChannel(){return H;};
	//! Get TX antenna numbers
	inline int getNumTx(){ return numTx; };
	//! Get RX antenna numbers
	inline int getNumRx(){ return numRx; };
	//! Get number of frequency bins
	inline int getNumFreqBins(){ return numFreqBins; };


protected:
	vector< InPort<cvec>* > vpDataInPort; 		//!< input for transmitted symbol streams
	vector< OutPort<cvec>* > vpDataOutPort; 	//!< output streams from channel
	InPort<int> controlPort; 					//!< control port indicate when to create new Channel

	string fileName;							//!< path of configuration file

	int idxBS;									//!< BS index of the simulation test
	int idxMS;									//!< MS index of the simulation test
	int numTx;									//!< TX antenna numbers
	int numRx;									//!< RX antenna numbers
	int numFreqBins;							//!< number of frequency bins
	int numTime;								//!< number of time snapshots
	int currentSnapshot;						//!< current snapshot counter
	double SNRdB;								//!< SNR in dB scale

	Channel_Specification channel;				//!< COST2100 channel model specification
	vector<Array<cmat> > transfer;				//!< loaded transfer functions
	Array<cmat> H;								//!< channel matrix at one snapshot
};

}
#endif /* COST2100_CHANNEL_H_ */
