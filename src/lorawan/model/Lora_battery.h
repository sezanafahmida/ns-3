#include "ns3/end-device-lora-phy.h"
#include "ns3/gateway-lora-phy.h"
#include "ns3/class-a-end-device-lorawan-mac.h"
#include "ns3/gateway-lorawan-mac.h"
#include "ns3/simulator.h"
#include "ns3/log.h"
#include "ns3/pointer.h"
#include "ns3/constant-position-mobility-model.h"
#include "ns3/lora-helper.h"
#include "ns3/node-container.h"
#include "ns3/mobility-helper.h"
#include "ns3/position-allocator.h"
#include "ns3/double.h"
#include "ns3/random-variable-stream.h"
#include "ns3/periodic-sender-helper.h"
#include "ns3/command-line.h"
#include "ns3/network-server-helper.h"
#include "ns3/correlated-shadowing-propagation-loss-model.h"
#include "ns3/building-penetration-loss.h"
#include "ns3/building-allocator.h"
#include "ns3/buildings-helper.h"
#include "ns3/forwarder-helper.h"
#include <algorithm>
#include <ctime>
#include "ns3/node-list.h"
#include <numeric>
#include <math.h>
#include <jsoncpp/json/json.h> 
#include <jsoncpp/json/writer.h> 
#include "ns3/lorawan-mac-header.h"

namespace ns3 {
namespace lorawan {

/**
 * Class representing a battery and green energy source for a LoRaWAN device.
 */

class loraBattery
{
public:
int id;
std::vector<double> SOC; 
double bCap= 7;//0.5; //150.0; //battery capacity (J)
double curSOC=bCap;  //initial charge in battery in J
//double maxSOC=0;
//double minSOC=0;
//double avgSOC=0;
//int lastUpdate=0;
int curDay=1;
//int prevUpdate =0;
double budget =0;
double curE = 0; //energy consumed by transmission (including retransmission attempts) accessed by updateSOC() periodically
double energyTS =0; // energy consumed by transmission (including retransmission attempts) accessed by decideTS() only after a packet has been generated
float s =0; //sampling period
int curSlot =0;
int curP = 0;  //nodes current sampling period number
double lifeSpan =1; //battery lifeSpan
double normAge=1; //normalized age
/*variables from akshar*/ 
int T = 15;  //number of timeslots in a sampling period
int P = 48; //total number of period 
int R = 8;
int Ts = P * T + 1;  //total number of slots
int p = 0;
float Age = 0; // Battery age Get updated age from NonLinear Degradatyion Model
int Retrans = 1; // Get number of retransmissions in previous slot
float energyUsed = 0.0; // Get Energy consumed in previous slot
float estimatedEnergyRequired = energyUsed / Retrans ; // Energy required for transmission 
int result = 1; // 1 if transmission in this period was sucessful else 0
std::vector<float> SelectedSlots; // Increment by 1 when a timeslot t is selected in period p
std::vector<std::vector<float>> SlotSuccess; // Increment by 1 if the selected timeslot t in a period p was sucessfull
std::vector<std::vector<float>> Prob; // Probability for each retransmission in a timeslot
std::vector<float> EstHarvestedE; // Estimated Harvested Energy
std::vector<float> w;  // Weight for energy harvested and energy required in the period
std::vector<float> gamma; // cost function
std::vector<std::vector<int>> X; // Decision Matrix
std::vector<float> Chargelevel; // Energy Level of battery
std::vector<float> H; // Storing energy harvested in a each min.
std::vector<std::vector<std::vector<int>>> c_table ; //table to store CAD values 
std::vector <int> TSVariation; //to store how many times each timeslot was selected by this node;
//std::vector <int> receptionSlots; //to store which packet was received in which timeslot 
std::vector<float> u; //utility
int EstRetrans =1; //estimated retransmission number for this slot
std::vector<float> latency; //packet latencies 
float avgLat=0;
Time lastGeneratedPkt; //holds the time of the most recently generated packet
float maxLat=0;
float buffer = 1; //used to emulate piggyback, holds the latest age until ack is received 
double prr = 0; //prr of node
double avgR = 0;///average retx attempts per node
double txE = 0; //tx energy consumption of node
double utilSum =0;//total gained utility
int success=0; //total successful acks
float avgErr =0; //average error
double totG =0; //total green energy generated
float theta=1;
int h = 0;
bool isLading= false;
template<typename L>
std::vector<float> linspace(L start_in, float end_in, int num_in);

template <typename S>
std::vector<size_t> sort_indexes(const std::vector<S>& v);
float  initialCharge();
float rtx_thr = 0.9;
int  getRetrans(int t);
void SortGamma(int p);
void ProbabilityUpdate(int t);
void initialize(int tsLen,std::vector <double> greenSource,int gsLen,float t,int heuristic) ;
void decideTS2();
void decideTS(); 
double GetTxEnergy(int dr);
double GetToA(int dr);
std::vector <double> ToA={1810.0,987.0,493.0,267.0,154.0,82.0};  //values corresponding to GetToA() function for a packet size 30 bytes(without headers)
};
}
}
