/*
 * This script simulates a complex scenario with multiple gateways and end
 * devices. The metric of interest for this script is the throughput of the
 * network.
 */

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


//#define T 15
//#define P 48
//#define R 8
//#define Ts P*T+1


using namespace ns3;
using namespace lorawan;

NS_LOG_COMPONENT_DEFINE ("ComplexLorawanNetworkExample");

// Network settings
int nDevices = 100;
int nGateways = 1;
double radius = 5000;
double simulationTime = 24; //In hours
int tsLen =120;
int gsLen =  525600;

// Channel model
bool realisticChannelModel = false;

int appPeriodSeconds = 1800;


// Output control
bool print = true;
int overall_attempt =0;
int overall_success =0;
int overall_failed =0;
double total_energy = 0;
int total_packet=0;
std::vector <double> greenSource;

std::vector <double> AggError; 
std::vector <double> AggError2; 
//std::vector<float> EstHarvestedE; // Daily Estimated Harvested Energy
int h =0; //enable or disabe heuristic
double networkLife = 1440; //network lifetimef
int day=1;
float alpha =0.7;
float networkLat = 0; 
/*helper function to retrieve the time on air*/
float theta =1;
double GetToA(int dr) {

std::vector <double> ToA={1810.0,987.0,493.0,267.0,154.0,82.0};  //values corresponding to GetToA() function for a packet size 30 bytes(without headers) 
return ToA[dr]/1000;

}



//functions to calculate Tx energy consumption

double GetTxEnergy(int dr) {
double txCurrent = 29.0/1000;

double txEnergy = txCurrent* GetToA(dr)* 3.3;
return txEnergy ;
}


//CAD variables//

struct CADPacketStatus
{
  Ptr<Packet const> packet;
  uint32_t senderId;
  uint8_t sf ;
  double frequencyMhz;
  double txPower; 
  Time sendTime;
};

typedef std::map<Ptr<Packet const>, CADPacketStatus> CadPacketData;

CadPacketData c_packetTracker;


//latency calc variables

struct n_PacketStatus
{
  Ptr<Packet const> packet;
  uint32_t senderId;
  uint8_t sf ;
  double frequencyMhz;
  double txPower; 
  Time sendTime;
};

typedef std::map<Ptr<Packet const>, n_PacketStatus> n_PacketData;

n_PacketData n_packetTracker;


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
int curDay=0;
//int prevUpdate =0;
double budget =0;
double curE = 0; //energy consumed by transmission (including retransmission attempts) accessed by updateSOC() periodically
double energyTS =0; // energy consumed by transmission (including retransmission attempts) accessed by decideTS() only after a packet has been generated
int s =0; //sampling period
int curSlot =0;
int curP = 0;  //nodes current sampling period
double lifeSpan =1; //battery lifeSpan
double normAge=0; //normalized age
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
std::vector<float> u; //utility
int EstRetrans =1; //estimated retransmission number for this slot
std::vector<float> latency; //packet latencies 
float avgLat=0;
Time lastGeneratedPkt; //holds the time of the most recently generated packet
float maxLat=0;

/********************************UPPER LAYER ******************/ 


template<typename L>

std::vector<float> linspace(L start_in, float end_in, int num_in)
{
	std::vector<float> linspaced;

	float start = static_cast<float>(start_in);
	float end = static_cast<float>(end_in);
	float num = static_cast<float>(num_in);

	if (num == 0) { return linspaced; }
	if (num == 1)
	{
		linspaced.push_back(start);
		return linspaced;
	}

	float delta = (end - start) / (num - 1);

	for (int i = 0; i < num - 1; ++i)
	{
		linspaced.push_back(start + delta * i);
	}
	linspaced.push_back(end); // I want to ensure that start and end
							  // are exactly the same as the input
	return linspaced;
}

template <typename S>
std::vector<size_t> sort_indexes(const std::vector<S>& v) {
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);
	stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });

	return idx;
}

float  initialCharge() {
	float initialEnergy = Chargelevel.back();
		return initialEnergy;
}
	

float rtx_thr = 0.9;
int  getRetrans(int t) {
	/*std::vector<float> Vmaxr; // { 10, 5, 50, 333, 20, 30 };

	for (int r = 0; r < R; r++) {
		Vmaxr.push_back(Prob[t][r]); // Prob[p][t][r]
		//cout << "element: " << Vmaxr[r] << endl;
	}
	//std:: cout << " Max-element: " << std::max_element(Vmaxr.begin(), Vmaxr.end()) - Vmaxr.begin() << std::endl;
	//cout << Vmaxr << endl;
	return max_element(Vmaxr.begin(), Vmaxr.end()) - Vmaxr.begin();*/ 
        for (int r = 0; r < R; r++) {
            if(Prob[t][r]>0.9) return r;
         }
         return 9;
        
}

void SortGamma(int p) {
	std::vector<float> gammaMax;
	for (int t = 0; t < T; t++) {
		gammaMax.push_back(gamma[t]);
	}
	
}

void ProbabilityUpdate(int t) {

                 
                 Prob[t][0] = SlotSuccess[t][0]/ SelectedSlots[t];
                 for (int r = 1; r < R+1; r++){
                 Prob[t][r] =  ((Prob[t][r - 1]) + ((SlotSuccess[t][r] / SelectedSlots[t])));}
                 
		/*for (int t = 0; t < T; t++) {
			for (int r = 0; r < R; r++) {
                               // Prob[t][r] = Prob[t][r] * (1 - Prob[t][i]);
			        // Prob[t][r] = float(SlotSuccess[t][r] / SelectedSlots[t]);
                               //  Prob[t][r] = 1;
			      /*   for(int i = 0; i < r; i++) 
                                  {
				    Prob[t][r] = Prob[t][r] * (1 - Prob[t][i]);
			          }
			          Prob[t][r] = Prob[t][r] * (1- (SlotSuccess[t][r] / SelectedSlots[t]));
                               //std::cout << "Node: "<< id <<  " t " << t <<  " r " << r << " Prob " << Prob[t][r] << "SlotSuccess" << SlotSuccess[t][r] << " SelectedSlots " <<SelectedSlots[t] << "\n";
			}

		//}
                //for (int t = 0; t < T; t++) {
		for (int r = 0; r < R; r++) {
			Prob[t][r] = 1;
			for(int i = 0; i < r; i++) {
				Prob[t][r] = Prob[t][r] * (1 - Prob[t][i]);
			}
			Prob[t][r] = Prob[t][r] * (SlotSuccess[t][r] / SelectedSlots[t]);
                        //std::cout << "Node: "<< id <<  " t " << t <<  " r " << r << " Prob " << Prob[t][r] << "SlotSuccess" << SlotSuccess[t][r] << " SelectedSlots " <<SelectedSlots[t] << "\n";
		//}
	}*/ 
       
     
}

void initialize() 
{
	EstHarvestedE.assign(T*P+T, 0.0);
        
        Chargelevel.assign(T*P,0.0);
	/*for (int i = 0; i < P; i++)
	{
		std::vector<std::vector<float>> v2d;
		std::vector<std::vector<float>> v2df;
		
		std::vector<float> s1d;
		std::vector<float> g1df;*/
                   
		for (int j = 0; j < T; j++)
		{
			std::vector<float> v1d;
			std::vector<float> v1df;
			//g1d.push_back(0);
			//s1d.push_back(1);
			//g1df.push_back(0.0);
                        float temp =5.0;
                        SelectedSlots.push_back(temp);
			for (int k = 0; k < temp; k++)
			{
				v1d.push_back(1);
				v1df.push_back(1/temp);
			}
                        for (int l = temp+1; l < R+1; l++)
			{
				v1d.push_back(0);
				v1df.push_back(0);
			}
                        
			//v2d.push_back(v1d);
			//v2df.push_back(v1df);
			//EstHarvestedE.push_back(0.0);
			//Chargelevel.push_back(0.0);
                        SlotSuccess.push_back(v1d);
                        Prob.push_back(v1df);
		}

		
		//X.push_back(g1d); // g1df
		
		
	//}

        for (int p=0;p<P;p++) 
        { 
          std::vector<int> g1d;
          for(int t=0;t<T;t++)
          {
           
            g1d.push_back(0);
            
          
          }
          X.push_back(g1d);
        }
	std::string value;
	int day = curDay% 365;
	int Cmin = day * 1440;
        
	for (int i = 0; i < P*T; i++)
	{              
		double greenEnergy = 0.0;
		for (int q = 0; q < (tsLen/60) ; q++)
		{	
			greenEnergy+= greenSource.at(Cmin%gsLen); //GetTxEnergy(0) * 3	
			Cmin++;  
		}
		EstHarvestedE[i] = greenEnergy;
              //  std::cout<< "Energy harvested in Decide TS() for timeslot " << i << " " << EstHarvestedE[i] << "\n";
	}
	//curDay++;
       // std::cout<< " Node " << id << " day " << curDay << std::endl;
       for (int i =0;i<T;i++) 
       {
          
         TSVariation.push_back(0);
       
       }

      u = linspace(1, 0.1, T);  //std::vector<float>(T, 1.0);
}







void decideTS()
{
	float Ebat = bCap*theta; //0.008; // Wh**** product of e_required
	//float Edod = 0.2 * Ebat;// Depth of Discharege
	float Emax = 1 * Ebat; // Max charge to store in battery
	float Emin = 0 * Ebat;
        int p = curP-1;   
	//normAge=1;
	
	int retransNo;


	Chargelevel[0] = curSOC;//0.003; // initialCharge();
        Chargelevel[p*T] = curSOC;
         
        w.assign(T, 0.0);
	gamma.assign(T, 0.0);
	for (int t = 0; t < T; t++) {
		retransNo = getRetrans(t)+1;
                if(retransNo>8) retransNo = 1000; //in this case the packet failed
		w[t] = (EstHarvestedE[p*T + t] - std::max(EstHarvestedE[p*T + t] , estimatedEnergyRequired * retransNo)) /  GetTxEnergy(0);

                //(*std::max_element(greenSource.begin(),greenSource.end())*(tsLen/60));
		if(h == 1) gamma[t] =  (u[t] +  ( normAge* w[t]));
                else if(h==2) gamma[t] = normAge* w[t];
                //  std::cout <<"t " << t << " U " << u[t] << " W " << w[t] << " E_req " << estimatedEnergyRequired * retransNo << " E_h " << EstHarvestedE[p*T + t]  << " Gamma " << gamma[t] << std::endl;
	}


	std::vector<float> v;
	for (int t = 0; t < T; t++) {
		v.push_back(gamma[t]);
	}
        std::vector<long unsigned int>temp_vec = sort_indexes(v);
	for (int i : temp_vec) {
		int t = 0;
		bool stop = false;
		while (t < T && stop == false) {
			retransNo = getRetrans(t) + 1; // as r starts from 0
                        if(retransNo>8) retransNo = 1000;
		        budget = fabs(EstHarvestedE[p*T + t] + (estimatedEnergyRequired * retransNo));
                       // std::cout << "EstE "<< EstHarvestedE[T * p + t] << "EstR "<< estimatedEnergyRequired << "Retr " << retransNo << "\n";
			if (i == t) {
				if ((Chargelevel[T * p + t] - (estimatedEnergyRequired * retransNo)) >= Emin) {
					X[p][t] = 1;
                                      //  std:: cout <<" Node "<< id << " period " << p << "," << " timeslot "<< t << "--" << Chargelevel[T * p + t] << "     X-- " << X[p][t] << "  -- REQ: " << estimatedEnergyRequired * retransNo  << " || E: " << EstHarvestedE[T * p + t] <<  " || retransNo: " << retransNo << " --  w " << w[t] << " || Gamma: "<< gamma[t] << std::endl;
					if (Chargelevel[T * p + t] + EstHarvestedE[p*T + t] - estimatedEnergyRequired * retransNo < Emax) {
						Chargelevel[T * p + t + 1] = Chargelevel[T * p + t] + EstHarvestedE[p*T + t] - estimatedEnergyRequired * retransNo;
					}
					else {
						Chargelevel[T * p + t + 1] = Emax;
					}
				}
				else {
					stop = true;
				}
			}
			else {
				if (Chargelevel[T * p + t] + EstHarvestedE[p*T + t] - estimatedEnergyRequired * retransNo < Emax) {
					Chargelevel[T * p + t + 1] = Chargelevel[T * p + t] + EstHarvestedE[p*T + t];
				}
				else {
					Chargelevel[T * p + t + 1] = Emax;
				}
			}
			//std:: cout <<" Node "<< id << " period " << p << "," << " timeslot "<< t << "--" << Chargelevel[T * p + t] << "     X-- " << X[p][t] << "  -- REQ: " << estimatedEnergyRequired * retransNo  << " || E: " << EstHarvestedE[T * p + t] <<  " || retransNo: " << retransNo << " --  w " << w[t] << " || Gamma: "<< gamma[t] << std::endl;
			t++;
		}
		if (X[p][i] == 1) {
                        curSlot = i;
                        TSVariation[i]+=1;
                        EstRetrans = getRetrans(curSlot)+1;
                        //if(i==0) std::cout<< "Node " << id << "slot" << curSlot << " Est retrans" << EstRetrans << std::endl;
			break;
		}
		//std::cout << std::endl << std::endl;
	}
	//system("pause>0");
}
}; //end of classf

//function to calculate max,min and avg SOC of all batteries


/*void getSummary(loraBattery bList[])
{
   std::ofstream bTrace;
   bTrace.open("battery_trace.csv");
   
 for(int i=0;i<nDevices;i++)
 {
  double max=-99999;
  double min = 99999;
  double sum =0;
  for(int j =0;j<bList[i].SOC.size();j++) 
  {
    if(bList[i].SOC.at(j) >max) max= bList[i].SOC.at(j);
    if(bList[i].SOC.at(j)< min) min =  bList[i].SOC.at(j);
    sum+= bList[i].SOC.at(j);
  }
  bList[i].maxSOC=max;
  bList[i].minSOC=min;
  bList[i].avgSOC = sum/bList[i].SOC.size();
  //std::cout<< i << " " << bList[i].maxSOC << " " << bList[i].minSOC <<" "<<bList[i].avgSOC << std::endl;
  bTrace << bList[i].maxSOC << "," << bList[i].minSOC <<" , "<<bList[i].avgSOC << std::endl;
 }

}*/


//function to count total successful packets and transmission attempts
void packetCounter (uint32_t id,uint8_t reqTx, bool success, Time firstAttempt, Ptr<Packet> packet)
{ 

if(success)
{   
    overall_attempt+= reqTx;
    overall_success++;
}
else
{
    
    overall_attempt+=reqTx;
    overall_failed++;
    
}
}


//trace callback function to log packets into CAD data structure once it has been sent
void logCAD(Ptr<Packet const> packet, uint32_t edId , uint8_t sf, double frequencyMhz, double txPower)
{
   //std::cout << "packet sent by device " << packet << std::endl;   
      CADPacketStatus data;
      data.packet = packet;
      data.sendTime = Simulator::Now ();
      data.senderId = edId;
      data.sf = sf;
      data.frequencyMhz= frequencyMhz;
      data.txPower = txPower;
      c_packetTracker.insert (std::pair<Ptr<Packet const>, CADPacketStatus> (packet, data));

}

//trace callback function to log packets into data structure once it has been sent
void logPacket(Ptr<Packet const> packet, uint32_t edId )
{
   //std::cout << "packet sent by device " << packet << std::endl;   
      n_PacketStatus data;
      data.packet = packet;
      data.sendTime = Simulator::Now ();
      data.senderId = edId;
      n_packetTracker.insert (std::pair<Ptr<Packet const>, n_PacketStatus> (packet, data));
      

}



void
GWPacketReceptionCallback (loraBattery bList[], Ptr<Packet const> packet, uint32_t gwId)
{
 // if (IsUplink (packet))
   // {
      // Remove the successfully received packet from the list of sent ones
    //  NS_LOG_INFO ("PHY packet " << packet<< " was successfully received at gateway "<< gwId);

      std::map<Ptr<Packet const>, n_PacketStatus>::iterator it = n_packetTracker.find (packet);
      if(it!=n_packetTracker.end())
      {
       Time delay = (Simulator::Now() - (*it).second.sendTime);
       int nodeId = (*it).second.senderId;
       bList[nodeId].latency.push_back(delay.GetSeconds());
       n_packetTracker.erase(it);
      }
//(*it).second.outcomes.insert (std::pair<int, enum PhyPacketOutcome> (gwId,
                                                                          // RECEIVED));
   // }
}

void FailedPacketDelay(loraBattery bList[])
{

for (auto it = n_packetTracker.begin ();it != n_packetTracker.end ();++it)
{
int delay = 1000;
int nodeId = (*it).second.senderId;
bList[nodeId].latency.push_back(delay);
}


}

//trace callback to remove packets from CAD DS once it has been received at the gateway (regardless of succesful or not)
void removePacket(Ptr<Packet const> packet,uint32_t gwId )
{

//std::cout << "packet received by gateway " << packet << std::endl;
     std::map<Ptr<Packet const>, CADPacketStatus>::iterator it = c_packetTracker.find (packet);
      if (it != c_packetTracker.end ())
       {
          //std::cout<< "Found the packet: removing" <<std::endl;
          c_packetTracker.erase(it);
        }
     else
        {
          std::cout<< "Packet not found in tracker" <<std::endl;
        }

}


void calcLatency(loraBattery bList[],uint8_t reqTx, bool success, Time firstAttempt, Ptr<Packet> packet)
{

uint32_t id = Simulator::GetContext();
float lat = (Simulator::Now().GetSeconds() - bList[id].lastGeneratedPkt.GetSeconds());
bList[id].latency.push_back(lat);
}



//function to provide feedback to upper layer
void feedback (loraBattery bList[],uint8_t reqTx, bool success, Time firstAttempt, Ptr<Packet> packet)
{ 
uint32_t id = Simulator::GetContext();
int t = bList[id].curSlot;

//float lat = (Simulator::Now().GetSeconds() - bList[id].lastGeneratedPkt.GetSeconds());
//bList[id].latency.push_back(lat);
bList[id].Retrans = int(reqTx);
//bList[id].energyUsed = bList[id].curE;
bList[id].estimatedEnergyRequired = alpha* (bList[id].energyUsed / bList[id].Retrans) + (1- alpha)*bList[id].estimatedEnergyRequired  ;
//if(id==2) std::cout << "Node " << id <<" energy used " << bList[id].energyUsed << " R " << bList[id].Retrans << " in period " << p <<  " timeslot "<< t << " outcome " << success << " Time " << Simulator::Now().GetMinutes() <<"\n";
//std::cout << " Node " << id << " energy used " << bList[id].energyUsed <<" budget " <<  bList[id].budget << "Tx attempt " << int(reqTx) << " outcome " << success << std::endl;
bList[id].SelectedSlots[t]+=1;

if(success)
{  

  bList[id].SlotSuccess[t][int(reqTx-1)] +=1;  
}

/*report the failure*/ 
else  
{
    
 bList[id].SlotSuccess[t][bList[id].R] += 1;  
  
}
bList[id].ProbabilityUpdate(t);
}

//prints the summary output file
void printToFile(LoraPacketTracker tracker, Time appStopTime)
{

/*create the output file*/

   std::ofstream logfile;
   logfile.open ("sim_results_0122.csv",std::ios_base::app);

logfile << h << "," << nDevices << "," << radius << "," << tracker.CountMacPacketsGlobally (Seconds (0), appStopTime + Hours (1)) << " ," << overall_attempt << " , " << overall_success << " , " << overall_failed << ", "<< total_energy <<" ," << total_packet<< ", " << double(overall_attempt)/double(total_packet) << "," <<networkLife <<  "," << networkLat << "," << day <<"\n";
}

//prints the timeslot variation of the network
void printTSVariation(loraBattery bList[])
{

std::ofstream tsfile;
tsfile.open("timeslot.csv");
 for(int i =0;i<nDevices;i++)  
 {  
    //tsfile<< " Node " << i << " Timeslot ";
    for(int j=0;j<bList[i].T ;j++)
    {
     tsfile << bList[i].SelectedSlots[j] << "," ;
    }
    tsfile << std::endl;
 }
}

//prints the probability distribution of retx

void printProb(int id, loraBattery bList[])

{
std::stringstream ss;
ss<<"Prob/prob_"<< id << ".csv";
std::ofstream probfile;
probfile.open(ss.str());

for (int t=0;t<bList[id].T;t++) 
{

 for (int r=0;r<bList[id].R;r++)
 {

  probfile<< bList[id].Prob[t][r] << " ";
  
 }
 probfile << std::endl;
 
}


}

void getProb(int id, loraBattery bList[])
{

std::stringstream ss;
ss<<"Prob/prob_"<< id << ".csv";

std::ifstream input;
input.open(ss.str());


if(!input.is_open())
{
std::cout<< "prob.csv not found, prob starting from 100%, time starting from default \n";
for(int i =0;i<nDevices;i++) 
{
bList[i].curSOC = bList[i].bCap;
}
return;
}


int index=0;
double SOC=0;
int time=0;
while(input>>SOC>>time)
{
if(index<nDevices) 
{
//std::cout<<index<<std::endl;
bList[index].curSOC = SOC;
//bList[index].prevUpdate = time;
//bList[index].lastUpdate = bList[index].prevUpdate;

}
index++;
}


}


//prints the SOC trace (in perc of original Bcap) for a battery
void printTrace(int id, loraBattery Blist[])
{

std::stringstream ss;
ss<<"Trace/trace_"<< id << ".csv";
std::ofstream tracefile;
tracefile.open(ss.str());

std::vector<double> tempSOC = Blist[id].SOC;

for(double val:tempSOC)
{
tracefile<< val/Blist[id].bCap << std::endl;
}

}

void printLastSOC(loraBattery Blist[])
{

std::ofstream out;
out.open("SOC.csv");

for(int i=0;i<nDevices;i++) {

//std::vector<double> tempSOC = Blist[i].SOC;

//int n = tempSOC.size();

out<< Blist[i].curSOC <<std::endl;

}
out.close();

}


void printAvgError()
{
std::ofstream out;
out.open("AvgError.csv");

for(int k=0;k<AggError.size();k++)
    {
    out<< AggError[k] <<","<<AggError2[k] << "\n";  
    }

out.close();
}

//reads the green energy trace from file and returns a vector

std::vector <double> readFile()
{
std::ifstream input;
std::vector<double> temp_trace;
double temp=0;
input.open("OneMin_EnergyHarvested.csv");
if(!input.is_open())
{
std::cout<< "OneMin_EnergyHarvested.csv not found. Green energy source not initialized\n";
std::vector<double> temp(gsLen, 0.0);
return temp;

}

while(input>>temp)
{
//std::cout<<temp<< std::endl;
temp_trace.push_back(temp*7);
}
input.close();
return temp_trace;
}


//reads the battery degradation from file and updates battery_age

void updateLifeSpan(loraBattery bList[])
{

std::ifstream input;
input.open("output.csv");

if(!input.is_open())
{
std::cout<< "output.csv not found, battery lifespan starting from 1\n";
return;
}
double d=0;
int index=0;
while((input>>d)&&(index<nDevices))
{
//bList[index].lifeSpan -= d;
bList[index].Age += d;
//std::cout<<"i:"<<index<< " "<<d<<std::endl;
//if(index==b->id) break;
index++;
}
input.close();
}

/*normalize the life span of the nodes by the maximum lifespan*/
void normAge(loraBattery bList[])
{

//double max = -1000000;
//double min =100000;
std::vector <double> temp;
for (int i =0;i<nDevices;i++)
{
//if(bList[i].Age > max) max = bList[i].Age; 
temp.push_back(bList[i].Age);
}

std::sort(temp.begin(),temp.end());

double min = temp.at(0);
double max = temp.at(nDevices-1);

for (int i =0;i<nDevices;i++)
{
if(max>0) bList[i].normAge = (bList[i].Age - min) /(max -min) ;
else bList[i].normAge =1;
std::cout << "Node " << i << " Age " <<  bList[i].Age << " " <<max << " " << bList[i].normAge << "\n";
}

}


void getSOC(loraBattery bList[])
{
std::ifstream input;
input.open("SOC.csv");

if(!input.is_open())
{
std::cout<< "SOC.csv not found, SOC starting from 100%, time starting from zero \n";
for(int i =0;i<nDevices;i++) 
{
bList[i].curSOC = bList[i].bCap;
}
return;
}


int index=0;
double SOC=0;
int time=0;
while(input>>SOC>>time)
{
if(index<nDevices) 
{
//std::cout<<index<<std::endl;
bList[index].curSOC = SOC;
//bList[index].prevUpdate = time;
//bList[index].lastUpdate = bList[index].prevUpdate;

}
index++;
}
}

void getDay(loraBattery blist[])
{
	std::ifstream input;
	input.open("Day.csv");

	if(!input.is_open())
	{
		std::cout << "Day.csv not found, Day starting from 0 \n";
		return;
	}

	int Day = 0;

        int index=0;
	while(input >> Day)
	{
		if (index < nDevices)
		{
			blist[index].curDay = Day;
			}
	index++;
	}
}

void updateDay(loraBattery Blist[])
{

std::ofstream out;
out.open("Day.csv");

for(int i=0;i<nDevices;i++) {

out<< Blist[i].curDay <<std::endl;

}
out.close();

}

//counts new packets
void countNewPacket(loraBattery bList[], Ptr<Packet const> packet){

  total_packet+=1;
  int id = Simulator::GetContext(); 
 // if(id == 2)  std::cout << "Count new packet called "<< std::endl;
  //error estimation
  if((bList[id].curP >0) && (id ==7)) { //std::cout <<" id " << id <<" "<< bList[id].EstRetrans << " " << bList[id].Retrans << "\n";
  float diff = bList[id].EstRetrans - bList[id].Retrans;
  AggError[bList[id].curP-1]= diff;}
  if((bList[id].curP >0) && (id==23)) { //std::cout <<" id " << id <<" "<< bList[id].EstRetrans << " " << bList[id].Retrans << "\n";
  float diff = bList[id].EstRetrans - bList[id].Retrans;
  AggError2[bList[id].curP-1]= diff;}
  bList[id].lastGeneratedPkt = Simulator::Now();
  bList[id].curP += 1;
  bList[id].curE =0;
  bList[id].energyUsed = 0;
 

}


//sets timeslot
  void getTs(loraBattery bList[],Ptr<Packet const> packet){

  
//  Ptr<RandomVariableStream> dr_high = CreateObjectWithAttributes<UniformRandomVariable> ( "Min", DoubleValue (2), "Max", DoubleValue (4)); 
//  Ptr<RandomVariableStream> dr_low = CreateObjectWithAttributes<UniformRandomVariable> ( "Min", DoubleValue (0), "Max", DoubleValue (1)); //max  included
  int id = Simulator::GetContext();
  if(bList[id].curP>=bList[id].P) { return;}
  Ptr<Node> node = NodeList::GetNode(id);
  Ptr<ns3::lorawan::LorawanMac> edMac = node->GetDevice (0)->GetObject<ns3::lorawan::LoraNetDevice> ()->GetMac ();
  Ptr<ns3::lorawan::EndDeviceLorawanMac> edLorawanMac = edMac->GetObject<ns3::lorawan::EndDeviceLorawanMac> ();
  
  bList[id].decideTS();

 // if(id == 2) std::cout << "Calling Decide TS for period " << bList[id].curP  << " slot " << bList[id].curSlot<< std::endl;
  edLorawanMac->ts = Minutes(bList[id].curSlot * (tsLen/60));
  /*int dr =0;
  if(bList[id].lifeSpan < 0.8) 
  {
  dr = dr_low->GetInteger();
  }
  else
  {
  dr = dr_high->GetInteger();
  }
 // edLorawanMac->SetDataRate(dr);*/
 // std::cout << "Node " << id << " transmitting in timeslot " << bList[id].curSlot << " Min " << bList[id].curSlot * (tsLen/60) << " dr " << dr << std::endl;

}

   //emulates cad
  void doCad(loraBattery bList[],Ptr<Packet const> packet, double freq, uint8_t sf){
  //std::cout<< "CAD "<< std::endl;
  Ptr<RandomVariableStream> backoff_low = CreateObjectWithAttributes<UniformRandomVariable> ( "Min", DoubleValue (0), "Max", DoubleValue (5));
  Ptr<RandomVariableStream> backoff_high = CreateObjectWithAttributes<UniformRandomVariable> ( "Min", DoubleValue (5), "Max", DoubleValue (10));
  Ptr<RandomVariableStream> random_sf = CreateObjectWithAttributes<UniformRandomVariable> ( "Min", DoubleValue (0), "Max", DoubleValue (5));
  int id = Simulator::GetContext();
  Ptr<Node> node = NodeList::GetNode(id);
  Ptr<ns3::lorawan::LorawanMac> edMac = node->GetDevice (0)->GetObject<ns3::lorawan::LoraNetDevice> ()->GetMac ();
  Ptr<ns3::lorawan::EndDeviceLorawanMac> edLorawanMac = edMac->GetObject<ns3::lorawan::EndDeviceLorawanMac> (); 

  
 bool cadF = false;
 for (auto itCad = c_packetTracker.begin (); itCad != c_packetTracker.end (); ++itCad)
  {
  if(((*itCad).second.packet!=packet) && ((*itCad).second.sf == sf) && ((*itCad).second.frequencyMhz == freq))
   {
     cadF=true;
     break;
   }
  } 
  double bo =0;
  int dr= random_sf->GetInteger();
  if(cadF){
   
     if(bList[id].lifeSpan < 0.8 )  
     {
      bo = backoff_high->GetValue();
      edLorawanMac->SetDataRate(random_sf->GetInteger());
     }
     else 
     {
     bo = backoff_low->GetValue();
     edLorawanMac->SetDataRate(dr);
     }
    // std::cout << "packet found in channel " << freq << " sf " << int(sf) <<" node "<< id <<  " backing off for "  << bo << " Seconds "<< " SF " << dr << std::endl;
     edLorawanMac->cadBo = Seconds(bo);
  }
  else 
  {
    edLorawanMac->cadBo = Seconds(0);
  }
}


//checks budget
void checkBudget(loraBattery bList[],Ptr<Packet const> packet)
{

  int id = Simulator::GetContext();
  Ptr<Node> node = NodeList::GetNode(id);
  Ptr<ns3::lorawan::LorawanMac> edMac = node->GetDevice (0)->GetObject<ns3::lorawan::LoraNetDevice> ()->GetMac ();
  Ptr<ns3::lorawan::ClassAEndDeviceLorawanMac> edLorawanMac = edMac->GetObject<ns3::lorawan::ClassAEndDeviceLorawanMac> ();
  if( bList[id].energyUsed > bList[id].budget) 
  {
   edLorawanMac->haveBudget = false;
  // std::cout << " Node " << id << " crossed budget " << "E " <<  bList[id].energyUsed << " budget " << bList[id].budget << std::endl;
  }

}

//changes theta
void changeTheta(double v)
{
theta =v;

}


//calculates energy consumption  called everytime the node sends a packet regardless of retransmission or new packet

void energyNewPkt(loraBattery bList[], Ptr<Packet const> packet, uint32_t id ) {
  
//  int id = Simulator::GetContext();
  Ptr<Node> node = NodeList::GetNode(id);
  Ptr<ns3::lorawan::LorawanMac> edMac = node->GetDevice (0)->GetObject<ns3::lorawan::LoraNetDevice> ()->GetMac ();
  Ptr<ns3::lorawan::EndDeviceLorawanMac> edLorawanMac = edMac->GetObject<ns3::lorawan::EndDeviceLorawanMac> ();
  int dr = edLorawanMac->GetDataRate();
  total_energy += GetTxEnergy(dr);  
  double curEnergy= GetTxEnergy(dr);
  bList[id].curE += curEnergy;
  
  bList[id].energyUsed +=curEnergy; 
 // if(id ==2) std::cout << " Energy used updated " << bList[id].energyUsed << " Time " << Simulator::Now().GetMinutes() <<"\n" ;
 /* int curMin = int(Simulator::Now().GetMinutes());
 // std::cout<<curMin<<std::endl;
  double greenEnergy =0;
  for (int i=bList[id].lastUpdate;i<curMin && curMin<greenSource.size();i++) 
  {
   greenEnergy+= greenSource.at(i);
  }
  bList[id].curSOC = bList[id].curSOC + greenEnergy - curEnergy;  
  bList[id].SOC.push_back(bList[id].curSOC);
  bList[id].lastUpdate=curMin;*/
}



//updates the SOC trace at regular intervals
void updateSOC(loraBattery* b) {
  
  int curMin = int(Simulator::Now().GetMinutes());  // b->prevUpdate + int(Simulator::Now().GetMinutes());
  int tsMin = tsLen/60;
  int tsIndex = floor(curMin/tsMin) ;
 // std::cout<<curMin<<std::endl;
 // b->lastUpdate = b->prevUpdate+curMin;
/*  double greenEnergy =0;
  for (int i=b->lastUpdate;i<curMin && curMin<greenSource.size();i++) 
  {
  greenEnergy+= greenSource.at(i%gsLen);
 // std::cout<<"green energy " <<  i%gsLen << " " << greenSource.at(i%gsLen) << std::endl;
   
  
  } */
  double greenEnergy = b->EstHarvestedE[tsIndex];
 // std::cout<< greenEnergy/3000 <<std::endl;
// std::cout<< "Energy harvested for node " << b->id << " in Update SOC() for timeslot " << tsIndex << " " << greenEnergy << "\n";
  Ptr<Node> node = NodeList::GetNode(b->id);
  Ptr <ns3::lorawan::PeriodicSender> ps = node->GetApplication(0)->GetObject<ns3::lorawan::PeriodicSender>();
  double curCap = std::min(b->bCap*theta, b->bCap* (1-b->Age));
  b->curSOC = std::min((b->curSOC + greenEnergy - b->curE),curCap ); 
//  std::cout<< "CurE for node " << b->id <<  " in Update SOC() for timeslot " << tsIndex << " " << b->curE << "\n"; 
//  std::cout<< "CurSOC for node " << b->id <<  " in Update SOC() for timeslot " << tsIndex << " " << b->curSOC << "\n";
  if(b->curSOC <=0) 
  {
  
  ps->alive = false;
  b->curSOC =0;
  if(Simulator::Now().GetMinutes() < networkLife) networkLife = Simulator::Now().GetMinutes();
 // std::cout << "Battery " << b->id << "died on min " <<curMin << " Network life " << networkLife << std::endl; 
  
  }
  else ps->alive =true;
  
  b->SOC.push_back(b->curSOC);
//  b->lastUpdate = curMin;
  b->curE =0;
  Simulator::Schedule(Seconds(tsLen), &updateSOC, b);
}

void writeToJson(loraBattery* bList)
{
std::ofstream file;
file.open("node_info.json");

Json::Value node_vec(Json::arrayValue);

for(int i=0;i<nDevices;i++)
{

Json::Value temp;
Json::Value vec(Json::arrayValue);

for(float s: bList[i].SelectedSlots)
{
vec.append(Json::Value(s));
}

for (int t=0;t<bList[i].T;t++)
{
Json::Value vec2(Json::arrayValue);
for(float ss: bList[i].SlotSuccess[t])
{
vec2.append(Json::Value(ss));
}
std::stringstream ss;
ss<<"slot"<< t;
temp[ss.str()] = vec2;
}

temp["id"] = bList[i].id;
temp["soc"] = bList[i].curSOC;
temp["day"] =bList[i].curDay;
temp["selectedSlots"] = vec;
temp["avgLatency"] = bList[i].avgLat;
temp["maxLatency"] = bList[i].maxLat;
temp["estE"] = bList[i].estimatedEnergyRequired;
//temp["success"] = vec2;
node_vec.append(temp);
}
file<< node_vec << std::endl;
file.close();
}



void readFromJson(loraBattery* bList)
{
std::ifstream infoFile("node_info.json",std::ifstream::binary);
if(!infoFile.is_open())
{
std::cout<< "node_info.json not found, starting with default values\n";
return;
}
Json::Reader reader;
Json::Value node_info;
bool parseSuccess = reader.parse(infoFile, node_info,false);

if(parseSuccess){

for(Json::Value::ArrayIndex i =0;i!=node_info.size() && i<nDevices ;i++)
{

Json::Value temp = node_info[i];
bList[i].curSOC = temp["soc"].asFloat();
bList[i].curDay = temp["day"].asInt()+1;
bList[i].avgLat = temp["avgLatency"].asFloat();
bList[i].maxLat = temp["maxLatency"].asFloat();
bList[i].estimatedEnergyRequired = temp["estE"].asFloat();
Json::Value vec = temp["selectedSlots"];

for(Json::Value::ArrayIndex j =0;j!=vec.size();j++)
{
bList[i].SelectedSlots[j] = vec[j].asInt();
}

for (int t =0;t<bList[i].T;t++)
{
std::stringstream ss;
ss<< "slot"<< t;
Json::Value vec2 = temp[ss.str()];

for(Json::Value::ArrayIndex r =0;r!=vec2.size();r++)
{
bList[i].SlotSuccess[t][r] = vec2[r].asInt();
}

}

}


}

infoFile.close();
}


void printAvgLatency(loraBattery bList[])
{
std::ofstream out;
out.open("AvgLatency.csv");
float networkSum =0;
float networkMax =-9999;
for (int i =0;i<nDevices;i++)
{
float sum=0;
float max=-99999;
float avg = 1;
for( int k=0;k<bList[i].latency.size();k++)
{
sum+= bList[i].latency[k];
if(bList[i].latency[k]>max) {max = bList[i].latency[k];}
}
sum+= bList[i].avgLat*bList[i].latency.size()*bList[i].curDay; //taking into account the previous day's latency
avg = sum/(bList[i].latency.size()*(bList[i].curDay+1));
bList[i].avgLat = avg;
if(max>bList[i].maxLat) bList[i].maxLat = max;
out << avg <<","<<bList[i].maxLat << std::endl;
networkSum+=sum;
}
networkLat = networkSum/(total_packet*(day));
out.close();
}

int
main (int argc, char *argv[])
{

  CommandLine cmd;
  cmd.AddValue ("nDevices", "Number of end devices to include in the simulation", nDevices);
  cmd.AddValue ("radius", "The radius of the area to simulate", radius);
  cmd.AddValue ("simulationTime", "The time for which to simulate", simulationTime);
  cmd.AddValue ("appPeriod",
                "The period in seconds to be used by periodically transmitting applications",
                appPeriodSeconds);
  cmd.AddValue ("print", "Whether or not to print various informations", print);
 // cmd.AddValue ("ts", "chosen timeslot for tx", chosenTs);
  cmd.AddValue ("l", "timeslot length for tx ", tsLen);
  cmd.AddValue ("h", "enable heuristic ", h);
   cmd.AddValue ("day", "day of the simulation", day);
  cmd.AddValue("theta","SOC cap", theta);
  cmd.Parse (argc, argv);

  // Set up logging
  LogComponentEnable ("ComplexLorawanNetworkExample", LOG_LEVEL_ALL);
  // LogComponentEnable("LoraChannel", LOG_LEVEL_INFO);
 // LogComponentEnable("LoraPhy", LOG_LEVEL_ALL);
  // LogComponentEnable("EndDeviceLoraPhy", LOG_LEVEL_ALL);
  // LogComponentEnable("GatewayLoraPhy", LOG_LEVEL_ALL);
  // LogComponentEnable("LoraInterferenceHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("LorawanMac", LOG_LEVEL_ALL);
  // LogComponentEnable("EndDeviceLorawanMac", LOG_LEVEL_ALL);
  // LogComponentEnable("ClassAEndDeviceLorawanMac", LOG_LEVEL_ALL);
  // LogComponentEnable("GatewayLorawanMac", LOG_LEVEL_ALL);
  // LogComponentEnable("LogicalLoraChannelHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("LogicalLoraChannel", LOG_LEVEL_ALL);
  // LogComponentEnable("LoraHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("LoraPhyHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("LorawanMacHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("PeriodicSenderHelper", LOG_LEVEL_ALL);
  // LogComponentEnable("PeriodicSender", LOG_LEVEL_ALL);
  // LogComponentEnable("LorawanMacHeader", LOG_LEVEL_ALL);
  // LogComponentEnable("LoraFrameHeader", LOG_LEVEL_ALL);
  // LogComponentEnable("NetworkScheduler", LOG_LEVEL_ALL);
  // LogComponentEnable("NetworkServer", LOG_LEVEL_ALL);
  // LogComponentEnable("NetworkStatus", LOG_LEVEL_ALL);
  // LogComponentEnable("NetworkController", LOG_LEVEL_ALL);

   greenSource = readFile();
   loraBattery bList[nDevices];
  
   for (int i=0;i<nDevices;i++)
   {
    bList[i].id = i;
    
   }
   
   updateLifeSpan(bList);  //update the degradation from file

  // getDay(bList); 		//get the day of the year
   normAge(bList);    //normalize
    
  /***********
   *  Setup  *
   ***********/
  
  // Create the time value from the period
  Time appPeriod = Seconds (appPeriodSeconds);

  // Mobility
  MobilityHelper mobility;
  mobility.SetPositionAllocator ("ns3::UniformDiscPositionAllocator", "rho", DoubleValue (radius),
                                 "X", DoubleValue (0.0), "Y", DoubleValue (0.0));
  mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");

  /************************
   *  Create the channel  *
   ************************/

  // Create the lora channel object
  Ptr<LogDistancePropagationLossModel> loss = CreateObject<LogDistancePropagationLossModel> ();
  loss->SetPathLossExponent (3.76);
  loss->SetReference (1, 7.7);

  if (realisticChannelModel)
    {
      // Create the correlated shadowing component
      Ptr<CorrelatedShadowingPropagationLossModel> shadowing =
          CreateObject<CorrelatedShadowingPropagationLossModel> ();

      // Aggregate shadowing to the logdistance loss
      loss->SetNext (shadowing);

      // Add the effect to the channel propagation loss
      Ptr<BuildingPenetrationLoss> buildingLoss = CreateObject<BuildingPenetrationLoss> ();

      shadowing->SetNext (buildingLoss);
    }

  Ptr<PropagationDelayModel> delay = CreateObject<ConstantSpeedPropagationDelayModel> ();

  Ptr<LoraChannel> channel = CreateObject<LoraChannel> (loss, delay);

  /************************
   *  Create the helpers  *
   ************************/

  // Create the LoraPhyHelper
  LoraPhyHelper phyHelper = LoraPhyHelper ();
  phyHelper.SetChannel (channel);

  // Create the LorawanMacHelper
  LorawanMacHelper macHelper = LorawanMacHelper ();

  // Create the LoraHelper
  LoraHelper helper = LoraHelper ();
  helper.EnablePacketTracking (); // Output filename
  // helper.EnableSimulationTimePrinting ();

  //Create the NetworkServerHelper
  NetworkServerHelper nsHelper = NetworkServerHelper ();

  //Create the ForwarderHelper
  ForwarderHelper forHelper = ForwarderHelper ();

  /************************
   *  Create End Devices  *
   ************************/

  // Create a set of nodes
  NodeContainer endDevices;
  endDevices.Create (nDevices);

  // Assign a mobility model to each node
  mobility.Install (endDevices);

  // Make it so that nodes are at a certain height > 0
  for (NodeContainer::Iterator j = endDevices.Begin (); j != endDevices.End (); ++j)
    {
      Ptr<MobilityModel> mobility = (*j)->GetObject<MobilityModel> ();
      Vector position = mobility->GetPosition ();
      position.z = 1.2;
      mobility->SetPosition (position);
    }

  // Create the LoraNetDevices of the end devices
  uint8_t nwkId = 54;
  uint32_t nwkAddr = 1864;
  Ptr<LoraDeviceAddressGenerator> addrGen =
      CreateObject<LoraDeviceAddressGenerator> (nwkId, nwkAddr);

  // Create the LoraNetDevices of the end devices
  macHelper.SetAddressGenerator (addrGen);
  phyHelper.SetDeviceType (LoraPhyHelper::ED);
  macHelper.SetDeviceType (LorawanMacHelper::ED_A);
  helper.Install (phyHelper, macHelper, endDevices);

  // Now end devices are connected to the channel

  // Connect trace sources (if any) , enable confirmed traffic
 
  for (NodeContainer::Iterator j = endDevices.Begin (); j != endDevices.End (); ++j)
    {
      Ptr<Node> node = *j;
      int devID = node->GetId();
      Ptr<LorawanMac> edMac1 = node->GetDevice (0)->GetObject<LoraNetDevice> ()->GetMac ();
      Ptr<ClassAEndDeviceLorawanMac> edLorawanMac1 = edMac1->GetObject<ClassAEndDeviceLorawanMac> ();
      Ptr<LoraPhy> edPhy = node->GetDevice(0)->GetObject<LoraNetDevice> ()->GetPhy();
      edPhy->TraceConnectWithoutContext("StartSending", MakeBoundCallback(&energyNewPkt,bList));
      edLorawanMac1->SetMType (LorawanMacHeader::CONFIRMED_DATA_UP);
      edLorawanMac1->TraceConnectWithoutContext("RequiredTransmissions", MakeBoundCallback(&packetCounter, devID));
      edLorawanMac1->TraceConnectWithoutContext("RequiredTransmissions", MakeBoundCallback(&calcLatency, bList));
     //s edLorawanMac1->TraceConnectWithoutContext("SentNewPacket", MakeBoundCallback(&countNewPacket,bList));
     // edPhy->TraceConnectWithoutContext("StartSending", MakeCallback(&logPacket));
      if(h) 
      {
      edLorawanMac1->h = true;
    //  edLorawanMac1->TraceConnectWithoutContext("SentNewPacket", MakeBoundCallback(&getTs,bList));
    //  edLorawanMac1->TraceConnectWithoutContext("Dosend", MakeBoundCallback(&doCad,bList));
     // edLorawanMac1->TraceConnectWithoutContext("CalcBudget", MakeBoundCallback(&checkBudget,bList));
      edLorawanMac1->TraceConnectWithoutContext("RequiredTransmissions", MakeBoundCallback(&feedback, bList));
    //  edPhy->TraceConnectWithoutContext("SendingPacket", MakeCallback(&logCAD));
      }
      
    }

  /*********************
   *  Create Gateways  *
   *********************/

  // Create the gateway nodes (allocate them uniformely on the disc)
  NodeContainer gateways;
  gateways.Create (nGateways);





  Ptr<ListPositionAllocator> allocator = CreateObject<ListPositionAllocator> ();
  // Make it so that nodes are at a certain height > 0
  allocator->Add (Vector (0.0, 0.0, 15.0));
  mobility.SetPositionAllocator (allocator);
  mobility.Install (gateways);

  // Create a netdevice for each gateway
  phyHelper.SetDeviceType (LoraPhyHelper::GW);
  macHelper.SetDeviceType (LorawanMacHelper::GW);
  helper.Install (phyHelper, macHelper, gateways);

  /**********************
   *  Handle buildings  *
   **********************/

  double xLength = 130;
  double deltaX = 32;
  double yLength = 64;
  double deltaY = 17;
  int gridWidth = 2 * radius / (xLength + deltaX);
  int gridHeight = 2 * radius / (yLength + deltaY);
  if (realisticChannelModel == false)
    {
      gridWidth = 0;
      gridHeight = 0;
    }
  Ptr<GridBuildingAllocator> gridBuildingAllocator;
  gridBuildingAllocator = CreateObject<GridBuildingAllocator> ();
  gridBuildingAllocator->SetAttribute ("GridWidth", UintegerValue (gridWidth));
  gridBuildingAllocator->SetAttribute ("LengthX", DoubleValue (xLength));
  gridBuildingAllocator->SetAttribute ("LengthY", DoubleValue (yLength));
  gridBuildingAllocator->SetAttribute ("DeltaX", DoubleValue (deltaX));
  gridBuildingAllocator->SetAttribute ("DeltaY", DoubleValue (deltaY));
  gridBuildingAllocator->SetAttribute ("Height", DoubleValue (6));
  gridBuildingAllocator->SetBuildingAttribute ("NRoomsX", UintegerValue (2));
  gridBuildingAllocator->SetBuildingAttribute ("NRoomsY", UintegerValue (4));
  gridBuildingAllocator->SetBuildingAttribute ("NFloors", UintegerValue (2));
  gridBuildingAllocator->SetAttribute (
      "MinX", DoubleValue (-gridWidth * (xLength + deltaX) / 2 + deltaX / 2));
  gridBuildingAllocator->SetAttribute (
      "MinY", DoubleValue (-gridHeight * (yLength + deltaY) / 2 + deltaY / 2));
  BuildingContainer bContainer = gridBuildingAllocator->Create (gridWidth * gridHeight);

  BuildingsHelper::Install (endDevices);
  BuildingsHelper::Install (gateways);

  // Print the buildings
  if (print)
    {
      std::ofstream myfile;
      myfile.open ("buildings.txt");
      std::vector<Ptr<Building>>::const_iterator it;
      int j = 1;
      for (it = bContainer.Begin (); it != bContainer.End (); ++it, ++j)
        {
          Box boundaries = (*it)->GetBoundaries ();
          myfile << "set object " << j << " rect from " << boundaries.xMin << "," << boundaries.yMin
                 << " to " << boundaries.xMax << "," << boundaries.yMax << std::endl;
        }
      myfile.close ();
    }

  /**********************************************
   *  Set up the end device's spreading factor  *
   **********************************************/

  macHelper.SetSpreadingFactorsUp (endDevices, gateways, channel);

  NS_LOG_DEBUG ("Completed configuration");

  /*********************************************
   *  Install applications on the end devices  *
   *********************************************/

  Time appStopTime = Hours (simulationTime);
  PeriodicSenderHelper appHelper = PeriodicSenderHelper ();
  appHelper.SetPacketSize (30);
  Ptr<RandomVariableStream> rv = CreateObjectWithAttributes<UniformRandomVariable> ( "Min", DoubleValue (8), "Max", DoubleValue (15));
  ApplicationContainer appContainer;
  
 
  for (NodeContainer::Iterator j = endDevices.Begin (); j != endDevices.End (); ++j)
    { 
      /*provide random period to applications*/
      int randomValue = 30;//rv->GetInteger()*2; 
      Time random_time = Minutes(randomValue);
      appHelper.SetPeriod(random_time);
      appHelper.initDelayValue= 0; //chosenTs*tsLen;
      appContainer.Add(appHelper.Install(*j)); //attach application to nodes 

      /*initialize the bList array*/
      Ptr <Node> node = (*j);
      int id = node->GetId();
      bList[id].s = randomValue;
      bList[id].P = std::ceil((appStopTime.GetMinutes()/bList[id].s));
      bList[id].bCap = bList[id].P * GetTxEnergy(0)*2;
      bList[id].curSOC = bList[id].bCap;
      int tsMin =  int(tsLen/60);
      bList[id].T = int(bList[id].s/ tsMin);
      bList[id].Ts = bList[id].P*bList[id].T+1;
      bList[id].initialize();
      
     // 
      
    }
   //  getSOC(bList);           //update the state of charge from file 
     readFromJson(bList);
     
    
    for (int i=0;i<nDevices;i++)
   {
    std::cout << "Node " << i << " sampling period " << bList[i].s << " battery capacity "  << "P" << bList[i].P << " "<< bList[i].bCap << " SOC " <<  bList[i].curSOC <<std::endl;
    for (int t=0;t< bList[i].T ;t++)
    {
     bList[i].ProbabilityUpdate(t);
    }

   } 
   
   for (int p=0;p<bList[7].P; p++)
   {
     AggError.push_back(0);
     AggError2.push_back(0); 
   }
//   std::cout << AggError.size() << " LLL \n " ;*/
  if(h)  ///add gw traces
  {
  for (NodeContainer::Iterator j = gateways.Begin (); j != gateways.End (); ++j)
  { 
    Ptr<Node> gw = *j;
    Ptr<LoraPhy> gwPhy = gw->GetDevice(0)->GetObject<LoraNetDevice> ()->GetPhy();
   // gwPhy->TraceConnectWithoutContext("ReceivedPacket", MakeCallback(&removePacket));
  // gwPhy->TraceConnectWithoutContext("ReceivedPacket",MakeBoundCallback(&GWPacketReceptionCallback,bList));
  }
  }  
  
  for (ApplicationContainer::Iterator j = appContainer.Begin ();j != appContainer.End (); ++j)
    { 
      
      Ptr<Node> node = (*j)->GetNode();
      uint32_t devID = node->GetId();
      (*j)->TraceConnectWithoutContext ("GeneratedNewPacket", MakeBoundCallback(&countNewPacket,bList));
      if(h)(*j)->TraceConnectWithoutContext ("GeneratedNewPacket", MakeBoundCallback(&getTs,bList));
    }


for(int i=0;i< nDevices;i++)
{
//std::cout << bList[i].s << " " << bList[i].P << " " << bList[i].T << std::endl;
}
  
  appContainer.Start (Seconds (0));
  appContainer.Stop (appStopTime);

  /**************************
   *  Create Network Server  *
   ***************************/

  // Create the NS node
  NodeContainer networkServer;
  networkServer.Create (1);

  // Create a NS for the network
  nsHelper.SetEndDevices (endDevices);
  nsHelper.SetGateways (gateways);
  nsHelper.Install (networkServer);

  //Create a forwarder for each gateway
  forHelper.Install (gateways);

  ////////////////
  // Simulation //
  ////////////////
 //Simulator::Schedule(Hours(24), &updateLifeSpan);
//  Simulator::Schedule(Hours(6), &changeTheta, 0.8);  //cap battery charge to 80% 
//  Simulator::Schedule(Hours(16), &changeTheta, 1);  //charge battery to full while there is still sunlight
  for (int i=0;i<nDevices;i++) {
  Simulator::Schedule(Seconds(tsLen), &updateSOC, &bList[i]);
  //Simulator::Schedule(Hours(24), &updateLifeSpan, &bList[i]);
  } 
  Simulator::Stop (appStopTime);

  NS_LOG_INFO ("Running simulation...");
  Simulator::Run ();

  Simulator::Destroy ();

  ///////////////////////////
  // Print results to file //
  ///////////////////////////

  
  NS_LOG_INFO ("Computing performance metrics...");
//  printLastSOC(bList);
//  updateDay(bList);
  LoraPacketTracker &tracker = helper.GetPacketTracker ();
  std::cout << tracker.CountMacPacketsGlobally (Seconds (0), appStopTime + Hours (1)) << " " << tracker.CountMacPacketsGloballyCpsr (Seconds (0), appStopTime + Hours (1)) << std::endl;
  std::cout << overall_attempt << "  " << overall_success << "  " << overall_failed << " "<< total_energy <<" " << total_packet<< " " << double(overall_attempt)/double(total_packet) << "\n";
   

/*for (auto itPhy = c_packetTracker.begin ();
       itPhy != c_packetTracker.end ();
       ++itPhy)
 {

//std::cout<<(*itPhy).second.packet <<" " << int((*itPhy).second.sf) << " "<<  (*itPhy).second.frequencyMhz << std::endl;

 }*/

 for(int i =0;i<nDevices;i++)  
 {  
    std::cout<< " Node " << i << " Timeslot ";
    for(int j=0;j<bList[i].T ;j++)
    {
     std::cout << bList[i].SelectedSlots[j] << " " ;
    }
    std::cout << std::endl;
    
    std::cout << std::endl;
 }

   // std::cout<< " Error " << " Size " << AggError.size();
 /*   for(int k=0;k<AggError.size();k++)
    {
    std::cout<< AggError[k] << " ";  
    }*/
  
   printTSVariation(bList);
   printAvgError();
   printAvgLatency(bList);

   
   writeToJson(bList);

   printToFile(tracker, appStopTime);
   
   for(int i=0;i<nDevices;i++)
   {
   printTrace(i, bList); 
   
   printProb(i,bList);
   }
  std::cout<< "h\n";
   //getSummary(bList);
  return 0;
}
