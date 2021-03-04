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
#include <vector>
#include <cmath>


using namespace ns3;
using namespace lorawan;

NS_LOG_COMPONENT_DEFINE ("ComplexLorawanNetworkExample");

// Network settings
int nDevices = 200;  //number of nodes in the network
int nGateways = 1;
double radius = 5000;
double simulationTime = 600;
int maxPeriod = 20;
int minPeriod = 10;
int CHANNEL_NUM = 48;
int SF_NUM=6;
int algoNum=1;  // 0 for bestfit-decresing, 6 for worstfit-decreasing, 2 for firstfit-decreasing, 3 for randomfit, 4 BFMinU , 5 WFMinU , 1 FFMinU
int schedulable=1; //0 for not schedulable, 1 for schedulable
double U_total =0;  //total network utilization (after assignment)
double U_est =0;  //total estimated network utilization 
double U_target= 100;
int total_assigned=0;
double total_txEnergy=0;
double MAX_TIMESLOT = 1.810;
// Channel model
bool realisticChannelModel = false;

int appPeriodSeconds = 600;

double hyperPeriod=1;

// Output control
bool print = true;


//number of retransmissions and timeslots

int maxRtx=2;
int slotNum=2;

/*helper class representing a real-time LoRa node, used in this script only*/
class LoRaNode
{
public:
    int id;
    double u;
    double p;
    double wcet;
    double u_est;
    int maxDR;
    double distance;
    double minU;

};



/*Class for virtual channel*/
class Vchnl
{
public:
    int id;
    std:: vector <LoRaNode> assignedNodes;  //vector with assigned Nodes in this channel
    double c=1;                   //Channel Capacity
    double timeslot ;                //timeslot length
    int numRtx=maxRtx; // Retransmission number for this channel, set to maxRtx for now;
    int hyperperiod=1;
    double txEnergy=0;
};

  


/*comparision function to sort the Nodes based on non-increasing (highest one first) of maximum utilization*/
bool compBymaxU(const LoRaNode& n1,const LoRaNode& n2) {
    
    return n1.u_est > n2.u_est;
}

/*comparision function to sort the Nodes based on non-decreasing order (lowest one first) of periods (=relative deadline in this case)*/
bool compByP(const LoRaNode& n1,const LoRaNode& n2) {
    
    return n1.p < n2.p;
}

/*comparision function to sort the Nodes based on non-increasing order (highest one first) of distance*/
bool compByDistance(const LoRaNode& n1,const LoRaNode& n2) {
    
    return n1.distance > n2.distance;
}

/*comparision function to sort the Nodes based on non-increasing order (highest one first) of minimum utilization (higher datarate->lower SF, this gives a minimum SF threshold for a node which in turn gives the minimum utilization)*/
bool compByminU(const LoRaNode& n1,const LoRaNode& n2) {
    
    return n1.minU > n2.minU;
}


bool compBymaxDR(const LoRaNode& n1,const LoRaNode& n2) {
    
    
  if(n1.maxDR!= n2.maxDR) 
    return n1.maxDR < n2.maxDR;
  else 
    return n1.minU > n2.minU;
}

// Recursive function to return gcd of a and b  
long long gcd(long long int a, long long int b) 
{ 
  if (b == 0) 
    return a; 
  return gcd(b, a % b); 
} 
  
// Function to return LCM of two numbers  
long long lcm(int a, int b) 
{ 
    return (a / gcd(a, b)) * b; 
} 


/*function to print the channels and assigned nodes in each channel*/
void printChannels(Vchnl chnlList[])
{   
    for (int i=0;i<CHANNEL_NUM;i++)
    {
        Vchnl currentCh = chnlList[i];
        U_total += (1-chnlList[i].c);
        std::cout << " Channel " << i << " timeslot length "<< chnlList[i].timeslot<<  " Remaining Capacity " << chnlList[i].c << " Channel Utilization " << 1 - chnlList[i].c << std:: endl;
        
        for (uint32_t j =0;j<currentCh.assignedNodes.size();j++){
            
            LoRaNode currentNode= currentCh.assignedNodes.at(j);
            std::cout << currentNode.id << " ";
            currentCh.hyperperiod= lcm(currentCh.hyperperiod,currentNode.p*1000);
           // total_txEnergy += (0.029*currentNode.wcet*currentCh.timeslot*3.3*(simulationTime/currentNode.p));
        }
         

        for (uint32_t j =0;j<currentCh.assignedNodes.size();j++){
        LoRaNode currentNode= currentCh.assignedNodes.at(j);
        currentCh.txEnergy += (0.029*currentNode.wcet*currentCh.timeslot*3.3*(currentCh.hyperperiod/(currentNode.p*1000)));
        
 	}
        std:: cout << "txEnergy " <<currentCh.txEnergy << std:: endl;
        total_txEnergy += currentCh.txEnergy;
        
    }
    std::cout << " Total Network Utilization " << U_total << std:: endl;
  
}

std::vector<double> readFromFile(int nDevices,double U_target)
{

std::stringstream ss;
std::ifstream input;
std::vector<double> u_set;
double U=0;
ss<<"input/"<<nDevices<<"_"<<U_target<< ".txt";
std:: cout<<ss.str()<< "\n";
input.open(ss.str());
if(!input.is_open())
{
std::cout<< "ERROR\n";
}
std::string line;
while(input>>U)
{
//std::cout<<U<< std::endl;
u_set.push_back(U);
}
return u_set;
}


/*function to print output to the file*/

void printToFile(Vchnl chnlList[], LoRaNode nodes[])
{

/*create the output file*/

   std::ofstream logfile;
   logfile.open ("sim_results_0304.csv",std::ios_base::app);

   logfile << algoNum << "," << nDevices << "," << radius << "," << CHANNEL_NUM << "," << U_total <<"," <<U_est << "," << schedulable << ","<< double(total_assigned)/double(nDevices) << ","  << total_txEnergy <<"\n";

 
}

double DBF_func(LoRaNode node_i, double t, double act_wcet)

{
if(t<node_i.p) return 0;
else 
{
double act_u = act_wcet/node_i.p;
return  act_wcet + act_u*(t-node_i.p);
}

}


/*function to check the DBF condition in PARTITION algorithm by baruah*/

bool Checkcondition(std::vector<LoRaNode> assignedNodes, LoRaNode newNode, double act_wcet )
{

double d_i = newNode.p; 
double dbf=0;
for (int j=0;j<assignedNodes.size();j++)
{

dbf += DBF_func(assignedNodes.at(j), d_i, act_wcet);

}
return ((d_i - dbf)>= act_wcet); 
}

/*function to print the nodes in the network*/

void printNodes(LoRaNode nodes[])
{
    for (int i=0;i<nDevices;i++){
        U_est += nodes[i].u_est;
        std::cout <<"Node " << nodes[i].id << " period " << nodes[i].p << " wcet "<< nodes[i].wcet << " Estimated Utilization " << nodes[i].u_est << " Actual Utilization " << nodes[i].u <<" distance from gateway " << nodes[i].distance <<" Maximum Datarate " << nodes[i].maxDR  <<" Minimum Utilization " << nodes[i].minU << std:: endl;
    }
 
   std::cout << " Total estimated Network Utilization " << U_est << " Hyper Period " << hyperPeriod << std:: endl;
    
}

/*UUnifast algorithm by Enrico Bini*/
std::vector <double> UUnifast(int nDevices, double U_target)
{
 std::vector <double> u_set;
 
 double sum=U_target;
 double nextSum=0;
 bool done = false;
  int iter_num=0;
 Ptr <RandomVariableStream> rv = CreateObjectWithAttributes<UniformRandomVariable> ("Min", DoubleValue (0), "Max", DoubleValue (1));
 while(!done && iter_num<1000){
    u_set.clear();
     double sum=U_target;
     double nextSum=0;
    for (int i=0; i< nDevices-1 ;i++) {

     double exp = 1/(double(nDevices)-double(i));
     nextSum = sum* pow(rv->GetValue(), exp);
     u_set.push_back(sum-nextSum);
     sum=nextSum;
     } 
     u_set.push_back(sum);
     done=true;
    for (int i=0; i< nDevices ;i++)
    {
     if(u_set.at(i)>1) 
     {
     done = false;
     iter_num++;
     break;
     }
    }
 
 }
 if(done) return u_set;
 else 
 { 
 std::cout<< "Task Generation Failed"<<std::endl;
 return u_set;
 }
}

/*Randomly Assigns Nodes to the channels */

void RandomFit(LoRaNode nodes[], int size, Vchnl chnlList[]){
   
    std::sort(nodes,nodes+size,compByminU);  //sort the nodes in decreasing order of min utilization
    Ptr <RandomVariableStream> rv_ch = CreateObjectWithAttributes<UniformRandomVariable> ("Min", DoubleValue (0), "Max", DoubleValue (CHANNEL_NUM-1));
    for(int i=0;i<size;i++){
       LoRaNode currentNode = nodes[i];
       int j = rv_ch->GetInteger();   
       int vChDR = j%SF_NUM;                   
       //calculate actual execution time and utilization based on timeslot lenght of virtual channel
       double act_time = (chnlList[j].timeslot/MAX_TIMESLOT)*currentNode.wcet;
       double actU = act_time/currentNode.p;
          
       if((chnlList[j].c >= actU) && (vChDR<=currentNode.maxDR)){
        currentNode.u=actU;
        nodes[i].u=actU;
        chnlList[j].assignedNodes.push_back(currentNode);
        chnlList[j].c -=actU;
        total_assigned+=1;
        }
        else    //looped through all channels, couldnt find a channel to match the utilization of some node
        {
        std::cout<<" Not feasible to assign these nodes to " << CHANNEL_NUM << " channels "<< std:: endl;
        schedulable=0;
        }
    }
}



/*Assigns Nodes to the channels based on first fit decreasing algorithm*/

void FirstFit(LoRaNode nodes[], int size, Vchnl chnlList[]){
   
    std::sort(nodes,nodes+size,compBymaxU);  //sort the nodes in decreasing order of maximum utilization
    
    for(int i=0;i<size;i++){
        LoRaNode currentNode = nodes[i];
        int j;                      //need this index later thats why declared here instead of inline declaration
        int vChDR;
        for( j=0;j<CHANNEL_NUM;j++)
        {
            //calculate actual execution time and utilization based on timeslot lenght of virtual channel
            vChDR = j%SF_NUM;
            double act_time = (chnlList[j].timeslot/MAX_TIMESLOT)*currentNode.wcet;
            double actU = act_time/currentNode.p;
           
            if((chnlList[j].c >= actU) && (vChDR<=currentNode.maxDR)){
            //std::cout<< chnlList[j].id << " has been assigned Node "<<currentNode.id << endl;
                currentNode.u=actU;
                nodes[i].u = actU;
                chnlList[j].assignedNodes.push_back(currentNode);
                chnlList[j].c -=actU;
                total_assigned+=1;
             //   std::cout<< chnlList[j].id << " has capacity "<<chnlList[j].c << endl;
                break;
            }
        }
        
        if(j==CHANNEL_NUM)     //looped through all channels, couldnt find a channel to match the utilization of some node
        {
        
        std::cout<<" Not feasible to assign these nodes to " << CHANNEL_NUM << " channels "<< std:: endl;
        schedulable=0;
        }
    }
}



/*Assigns Nodes to the channels based on first fit min util decreasing algorithm*/

void FirstFitMinUtil(LoRaNode nodes[], int size, Vchnl chnlList[]){
   
    std::sort(nodes,nodes+size,compBymaxDR);  //sort the nodes in decreasing order of maximum utilization
    
    for(int i=0;i<size;i++){
        LoRaNode currentNode = nodes[i];
        int j;                      //need this index later thats why declared here instead of inline declaration
        int vChDR;
        for( j=0;j<CHANNEL_NUM;j++)
        {
            //calculate actual execution time and utilization based on timeslot lenght of virtual channel
            vChDR = j%SF_NUM;
            double act_time = (chnlList[j].timeslot/MAX_TIMESLOT)*currentNode.wcet;
            double actU = act_time/currentNode.p;
           
            if((chnlList[j].c >= actU) && (vChDR<=currentNode.maxDR)){
            //std::cout<< chnlList[j].id << " has been assigned Node "<<currentNode.id << endl;
                currentNode.u=actU;
                nodes[i].u = actU;
                chnlList[j].assignedNodes.push_back(currentNode);
                chnlList[j].c -=actU;
                total_assigned+=1;
             //   std::cout<< chnlList[j].id << " has capacity "<<chnlList[j].c << endl;
                break;
            }
        }
        
        if(j==CHANNEL_NUM)     //looped through all channels, couldnt find a channel to match the utilization of some node
        {
        
        std::cout<<" Not feasible to assign these nodes to " << CHANNEL_NUM << " channels "<< std:: endl;
        schedulable=0;
        }
    }
}



/*Assigns Nodes to the channels based on PARTITION algorithm by Baruah et.al.*/

void Partition(LoRaNode nodes[], int size, Vchnl chnlList[]){
   
    std::sort(nodes,nodes+size,compByP);  //sort the nodes in increasing order of relative deadline
    
    for(int i=0;i<size;i++){
        LoRaNode currentNode = nodes[i];
        int j;                      //need this index later thats why declared here instead of inline declaration
        int vChDR;
        for( j=0;j<CHANNEL_NUM;j++)
        {
            //calculate actual execution time and utilization based on timeslot lenght of virtual channel
            vChDR = j%SF_NUM;
            double act_time = (chnlList[j].timeslot/MAX_TIMESLOT)*currentNode.wcet;
            double actU = act_time/currentNode.p;
           
            if((chnlList[j].c >= actU) && (Checkcondition(chnlList[j].assignedNodes, currentNode, act_time)) && (vChDR<=currentNode.maxDR) ){
            //std::cout<< chnlList[j].id << " has been assigned Node "<<currentNode.id << endl;
                currentNode.u=actU;
                nodes[i].u = actU;
                chnlList[j].assignedNodes.push_back(currentNode);
                chnlList[j].c -=actU;
                total_assigned+=1;
             //   std::cout<< chnlList[j].id << " has capacity "<<chnlList[j].c << endl;
                break;
            }
        }
        
        if(j==CHANNEL_NUM)     //looped through all channels, couldnt find a channel to match the utilization of some node
        {
        
        std::cout<<" Not feasible to assign these nodes to " << CHANNEL_NUM << " channels "<< std:: endl;
        schedulable=0;
        }
    }
}





/*Assigns Nodes to the channels based on best fit decreasing algorithm*/

void BestFit(LoRaNode nodes[], int size, Vchnl chnlList[]){
    
    std::sort(nodes,nodes+size,compBymaxU);  //sort the nodes in decreasing order of maximum util (highest one first)
   
    for(int i=0;i<size;i++){
        LoRaNode currentNode = nodes[i];
        int j;
        double min = 1000;
        int best=0;    //store the best index
        double bestU=0;
        double act_time, actU=0;
        int vChDR=0;
        for( j=0;j<CHANNEL_NUM;j++)
        {
            //calculate actual execution time and utilization based on timeslot lenght of virtual channel
            act_time = (chnlList[j].timeslot/MAX_TIMESLOT)*currentNode.wcet;
            //std::cout<< " actual time for node : "<< i << " " << chnlList[j].timeslot << "  " <<  MAX_TIMESLOT<< " in channel " << j  << std:: endl;
            actU = act_time/currentNode.p;
            vChDR = j%SF_NUM;
        //    std::cout<< "Act u for channel " << j << ":" << actU <<endl;
            if((chnlList[j].c >= actU) && ((chnlList[j].c - actU)< min ) && (vChDR<=currentNode.maxDR)){
                
                best = j;
                bestU = actU;
                min = chnlList[j].c - actU;
                
            }
            
        }
        
        
        if(min==1000)
        {
        std::cout<<" Not feasible to assign this node to any of the " << CHANNEL_NUM << " channels "<< std:: endl;
        schedulable=0;
        std::cout<< "Node " << currentNode.id << " utization " << currentNode.u_est << " min SF " << currentNode.maxDR << std:: endl;
         }
        else {
          //  std::cout<< best << " has been assigned Node "<<currentNode.id << " utilization " << bestU << std:: endl;
            currentNode.u=bestU;
            nodes[i].u = bestU;
            chnlList[best].assignedNodes.push_back(currentNode);
            chnlList[best].c -=bestU;
          //  std::cout<< best << " has capacity "<<chnlList[best].c << endl;
            
            total_assigned+=1;
        }
    }
}

/*Assigns Nodes to the channels based on best fit min util decreasing algorithm*/

void BestFitMinUtil(LoRaNode nodes[], int size, Vchnl chnlList[]){
    
    std::sort(nodes,nodes+size,compBymaxDR);  //sort the nodes in decreasing order of min util (highest one first)
   
    for(int i=0;i<size;i++){
        LoRaNode currentNode = nodes[i];
        int j;
        double min = 1000;
        int best=0;    //store the best index
        double bestU=0;
        double act_time, actU=0;
        int vChDR=0;
        for( j=0;j<CHANNEL_NUM;j++)
        {
            //calculate actual execution time and utilization based on timeslot lenght of virtual channel
            act_time = (chnlList[j].timeslot/MAX_TIMESLOT)*currentNode.wcet;
            //std::cout<< " actual time for node : "<< i << " " << chnlList[j].timeslot << "  " <<  MAX_TIMESLOT<< " in channel " << j  << std:: endl;
            actU = act_time/currentNode.p;
            vChDR = j%SF_NUM;
        //    std::cout<< "Act u for channel " << j << ":" << actU <<endl;
            if((chnlList[j].c >= actU) && ((chnlList[j].c - actU)< min ) && (vChDR<=currentNode.maxDR)){
                
                best = j;
                bestU = actU;
                min = chnlList[j].c - actU;
                
            }
            
        }
        
        
        if(min==1000)
        {
        std::cout<<" Not feasible to assign this node to any of the " << CHANNEL_NUM << " channels "<< std:: endl;
        schedulable=0;
        std::cout<< "Node " << currentNode.id << " utization " << currentNode.u_est << " min SF " << currentNode.maxDR << std:: endl;
         }
        else {
          //  std::cout<< best << " has been assigned Node "<<currentNode.id << " utilization " << bestU << std:: endl;
            currentNode.u=bestU;
            nodes[i].u = bestU;
            chnlList[best].assignedNodes.push_back(currentNode);
            chnlList[best].c -=bestU;
          //  std::cout<< best << " has capacity "<<chnlList[best].c << endl;
            
            total_assigned+=1;
        }
    }
}




/*Assigns Nodes to the channels based on worst fit decreasing algorithm*/

void WorstFit(LoRaNode nodes[], int size , Vchnl chnlList[]){
   
    std::sort(nodes,nodes+size,compBymaxU);
    for(int i=0;i<size;i++){
        LoRaNode currentNode = nodes[i];
        int j;
        double max = -999;
        int best=0;
        double bestU=0;
        double act_time=0, actU=0;
        int vChDR=0;
        for( j=0;j<CHANNEL_NUM;j++)
        {   act_time = (chnlList[j].timeslot/MAX_TIMESLOT)*currentNode.wcet;
            actU = act_time/currentNode.p;
            vChDR = j%SF_NUM;
            if((chnlList[j].c >= actU) && ((chnlList[j].c - actU)> max ) && (vChDR<=currentNode.maxDR)){
                
                best = j;
                bestU=actU;
                max = chnlList[j].c - actU;
               
            }
            
        }
        
        
        if(max==-999)
        {
        std::cout<<" Not feasible to assign these nodes to " << CHANNEL_NUM << " channels "<< std:: endl;
        schedulable=0;
        }
        else {
            currentNode.u = bestU;
            nodes[i].u= bestU;
            chnlList[best].assignedNodes.push_back(currentNode);
            chnlList[best].c -=bestU;
             total_assigned+=1;

        }
    }
}


/*Assigns Nodes to the channels based on worst fit decreasing min Util algorithm*/

void WorstFitMinUtil(LoRaNode nodes[], int size , Vchnl chnlList[]){
   
    std::sort(nodes,nodes+size,compBymaxDR);    //sort the nodes in decreasing order of min util (highest one first)
    for(int i=0;i<size;i++){
        LoRaNode currentNode = nodes[i];
        int j;
        double max = -999;
        int best=0;
        double bestU=0;
        double act_time=0, actU=0;
        int vChDR=0;
        for( j=0;j<CHANNEL_NUM;j++)
        {   act_time = (chnlList[j].timeslot/MAX_TIMESLOT)*currentNode.wcet;
            actU = act_time/currentNode.p;
            vChDR = j%SF_NUM;
            if((chnlList[j].c >= actU) && ((chnlList[j].c - actU)> max ) && (vChDR<=currentNode.maxDR)){
                
                best = j;
                bestU=actU;
                max = chnlList[j].c - actU;
               
            }
            
        }
        
        
        if(max==-999)
        {
        std::cout<<" Not feasible to assign these nodes to " << CHANNEL_NUM << " channels "<< std:: endl;
        schedulable=0;
        }
        else {
            currentNode.u = bestU;
            nodes[i].u= bestU;
            chnlList[best].assignedNodes.push_back(currentNode);
            chnlList[best].c -=bestU;
             total_assigned+=1;

        }
    }
}

/*calculates the minimum Spreading Factor/maximum DR required for each node according to its relative position from the gateway and saves it in the <LoRaNode> array*/

void SetmaxDR(NodeContainer endDevices, NodeContainer gateways, Ptr<LoraChannel> channel , LoRaNode nodeList[])
{

for (NodeContainer::Iterator j = endDevices.Begin (); j != endDevices.End (); ++j)
    {
      Ptr<Node> object = *j;
      int nodeId = object->GetId();
      Ptr<MobilityModel> position = object->GetObject<MobilityModel> ();
      NS_ASSERT (position != 0);
      Ptr<NetDevice> netDevice = object->GetDevice (0);
      Ptr<LoraNetDevice> loraNetDevice = netDevice->GetObject<LoraNetDevice> ();
      NS_ASSERT (loraNetDevice != 0);
      Ptr<ClassAEndDeviceLorawanMac> mac =
          loraNetDevice->GetMac ()->GetObject<ClassAEndDeviceLorawanMac> ();
      NS_ASSERT (mac != 0);

      // Try computing the distance from each gateway and find the best one
      Ptr<Node> bestGateway = gateways.Get (0);
      Ptr<MobilityModel> bestGatewayPosition = bestGateway->GetObject<MobilityModel> ();
       double distance = position->GetDistanceFrom(bestGatewayPosition);
        //  std::cout << " distance " << distance << "\n";
          nodeList[nodeId].distance = distance;
      // Assume devices transmit at 14 dBm
      double highestRxPower = channel->GetRxPower (14, position, bestGatewayPosition);

      for (NodeContainer::Iterator currentGw = gateways.Begin () + 1; currentGw != gateways.End ();
           ++currentGw)
        {
          // Compute the power received from the current gateway
          Ptr<Node> curr = *currentGw;
          Ptr<MobilityModel> currPosition = curr->GetObject<MobilityModel> ();
          double currentRxPower = channel->GetRxPower (14, position, currPosition); // dBm
         
          if (currentRxPower > highestRxPower)
            {
              bestGateway = curr;
              bestGatewayPosition = curr->GetObject<MobilityModel> ();
              highestRxPower = currentRxPower;
            }
        }

      // NS_LOG_DEBUG ("Rx Power: " << highestRxPower);
      double rxPower = highestRxPower;

      // Get the ED sensitivity
      Ptr<EndDeviceLoraPhy> edPhy = loraNetDevice->GetPhy ()->GetObject<EndDeviceLoraPhy> ();
      const double *edSensitivity = edPhy->sensitivity;

      if (rxPower > *edSensitivity)
        {
          nodeList[nodeId].maxDR = 5;
        }
      else if (rxPower > *(edSensitivity + 1))
        {
          nodeList[nodeId].maxDR=4;
        }
      else if (rxPower > *(edSensitivity + 2))
        {
          nodeList[nodeId].maxDR =3;
        }
      else if (rxPower > *(edSensitivity + 3))
        {
          nodeList[nodeId].maxDR =2;
        }
      else if (rxPower > *(edSensitivity + 4))
        {
          nodeList[nodeId].maxDR =1;
        }
      else if (rxPower > *(edSensitivity + 5))
        {
          nodeList[nodeId].maxDR =0;
        }
      else // Device is out of range. Assign SF12.
        {
          // NS_LOG_DEBUG ("Device out of range");
          nodeList[nodeId].maxDR=0;
          // NS_LOG_DEBUG ("sfQuantity[6] = " << sfQuantity[6]);
        }
    }
}


void SetmaxDRGivenDistribution (NodeContainer endDevices,
                                                        NodeContainer gateways,
                                                        std::vector<double> distribution,LoRaNode nodeList[])
{
  NS_LOG_FUNCTION_NOARGS ();

  std::vector<int> sfQuantity (7, 0);
  Ptr<UniformRandomVariable> uniformRV = CreateObject<UniformRandomVariable> ();
  std::vector<double> cumdistr (6);
  cumdistr[0] = distribution[0];
  for (int i = 1; i < 7; ++i)
    {
      cumdistr[i] = distribution[i] + cumdistr[i - 1];
    }

  NS_LOG_DEBUG ("Distribution: " << distribution[0] << " " << distribution[1] << " "
                                 << distribution[2] << " " << distribution[3] << " "
                                 << distribution[4] << " " << distribution[5]);
  NS_LOG_DEBUG ("Cumulative distribution: " << cumdistr[0] << " " << cumdistr[1] << " "
                                            << cumdistr[2] << " " << cumdistr[3] << " "
                                            << cumdistr[4] << " " << cumdistr[5]);

  for (NodeContainer::Iterator j = endDevices.Begin (); j != endDevices.End (); ++j)
    {
      Ptr<Node> object = *j;
      int nodeId=object->GetId();
   

      double prob = uniformRV->GetValue (0, 1);

      // NS_LOG_DEBUG ("Probability: " << prob);
      if (prob < cumdistr[0])
        {
          nodeList[nodeId].maxDR=5;
        }
      else if (prob > cumdistr[0] && prob < cumdistr[1])
        {
          nodeList[nodeId].maxDR=4;
        }
      else if (prob > cumdistr[1] && prob < cumdistr[2])
        {
          nodeList[nodeId].maxDR =3;
        }
      else if (prob > cumdistr[2] && prob < cumdistr[3])
        {
          nodeList[nodeId].maxDR=2;
        }
      else if (prob > cumdistr[3] && prob < cumdistr[4])
        {
          nodeList[nodeId].maxDR=1;
        }
      else
        {
          nodeList[nodeId].maxDR=0;
        }

    } // end loop on nodes



} //  end function


/*calculates the minimum utilization for each node using the maximum DR that it can use*/
void SetMinU(LoRaNode nodeList[],  Vchnl chnlList[] ,int size )

{
for(int i=0;i<size;i++)
{
int maxDR = nodeList[i].maxDR;
double min_time = (chnlList[maxDR].timeslot/MAX_TIMESLOT)*nodeList[i].wcet; 
nodeList[i].minU =  min_time/nodeList[i].p;
}

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
  cmd.AddValue ("Max","Max period ",maxPeriod);
  cmd.AddValue ("Min","Min period",minPeriod);
  cmd.AddValue ("numCh","Number of channels ",CHANNEL_NUM);
  cmd.AddValue ("algo",  "algo to use", algoNum);
  cmd.AddValue("U", "total utilization before assignment" , U_target);

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

  /***********
   *  Setup  *
   ***********/

/*Real-time LoRa arrays*/

   LoRaNode nodeList[nDevices];
   Vchnl chnlList[CHANNEL_NUM]; 

 
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

  // Connect trace sources
  for (NodeContainer::Iterator j = endDevices.Begin (); j != endDevices.End (); ++j)
    {
      Ptr<Node> node = *j;
      Ptr<LoraNetDevice> loraNetDevice = node->GetDevice (0)->GetObject<LoraNetDevice> ();
      Ptr<LoraPhy> phy = loraNetDevice->GetPhy ();
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


  /**********************************************
   *  Set up the end device's spreading factor  *
   **********************************************/

  macHelper.SetSpreadingFactorsUp (endDevices, gateways, channel);

  NS_LOG_DEBUG ("Completed configuration");

  /*********************************************
   *  Install applications on the end devices  *
   *********************************************/

  Time appStopTime = Seconds (simulationTime);
  PeriodicSenderHelper appHelper = PeriodicSenderHelper ();
  Ptr <RandomVariableStream> rv = CreateObjectWithAttributes<UniformRandomVariable> ("Min", DoubleValue (minPeriod), "Max", DoubleValue (maxPeriod));
  Ptr <RandomVariableStream> rv_harmonic = CreateObjectWithAttributes<UniformRandomVariable> ("Min", DoubleValue (13), "Max", DoubleValue (20)); //random exponent
  rv->SetStream(1);
  rv_harmonic->SetStream(1);
  ApplicationContainer appContainer;
  
  /*initialize LoRaNodes*/
 // std::vector<double> u_set= UUnifast(nDevices,U_target);
  std::vector<double> u_set= readFromFile(nDevices,U_target);
  for (NodeContainer::Iterator j = endDevices.Begin (); j != endDevices.End (); ++j)
    { 

      Ptr<Node> node = *j;
      int random_exponent = rv_harmonic->GetInteger();
     // std::cout<< "random exponent " << random_exponent <<std::endl;
      LoRaNode temp;
      temp.id = node->GetId();
      temp.p = (pow(2,random_exponent))/1000; //s
      temp.wcet = MAX_TIMESLOT*slotNum*maxRtx;
   //   temp.u_est= u_set.at(temp.id);//(double) temp.wcet/temp.p; 
      temp.u_est = temp.wcet/temp.p;
      temp.p = temp.wcet/temp.u_est;
      nodeList[temp.id]=temp;
      
      Time random_time = Seconds(temp.p);//Seconds(random_period);
      appHelper.SetPeriod(random_time);
      appContainer.Add(appHelper.Install(*j)); 

    }
//   printNodes(nodeList);
  appHelper.SetPacketSize (23);
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


  /**************************
   *  Real-time LoRa setup  *
   ***************************/

  

    
    //initialize channels
    for(int i=0;i<CHANNEL_NUM;i++)
    {
        chnlList[i].id=i;
       // std::cout << " SF MOD "<<i%SF_NUM << "\n";
        switch ((i%SF_NUM)) {
            case 0:
                chnlList[i].timeslot = 1.810;
                break;
            case 1:
                chnlList[i].timeslot=0.987;
                break;
            case 2:
                chnlList[i].timeslot = 0.493;
                break;
            case 3:
                chnlList[i].timeslot = 0.267;
                break;
            case 4:
                chnlList[i].timeslot = 0.154;
                break;
            case 5:
                chnlList[i].timeslot = 0.082;
                break;
            default:
                break;
        }
        
    }
    
    
 
   SetmaxDR(endDevices, gateways, channel, nodeList); 
   //std::vector<double> distribution = {0.1,0.1,0.1,0.1,0.1,0.5};
   //SetmaxDRGivenDistribution(endDevices, gateways, distribution, nodeList);
   SetMinU(nodeList, chnlList ,nDevices);
   switch(algoNum)
   {
     case 0: 
          BestFit(nodeList,nDevices,chnlList);
          break;
     case 1:
          FirstFitMinUtil(nodeList,nDevices,chnlList);
          break;
     case 2:
          FirstFit(nodeList,nDevices,chnlList);
          break;
     case 3:
          RandomFit(nodeList,nDevices,chnlList); 
          break;
     case 4:
          BestFitMinUtil(nodeList,nDevices,chnlList);
          break;
     case 5:
          WorstFitMinUtil(nodeList,nDevices,chnlList);
          break;
     case 6:
          WorstFit(nodeList,nDevices,chnlList);
          break;
     default:
          break;  
   }
    printNodes(nodeList);
    printChannels(chnlList);

  
  ////////////////
  // Simulation //
  ////////////////

  Simulator::Stop (appStopTime + Hours (1));

  NS_LOG_INFO ("Running simulation...");
  Simulator::Run ();

  Simulator::Destroy ();

  ///////////////////////////
  // Print results to file //
  ///////////////////////////
  
  NS_LOG_INFO ("Computing performance metrics...");

  LoraPacketTracker &tracker = helper.GetPacketTracker ();
  std::cout << tracker.CountMacPacketsGlobally (Seconds (0), appStopTime + Hours (1)) << std::endl;
  printToFile(chnlList,nodeList);

  
  return 0;
}
