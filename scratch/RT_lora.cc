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
#define CHANNEL_NUM 8
#define NODE_NUM 100

using namespace ns3;
using namespace lorawan;

NS_LOG_COMPONENT_DEFINE ("ComplexLorawanNetworkExample");

// Network settings
int nDevices = 200;
int nGateways = 1;
double radius = 7500;
double simulationTime = 600;
int maxPeriod = 50;
int minPeriod = 25;
// Channel model
bool realisticChannelModel = false;

int appPeriodSeconds = 600;

// Output control
bool print = true;



int numRtx=2;
class LoRaNode
{
public:
    int id;
    double u;
    double p;
    double wcet;
    double u_est;
    int minSF;
};



/*Class for virtual channel*/
class Vchnl
{
public:
    int id;
    std:: vector <LoRaNode> assignedNodes;  //vector with assigned Nodes in this channel
    double c=1;                   //Channel Capacity
    int timeslot ;                //timeslot length
};

Vchnl chnlList[CHANNEL_NUM];     //Global array for the channels in the network (may change to local later)


/*comparision function to sort the Nodes based on utilization*/
bool mycomp(const LoRaNode& n1,const LoRaNode& n2) {
    
    return n1.u_est > n2.u_est;
}

/*function to print the channels and assigned nodes in each channel*/
void printChannels()
{   double U_total =0;
    for (int i=0;i<CHANNEL_NUM;i++)
    {
        Vchnl currentCh = chnlList[i];
        U_total += (1-chnlList[i].c);
        std::cout << " Channel " << i << " Remaining Capacity " << chnlList[i].c << " Channel Utilization " << 1 - chnlList[i].c << std:: endl;
        
        for (uint32_t j =0;j<currentCh.assignedNodes.size();j++){
            std::cout << currentCh.assignedNodes.at(j).id << " ";
        }
        std:: cout<< std:: endl;
    
        
    }
    std::cout << " Total Network Utilization " << U_total << std:: endl;
}

/*function to print the nodes in the network*/

void printNodes(LoRaNode nodes[])
{
    for (int i=0;i<NODE_NUM;i++){
        
        std::cout <<"Node " << nodes[i].id << " period " << nodes[i].p << " wcet "<< nodes[i].wcet << " Estimated Utilization " << nodes[i].u_est << " Actual Utilization " << nodes[i].u << std:: endl;
    }
    
}


/*Assigns Nodes to the channels based on first fit decreasing algorithm*/

void FirstFit(LoRaNode nodes[], int size){
   
    std::sort(nodes,nodes+size,mycomp);  //sort the nodes in decreasing order of utilization
    
    for(int i=0;i<size;i++){
        LoRaNode currentNode = nodes[i];
        int j;                      //need this index later thats why declared here instead of inline declaration
        for( j=0;j<CHANNEL_NUM;j++)
        {
            //calculate actual execution time and utilization based on timeslot lenght of virtual channel
            double act_time = (chnlList[j].timeslot*currentNode.wcet)/1000;
            double actU = act_time/currentNode.p;
            
            if((chnlList[j].c >= actU) && (j>=currentNode.minSF)){
            //std::cout<< chnlList[j].id << " has been assigned Node "<<currentNode.id << endl;
                currentNode.u=actU;
                chnlList[j].assignedNodes.push_back(currentNode);
                chnlList[j].c -=actU;
             //   std::cout<< chnlList[j].id << " has capacity "<<chnlList[j].c << endl;
                break;
            }
        }
        
        if(j==CHANNEL_NUM)     //looped through all channels, couldnt find a channel to match the utilization of some node
        {
        std::cout<<" Not feasible to assign these nodes to " << CHANNEL_NUM << " channels "<< std:: endl;
        }
    }
}


/*Assigns Nodes to the channels based on best fit decreasing algorithm*/

void BestFit(LoRaNode nodes[], int size){
    
    std::sort(nodes,nodes+size,mycomp);  //sort the nodes in decreasing order of utilization
   
    for(int i=0;i<size;i++){
        LoRaNode currentNode = nodes[i];
        int j;
        double min = 1000;
        int best=0;    //store the best index
        double bestU=0;
        double act_time, actU=0;
        for( j=0;j<CHANNEL_NUM;j++)
        {
            //calculate actual execution time and utilization based on timeslot lenght of virtual channel
            act_time = (chnlList[j].timeslot*currentNode.wcet)/1000;
            actU = act_time/currentNode.p;
        //    std::cout<< "Act u for channel " << j << ":" << actU <<endl;
            if((chnlList[j].c >= actU) && ((chnlList[j].c - actU)< min ) && (j>=currentNode.minSF)){
                
                best = j;
                bestU = actU;
                min = chnlList[j].c - actU;
            }
            
        }
        
        
        if(min==1000)
        {
        std::cout<<" Not feasible to assign these nodes to " << CHANNEL_NUM << " channels "<< std:: endl;
        }
        else {
          //  std::cout<< best << " has been assigned Node "<<currentNode.id << endl;
            currentNode.u=bestU;
            chnlList[best].assignedNodes.push_back(currentNode);
            chnlList[best].c -=bestU;
          //  std::cout<< best << " has capacity "<<chnlList[best].c << endl;
        }
    }
}


/*Assigns Nodes to the channels based on worst fit decreasing algorithm*/

void WorstFit(LoRaNode nodes[], int size){
   
    std::sort(nodes,nodes+size,mycomp);
    for(int i=0;i<size;i++){
        LoRaNode currentNode = nodes[i];
        int j;
        double max = -999;
        int best=0;
        double bestU=0;
        double act_time=0, actU=0;
        for( j=0;j<CHANNEL_NUM;j++)
        {   act_time = (chnlList[j].timeslot*currentNode.wcet)/1000;
            actU = act_time/currentNode.p;
            if((chnlList[j].c >= actU) && ((chnlList[j].c - actU)> max ) && (j>=currentNode.minSF)){
                
                best = j;
                bestU=actU;
                max = chnlList[j].c - actU;
            }
            
        }
        
        
        if(max==-999)
        {
        std::cout<<" Not feasible to assign these nodes to " << CHANNEL_NUM << " channels "<< std:: endl;
        }
        else {
            currentNode.u = bestU;
            chnlList[best].assignedNodes.push_back(currentNode);
            chnlList[best].c -=bestU;

        }
    }
}

/*calculate the minimum Spreading Factor required for each node according to its relative position from the gateway and saves it in the <LoRaNode> array*/

void SetMinSF(NodeContainer endDevices, NodeContainer gateways, Ptr<LoraChannel> channel , LoRaNode nodeList[])
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
          nodeList[nodeId].minSF = 5;
        }
      else if (rxPower > *(edSensitivity + 1))
        {
          nodeList[nodeId].minSF=4;
        }
      else if (rxPower > *(edSensitivity + 2))
        {
          nodeList[nodeId].minSF =3;
        }
      else if (rxPower > *(edSensitivity + 3))
        {
          nodeList[nodeId].minSF =2;
        }
      else if (rxPower > *(edSensitivity + 4))
        {
          nodeList[nodeId].minSF =1;
        }
      else if (rxPower > *(edSensitivity + 5))
        {
          nodeList[nodeId].minSF =0;
        }
      else // Device is out of range. Assign SF12.
        {
          // NS_LOG_DEBUG ("Device out of range");
          nodeList[nodeId].minSF=0;
          // NS_LOG_DEBUG ("sfQuantity[6] = " << sfQuantity[6]);
        }
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
   LoRaNode nodeList[nDevices];
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
  rv->SetStream(1);
  ApplicationContainer appContainer;
  
  
  for (NodeContainer::Iterator j = endDevices.Begin (); j != endDevices.End (); ++j)
    { 

      Ptr<Node> node = *j;
      int random_period = rv->GetInteger();
      //std::cout << "node " << node->GetId() << " period " << random_period << std:: endl; 
      LoRaNode temp;
      temp.id = node->GetId();
      temp.p = random_period;
      temp.wcet = 1*numRtx;
      temp.u_est=  (double) temp.wcet/temp.p; //((double) rand() / (RAND_MAX));
      nodeList[temp.id]=temp;
      
      Time random_time = Seconds(random_period);
      appHelper.SetPeriod(random_time);
      appContainer.Add(appHelper.Install(*j)); 

    }
   printNodes(nodeList);
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
        switch (i%6) {
            case 0:
                chnlList[i].timeslot = 1810;
                break;
            case 1:
                chnlList[i].timeslot=987;
            case 2:
                chnlList[i].timeslot = 493;
            case 3:
                chnlList[i].timeslot = 267;
            case 4:
                chnlList[i].timeslot = 154;
            case 5:
                chnlList[i].timeslot = 82;
            default:
                break;
        }
        
    }
    
    
 
    

   SetMinSF(endDevices, gateways, channel, nodeList);    
   // FirstFit(node_array,NODE_NUM);
    BestFit(nodeList, NODE_NUM);
  //  WorstFit(nodeList, NODE_NUM);
    printChannels();

  
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


  
  return 0;
}
