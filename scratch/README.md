# RTLoRa
This is a custom-built NS-3 Simulation script for Real-Time LoRa on top of the lorawan ns3 module from https://github.com/signetlabdei/lorawan. 


#To use 
#install ns3 
#install lorawan module (instruction in https://github.com/signetlabdei/lorawan) 
#clone the scratch directory in your local scratch directory
#run the following from the ns3 root directory

./waf --run RT_lora

You can set the simulation parameters (simulationtime, number of channels, number of nodes etc) via cmd arguments
run with --help to see the available cmd arguments. s
