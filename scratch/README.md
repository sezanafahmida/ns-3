# RTLoRa
This is a custom-built NS-3 Simulation script for Real-Time LoRa on top of the lorawan ns3 module from https://github.com/signetlabdei/lorawan. 


To use:


install ns3 


install lorawan module (instruction in https://github.com/signetlabdei/lorawan) 


clone this scratch directory in your local scratch directory.


run the following from your ns3 root directory

```bash
./waf --run RT_lora
```

You can set the simulation parameters (simulationtime, number of channels, number of nodes etc) via cmd arguments


to see avaialble cmd parameters execute: 

```bash
./waf --run "RT_lora --help" 
```


