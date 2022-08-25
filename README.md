# Battery Lifespan

NS3 simulation script for battery lifespan paper.

- execute ./runH.sh for running proposed heuristic
- Aggregate Results stored in folder "Result_h"
- Detailed lifespan results stored in folder "Lifespan_heuristic"

- execute ./runLora.sh for running baseline LoRaWAN
- Aggregate Results stored in folder "Result_Lora"
- Detailed lifespan results stored in folder "Lifespan_LoRa" 
parameters
node_num=100
l=120 #timeslot length
theta=0.5 #Soc Threshold
h=1 #enable/disable heuristic
initD=0 #initial battery degradation
simDay #simulation time in days.

Bash script->ns3 script for one day->writes battery soc in file->python script to generate battery degradation->writes battery degradation to file->ns3 reads this file again at the beginning of each iteration. 
