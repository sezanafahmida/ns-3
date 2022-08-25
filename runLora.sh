node_num=100
l=120 #timeslot length
l_min=$((l/60))
initD=0
h=0
theta=0.5
simDay=365
rm /home/gp7532/ns-3/SOC.csv
rm /home/gp7532/ns-3/day.csv
rm /home/gp7532/ns-3/Mean/mean*
rm /home/gp7532/ns-3/Trace/trace*
rm /home/gp7532/ns-3/Prob/prob*
rm /home/gp7532/ns-3/node_info.json

for ((day=1;day<=$simDay;day++))
  do
  echo "$node_num, $day, $h"
 ./waf --run "complete-network-example --nDevices=$node_num --l=$l --h=$h --day=$day --theta=$theta" >sim_run.txt
  python3 Battery_degradation.py $l_min $day $node_num $initD 
  ret=$?
  if [ $ret -ne 0 ]; then
     break
  fi
  if [ $((day%30)) -eq 0 ] 
  then
     #cp /home/gp7532/ns-3/output.csv /home/gp7532/ns-3/Lifespan_LoRa/output_lora_day_${day}_${node_num}.csv
     cp /home/gp7532/ns-3/d_out.csv /home/gp7532/ns-3/Lifespan_LoRa/dout_lora_d${day}_n${node_num}_t${theta}_i${initD}.csv
     cp /home/gp7532/ns-3/GE.csv /home/gp7532/ns-3/Lifespan_LoRa/GE_lora_d${day}_n${node_num}_t${theta}_h${h}_i${initD}.csv
     cp /home/gp7532/ns-3/AvgError.csv /home/gp7532/ns-3/Lifespan_LoRa/AvgError_lora_d${day}_n${node_num}_t${theta}_h${h}_i${initD}.csv
  fi
  rm /home/gp7532/ns-3/Trace/trace*
  done
mv /home/gp7532/ns-3/output.csv /home/gp7532/ns-3/Result_Lora/output_lora_d${day}_n${node_num}_t${theta}_i${initD}.csv
mv /home/gp7532/ns-3/d_out.csv /home/gp7532/ns-3/Result_Lora/dout_lora_d${day}_n${node_num}_t${theta}_i${initD}.csv
mv /home/gp7532/ns-3/AvgLatency.csv /home/gp7532/ns-3/Result_Lora/AvgLatency_lora_d${day}_n${node_num}_t${theta}_i${initD}.csv
mv /home/gp7532/ns-3/AvgError.csv /home/gp7532/ns-3/Result_Lora/AvgError_lora_d${day}_n${node_num}_t${theta}_i${initD}.csv
mv /home/gp7532/ns-3/sim_results_0122.csv /home/gp7532/ns-3/Result_Lora/sim_result_lora_d${day}_n${node_num}_t${theta}_i${initD}.csv
mv /home/gp7532/ns-3/timeslot.csv /home/gp7532/ns-3/Result_Lora/timeslot_lora_d${day}_n${node_num}_t${theta}_i${initD}.csv
mv /home/gp7532/ns-3/receptionSlot.csv /home/gp7532/ns-3/Result_Lora/receptionSlot_lora_d${day}_n${node_num}_t${theta}_i${initD}.csv
mv /home/gp7532/ns-3/prr.csv /home/gp7532/ns-3/Result_Lora/prr_lora_d${day}_n${node_num}_t${theta}_i${initD}.csv
mv /home/gp7532/ns-3/avgR.csv /home/gp7532/ns-3/Result_Lora/avgR_lora_d${day}_n${node_num}_t${theta}_i${initD}.csv
mv /home/gp7532/ns-3/txE.csv /home/gp7532/ns-3/Result_Lora/txE_lora_d${day}_n${node_num}_t${theta}_i${initD}.csv
mv /home/gp7532/ns-3/GE.csv /home/gp7532/ns-3/Result_Lora/GE_lora_d${day}_n${node_num}_t${theta}_h${h}_i${initD}.csv
