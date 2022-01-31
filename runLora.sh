node_num=100
l=120 #timeslot length
l_min=$((l/60))
rm /home/gp7532/ns-3/SOC.csv
rm /home/gp7532/ns-3/day.csv
rm /home/gp7532/ns-3/Mean/mean*
rm /home/gp7532/ns-3/Trace/trace*
rm /home/gp7532/ns-3/Prob/prob*
rm /home/gp7532/ns-3/node_info.json

for day in {1..3650..1}
  do
  echo "$node_num, $l , $l_min, $day"
 ./waf --run "complete-network-example --nDevices=$node_num --l=$l --h=0 --day=$day" >sim_run.txt
  python3 /home/gp7532/Battery_degradation.py $l_min $day $node_num
  ret=$?
  if [ $ret -ne 0 ]; then
     break
  fi
  if [ $((day%30)) -eq 0 ] 
  then
     #cp /home/gp7532/ns-3/output.csv /home/gp7532/ns-3/Lifespan_LoRa/output_lora_day_${day}_${node_num}.csv
     cp /home/gp7532/ns-3/d_out.csv /home/gp7532/ns-3/Lifespan_LoRa/dout_lora_day_${day}_${node_num}.csv
  fi
  rm /home/gp7532/ns-3/Trace/trace*
  done
mv /home/gp7532/ns-3/output.csv /home/gp7532/ns-3/Result_Lora/output_lora_day_${day}_${node_num}.csv
mv /home/gp7532/ns-3/d_out.csv /home/gp7532/ns-3/Result_Lora/dout_lora_day_${day}_${node_num}.csv
mv /home/gp7532/ns-3/AvgLatency.csv /home/gp7532/ns-3/Result_Lora/AvgLatency_lora_day_${day}_${node_num}.csv
mv /home/gp7532/ns-3/AvgError.csv /home/gp7532/ns-3/Result_Lora/AvgError_lora_day_${day}_${node_num}.csv
mv /home/gp7532/ns-3/sim_results_0122.csv /home/gp7532/ns-3/Result_Lora/sim_result_lora_${day}_${node_num}.csv
