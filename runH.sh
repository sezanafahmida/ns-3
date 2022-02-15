node_num=100
l=120 #timeslot length
l_min=$((l/60))
theta=0.8
h=1
##for theta in $(seq 0.4 0.1 0.8)
for node_num in {100..100..100}
do
rm /home/gp7532/ns-3/SOC.csv
rm /home/gp7532/ns-3/day.csv
rm /home/gp7532/ns-3/Mean/mean*
rm /home/gp7532/ns-3/Trace/trace*
rm /home/gp7532/ns-3/Prob/prob*
rm /home/gp7532/ns-3/node_info.json
rm /home/gp7532/ns-3/AvgError.csv
#rm /home/gp7532/ns-3/Lifespan_heuristic/*.csv


for day in {1..1..1}
  do

  echo "$node_num, $day , $theta ,$h"
 ./waf --run "complete-network-example --nDevices=$node_num --l=$l --h=$h --day=$day --theta=$theta --b=0" --gdb # >sim_run.txt
  python3 Battery_degradation.py $l_min $day $node_num
  ret=$?
  if [ $ret -ne 0 ]; then
     break
  fi
  if [ $((day%30)) -eq 0 ] 
  then
     #cp /home/gp7532/ns-3/output.csv /home/gp7532/ns-3/Lifespan_heuristic/output_h_day_${day}_${node_num}.csv
     cp /home/gp7532/ns-3/d_out.csv /home/gp7532/ns-3/Lifespan_heuristic/dout_h_${day}_${node_num}_${theta}_$h.csv
  
  fi
  rm /home/gp7532/ns-3/Trace/trace*
  done

mv /home/gp7532/ns-3/output.csv /home/gp7532/ns-3/Result_h/output_h_${day}_${node_num}_${theta}_$h.csv
mv /home/gp7532/ns-3/d_out.csv /home/gp7532/ns-3/Result_h/dout_h_${day}_${node_num}_${theta}_$h.csv
mv /home/gp7532/ns-3/AvgLatency.csv /home/gp7532/ns-3/Result_h/AvgLatency_h_${day}_${node_num}_${theta}_$h.csv
mv /home/gp7532/ns-3/AvgError.csv /home/gp7532/ns-3/Result_h/AvgError_h_${day}_${node_num}_${theta}_$h.csv
mv /home/gp7532/ns-3/sim_results_0122.csv /home/gp7532/ns-3/Result_h/sim_result_h_${day}_${node_num}_${theta}_$h.csv
mv /home/gp7532/ns-3/timeslot.csv /home/gp7532/ns-3/Result_h/timeslot_h_${day}_${node_num}_${theta}_$h.csv
mv /home/gp7532/ns-3/receptionSlot.csv /home/gp7532/ns-3/Result_h/receptionSlot_h_${day}_${node_num}_${theta}_$h.csv
done
