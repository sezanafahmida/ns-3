#!/bin/bash
node_num=5
l=60 #timeslot length
l_min=$((l/60))
theta=0.5
h=1
initD=0
m=1;
day=2
nw=1 #new weight
incr=0 #new nodes
for theta in {0.05,0.5,1} # $(seq 0.05 0.45 0.5)
do
#for node_num in {100,250,500}
#do
#for nw in {0.1,0.4,0.8}
#do
rm /home/gp7532/ns-3/SOC.csv
rm /home/gp7532/ns-3/day.csv
rm /home/gp7532/ns-3/Mean/mean*
rm /home/gp7532/ns-3/Trace/trace*
rm /home/gp7532/ns-3/Prob/prob*
rm /home/gp7532/ns-3/node_info.json
rm /home/gp7532/ns-3/AvgError.csv
#rm /home/gp7532/ns-3/Lifespan_heuristic/*.csv

#for h in {1}
#do((day=1;day<=$simDay;day++))
#for ((day=1;day<=$simDay;day++))

for day in $(seq 1 1 2)
  do
  #if [ $((day%2)) -eq 0 ] ;then
  #   node_num=$((node_num+100))
  #   incr=100
 # fi
  echo "$node_num, $day , $theta ,$h"
 ./waf --run "complete-network-example --nDevices=$node_num --l=$l --h=$h --day=$day --theta=$theta --b=0 --nw=$nw --incr=$incr">sim_run_2.txt
  python3 Battery_degradation.py $l_min $day $node_num $initD
  ret=$?
  if [ $ret -ne 0 ]; then
     break
  fi
  if [ $((day%30)) -eq 0 ] 
  then 
     #cp /home/gp7532/ns-3/output.csv /home/gp7532/ns-3/Lifespan_heuristic/output_h_day_${day}_${node_num}.csv
     cp /home/gp7532/ns-3/d_out.csv /home/gp7532/ns-3/Lifespan_heuristic/dout_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}.csv
     cp /home/gp7532/ns-3/GE.csv /home/gp7532/ns-3/Lifespan_heuristic/GE_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}.csv
    # cp /home/gp7532/ns-3/AvgError.csv /home/gp7532/ns-3/Lifespan_heuristic/AvgError_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}.csv
  
  fi
  if [ $day -eq 1 ];then
     cp /home/gp7532/ns-3/d_out.csv /home/gp7532/ns-3/Lifespan_heuristic/dout_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}.csv
  fi

  
  rm /home/gp7532/ns-3/Trace/trace*
  done

mv /home/gp7532/ns-3/output.csv /home/gp7532/ns-3/Result_h/output_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}_nw_${nw}.csv
mv /home/gp7532/ns-3/d_out.csv /home/gp7532/ns-3/Result_h/dout_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}_nw_${nw}.csv
mv /home/gp7532/ns-3/AvgLatency.csv /home/gp7532/ns-3/Result_h/AvgLatency_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}_nw_${nw}.csv
#mv /home/gp7532/ns-3/AvgError.csv /home/gp7532/ns-3/Result_h/AvgError_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}.csv
mv /home/gp7532/ns-3/sim_results_0122.csv /home/gp7532/ns-3/Result_h/sim_result_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}_nw_${nw}.csv
mv /home/gp7532/ns-3/timeslot.csv /home/gp7532/ns-3/Result_h/timeslot_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}_nw_${nw}.csv
mv /home/gp7532/ns-3/receptionSlot.csv /home/gp7532/ns-3/Result_h/receptionSlot_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}_nw_${nw}.csv
mv /home/gp7532/ns-3/prr.csv /home/gp7532/ns-3/Result_h/prr_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}_nw_${nw}.csv
mv /home/gp7532/ns-3/avgR.csv /home/gp7532/ns-3/Result_h/avgR_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}_nw_${nw}.csv
mv /home/gp7532/ns-3/txE.csv /home/gp7532/ns-3/Result_h/txE_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}_nw_${nw}.csv
mv /home/gp7532/ns-3/GE.csv /home/gp7532/ns-3/Result_h/GE_h_d${day}_n${node_num}_t${theta}_h${h}_i${initD}_nw_${nw}.csv
done
#done
#done
