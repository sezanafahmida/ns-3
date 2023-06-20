import numpy as np

import pandas as pd
import rainflow  
import statistics
from math import  exp
import os
import sys
import random
from os import path

random.seed(0);

TimeSlots_Length = int(sys.argv[1])  #2 #min
day = int(sys.argv[2])
N = int(sys.argv[3])
initD = int(sys.argv[4])
def Linear_degradation(soc,days,prev_cyl,id,day,prevAvg,l,T, prev_cal):
    
    n = 0
    t = TimeSlots_Length * 60
   # T = 25  #temperature
    t_tot = (3600 * 24)  #* (days)
    
    dischargeDepth = []
    meanSoc = []
    cycleCount = []
    cyclePeriod = []
    Crate = []
    d_SoC = []
    d_DoD = []
    
    for rng, mean, count, i_start, i_end in rainflow.extract_cycles(soc): 
        dischargeDepth.append(rng)
        meanSoc.append(mean)
        cycleCount.append(count)
        cyclePeriod.append(i_end - i_start)
        Crate.append(rng * 2 * count / (t/3600))


    dischargeDepth.pop(0)
    meanSoc.pop(0)
    cycleCount.pop(0)

    # SoC stress model
    k_soc = 1.039
    SoC_ref = .6
    
    # DoD stress model
    k_DoD2 = 2.03 # DoD stress model nonlinear coefficient
    k_DoD1 = .2/(3000*.8 ** k_DoD2) # 3000 cycles @ 80% DoD till 80% end of life
    
    #cell temperature effect
    k_t = 0.0693
    T_ref = 25
    d_temp = exp(k_t * (T-T_ref) * (273 + T_ref)/(273 + T))
    
    # C-rate effect
    d_Crate = 1

    for i in range(0, len(dischargeDepth)):
        d_SoC.append(exp(k_soc * (meanSoc[i] - SoC_ref))) # SoC stress model
        d_DoD.append(k_DoD1 * dischargeDepth[i] ** k_DoD2) # DoD stress model
    
    # calender ageing
    k_cal = 3.31e-9/8
    NewAvg =  statistics.mean(meanSoc) 
    SoC_avg  = (prevAvg*l + NewAvg*len(meanSoc) ) /(l+len(meanSoc))
    l=l+len(meanSoc)  #no of charge-discharge cycle
    test = sum(d_DoD[i] * d_SoC[i] * d_Crate * d_temp * cycleCount[i] for i in range(0, len(dischargeDepth))) 
    
    d_cycle = prev_cyl + sum(d_DoD[i] * d_SoC[i] * d_Crate * d_temp * cycleCount[i] for i in range(0, len(dischargeDepth))) # total cycle ageing
    

    d_cal = prev_cal + (k_cal * t_tot * d_temp * exp(k_soc * (SoC_avg - SoC_ref)))   #0.5 -1.5 
    #if id == 1:
      #  print("dcal ", d_cal )
    d = d_cycle + d_cal
    #if id==80 and days>=100:
       # print( "dcal updated", d_cal, " days ", days, " t_tot ", t_tot, " k_cal ", k_cal, "d_temp", d_temp , "K_soc" , k_soc, " SOC avg ", SoC_avg , " ref ", SoC_ref)
    return round(d,10),round(d_cal,10),round(d_cycle,10),round(SoC_avg,10), l, round(d_temp,10)


def Nonlinear_degradation(d):
    a = .0575
    b = 121  
    L = 1 - a*exp(-b*d) - (1-a)*exp(-d)
    return L               

#print(sys.argv)
path1 = '/home/gp7532/ns-3/Trace'
#files = os.listdir(path1)
dFilename = '/home/gp7532/ns-3/d_out.csv'
#output = np.zeros(N);
fout = open("/home/gp7532/ns-3/output.csv","w")

#sorted_files = sorted(files, key=lambda x: int(x.split('trace_')[1].split(".")[0])) #sort files numerically
id = 0
if day>1:
    data=pd.read_csv(dFilename,names=["d","d_cycle","d_cal","d_nl","SoC_avg","len","d_temp"])
dout =  open("/home/gp7532/ns-3/d_out.csv","w")
tempData = pd.read_csv('/home/gp7532/ns-3/mi_daily_temp.csv')
#T = random.gauss(tempData.loc[tempData['newDay']==day%365]['Temperature'], 5) 
#print(T)
for id in range(0,N): 
    filename = path1+ "/" + "trace_" +str(id)+".csv"
    #print(filename) 
    if(initD ==1):
        dayOffset = random.randint(1,100);
    else:
        dayOffset = 0;#random.randint(1,100);
    #print(dayOffset);
    soc = np.loadtxt(filename, dtype =float)   #get the generated SOC trace
   # print(soc)
    if day>1:
        #print(data)
        prev_cyl = data.iloc[id,1]

        prev_cal = data.iloc[id,2]
        avg = data.iloc[id,4]
        l = data.iloc[id,5]
    else: 
        prev_cyl  = 0
        prev_cal  = 0
        avg = 0
        l=0
   #print(mean)
    T = random.gauss(tempData.loc[tempData['newDay']==day%365]['Temperature'], 5) 
    print(T)
    d,d_cal,d_cycle,SoC_avg,l,d_temp = Linear_degradation(soc,day+dayOffset,prev_cyl,id,day,avg,l,T,prev_cal) #get the new degradation, D_cycle
    z = Nonlinear_degradation(d)  #get non-linear Degradation
    if z>=0.2:
        print("battery died on day: ", day)
        sys.exit(2)
  #  print(z)
    print(z ,file=fout)           #store degradation in the output file
    print(d, ",", d_cycle,",",d_cal, ",", z,",",SoC_avg, ",",l, ",",d_temp, file=dout) 
    id+=1                                   #increment nodeID counter 

