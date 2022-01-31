import numpy as np

import pandas as pd
import rainflow  
import statistics
from math import  exp
import os
import sys

TimeSlots_Length = int(sys.argv[1])  #2 #min
day = int(sys.argv[2])
N = int(sys.argv[3])
def Linear_degradation(soc,days,d_cycle):
    
    n = 0
    t = TimeSlots_Length * 60
    T = 25
    t_tot = (3600 * 24)  * (days)
    
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
    SoC_avg =  statistics.mean(meanSoc)  
    d_cycle = d_cycle + sum(d_DoD[i] * d_SoC[i] * d_Crate * d_temp * cycleCount[i] for i in range(0, len(dischargeDepth))) # total cycle ageing
    
    d_cal = k_cal * t_tot * d_temp * exp(k_soc * (SoC_avg - SoC_ref))   
    d = d_cycle + d_cal
    #print(d)
    return d,d_cal,d_cycle


def Nonlinear_degradation(d):
    a = .0575
    b = 121  
    L = 1 - a*exp(-b*d) - (1-a)*exp(-d)
    return L               

#print(sys.argv)
path = '/home/gp7532/ns-3/Trace'
files = os.listdir(path)
#output = np.zeros(N);
fout = open("/home/gp7532/ns-3/output.csv","w")
dout =  open("/home/gp7532/ns-3/d_out.csv","w")
sorted_files = sorted(files, key=lambda x: int(x.split('trace_')[1].split(".")[0])) #sort files numerically
id = 0

for file in sorted_files: 
    filename = path+ "/" + file
    #print(filename) 
    meanFile = "/home/gp7532/ns-3/Mean/mean_"+ str(id)

    
    soc = np.loadtxt(filename, dtype =float)   #get the generated SOC trace
   # print(soc)
    if os.path.isfile(meanFile):
        prev = np.loadtxt(meanFile, dtype = float) #get the stored d_cycle
    else: 

        prev  =0
   #print(mean)
    d,d_cal,d_cycle = Linear_degradation(soc,day,prev) #get the new degradation, D_cycle
    z = Nonlinear_degradation(d)  #get non-linear Degradation
    if z>=0.2:
        print("battery died on day: ", day)
        sys.exit(2)
  #  print(z)
    print(z ,file=fout)           #store degradation in the output file
    print(d, ",", d_cycle,",",d_cal, ",", z, file=dout)
    meanOut = open(meanFile,"w")  
    print(d_cycle,file = meanOut)  #store updated d_cycle 
    #print(SoC_avg,d_cycle) 
    id+=1                                   #increment nodeID counter 
