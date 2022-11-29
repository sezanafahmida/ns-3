#include "Lora_battery.h"
namespace ns3 {
namespace lorawan {



template<typename L>
std::vector<float> 
loraBattery::linspace(L start_in, float end_in, int num_in)
{
	std::vector<float> linspaced;

	float start = static_cast<float>(start_in);
	float end = static_cast<float>(end_in);
	float num = static_cast<float>(num_in);

	if (num == 0) { return linspaced; }
	if (num == 1)
	{
		linspaced.push_back(start);
		return linspaced;
	}

	float delta = (end - start) / (num - 1);

	for (int i = 0; i < num - 1; ++i)
	{
		linspaced.push_back(start + delta * i);
	}
	linspaced.push_back(end); // I want to ensure that start and end
							  // are exactly the same as the input
	return linspaced;
}
template<typename S>
std::vector<size_t>
loraBattery::sort_indexes(const std::vector<S>& v) {
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);
	stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });

	return idx;
}

float
loraBattery::initialCharge() {
	float initialEnergy = Chargelevel.back();
		return initialEnergy;
}

int 
loraBattery::getRetrans(int t) {
	/*std::vector<float> Vmaxr; // { 10, 5, 50, 333, 20, 30 };

	for (int r = 0; r < R; r++) {
		Vmaxr.push_back(Prob[t][r]); // Prob[p][t][r]
		//cout << "element: " << Vmaxr[r] << endl;
	}
	//std:: cout << " Max-element: " << std::max_element(Vmaxr.begin(), Vmaxr.end()) - Vmaxr.begin() << std::endl;
	//cout << Vmaxr << endl;
	return max_element(Vmaxr.begin(), Vmaxr.end()) - Vmaxr.begin();*/ 
        for (int r = 0; r < R; r++) {
            if(Prob[t][r]>0.9) return r;
         }
         return 9;      
}

void
loraBattery::SortGamma(int p) {
	std::vector<float> gammaMax;
	for (int t = 0; t < T; t++) {
		gammaMax.push_back(gamma[t]);
	}
	
}

void
loraBattery::ProbabilityUpdate(int t) {

                 
                 Prob[t][0] = SlotSuccess[t][0]/ SelectedSlots[t];
                 for (int r = 1; r < R+1; r++){
                 Prob[t][r] =  ((Prob[t][r - 1]) + ((SlotSuccess[t][r] / SelectedSlots[t])));}
                 
		/*for (int t = 0; t < T; t++) {
			for (int r = 0; r < R; r++) {
                               // Prob[t][r] = Prob[t][r] * (1 - Prob[t][i]);
			        // Prob[t][r] = float(SlotSuccess[t][r] / SelectedSlots[t]);
                               //  Prob[t][r] = 1;
			      /*   for(int i = 0; i < r; i++) 
                                  {
				    Prob[t][r] = Prob[t][r] * (1 - Prob[t][i]);
			          }
			          Prob[t][r] = Prob[t][r] * (1- (SlotSuccess[t][r] / SelectedSlots[t]));
                               //std::cout << "Node: "<< id <<  " t " << t <<  " r " << r << " Prob " << Prob[t][r] << "SlotSuccess" << SlotSuccess[t][r] << " SelectedSlots " <<SelectedSlots[t] << "\n";
			}

		//}
                //for (int t = 0; t < T; t++) {
		for (int r = 0; r < R; r++) {
			Prob[t][r] = 1;
			for(int i = 0; i < r; i++) {
				Prob[t][r] = Prob[t][r] * (1 - Prob[t][i]);
			}
			Prob[t][r] = Prob[t][r] * (SlotSuccess[t][r] / SelectedSlots[t]);
                        //std::cout << "Node: "<< id <<  " t " << t <<  " r " << r << " Prob " << Prob[t][r] << "SlotSuccess" << SlotSuccess[t][r] << " SelectedSlots " <<SelectedSlots[t] << "\n";
		//}
	}*/ 
       
     
}

void 
loraBattery::initialize(int tsLen,std::vector <double> greenSource,int gsLen,float t, int heuristic) 
{
	theta = t;
        h= heuristic;
        EstHarvestedE.assign(T*P+T, 0.0);
        
        Chargelevel.assign(T*P,0.0);
	/*for (int i = 0; i < P; i++)
	{
		std::vector<std::vector<float>> v2d;
		std::vector<std::vector<float>> v2df;
		
		std::vector<float> s1d;
		std::vector<float> g1df;*/
                   
		for (int j = 0; j < T; j++)
		{
			std::vector<float> v1d;
			std::vector<float> v1df;
			//g1d.push_back(0);
			//s1d.push_back(1);
			//g1df.push_back(0.0);
                        float temp =5.0;
                        SelectedSlots.push_back(temp);
			for (int k = 0; k < temp; k++)
			{
				v1d.push_back(1);
				v1df.push_back(1/temp);
			}
                        for (int l = temp+1; l < R+1; l++)
			{
				v1d.push_back(0);
				v1df.push_back(0);
			}
                        
			//v2d.push_back(v1d);
			//v2df.push_back(v1df);
			//EstHarvestedE.push_back(0.0);
			//Chargelevel.push_back(0.0);
                        SlotSuccess.push_back(v1d);
                        Prob.push_back(v1df);
		}

		
		//X.push_back(g1d); // g1df
		
		
	//}

        for (int p=0;p<P;p++) 
        { 
          std::vector<int> g1d;
          for(int t=0;t<T;t++)
          {
           
            g1d.push_back(0);
            
          
          }
          X.push_back(g1d);
        }
	std::string value;
	int day1 = curDay% 365;
        
	int Cmin = day1 * 1440;
        
	for (int i = 0; i < P*T; i++)
	{              
		double greenEnergy = 0.0;
		for (int q = 0; q < (tsLen/60) ; q++)
		{	
                        if(isLading) greenEnergy+= greenSource.at(Cmin%gsLen)*3; //GetTxEnergy(0) * 3	
			else greenEnergy+= greenSource.at(Cmin%gsLen); //GetTxEnergy(0) * 3	
			Cmin++;  
                       // if(id<2) std::cout<< "Energy harvested for min " << day1 << " " << Cmin << " " << greenEnergy << "\n";
		}
		EstHarvestedE[i] = greenEnergy;
              //  std::cout<< "Energy harvested in Decide TS() for timeslot " << i << " " << EstHarvestedE[i] << "\n";
	}
	//curDay++;
       // std::cout<< " Node " << id << " day " << curDay << std::endl;
       for (int i =0;i<T;i++) 
       {
          
         TSVariation.push_back(0);
         
       
       }
       
       

      u = linspace(1, 0.1, T);  //std::vector<float>(T, 1.0); //utility 
      
      
}

void 
loraBattery::decideTS2()
{
        float Ebat = bCap*theta; //0.008; // Wh**** product of e_required
	//float Edod = 0.2 * Ebat;// Depth of Discharege
	float Emax = 1 * Ebat; // Max charge to store in battery
	float Emin = 0 * Ebat;
        int p = curP-1;   
	//normAge=1;
	
	int retransNo;


	Chargelevel[0] = curSOC;//0.003; // initialCharge();
        Chargelevel[p*T] = curSOC;
        
        w.assign(T, 0.0);
	gamma.assign(T, 0.0);
        
        for (int t = 0; t < T; t++) {
		retransNo = getRetrans(t)+1;
                if(retransNo>8) retransNo = 1000; //in this case the packet failed
		w[t] = (EstHarvestedE[p*T + t] - std::max(EstHarvestedE[p*T + t] , estimatedEnergyRequired * retransNo)) /  GetTxEnergy(0);

                //(*std::max_element(greenSource.begin(),greenSource.end())*(tsLen/60));
		gamma[t] =  (u[t] +  ( normAge* w[t]));
               
	}

       for (int t = 0; t < T; t++) {
		if (Chargelevel[T * p + t] + EstHarvestedE[p*T + t] < Emax) 
                {
		 Chargelevel[T * p + t + 1] = Chargelevel[T * p + t] + EstHarvestedE[p*T + t];
		}
		else {
		Chargelevel[T * p + t + 1] = Emax;
		}
	}

        std::vector<float> v;
	for (int t = 0; t < T; t++) {
		v.push_back(gamma[t]);
	}
        std::vector<long unsigned int>temp_vec = sort_indexes(v);
	for (int i : temp_vec) {
        retransNo = getRetrans(i)+1;
        if ((Chargelevel[T * p + i] - (estimatedEnergyRequired * retransNo)) >= Emin) {
					X[p][i] = 1;
                                        curSlot = i;
                                        TSVariation[i]+=1;
                                        EstRetrans = getRetrans(curSlot)+1;
                        //if(i==0) std::cout<< "Node " << id << "slot" << curSlot << " Est retrans" << EstRetrans << std::endl;
			                break;
                                 }
           }

}

void 
loraBattery::decideTS()
{
	float Ebat = bCap*theta; //0.008; // Wh**** product of e_required
	//float Edod = 0.2 * Ebat;// Depth of Discharege
	float Emax = 1 * Ebat; // Max charge to store in battery
	float Emin = 0 * Ebat;
        int p = curP-1;   
	//normAge=1;
	
	int retransNo;


	Chargelevel[0] = curSOC;//0.003; // initialCharge();
        Chargelevel[p*T] = curSOC;
         
        w.assign(T, 0.0);
	gamma.assign(T, 0.0);
	for (int t = 0; t < T; t++) {
		retransNo = getRetrans(t)+1;
                if(retransNo>8) retransNo = 1000; //in this case the packet failed
		w[t] = (EstHarvestedE[p*T + t] - std::max(EstHarvestedE[p*T + t] , estimatedEnergyRequired * retransNo)) /  GetTxEnergy(0);

                //(*std::max_element(greenSource.begin(),greenSource.end())*(tsLen/60));
		if(h == 1) gamma[t] =  (u[t] +  ( normAge* w[t]));
                //else if(h==2) gamma[t] = normAge* w[t];
                //  std::cout <<"t " << t << " U " << u[t] << " W " << w[t] << " E_req " << estimatedEnergyRequired * retransNo << " E_h " << EstHarvestedE[p*T + t]  << " Gamma " << gamma[t] << std::endl;
	}


	std::vector<float> v;
	for (int t = 0; t < T; t++) {
		v.push_back(gamma[t]);
	}
        std::vector<long unsigned int>temp_vec = sort_indexes(v);
	for (int i : temp_vec) {
		int t = 0;
		bool stop = false;
		while (t < T && stop == false) {
			retransNo = getRetrans(t) + 1; // as r starts from 0
                        if(retransNo>8) retransNo = 1000;
		        budget = fabs(EstHarvestedE[p*T + t] + (estimatedEnergyRequired * retransNo));
                       // std::cout << "EstE "<< EstHarvestedE[T * p + t] << "EstR "<< estimatedEnergyRequired << "Retr " << retransNo << "\n";
			if (i == t) {
				if ((Chargelevel[T * p + t] - (estimatedEnergyRequired * retransNo)) >= Emin) {
					X[p][t] = 1;
                                      //  std:: cout <<" Node "<< id << " period " << p << "," << " timeslot "<< t << "--" << Chargelevel[T * p + t] << "     X-- " << X[p][t] << "  -- REQ: " << estimatedEnergyRequired * retransNo  << " || E: " << EstHarvestedE[T * p + t] <<  " || retransNo: " << retransNo << " --  w " << w[t] << " || Gamma: "<< gamma[t] << std::endl;
					if (Chargelevel[T * p + t] + EstHarvestedE[p*T + t] - estimatedEnergyRequired * retransNo < Emax) {
						Chargelevel[T * p + t + 1] = Chargelevel[T * p + t] + EstHarvestedE[p*T + t] - estimatedEnergyRequired * retransNo;
					}
					else {
						Chargelevel[T * p + t + 1] = Emax;
					}
				}
				else {
					stop = true;
				}
			}
			else {
				if (Chargelevel[T * p + t] + EstHarvestedE[p*T + t] - estimatedEnergyRequired * retransNo < Emax) {
					Chargelevel[T * p + t + 1] = Chargelevel[T * p + t] + EstHarvestedE[p*T + t];
				}
				else {
					Chargelevel[T * p + t + 1] = Emax;
				}
			}
			//std:: cout <<" Node "<< id << " period " << p << "," << " timeslot "<< t << "--" << Chargelevel[T * p + t] << "     X-- " << X[p][t] << "  -- REQ: " << estimatedEnergyRequired * retransNo  << " || E: " << EstHarvestedE[T * p + t] <<  " || retransNo: " << retransNo << " --  w " << w[t] << " || Gamma: "<< gamma[t] << std::endl;
			t++;
		}
		if (X[p][i] == 1) {
                        curSlot = i;
                        TSVariation[i]+=1;
                        EstRetrans = getRetrans(curSlot)+1;
                        //if(i==0) std::cout<< "Node " << id << "slot" << curSlot << " Est retrans" << EstRetrans << std::endl;
			break;
		}
		//std::cout << std::endl << std::endl;
	}
	//system("pause>0");
}

double 
loraBattery::GetToA(int dr) {
return ToA[dr]/1000;

}



double 
loraBattery:: GetTxEnergy(int dr) {
double txCurrent = 29.0/1000;
double txEnergy = txCurrent* GetToA(dr)* 3.3;
return txEnergy ;
}

}
}

