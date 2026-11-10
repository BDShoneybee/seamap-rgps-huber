// Portable research reconstruction of the deployed Huber update equations.
// Input CSV: time,value (nan denotes an unavailable observation).
// Six-observation warmup; chronological input; missing queries do not mutate state.
// This is not the Qt/Eigen production application. Build: g++ -std=c++14 -O2.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
double median(std::vector<double> v) {
    std::sort(v.begin(),v.end());auto n=v.size();
    return n%2?v[n/2]:(v[n/2-1]+v[n/2])/2;
}
double mad(std::vector<double> v) {
    double m=median(v);for(auto& z:v)z=std::abs(z-m);return median(v);
}
int main() {
    std::vector<double> bt,by; bool ready=false,hasdy=false;
    double x=0,v=0,s=0,s0=0,a=0,last=0,prev=0,pdy=0,prior_time=-1e300;
    std::string line;std::cout<<"time,estimate,rejected\n"<<std::setprecision(17);
    while(std::getline(std::cin,line)) {
        if(line.empty()||line.rfind("time",0)==0)continue;
        auto c=line.find(',');double t=std::stod(line.substr(0,c)),z=std::stod(line.substr(c+1));
        if(t<=prior_time)throw std::runtime_error("Input must have increasing timestamps");
        prior_time=t;double out=z;bool rejected=false;
        if(!ready) {
            if(std::isfinite(z)){bt.push_back(t);by.push_back(z);}
            if(by.size()==6) {
                std::vector<double> dy,sl;
                for(unsigned i=1;i<6;i++){dy.push_back(by[i]-by[i-1]);sl.push_back(dy.back()/(bt[i]-bt[i-1]));}
                s0=std::max(.35*1.4826*mad(dy),1e-4);s=s0;v=median(sl);x=z;last=t;prev=z;ready=true;
            }
        } else {
            double dt=std::max(t-last,1e-3),vp=v*std::exp(-dt/12),xp=x+vp*dt;
            if(!std::isfinite(z))out=dt>100?NAN:(dt>8?x:xp);
            else {
                double nu=z-xp,dy=z-prev;
                if(hasdy)a=dy*pdy<0?std::min(1.,.88*a+.12):std::max(0.,.90*a-.06);
                pdy=dy;hasdy=true;prev=z;
                double sn=std::max(.97*s+.03*std::abs(nu),.35*s0);
                rejected=std::abs(nu)>4*s;
                double w=std::min(1.,2.5*s/std::max(std::abs(nu),1e-12));
                double kx=(1-std::exp(-dt/10))*(1-.45*a);
                double kv=(1-std::exp(-dt/24))/dt*(1-.60*a);
                x=rejected?xp:xp+kx*w*nu;v=rejected?vp:vp+kv*w*nu;
                s=sn;last=t;out=x;
            }
        }
        std::cout<<t<<','<<out<<','<<rejected<<'\n';
    }
}
