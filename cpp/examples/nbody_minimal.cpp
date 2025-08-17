#include <iostream>
#include <vector>
#include <cmath>
#include "../include/info_gravity.hpp"


constexpr double Gnewt = 6.67430e-11;


static Vec3 accelNewtonian(const std::vector<Particle>& ps, size_t i) {
Vec3 a{0,0,0};
const auto& pi = ps[i];
for (size_t j = 0; j < ps.size(); ++j) if (j != i) {
auto d = ps[j].x - pi.x;
double r2 = d.x*d.x + d.y*d.y + d.z*d.z + 1e-12;
double r = std::sqrt(r2);
double invr3 = 1.0 / (r2 * r);
a += d * (Gnewt * ps[j].m * invr3);
}
return a;
}


int main(){
// Binario simples para demonstrar a inclusao de a_I
std::vector<Particle> ps(2);
ps[0] = Particle{Vec3{-1.0, 0.0, 0.0}, Vec3{0.0, 0.2, 0.0}, 1.0};
ps[1] = Particle{Vec3{ 1.0, 0.0, 0.0}, Vec3{0.0,-0.2, 0.0}, 1.0};



InfoParams P; P.alpha = 0.8; P.rho0 = 1.0; P.epsGrad = 1e-12;
auto rhoProv = [](const Vec3& x, double& rhoI, Vec3& gradRhoI){
const double sigma2 = 1.0;
double r2 = x.x*x.x + x.y*x.y + x.z*x.z;
rhoI = std::exp(-0.5 * r2 / sigma2);
gradRhoI = Vec3{-x.x * rhoI / sigma2, -x.y * rhoI / sigma2, -x.z * rhoI / sigma2};
};
InfoGravity info(P, rhoProv);


auto aInfo = [&info](const std::vector<Particle>& ps, size_t i)->Vec3{
(void)ps; 
return info.acceleration(ps[i].x);
};


NBodyStepper stepper(accelNewtonian, aInfo, /*dt=*/1e-3);


for (int k=0; k<10000; ++k) {
stepper.step(ps);
if (k % 1000 == 0) {
std::cout << "Step " << k << "\n";
for (size_t i=0;i<ps.size();++i){
std::cout << "p"<<i<<" x="<<ps[i].x.x<<","<<ps[i].x.y<<","<<ps[i].x.z
<<" v="<<ps[i].v.x<<","<<ps[i].v.y<<","<<ps[i].v.z<<"\n";
}
}
}


return 0;
}