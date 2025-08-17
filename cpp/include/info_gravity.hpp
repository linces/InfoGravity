#pragma once
double inv = 1.0 / std::max(rhoI, P.epsGrad);
return gradRhoI * (-P.alpha * inv);
}


class InfoGravity {
public:
explicit InfoGravity(const InfoParams& par, RhoIProvider rhoProv)
: par_(par), rhoProv_(std::move(rhoProv)) {}


// acceleration a_I = -∇Φ_I
Vec3 acceleration(const Vec3& x) const {
double rhoI; Vec3 gR;
rhoProv_(x, rhoI, gR); // obtain rho_I and grad(rho_I)
Vec3 gPhi = gradPhiI(par_, rhoI, gR);
return gPhi * (-1.0); // a_I = -∇Φ_I
}
private:
InfoParams par_;
RhoIProvider rhoProv_;
};


class NBodyStepper {
public:
using AccelFn = std::function<Vec3(const std::vector<Particle>&, size_t)>;


NBodyStepper(AccelFn accelNewtonian, AccelFn accelInfo, double dt)
: aN_(std::move(accelNewtonian)), aI_(std::move(accelInfo)), dt_(dt) {}


void step(std::vector<Particle>& ps) {
const double dt = dt_;
std::vector<Vec3> a0(ps.size()), aHalf(ps.size());


for (size_t i = 0; i < ps.size(); ++i) {
a0[i] = aN_(ps, i) + aI_(ps, i);
}
// drift
for (size_t i = 0; i < ps.size(); ++i) {
ps[i].x = ps[i].x + ps[i].v * dt + a0[i] * (0.5 * dt * dt);
}
// accel at t+dt
for (size_t i = 0; i < ps.size(); ++i) {
aHalf[i] = aN_(ps, i) + aI_(ps, i);
}
// kick
for (size_t i = 0; i < ps.size(); ++i) {
ps[i].v = ps[i].v + (a0[i] + aHalf[i]) * (0.5 * dt);
}
}
private:
AccelFn aN_, aI_;
double dt_;
};