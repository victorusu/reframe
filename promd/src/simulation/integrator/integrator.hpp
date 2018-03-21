#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#if defined(Thrust)
#include "integrator_gpu.hpp"
#else
// #include "integrator_cpu.hpp"
// #include "integrator_velocity_verlet_cpu.hpp"
#include "integrator_leap_frog_cpu.hpp"
#endif


#endif //INTEGRATOR_HPP