#ifndef NONBONDED_HPP
#define NONBONDED_HPP

#if defined(Thrust)
#include "nonbonded_gpu.hpp"
#else
// #include "nonbonded_cpu.hpp"
// #include "nonbonded_peratom_cpu.hpp"
#include "nonbonded_peratom_domain_cpu.hpp"
#endif


#endif //NONBONDED_HPP
