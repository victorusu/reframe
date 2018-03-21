#ifndef PAIRLIST_HPP
#define PAIRLIST_HPP

#if defined(Thrust)
#include "pairlist_gpu.hpp"
#else
// #include "pairlist_cpu.hpp"
// #include "pairlist_peratom_cpu.hpp"
#include "pairlist_peratom_domain_cpu.hpp"
#endif


#endif //PAIRLIST_HPP