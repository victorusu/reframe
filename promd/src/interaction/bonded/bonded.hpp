#ifndef BONDED_HPP
#define BONDED_HPP

#if defined(Thrust)
#include "bonded_gpu.hpp"
#else
// #include "bonded_cpu.hpp"
// #include "bonded_peratom_cpu.hpp"
// #include "bonded_cpu.hpp"
#include "bonded_cpu_parallel.hpp"
#endif


#endif //BONDED_HPP
