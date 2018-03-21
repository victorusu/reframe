-include Makefile.in

all:     noncuda
omp:     promd_omp
noncuda: promd promd_omp
cuda:    promd_cuda
mpi:     promd_mpi

#combined: promd_thrust_omp promd_thrust_tbb

OBJS_PROMD:= $(SRCS:.cu=.o)
OBJS_OMP  := $(SRCS:.cu=.omp.o)
OBJS_CUDA := $(SRCS:.cu=.cuda.o)
OBJS_MPI  := $(SRCS:.cu=.mpi.o)

$(OBJS_PROMD):
	$(CXX) $(CXXFLAGS) -c $(CXXLANGDETECTION) -o $(@) $(@:.o=.cu)

$(OBJS_OMP):
	$(CXX) $(CXXFLAGS) -c $(CXXLANGDETECTION) -fopenmp -o $(@) $(@:.omp.o=.cu)

$(OBJS_CUDA):
	$(NVCC) $(CXXFLAGS) -c $(CXXLANGDETECTION) $(THRUSTCUDA) -o $(@) $(@:.cuda.o=.cu)

$(OBJS_MPI):
	$(CXX) $(CXXFLAGS) -c $(CXXLANGDETECTION) -fopenmp -DUSE_MPI -o $(@) $(@:.mpi.o=.cu)

promd: $(OBJS_PROMD)
	@echo 'Linking: $@'
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(@) $(OBJS_PROMD)
	@echo 'Finished building: $@'
	@echo ' '

promd_mpi: $(OBJS_MPI)
	@echo 'Linking: $@'
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -fopenmp -DUSE_MPI -o $(@) $(OBJS_MPI)
	@echo 'Finished building: $@'
	@echo ' '

promd_omp: $(OBJS_OMP)
	@echo 'Linking: $@'
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -fopenmp -o $(@) $(OBJS_OMP)
	@echo 'Finished building: $@'
	@echo ' '

promd_cuda: $(OBJS_CUDA)
	@echo 'Linking: $@'
	$(NVCC) $(CXXFLAGS) $(LDFLAGS) $(THRUSTCUDA) -o $(@) $(OBJS_CUDA)
	@echo 'Finished building: $@'
	@echo ' '

# promd_thrust_omp: $(SRCS)
# 	@echo 'Linking: $@'
# 	$(CXX) $(CXXFLAGS) $(THRUSTINC) $(THRUSTOMP) $(CXXLANGDETECTION) -o $(@) $(<) -DPRECISION=$(PRECISION)
# 	@echo 'Finished building: $@'
# 	@echo ' '
# promd_tbb: $(SRCS)
# 	@echo 'Linking: $@'
# 	$(CXX) $(CXXFLAGS) $(THRUSTINC) $(THRUSTTBB) $(CXXLANGDETECTION) -o $(@) $(<) -DPRECISION
# 	@echo 'Finished building: $@'
# 	@echo ' '

# Cilk: *.cpp
# 	@echo 'Linking: $@'
# 	$(CXX) $(CXXFLAGS) -fcilkplus -lcilkrts -DCILK -o $(@) $(<:%.cu=%.cpp)
# 	@echo 'Finished building: $@'
# 	@echo ' '

# promd_cuda_omp: $(SRCS)
# 	@echo 'Linking: $@'
# 	$(NVCC) $(CXXFLAGS) -DThrust -DThrustOMP -o $(@) $(<:%.cu=%.cpp) -Xcompiler -fopenmp -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP -lgomp -DPRECISION
# 	@echo 'Finished building: $@'
# 	@echo ' '

# promd_cuda_tbb: $(SRCS)
# 	@echo 'Linking: $@'
# 	$(NVCC) $(CXXFLAGS) -DThrust -o $(@) $(<:%.cu=%.cpp) -Xcompiler -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_TBB -ltbb -DPRECISION
# 	@echo 'Finished building: $@'
# 	@echo ' '


# promd_vexcl: *.cpp
# 	@echo 'Linking: $@'
# 	#$(NVCC) $(CXXFLAGS) -std=c++11 -DVexCL -L/home/hvictor/bin/AMDAPPSDK-3.0/lib/x86_64 -lOpenCL -I/home/hvictor/bin/AMDAPPSDK-3.0/include -I/home/hvictor/Work/group-talk-progs/vexcl -lboost_system -o $(@) $(<)
# 	#$(CXX) $(CXXFLAGS) -std=c++11 -DVexCL -L/home/hvictor/bin/AMDAPPSDK-3.0/lib/x86_64 -lOpenCL -I/home/hvictor/bin/AMDAPPSDK-3.0/include -I/home/hvictor/Work/group-talk-progs/vexcl -lboost_system -o $(@) $(<)
# 	$(NVCC) $(CXXFLAGS) -std=c++11 -DVexCL -lOpenCL -I /home/hvictor/Work/group-talk-progs/vexcl  -lboost_system -o $(@) $(<)
# 	@echo 'Finished building: $@'
# 	@echo ' '

links: $(SRCS)
	@echo 'Making links: $(<:%.cu=%.cpp)'
	ln -f -s $(<) $(<:%.cu=%.cpp)
	@echo 'Finished making links'
	@echo ' '

clean:
	@echo 'Removing objects an links'
	-$(RM) $(OBJS_PROMD) $(OBJS_OMP) $(OBJS_CUDA) $(OBJS_MPI)
	-@echo ' '

distclean: clean
	@echo 'Removing program'
	-$(RM) promd promd_omp promd_mpi promd_tbb promd_thrust_omp promd_thrust_tbb promd_thrust_cuda promd_cuda
	-@echo ' '

.PHONY : clean
