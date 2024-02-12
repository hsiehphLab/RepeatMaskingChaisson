all: comask hardmask bemask mask2bed nl


hardmask: HardMask.cpp
	g++ -O2 $< -o $@ -I /panfs/jay/groups/7/hsiehph/gordo893/packages/htslib -I $(CONDA_PREFIX)/include -L/home/hsiehph/gordo893/packages/htslib -L$(CONDA_PREFIX)/lib -lhts -Wl,-rpath,$(CONDA_PREFIX)/lib

comask: CombineMask.cpp
	g++ -O2  $< -o $@  -I /panfs/jay/groups/7/hsiehph/gordo893/packages/htslib -I $(CONDA_PREFIX)/include -L/home/hsiehph/gordo893/packages/htslib -L$(CONDA_PREFIX)/lib -lhts -Wl,-rpath,$(CONDA_PREFIX)/lib  -lhts -lz -lpthread

bemask: MaskBed.cpp
	g++ -O2  $< -o $@  -I /panfs/jay/groups/7/hsiehph/gordo893/packages/htslib -I $(CONDA_PREFIX)/include -L/home/hsiehph/gordo893/packages/htslib -L$(CONDA_PREFIX)/lib -lhts -Wl,-rpath,$(CONDA_PREFIX)/htslib/lib  -lhts -lz -lpthread

mask2bed: MaskedToBed.cpp 
	g++ -O2  $< -o $@  -I /panfs/jay/groups/7/hsiehph/gordo893/packages/htslib -I $(CONDA_PREFIX)/include -L/home/hsiehph/gordo893/packages/htslib -L$(CONDA_PREFIX)/lib -lhts -Wl,-rpath,$(PWD)/htslib/lib  -lhts -lz -lpthread

countn: CountN.cpp
	g++ -O2 CountN.cpp -o countn

nl: CountRep.cpp
	g++ -O2 CountRep.cpp -o nl

bamToFreq: BamToFreq.cpp $(CONDA_PREFIX)/lib/libhts.so
	g++ -O2 $< -o $@ -I /panfs/jay/groups/7/hsiehph/gordo893/packages/htslib -I $(CONDA_PREFIX)/include -L/home/hsiehph/gordo893/packages/htslib -L $(CONDA_PREFIX)/lib -lhts -lpthread -lz -Wl,-rpath,$(CONDA_PREFIX)/lib

