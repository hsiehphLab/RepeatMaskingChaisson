all: comask hardmask bemask mask2bed nl


hardmask: HardMask.cpp
	g++ -O2 $< -o $@ -I /projects/standard/hsiehph/shared/software/packages/htslib/install/include -I $(CONDA_PREFIX)/include -L/projects/standard/hsiehph/shared/software/packages/htslib/install/lib -L$(CONDA_PREFIX)/lib -lhts -Wl,-rpath,$(CONDA_PREFIX)/lib

comask: CombineMask.cpp
	g++ -O2  $< -o $@  -I /projects/standard/hsiehph/shared/software/packages/htslib/install/include -I $(CONDA_PREFIX)/include -L/projects/standard/hsiehph/shared/software/packages/htslib/install/lib -L$(CONDA_PREFIX)/lib -lhts -Wl,-rpath,$(CONDA_PREFIX)/lib  -lhts -lz -lpthread

bemask: MaskBed.cpp
	g++ -O2  $< -o $@  -I /projects/standard/hsiehph/shared/software/packages/htslib/install/include -L/projects/standard/hsiehph/shared/software/packages/htslib/install/lib -lhts  -lhts -lz -lpthread

mask2bed: MaskedToBed.cpp 
	g++ -O2  $< -o $@  -I /projects/standard/hsiehph/shared/software/packages/htslib/install/include -I $(CONDA_PREFIX)/include -L/projects/standard/hsiehph/shared/software/packages/htslib/install/lib -L$(CONDA_PREFIX)/lib -lhts -Wl,-rpath,$(PWD)/htslib/lib  -lhts -lz -lpthread

countn: CountN.cpp
	g++ -O2 CountN.cpp -o countn

nl: CountRep.cpp
	g++ -O2 CountRep.cpp -o nl

bamToFreq: BamToFreq.cpp $(CONDA_PREFIX)/lib/libhts.so
	g++ -O2 $< -o $@ -I /projects/standard/hsiehph/shared/software/packages/htslib/install/include -I $(CONDA_PREFIX)/include -L/projects/standard/hsiehph/shared/software/packages/htslib/install/lib -L $(CONDA_PREFIX)/lib -lhts -lpthread -lz -Wl,-rpath,$(CONDA_PREFIX)/lib

