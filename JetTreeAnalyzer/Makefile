ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

MKDIR_BIN=mkdir -p $(PWD)/bin

all: mkdirBin  bin/era_run3_2023_correction.exe

mkdirBin:
	$(MKDIR_BIN)

bin/era_run3_2023_correction.exe: src/era_run3_2023_correction.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/era_run3_2023_correction.exe src/era_run3_2023_correction.C

# Run 2 basic-analysis revisit; build explicitly with
# make bin/run2_default_data_vn_revisit.exe
bin/run2_default_data_vn_revisit.exe: src/run2_default_data_vn_revisit.C | mkdirBin
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o $@ $<


clean:
	rm -f $(PWD)/include/*~
	rm -f $(PWD)/src/*~
	rm -f $(PWD)/src/*.so
	rm -f $(PWD)/src/*.d
	rm -f $(PWD)/src/*.pcm
	rm -f $(PWD)/bin/*.exe
	rmdir bin
