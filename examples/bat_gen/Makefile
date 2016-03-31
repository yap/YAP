###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
###################################################################

OBJSDIR = obj
OUTPUTDIR = output

CXXSRCS := bat_gen.cxx
PRGSRCS := runBatGen.cxx

# compiler and flags
CXX       = g++
CXXFLAGS  = -std=c++11 -g -O2 -Wall -fPIC -Wno-deprecated  -Wno-missing-braces #-fopenmp
LD        = /usr/bin/ld -m elf_x86_64
LDFLAGS   = -g -O2  #-fopenmp

# ----------------------------------------------------------------------
# The following definitions rely on the script bat-config being
# available in $PATH. If BAT is not installed in the standard system
# directories, update $PATH accordingly.

CXXFLAGS += $(shell bat-config --cflags) -I$(YAPDIR)/include -I$(ROOTSYS)/include
LIBS := $(shell bat-config --libs) $(YAPDIR)/build/src/libYAP.a $(ROOTSYS)/lib/libPhysics.so

CXXOBJS := $(addprefix $(OBJSDIR)/,$(addsuffix .o,$(basename $(CXXSRCS))))
PRGOBJS := $(addprefix $(OBJSDIR)/,$(addsuffix .o,$(basename $(PRGSRCS))))
MYPROGS = $(basename $(PRGSRCS))

VPATH = $(OBJSDIR)

GARBAGE = $(CXXOBJS) $(PRGOBJS) link.d $(MYPROGS)

# targets
all : $(MYPROGS)

.PHONY : all clean print

$(OBJSDIR):
	@mkdir $(OBJSDIR)

$(OUTPUTDIR):
	@mkdir $(OUTPUTDIR)

link.d : $(addsuffix .h,$(basename $(CXXSRCS))) $(CXXSRCS) $(PRGSRCS)
	@echo " " > link.d
	@$(foreach src,$(filter-out %.h,$^), $(CXX) -MM -MT $(OBJSDIR)/$(notdir $(basename $(src))).o $(CXXFLAGS) $(src) >> link.d;)
	@$(foreach prog,$(MYPROGS), echo $(prog) : $(OBJSDIR)/$(prog).o >> link.d;)

-include link.d

$(CXXOBJS) $(PRGOBJS) : | $(OBJSDIR)
	$(CXX) $(CXXFLAGS) -c $(filter $(basename $(@F)).%,$(filter-out %.h,$^)) -o $@

$(MYPROGS) : $(CXXOBJS) | $(OUTPUTDIR)
	$(CXX) $(LDFLAGS) $^ $(LIBS) -o $@

clean :
	rm -f $(GARBAGE)
	rmdir $(OBJSDIR)
	rmdir --ignore-fail-on-non-empty $(OUTPUTDIR)

print :
	@echo compiler  : $(CXX)
	@echo c++ srcs  : $(CXXSRCS) $(PRGSRCS)
	@echo c++ objs  : $(CXXOBJS) $(PRGOBJS)
	@echo c++ flags : $(CXXFLAGS)
	@echo ld flags  : $(LDFLAGS)
	@echo libs      : $(LIBS)
