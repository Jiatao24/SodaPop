RM =rm
CP =cp
CC =g++

IDIR =./include
CXXFLAGS =-std=c++11 -Wall -O2 -I$(IDIR)
# CXXFLAGS =-std=c++11 -Wall -g -I$(IDIR)
LINK = $(CXX) $(CXXFLAGS)
COMPILE = $(CXX) $(LIBS) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -c

SODAPOP = sodapop
SNAP2ASCII = sodasnap
SUMM2SNAP = sodasumm

INSTALLDIR = $(HOME)/local/bin

all: $(SODAPOP) $(SNAP2ASCII) $(SUMM2SNAP)
install:
	@echo \#\#\# Installing binaries to $(INSTALLDIR)/...
	$(CP) $(SODAPOP) $(INSTALLDIR)/
	$(CP) $(SNAP2ASCII) $(INSTALLDIR)/
	$(CP) $(SUMM2SNAP) $(INSTALLDIR)/

uninstall:
	@echo \#\#\# Uninstalling binaries from $(INSTALLDIR)/...
	$(RM) -r $(INSTALLDIR)/$(SODAPOP)
	$(RM) -r $(INSTALLDIR)/$(SNAP2ASCII)
	$(RM) -r $(INSTALLDIR)/$(SUMM2SNAP)

$(SODAPOP): sodapop.o rng.o Cell.o
	$(LINK) -o sodapop sodapop.o rng.o Cell.o
$(SNAP2ASCII): snap2ascii.o rng.o Cell.o
	$(LINK) -o sodasnap snap2ascii.o rng.o Cell.o
$(SUMM2SNAP): summ2snap.o rng.o Cell.o
	$(LINK) -o sodasumm summ2snap.o rng.o Cell.o

rng.o: ./src/rng.cpp
	$(COMPILE) -o rng.o ./src/rng.cpp
Cell.o: ./src/Cell.cpp
	$(COMPILE) ./src/Cell.cpp
sodapop.o: ./src/evolve.cpp ./src/Gene.h ./src/Cell.h ./src/global.h
	$(COMPILE) -o sodapop.o ./src/evolve.cpp
snap2ascii.o: ./tools/snap2ascii.cpp ./src/global.h
	$(COMPILE) -o snap2ascii.o ./tools/snap2ascii.cpp
summ2snap.o: ./tools/summ2snap.cpp ./src/PolyCell.h
	$(COMPILE) -o summ2snap.o ./tools/summ2snap.cpp

clean:
	rm -f *.o
