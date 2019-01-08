#object_files (testEffi.o)
OBJECTS := $(wildcard *.o)

#root_stuff (root libraries and needed root options)
ROOTLIBS  := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags --libs) -lRooFit -lRooFitCore -lMathMore
ROOTCINT  := $(shell which rootcint)

#directories
SOURCEDIR  := ./src
INCLUDEDIR := ./interface

#exe_files
EXECUTABLE0 := fit_recoMC_projEff
EXECUTABLE1 := fit_genMC
CLASS1      := AngularRT
CLASS2      := AngularWT
CLASS3      := PdfRT
CLASS4      := PdfWT
CLASSDICT   := AngularDict.cc
CLASSDICT1  := AngDict.cc

#compiling options
DEBUGFLAGS := -O3 -Wall -std=c++11
CXXFLAGS := $(DEBUGFLAGS) 

#compile class
LIBS := $(SOURCEDIR)/$(CLASS1).cc $(SOURCEDIR)/$(CLASS2).cc $(SOURCEDIR)/$(CLASS3).cc $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT1)

$(CLASSDICT1): $(INCLUDEDIR)/$(CLASS1).h $(INCLUDEDIR)/$(CLASS2).h $(INCLUDEDIR)/$(CLASS3).h $(INCLUDEDIR)/$(CLASS4).h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@ -c $^

$(EXECUTABLE0): $(EXECUTABLE0).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE1): $(EXECUTABLE1).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)


#cleaning options
.PHONY: clean
clean:
	rm -f $(OBJECTS) && rm -f $(EXECUTABLE0) $(EXECUTABLE1)
