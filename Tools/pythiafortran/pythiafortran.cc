//STARTHEADER
// $Id: pythiafortran.cc 2015-02-11 07:55:00MEZ kesenheimer $
//
// Copyright (c) M. Kesenheimer
//
//----------------------------------------------------------------------
//
//----------------------------------------------------------------------
//ENDHEADER

#include <iostream>
#include <stdio.h>
#include <string>
#include "Pythia8/Pythia.h"

//using namespace Pythia8;
namespace py = Pythia8;
//using namespace std;

extern "C" {

// global pythia objects (todo: use singleton class)
py::Pythia pythia;
py::Event& event = pythia.event;

void pythia_read_slha_(const char * str) {
        std::cout<<"input file: "<<str<<std::endl;
        std::string strfilename = "SLHA:file = ";
        strfilename.append(str);
        pythia.readString(strfilename);
}

void pythia_init_() {
        // todo
        //pythia.readString("Beams:eCM = 8000.");
        //pythia.readString("HardQCD:all = on");
        //pythia.readString("PhaseSpace:pTHatMin = 20.");
        pythia.init();
}

void pythia_event_() {
    for (int iEvent = 0; iEvent < 100; ++iEvent) {
        if (!pythia.next()) continue;
        //do something here:
    }
}

//end exern "C"
}

