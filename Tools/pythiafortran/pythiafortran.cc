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
#include "Pythia8/Pythia.h"

//using namespace Pythia8;
namespace py = Pythia8;
//using namespace std;

extern "C"
{   
// f77 interface to 
//
// Corresponds to the following Fortran subroutine
// interface structure:
//
//   SUBROUTINE FASTJETPPSEQREC(P,NPART,R,PALG,F77JETS,NJETS)
//   DOUBLE PRECISION P(4,*), R, PALG, F, F77JETS(4,*)
//   INTEGER          NPART, NJETS
// 
// where on input
//
//   P        the input particle 4-momenta
//   NPART    the number of input momenta
//   R        the radius parameter
//   PALG     the power for the generalised kt alg 
//            (1.0=kt, 0.0=C/A,  -1.0 = anti-kt)
//
// and on output 
//
//   F77JETS  the output jet momenta (whose second dim should be >= NPART)
//            sorted in order of decreasing p_t.
//   NJETS    the number of output jets 
//.
//
//void pythia_init_(const double * p, const int & npart,                   
//                 const double & R, const double & palg, const double & ptmin,
//                 double * f77jets, int & njets, int * f77jetvec) {
void pythia_init_(const int & arg) {
        
        int a = arg;
    
        //Pythia test
        py::Pythia pythia;
    
        std::cout<<"Hello World! a = "<<a<<std::endl;
   
}
}

