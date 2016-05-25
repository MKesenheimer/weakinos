## POWHEG-BOX-V2 Project Weakino Pair Production

### Synopsis

More information on how to set up the project properly can be found at http://powhegbox.mib.infn.it/.

### Compiling
To get started, compile the static libraries

* libdhelas3.a

* liblooptools.a

* libSLHA.a

for the operating system you are using. To this end, call the configuration scripts in the main directory by simply typing

        $ ./configure.sh

Afterwards you can compile the libraries by typing

        $ make libdhelas3.a
        $ make liblooptools.a
        $ make libSLHA.a
        
or, in short, 

        $ make libs

If you want to use your own libraries, copy them into ./Tools/ or provide
paths to the libraries in the Makefile.

Afterwards, change into a desired process directory and type

        $ make clean-results && make -j4 do
        
to compile and run the program. Use the flag -jn, where n is an integer, to 
compile on multiple cores. This will speed up the compilation.

Important note for Mac OSX and some Linux users:
In order to link the object files properly with newer compiler versions
it might be advisable to recompile all own libraries using the -lstdc++ flag.


### Precompiler Flags
In the current code version several C preprocessor (cpp) flags are implemented.
The preprocessor runs in traditional mode for gfortran. Any restrictions of the 
file format, especially the limits on line length, apply for 
preprocessed output as well, so it might be advisable to use the 

        -ffree-line-length-none 
or 

        -ffixed-line-length-none

options (activated as default). If you want to change a preprocessor flag
it is imperative to run

        $ make clean

before recompiling the source code.
The flags

* FORM_BORN, MAD_BORN

* DR_I, DR_II, DSUB_I, DSUB_II

are mandatory, you should not change them.

Please refer to the Makefile for more details. 


### Running

The following commands are all executed in the desired process directory.

        $ make do

compiles the source and runs the program in ./testrun.

        $ make clean

removes all object files in ./build. This has no effect on the compiled program.

        $ make clean-results

removes the results in ./testrun.

        $ make clean-all

removes the results, the object files and the compiled programs.

To remove all object files and to clean the libraries type

        $ make clean

in the main directory.

All parameters are read from a single slha-file in ./testrun. Runtime variables, such as 
integration points, number of events to generate, etc. have to be specified in powheg.input.
If you want to change the Z-mass, Z-width or alpha_em you can do this in powheg.input, too.
Please refer to the provided manual.

### Scripts

We have added several helpful scripts to the whole package, which can be used to generate 
results or clean old runs. The most important one is ./Scripts/runparallel.sh, which is used 
to run the POWHEG-BOX-V2 executables fully automated in parallel mode. 
Type

        $ ./runparallel.sh -h

to get an overview of the functionality of the script. This script works even with the MOAB 
submitting system msub.


### References

When using the code, please cite the following publications:

J.Baglio, B.Jager, M.Kesenheimer, "Electroweakino pair production at the LHC: NLO SUSY-QCD corrections and parton-shower effects", arXiv:1605.06509 [hep-ph].

P.Nason, "A New method for combining NLO QCD with shower Monte Carlo algorithms", JHEP 0411 (2004) 040, hep-ph/0409146.  

S.Frixione, P.Nason and C.Oleari, "Matching NLO QCD computations with Parton Shower simulations: the POWHEG method", JHEP 0711 (2007) 070, arXiv:0709.2092. 

S.Alioli, P.Nason, C.Oleari and E.Re, "A general framework for implementing NLO calculations in shower Monte Carlo programs: the POWHEG BOX", JHEP 1006 (2010) 043, arXiv:1002.2581.

