"""
This module defines the class ParticlesParser. The ParticlesParser class
implements a particles.dat file parser required by the various part of the
MG/ME project which needs to access this file.
 
Diagnostic usage: python particles_parser.py [options]

Options:
  -f ...,--file=...       use the specified file as an input instead of
                          the local particles.dat
  -h,--help               display this help

This file is part of MadGraph/MadEvent, a free software to generate
matrix element amplitudes and parton level events.
Visit http://madgraph.phys.ucl.ac.be. 
"""

__author__="Michel Herquet, the MadTeam, last change by $Author: madgraph $"
__version__= "$Revision: 1.1 $"
__date__= "$Date: 2009/08/14 09:54:28 $"
__copyright__= "Copyright (c) 2008 MadGraph/MadEvent team"
__license__= "GNU"

###############################################################################
##
##  IMPORT & GLOBAL
##
###############################################################################

import sys
import getopt
import copy

from hepobjects import Particle

###############################################################################
##
##  CLASSES
##
###############################################################################

class ParticlesParser(list):
    """The ParticlesParser class implements a particles.dat file parser."""
    
    class ParticlesParserError(Exception):
        """Exception raised if an error occurs during particles.dat pasring."""
        pass
    
    def __init__(self):
        """Initialize an empty ParticlesParser object."""
        list.__init__(self)
    
    def read(self,file):
        """Read the content of file, parse it and add it to the current object"""
        try:
            fsock = open(file,"r")
        except IOError,(errno, strerror):
            print "I/O error(%s): %s" % (errno, strerror)
        else:   
            for line in fsock:
                mypart = Particle()
                
                if line.find("# MULTIPARTICLES") != -1:
                    break # stop scanning if old MULTIPARTICLES tag found
                
                line = line.split("#",2)[0] # remove any comment
                line = line.strip() # makes the string clean
                line = line.lower() # makes the string conv. lowercase
                
                if line != "":
                    values=line.split()
                    if len(values) != 9:   # Not enough tags on the line
                        raise self.ParticlesParserError, \
                        "Unvalid initialization string:"+line
                    else:
                        try:
                            mypart.set("name",{"part":values[0],"apart":values[1]})
                            mypart.set("spin",values[2])
                            mypart.set("line_type",values[3])
                            mypart.set("mass",values[4])
                            mypart.set("width",values[5])
                            mypart.set("color_rep",values[6])
                            mypart.set("draw_name",values[7])
                            mypart.set("pdg_code",int(values[8]))
                        except Particle.ParticleError,why:
                                print "Warning:",why,", particle ignored"
                        else:
                            self.add(mypart)
        fsock.close()
    
    def add(self,mypart):
        """Add a particle object to the ParticlesParser item"""
        
        if not isinstance(mypart,Particle):
            raise ValueError # Trying to compare apple and ...
        
        try:
            # check if the a particle with same name or PDG code is not already there
            for curr_part in self:
                if mypart.same_name(curr_part):
                    raise self.ParticlesParserError,"Particle %s is already there!" \
                    % repr(mypart["name"])
                if mypart["pdg_code"] == curr_part["pdg_code"]:
                    raise self.ParticlesParserError,"Particle %s has the same PDG code as %s !" \
                    % (repr(mypart["name"]),repr(curr_part["name"]))
                
            self.append(mypart)
            
        except self.ParticlesParserError,why:
                    print "Warning:",why,", particle ignored"
                    
    def compare(self,ref_particlesparser):
        """Compare the current parser with another one. Reports new, modified and
        removed particles as a dictionnary of three parser objects."""
        
        new_parts = ParticlesParser()
        mod_parts = ParticlesParser()
        rm_parts = ParticlesParser()
        
        for new_part in self:
            if not self._gen_or([new_part.is_identical(old_part) \
                                 for old_part in ref_particlesparser]):
                if not self._gen_or([new_part.same_name(old_part) \
                                     for old_part in ref_particlesparser]):
                    new_parts.add(new_part)
                else:
                    mod_parts.add(new_part)
        
        for old_part in ref_particlesparser:
            if not self._gen_or([old_part.same_name(new_part) \
                           for new_part in self]):
                rm_parts.add(old_part)
        
        return {"new":new_parts,"modified":mod_parts,"removed":rm_parts}
        
    def _gen_or(self,bool_list):
        "OR operator generalized to a list, i.e. True if at least one True."
        for bool in bool_list:
            if bool: return True
        return False
    
    def catch(self,str_name):
        """Try to find a Particle object with name str_name. Return the corresponding
        Particle if possible (with a correct kind switch if str_name is the anti particle
        name) or None if not."""
          
        for part in self:
            if part["name"]["part"] == str_name:
                return_part = copy.deepcopy(part)
                return return_part
            if part["name"]["apart"] == str_name:
                return_part = copy.deepcopy(part)
                return_part.set("kind","apart")
                return return_part
        return None
    
    def to_string(self):
        """Return a string containing the particles.dat file content corresponding to the
        Parser content."""
        
        file_str = ""
        
        for part in self:
            file_str = file_str + \
                       "   ".join([part["name"]["part"],part["name"]["apart"],
                                  part["spin"],part["line_type"],part["mass"],part["width"],
                                  part["color_rep"],part["draw_name"],repr(part["pdg_code"]),'\n'])
        return file_str
                
###############################################################################
##
##  DIAGNOSTIC
##
###############################################################################

def usage():
    """Print the script usage documentation.
    """
    print __doc__
    
def main(argv):
    """Main function, called if script is ran from the command line.
    """
    
    file = "particles.dat"
    
    try:
        options, discard = getopt.getopt(argv,"f:h",\
                                         ["file=","help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
    for option, value in options:
        if option in ("-h", "--help"):
            usage()
            sys.exit()                       
        elif option in ("-f", "--file"):
            file = value
    
    mypp = ParticlesParser()
    mypp.read(file)

    for mypart in mypp:
        print "Particle %s read: %s" % (mypart["name"]["part"],repr(mypart))
    
if __name__=="__main__":
    main(sys.argv[1:])