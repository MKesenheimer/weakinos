"""
This module defines the class InteractionsParser. The InteractionsParser class
implements a interactions.dat file parser required by the various part of the
MG/ME project which needs to access this file.
 
Diagnostic usage: python interactions_parser.py [options]

Options:
  -f ...,--file=...       use the specified file as an input instead of
                          the local interactions.dat
  -p ...,--part=...       use the specified file as an input instead of
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

from hepobjects import Particle,Interaction
from particles_parser import ParticlesParser

###############################################################################
##
##  CLASSES
##
###############################################################################

class InteractionsParser(list):
    """The InteractionsParser class implements an interactions.dat file parser."""
    
    ref_part_parser = ParticlesParser()
    
    class InteractionsParserError(Exception):
        """Exception raised if an error occurs during interactions.dat parsing."""
        pass
    
    def __init__(self):
        """Initialize an empty InteractionsParser object."""
        list.__init__(self)
    
    def associate_part_parser(self,part_parser):
        """Attach the current InteractionParser to a given ParticleParser, needed
        to parse the interaction.dat file in a smart way"""
        
        if not isinstance(part_parser,ParticlesParser):
            raise self.InteractionsParserError, \
                        "Unvalid reference ParticlesParser object!"
        
        self.ref_part_parser = part_parser
        
        
    def read(self,file):
        """Read the content of file, parse it and add it to the current object"""
        
        try:
            fsock = open(file,"r")
        except IOError,(errno, strerror):
            print "I/O error(%s): %s" % (errno, strerror)
        else:   
            for line in fsock:
                myinter = Interaction()
                
                line = line.split("#",2)[0] # remove any comment
                line = line.strip() # makes the string clean
                # line = line.lower() # makes the string conv. lowercase
                
                if line != "":
                    
                    # Create the list of particles
                    values = line.split()
                    part_list = []

                    try:
                        for str_name in values:
                            curr_part = self.ref_part_parser.catch(str_name.lower())
                            if isinstance(curr_part,Particle):
                                # Look at the total number of strings, stop if anyway not enough
                                # required if a variable name corresponds to a particle! (eg G)
                                if len(values) >= 2*len(part_list)+1:
                                    part_list.append(curr_part)
                                else: break
                            # also stops if string does not correspond to a particle name
                            else: break
                                

                    
                        if len(part_list) < 3:
                            raise Interaction.InteractionError,\
                                "Vertex with less than 3 known particles found."
                        # Initialize the Interaction object
                        myinter.initialize(part_list)
                        # Use the other strings to fill variable names and tags
                        myinter.set("variables",values[myinter["npart"]:\
                                                       myinter["npart"]+myinter["num_var"]])
                        myinter.set("tags",values[myinter["npart"]+myinter["num_var"]:])
                        
                        if len(myinter["tags"]) >= len(part_list) :
                            raise Interaction.InteractionError,\
                                "Vertex with anomalous number of tags found."
                            
                        self.add(myinter)
                        
                    except Interaction.InteractionError,why:
                        print "Warning:",why," Interaction ignored"
                        
                    # print myinter
        fsock.close()

    def add(self,inter):
        """Add an Interaction object to the InteractionsParser item"""
    
        if not isinstance(inter,Interaction):
            raise ValueError # Trying to compare apple and ...
        
        try:
            # Check if the interaction is not already there, ordering IS important (cfr mssm)
            for curr_inter in self:
                if inter["particles"] == curr_inter["particles"]:
                    raise self.InteractionsParserError,\
                    "Interaction %s is already there!" \
                    % ",".join([part["name"][part["kind"]] for part in inter["particles"]])
                
            self.append(inter)
            
        except self.InteractionsParserError,why:
                    print "Warning:",why," Interaction ignored"
            
    def compare(self,ref_interactionsparser):
        """Compare the current parser with another one. Reports new, modified and
        removed interactions as a dictionnary of three parser objects."""
        
        new_inters = InteractionsParser()
        mod_inters = InteractionsParser()
        rm_inters = InteractionsParser()
        
        #if self.ref_part_parser != ref_interactionsparser.ref_part_parser:
        #    raise self.InteractionsParserError,\
        #            "Reference ParticlesParser objects in compare are not the same!"
        
        for new_inter in self:
            if not self.ref_part_parser._gen_or([(new_inter == old_inter) \
                                 for old_inter in ref_interactionsparser]):
                if not self.ref_part_parser._gen_or([(new_inter["particles"] == old_inter["particles"]) \
                                     for old_inter in ref_interactionsparser]):
                    new_inters.add(new_inter)
                else:
                    mod_inters.add(new_inter)
        
        for old_inter in ref_interactionsparser:
            if not self.ref_part_parser._gen_or([(new_inter["particles"] == old_inter["particles"]) \
                           for new_inter in self]):
                rm_inters.add(old_inter)
        
        return {"new":new_inters,"modified":mod_inters,"removed":rm_inters}

    def to_string(self):
        """Return a string containing the interactions.dat file content corresponding to the
        Parser content."""
        
        file_str = ""
        
        for inter in self:
            part_str = ""
            for part in inter["particles"]:
                if part["kind"] == "part":
                    part_str = part_str + "%s   " % part["name"]["part"]
                elif part["kind"] == "apart":
                    part_str = part_str + "%s   " % part["name"]["apart"]
            
            file_str = file_str + part_str + "  ".join(inter["variables"]) \
                    + "  " + "  ".join(inter["tags"]) + '\n'
        
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
    
    file = "interactions.dat"
    part_file = "particles.dat"
    
    try:
        options, discard = getopt.getopt(argv,"f:p:h",\
                                         ["file=","part=","help"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
    for option, value in options:
        if option in ("-h", "--help"):
            usage()
            sys.exit()                       
        elif option in ("-f", "--file"):
            file = value
        elif option in ("-p", "--part"):
            part_file = value
    
    mypp = ParticlesParser()
    mypp.read(part_file)
        
    myip = InteractionsParser()
    myip.associate_part_parser(mypp)
    myip.read(file)
    
    for myinter in myip:
        print "Interaction "+",".join([part["name"][part["kind"]] for part in myinter["particles"]])\
        +" with type "+myinter["type"]+", variables "+repr(myinter["variables"])+" and tags "\
        +repr(myinter["tags"])+" found."
    
if __name__=="__main__":
    main(sys.argv[1:])