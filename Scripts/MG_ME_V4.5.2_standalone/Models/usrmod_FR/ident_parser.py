"""
This module defines the class IdentParser. The IdentParser class
implements a ident_card.dat file format parser required
by the various part of the MG/ME project which needs
to access this kind of file.
 
Diagnostic usage: python ident_parser.py [options]

Options:
  -f ...,--file=...       use the specified file as an input instead of
                          the local ident_card.dat
  -h,--help               display this help

This file is part of MadGraph/MadEvent, a free software to generate
matrix element amplitudes and parton level events.
Visit http://madgraph.phys.ucl.ac.be. 
"""

__author__="Michel Herquet, the MadTeam, last change by $Author: madgraph $"
__version__= "$Revision: 1.2 $"
__date__= "$Date: 2009/08/22 08:01:09 $"
__copyright__= "Copyright (c) 2008 MadGraph/MadEvent team"
__license__= "GNU"

###############################################################################
##
##  IMPORT & GLOBAL
##
###############################################################################

import sys
import getopt
from particles_parser import ParticlesParser

###############################################################################
##
##  CLASSES
##
###############################################################################

class IdentParser(list):
    """Implements a simple parser for the ident_card.dat file format"""
    
    ref_part_parser = ParticlesParser()
    
    class IdentParserError(Exception):
        """Exception raised if an error occurs during the parsinf of
        ident_card.dat"""
        pass
    
    def __init__(self):
        """Initialize an empty IdentParser object."""
        list.__init__(self)
    
    def read(self,file):
        """Read the content of file, parse it and add it to the current object"""

        try:
            fsock = open(file,"r")
        except IOError,(errno, strerror):
            print "I/O error(%s): %s" % (errno, strerror)
        else:
            for line in fsock:
                line = line.split("*")[0] # split comments
                line = line.strip() # makes the string clean
                line = line.lower() # makes the string conv. lowercase
        
                if line != "":
                    param_line={}
                    values = line.split()
                    try:
                        param_line["block"] = values[0]
                        param_line["var"] = values[-1]
                        param_line["ids"] = [int(x) for x in values[1:-1]]
                        
                        # If the param is a mass or a decay width, check validity against known particle list
                        if param_line["block"] == "mass":
                            found = False
                            for part in self.ref_part_parser:
                                if (part["pdg_code"] == param_line["ids"][0]) \
                                    and (part["mass"] == param_line["var"]):
                                    found = True
                            if not found:
                                raise self.IdentParserError,\
                                "Mass parameter %s has no associated particle or wrong PDG code." \
                                % param_line["var"]
                        
                        if param_line["block"] == "decay":
                            found = False
                            for part in self.ref_part_parser:
                                if (part["pdg_code"] == param_line["ids"][0]) \
                                    and (part["width"] == param_line["var"]):
                                    found = True
                            if not found:
                                raise self.IdentParserError,\
                                "Width parameter %s has no associated particle or wrong PDG code." \
                                % param_line["var"]
                
                        self.append(param_line)
                        
                    except self.IdentParserError,why:
                            print "Warning:",why," Parameter ignored"
                            
            # Check if all mass and decay parameters appears in the current parser
            var_list = [param["var"] for param in self]
            var_list.append("zero")
            
            for part in self.ref_part_parser:
                myparam = {}
                if part["mass"] not in var_list:
                    myparam["block"] = "mass"
                    myparam["ids"] = [part["pdg_code"]]
                    myparam["var"] = part["mass"]
                    print "Warning: Mass parameter %s of particle %s does not exist. Parameter added." \
                                % (part["mass"],part["name"]["part"])
                    self.append(myparam)
                    myparam = {}
                    
                if part["width"] not in var_list:
                    myparam["block"] = "decay"
                    myparam["ids"] = [part["pdg_code"]]
                    myparam["var"] = part["width"]
                    print "Warning: Width parameter %s of particle %s does not exist. Parameter added." \
                                % (part["width"],part["name"]["part"])
                    self.append(myparam)
                    myparam = {}
                
            fsock.close()
            
    def to_string(self):
        """Return a string containing the ident_card.dat file content corresponding to the
        Parser content."""
        
        file_str = ""
        
        for param in self:
            file_str = file_str + "%s %s %s\n" % (param["block"],\
                                                " ".join([repr(x) for x in param["ids"]]),\
                                                param["var"])
        
        return file_str
    
    def compare(self,ref_idp):
        """Compare the current IdentParser with a reference one. Reports new, modified and
        removed parameters as a dictionnary of three parser objects."""
        
        new_params = IdentParser()
        rm_params = IdentParser()
        mod_params = IdentParser()
        
        for new_param in self:
            same_block = filter(lambda x: x["block"] == new_param["block"], ref_idp)
            if not self._gen_or([(new_param == old_param) for old_param in same_block]):
                if not self._gen_or([(new_param["ids"] == old_param["ids"]) \
                                     for old_param in same_block]):
                    new_params.append(new_param)
                else:
                    mod_params.append(new_param)
        
        for old_param in ref_idp:
            same_block = filter(lambda x: x["block"] == old_param["block"], self)
            if not self._gen_or([(new_param["ids"] == old_param["ids"]) \
                           for new_param in same_block]):
                rm_params.append(old_param)
        
        return {"new":new_params,"modified":mod_params,"removed":rm_params}
        
    
    def _gen_or(self,bool_list):
        "OR operator generalized to a list, i.e. True if at least one True."
        for bool in bool_list:
            if bool: return True
        return False
    
    def associate_part_parser(self,part_parser):
        """Attach the current IdentParser to a given ParticleParser, needed
        to parse the mass and decay info file in a smart way"""
        
        if not isinstance(part_parser,ParticlesParser):
            raise self.InteractionsParserError, \
                        "Unvalid reference ParticlesParser object!"
        
        self.ref_part_parser = part_parser
    
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
    
    file = "ident_card.dat"
    part_file = "particles.dat"
    try:
        options, discard = getopt.getopt(argv,"f:hp:",\
                                         ["file=","help","part="])
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
    myidp = IdentParser()
    myidp.associate_part_parser(mypp)
    myidp.read(file)

    for myparam in myidp:
        print "Parameter: "+",".join([repr(x) for x in myparam["ids"]])\
        +", block: "+myparam["block"]+", variable: " + myparam["var"]
    
if __name__=="__main__":
    main(sys.argv[1:])
