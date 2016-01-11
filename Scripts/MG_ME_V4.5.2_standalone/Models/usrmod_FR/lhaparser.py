"""
This module defines the class LHAParser. The LHAParser class
implements a Les Houches Accord file format parser required
by the various part of the MG/ME project which needs
to access this kind of file.
 
Diagnostic usage: python lhaparser.py [options]

Options:
  -f ...,--file=...       use the specified file as an input instead of
                          the local param_card.dat
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


###############################################################################
##
##  CLASSES
##
###############################################################################

class LHAParameter(dict):
    """The LHAParameter class represents an abstract LHA parameter object with a block name, a
    list of integer id numbers, a numerical parameter value (stored as a string), a possible
    associated comment string and a type format, either standard, decay or br. For br type
    parameters, the block name corresponds to the pdg code of the mother particle."""
    
    class LHAParameterError(Exception):
        """Exception raised if an error occurs in the definition
        of a parameter"""
        pass
    
    def __init__(self,block="",ids=[],value=0.0,comment="",type="standard"):
        """Create a Parameter object. Empty by default."""
        dict.__init__(self)
        
        self.set("block",block)
        self.set("ids",ids)
        self.set("value",value)
        self.set("comment",comment)
        self.set("type",type)
    
    def set(self,prop_name,value):
        """Set the value of prop_name to value, with relevant checks."""
        
        if prop_name == "block":
            if not isinstance(value,str):
                raise self.LHAParameterError, repr(value)+\
                " is not a valid block string"
            self["block"] = value
            
        elif prop_name == "ids":
            if not isinstance(value,list):
                raise self.LHAParameterError, repr(value)+" is not a valid list of id values!"
            for id in value:
                try:
                    dummy = int(id)
                except:
                    raise self.LHAParameterError, \
                    repr(value)+" at position "+repr(values.index(id))\
                    +" has not an integer value!"
            self["ids"] = [int(x) for x in value]
            
        elif prop_name == "value":
            if self["block"] in ["decay","br"]:
                try:
                        dummy = float(value)
                except:
                    raise self.LHAParameterError, \
                    repr(value)+" is not a correct float (or int) value for parameter "\
                    +repr(self["ids"])+" of block"+repr(repr(self["block"]))
            self["value"] = value
            
        elif prop_name == "comment":
            if not isinstance(value,str):
                raise self.LHAParameterError, repr(value)+\
                " is not a valid comment string for parameter "+repr(self["ids"])\
                +" of block"+repr(repr(self["block"]))
            self["comment"] = value
            
        elif prop_name == "type":
            if value not in ["standard","decay","br","qnumbers"]:
                raise self.LHAParameterError, repr(value)+\
                " is not a valid type format for parameter "+repr(self["ids"])\
                +" of block"+repr(repr(self["block"]))+". Should be either standard, decay or br."
            self["type"] = value
            
        else:
            raise self.LHAParameterError, repr(prop_name)+\
                " is not a valid property name for parameter "+repr(self["ids"])\
                +" of block"+repr(repr(self["block"]))

        

class LHAParser(list):
    """The LHAParser class implements a LH accord file format parser."""
    
    class LHAParserError(Exception):
        """Exception raised if an error occurs during LHA file parsing."""
        pass
    
    def __init__(self):
        """Initialize an empty LHAParser object."""
        list.__init__(self)
    
    def read(self,file):
        """Read the content of file, parse it and add it to the current object"""
        try:
            fsock = open(file,"r")
        except IOError,(errno, strerror):
            print "I/O error(%s): %s" % (errno, strerror)
        else:
            curr_block = ""
            curr_format = "standard"
            comment = ""
            
            for line in fsock:
                split_line = line.split("#") # split comments
                
                line = split_line[0]
                line = line.strip() # makes the string clean
                line = line.lower() # makes the string conv. lowercase
                
                if len(split_line)>1:
                    comment = split_line[1]
                    comment = comment.strip() # makes the comment string clean

                
                if line != "":
                    values = line.split()
    
                    
                    if values[0] == "block":
                        if values[1] == "qnumbers":
                            curr_block = values[2]
                            curr_format = "qnumbers"
                            continue
                        else:
                            curr_block = values[1]
                            curr_format = "standard"
                            continue
                    
                    elif values[0] == "decay":
                        curr_block = "decay"
                        curr_format = "decay"

                    
                    if curr_format == "standard":
                        myparam = LHAParameter(curr_block,values[0:-1],values[-1],comment,curr_format)
                    elif curr_format == "decay":
                        myparam = LHAParameter(curr_block,values[1:-1],values[-1],comment,curr_format)
                    elif curr_format == "br":
                        myparam = LHAParameter(curr_block,values[1:],values[0],comment,curr_format)
                    elif curr_format == "qnumbers":
                        myparam = LHAParameter(curr_block,values[0:-1],values[-1],comment,curr_format)
                    
                    if curr_format == "decay":
                        curr_format = "br"
                        curr_block = repr(myparam["ids"][0])
                        
                    try:
                        self.add(myparam)
                    except self.LHAParserError,why:
                        print "Warning:",why,", parameter ignored"
                    
            fsock.close()
    
    def add(self,param):
        """Add a LHAParameter object to the LHAParser item"""
    
        if not isinstance(param,LHAParameter):
            raise ValueError # Trying to compare apple and ...
        
        for curr_param in filter(lambda x: x["block"] == param["block"], self):
            if curr_param["ids"] == param["ids"]:
                raise self.LHAParserError, "Parameter "+repr(param["ids"])\
                +" of block "+param["block"]+" is already there!"
        
        #if param["type"] == "br":
        #    total_br = 0.0
        #    for curr_param in filter(lambda x: x["block"] == param["block"], self):
        #        total_br = total_br + float(curr_param["value"])
        #    total_br = total_br + float(param["value"])
        #    if total_br - 1.0 > 1e-08:
        #        raise self.LHAParserError, "Total BR for "+param["block"]+\
        #        " is larger than 1.0 because of parameter "+repr(param["ids"])
            
        self.append(param)
    
    def compare(self,ref_lhap):
        """Compare the current LHAParser with a reference one. Reports new, modified and
        removed parameters as a dictionnary of three parser objects."""
        
        new_params = LHAParser()
        rm_params = LHAParser()
        mod_params = LHAParser()
        
        for new_param in self:
            same_block = filter(lambda x: x["block"] == new_param["block"], ref_lhap)
            if not self._gen_or([(new_param == old_param) for old_param in same_block]):
                if not self._gen_or([(new_param["ids"] == old_param["ids"]) \
                                     for old_param in same_block]):
                    new_params.add(new_param)
                else:
                    mod_params.add(new_param)
        
        for old_param in ref_lhap:
            same_block = filter(lambda x: x["block"] == old_param["block"], self)
            if not self._gen_or([(new_param["ids"] == old_param["ids"]) \
                           for new_param in same_block]):
                rm_params.add(old_param)
        
        return {"new":new_params,"modified":mod_params,"removed":rm_params}
        
    
    def _gen_or(self,bool_list):
        "OR operator generalized to a list, i.e. True if at least one True."
        for bool in bool_list:
            if bool: return True
        return False
    
    def to_string(self):
        """Return a string containing the lha file corresponding to the
        Parser content."""
        
        file_str = ""
        curr_block = ""
        first_qnum = True
        
        for param in self:
            
            if param["type"] == "standard":
                if curr_block != param["block"]:
                    file_str = file_str + "Block %s" % param["block"].upper() + '\n'
                    curr_block = param["block"]
                file_str = file_str + "    %s    %s  # %s\n" \
                           % ("  ".join([repr(s) for s in param["ids"]]),\
                              param["value"],param["comment"])
            
            if param["type"] == "decay":
                curr_block = repr(param["ids"][0])
                file_str = file_str + "Decay    %s    %s  # %s\n" \
                           % ("  ".join([repr(s) for s in param["ids"]]),\
                              param["value"],param["comment"])
            
            if param["type"] == "br":
                if curr_block != param["block"]:
                    raise self.LHAParserError, "BR "+repr(param["ids"])\
                            +" of block "+param["block"]+" is not in the right decay block!"
                file_str = file_str + "    %s    %s  # %s\n" \
                           % ("  ".join([repr(s) for s in param["ids"]]),\
                              param["value"],param["comment"])
            
            if param["type"] == "qnumbers":
                if curr_block != param["block"]:
                    first_qnum = True
                if first_qnum:
                    file_str = file_str + "Block QNUMBERS %s" % param["block"] + '\n'
                    first_qnum = False
                    curr_block = param["block"]
                file_str = file_str + "    %s    %s  # %s\n" \
                           % ("  ".join([repr(s) for s in param["ids"]]),\
                              param["value"],param["comment"])
                    
            
        
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
    
    file = "param_card.dat"
    
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
    
    mylhap = LHAParser()
    mylhap.read(file)

    for myparam in mylhap:
        print "Parameter: "+",".join([repr(x) for x in myparam["ids"]])\
        +", block: "+myparam["block"]+",value: "+myparam["value"]+", comment: "\
        +myparam["comment"]+", type: "+myparam["type"]
    
    
if __name__=="__main__":
    main(sys.argv[1:])
