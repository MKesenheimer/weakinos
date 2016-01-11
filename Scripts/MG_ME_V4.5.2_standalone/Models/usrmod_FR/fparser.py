"""
This module defines the class FortranParser and all the associated classes.
The FortranParser class implements a (partial) fortran file parser
required by the various part of the MG/ME project which needs
to access/write fortran files.
 
Diagnostic usage: python fparser.py [options]

Options:
  -f ...,--file=...       use the specified file as an input
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
import re

# FORTRAN PATTERNS (global to compile regexps only once)

# Comment
fpat_comment = re.compile(r'c\s*(?P<comment>.*)',re.IGNORECASE)

# Variable declaration
# Note: ?: is a group which is not stored in the result, ?P<...> is a named group
fpat_declaration = re.compile(r"""\s{6,} # starts with at least 6 spaces
             (?P<type>(?:character|logical|integer|real|complex|double\s*precision|double\s*complex) # type name
             (?:\*[0-9]+)?) # Optional star followed by numbers (also contained in the type name)
             \s+            # White spaces, at least one
             (?P<vars>.*)  # the rest
             """,re.IGNORECASE | re.VERBOSE)

# Common block declaration
fpat_common = re.compile(r"""\s{6,}  # starts at least with 6 spaces
             common\s*\/\s*(?P<name>\w*)\s*\/\s* # common/name/ with optionnal spaces
             (?P<vars>.*)  # the rest
             """,re.IGNORECASE | re.VERBOSE)

# Call routine
fpat_call = re.compile(r"""\s{6,}  # starts with at least 6 spaces
             call\s*(?P<name>\w*)\(      # call name( with optionnal spaces
             (?P<vars>.*)\)  # the rest)
             """,re.IGNORECASE | re.VERBOSE)

# Simple write statement
fpat_write = re.compile(r"""\s{6,}  # starts with at least 6 spaces
             write\s*(?P<format>\([^\)]*\))\s*      # write and format
             (?P<vars>.*)  # the rest
             """,re.IGNORECASE | re.VERBOSE)

# A simple variable expression (search all in a list)
# Any set of alpha num characters, opt. followed by (...) where ... is anything except )
fpat_variable = re.compile(r'\w+(?:\([^\)]*\))?',re.IGNORECASE)

# Same as before but also find numbers, "..." and '...'
fpat_var_string_num = re.compile(r"""[0-9]*(?:\.[0-9]*)(?:e(?:\+|\-)[0-9]*) # a number, with opt decimal and opt e...
                                |\"[^\"]*\" # a "" string
                                |\'[^\']*\' # a '' string
                                |\w+(?:\([^\)]*\))? # a variable
                                """,re.IGNORECASE | re.VERBOSE)

# Variable definition
fpat_definition = re.compile(r"""\s{6,} # starts with at least 6 spaces
             (?P<name>\w+(?:\([^\)]*\))?) # variable name
             \s*=\s*  # Optional white spaces and = sign
             (?P<expr>.*)  # the rest
             """,re.IGNORECASE | re.VERBOSE)

###############################################################################
##
##  CLASSES
##
###############################################################################

class FortranStatement(dict):
    """This class implements various kinds of Fortran statements"""
    
    class FortranStatement(Exception):
        """Exception raised if an error occurs during a fortran statement parsing."""
        pass
    
    def __init__(self,init_string=""):
        """Initialize a FortranStatement object."""
        
        dict.__init__(self)
        
        self["type"] = "unknown"
        self["init_string"] = init_string
        self["params"] = {}
    
        if self["init_string"] is not "":
            self.parse_init_string()
            
    def parse_init_string(self):
        """Initialize the type and params list according to the current init_string.
        Here comes the actual Fortran syntax parsing algorithm."""
        
        init_string = self["init_string"].split('!',2)[0]
        
        # Find comments
        m = fpat_comment.match(init_string)
        if m:
            self["type"] = "comment"
            self["params"]["comment_string"] = m.group('comment')
            return True
        
        # Find variable declarations
        m = fpat_declaration.match(init_string)
        if m is not None:
            self["type"] = "declaration"
            self["params"]["type"] = m.group('type').lower()
            self["params"]["vars"] = fpat_variable.findall(m.group('vars').lower())
            return True
        
        # Find common block declarations
        m = fpat_common.match(init_string)
        if m is not None:
            self["type"] = "common"
            self["params"]["name"] = m.group('name').lower()
            self["params"]["vars"] = fpat_variable.findall(m.group('vars').lower())
            return True
        
        # Find calls
        m = fpat_call.match(init_string)
        if m is not None:
            self["type"] = "call"
            self["params"]["func_name"] = m.group('name').lower()
            self["params"]["vars"] = fpat_var_string_num.findall(m.group('vars').lower())
            return True
        
        # Find writes
        m = fpat_write.match(init_string)
        if m is not None:
            self["type"] = "write"
            self["params"]["format"] = m.group('format').lower()
            self["params"]["vars"] = fpat_var_string_num.findall(m.group('vars').lower())
            return True
        
        # Find definitions
        m = fpat_definition.match(init_string)
        if m is not None:
            self["type"] = "definition"
            self["params"]["name"] = m.group('name').lower()
            self["params"]["expr"] = m.group('expr').lower()
            return True
        
        return False
    
class FortranParser(list):
    """The FortranParser class implements a Fortran format parser."""
    
    class FortranParserError(Exception):
        """Exception raised if an error occurs during Fortran file parsing."""
        pass
    
    def __init__(self):
        """Initialize an empty FortranParser object."""
        list.__init__(self)
    
    def read(self,file):
        """Read the content of file, parse it and add it to the current object"""
        try:
            fsock = open(file,"r")
        except IOError,(errno, strerror):
            print "I/O error(%s): %s" % (errno, strerror)
        else:
            real_line = ""
            for line in fsock:
                
                line = line.strip('\n')
                
                # Pass empty lines
                if len(line) == 0:
                    continue
                
                # If line starts with a mulitag symbol, add it to the current real line
                
                if (len(line) > 5) and (line[5] in ('+','$','&')):
                        real_line = real_line + line[6:]
                
                # Otherwise, parse the current real line then store the current line in it
                else:
                    if real_line is not "":
                        
                        # deal correctly with strings cut with a multiline
                        real_line = real_line.replace('\"//\"','')
                        
                        try:
                            my_fs = FortranStatement(real_line)
                            self.add(my_fs)
                        except self.FortranParserError,why:
                            print "Warning:",why," Fortran statement ignored"
                    real_line = line
            
            # Required if the last line is also the last line of a multiline
            if real_line is not "":
                real_line = real_line.replace('\"//\"','')
                try:
                    my_fs = FortranStatement(real_line)
                    self.add(my_fs)
                except self.FortranParserError,why:
                    print "Warning:",why," Fortran statement ignored"

    def add(self,fs):
        """Add a FortranStatement object to the FortranParser item"""
        if not isinstance(fs,FortranStatement):
            raise ValueError # Trying to compare apple and ...
        try:
            self.append(fs)
        except self.FortranParserError,why:
                    print "Warning:",why," Fortran statement ignored"
    
    def to_string(self):
        """Return a string containing the fortran code corresponding to the
        Parser content."""
        
        file_str = ""
        line_str = ""
        
        for fs in self:
            
            if fs["type"] == "unknown":
                line_str = fs["init_string"]
            elif fs["type"] == "comment":
                line_str = "C     %s" % (fs["params"]["comment_string"])
            elif fs["type"] == "declaration":
                line_str = "      %s %s" % (fs["params"]["type"],','.join(fs["params"]["vars"]))
            elif fs["type"] == "common":
                line_str = "      common /%s/ %s" % (fs["params"]["name"],','.join(fs["params"]["vars"]))
            elif fs["type"] == "call":
                line_str = "      call %s(%s)" % (fs["params"]["func_name"],','.join(fs["params"]["vars"]))
            elif fs["type"] == "write":
                line_str = "      write %s %s" % (fs["params"]["format"],','.join(fs["params"]["vars"]))
            elif fs["type"] == "definition":
                line_str = "      %s = %s" % (fs["params"]["name"],fs["params"]["expr"])
                
            if fs["type"] == "comment":
                file_str = file_str + line_str + '\n' # Dont hash comments
            else:
                file_str = file_str + self._hash_str(line_str)

        return file_str
    
    def _hash_str(self,str,limit=70,multiline_symbol='+',num_spaces=5,split_str=False):
        """Divide a string into lines with a maximal number of characters given
        by limit. Each new line starts with num_spaces spaces followed by
        multiline_symbol. If split_str is False, a fortran string is always closed
        and re-opened on the new-line."""
                
        return_str = ""
        in_string = False
        counter = num_spaces
        
        while len(str) > 0:
            # remove first char
            curr_char = str[:1]
            str = str[1:]
            counter = counter + 1
            # if ", remember it
            if curr_char == '\"': in_string = not in_string
            
            # add char
            return_str = return_str + curr_char
            
            # if reach the limit, deal with carrier return: case in a str
            if (counter == limit - 1) and in_string and (not split_str):
                return_str = return_str + '\"\n' \
                             + " " * num_spaces \
                             + multiline_symbol + '//\"'
                counter = num_spaces + len(multiline_symbol) + 3
            
            # if reach the limit, deal with carrier return: case in a str
            if (counter == limit) and (not in_string):
                return_str = return_str + '\n' \
                             + " " * num_spaces \
                             + multiline_symbol
                counter = num_spaces + len(multiline_symbol) 
        
        return return_str+'\n'

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
    
    file = ""
    
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
    
    if file == "":
        file = raw_input("Please enter the name of the Fortran file to read: ")
        
    my_fp = FortranParser()
    my_fp.read(file)
    
    for statement in my_fp:
        print "LINE:",statement["init_string"]
        print "TYPE:",statement["type"],"PARAMS:",statement["params"],"\n"
    
if __name__=="__main__":
    main(sys.argv[1:])