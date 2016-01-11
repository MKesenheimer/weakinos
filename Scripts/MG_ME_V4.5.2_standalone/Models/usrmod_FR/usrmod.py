#!/usr/bin/env python

"""
The usrmod v2 main script. This script can be ran from any FeynRules
compatible model directory to allow simple, user defined, modifications
of the particle, interaction and parameter content.
 
Usage: usrmod.py [options] save|create|restore 

Running modes:
    save              the first step for new model creation, saves all
                      the relevant information contained in the current
                      model files into a single archive file
    create            the second step for new model creation, compare existing
                      model files with an archive file to identify user
                      modifications and update all files to create a valid
                      model directory
    restore           restore all model files to a previous state, saved in an
                      archive file
    
Options:
  -f ...,--file=...       use a different name for the database file, useful
                          to work with various model versions in one directory
  -p ...,--path=...       specify a different working path than the current one
  -v,--verbose            switch the verbose mode on    
  -h,--help               display this help

This file is part of MadGraph/MadEvent, a free software to generate
matrix element amplitudes and parton level events.
Visit http://madgraph.phys.ucl.ac.be.
"""

__author__="Michel Herquet, the MadTeam, last change by $Author: madgraph $"
__version__= "$Revision: 1.5 $"
__date__= "2009-08-13"
__date__= "$Date: 2009/12/22 14:22:14 $"
__copyright__= "Copyright (c) 2009 MadGraph/MadEvent team"
__license__= "GNU"

###############################################################################
##
##  IMPORT & GLOBAL
##
###############################################################################

import os
import sys
import getopt
import pickle
import re

from particles_parser import ParticlesParser 
from interactions_parser import InteractionsParser
from lhaparser import LHAParser,LHAParameter
from ident_parser import IdentParser
from fparser import FortranParser, FortranStatement

model_files = []

verbose = False


###############################################################################
##
##  METHODS
##
###############################################################################

def usage():
    """Print the script usage documentation.
    """
    print __doc__

def header():
    """Print the program banner"""
    
    header = """
    ****************************************************************
    **                                                            **
    **                          USRMOD v2                         **
    **      User defined MadGraph/MadEvent models made easy...    **
    **                                                            **
    ****************************************************************
    """
    
    print header
    
def find_files():
    """Fill the global model_files list with all the relevant fortran
    file names in the current directory."""
    
    global model_files
    
    model_files.extend(["helas_couplings.inc",\
                        "coupl.inc",\
                        "coupl_write.inc",\
                        "input.inc",\
                        "intparam_definition.inc",\
                        "param_read.inc",\
                        "param_write.inc"
                        ])
    
    for file in os.listdir('.'):
        if re.match(r'couplings[0-9]+.f',file):
            model_files.append(file)
        
def ask_ok(prompt, retries=3, complaint='Please answer Yes or No.'):
    """Ask a question (prompt) and wait for a Yes or No answer.
       Return true if yes, false if no.
    """   
    while retries >= 0:
        ok = raw_input(prompt)
        if ok.lower() in ('y', 'ye', 'yes'): return True
        if ok.lower() in ('n', 'no', 'nop', 'nope'): return False
        retries = retries - 1
        print complaint
    raise IOError, 'No valid answer obtained'

def bytestr(size, precision=1):
    """Return a string representing the greek/metric suffix of a size"""
    
    _abbrevs = [
    (1<<50L, ' PB'),
    (1<<40L, ' TB'),
    (1<<30L, ' GB'),
    (1<<20L, ' MB'),
    (1<<10L, ' kB'),
    (1, ' bytes')
    ]
    
    if size==1:
        return '1 byte'
    for factor, suffix in _abbrevs:
        if size >= factor:
            break

    float_string_split = `size/float(factor)`.split('.')
    integer_part = float_string_split[0]
    decimal_part = float_string_split[1]
    if int(decimal_part[0:precision]):
        float_string = '.'.join([integer_part, decimal_part[0:precision]])
    else:
        float_string = integer_part
    return float_string + suffix

def read_files():
    """Read and parse all the existing model file information, then save it in a
    dictionnary structure"""
    
    global verbose
    global model_files
    
    save_dict = {}
    
    mypp = ParticlesParser()
    mypp.read("particles.dat")
    save_dict["particles.dat"] = mypp
    if verbose: print "File particles.dat successfully read, found %i particles" \
        % len(mypp)

    myip = InteractionsParser()
    myip.associate_part_parser(mypp)
    myip.read("interactions.dat")
    save_dict["interactions.dat"] = myip
    if verbose: print "File interactions.dat successfully read, found %i identified interactions" \
        % len(myip)
    
    mylhap = LHAParser()
    mylhap.read("param_card.dat")
    save_dict["param_card.dat"] = mylhap
    if verbose: print "File param_card.dat successfully read, found %i parameters" \
        % len(mylhap)
    
    myidp = IdentParser()
    myidp.associate_part_parser(mypp)
    myidp.read("ident_card.dat")
    save_dict["ident_card.dat"] = myidp
    if verbose: print "File ident_card.dat successfully read, found %i parameters" \
        % len(myidp)
    
    for file in model_files:
        myfp = FortranParser()
        myfp.read(file)
        save_dict[file] = myfp
        if verbose: print "Fortran file %s successfully read, found %i statements" \
            % (file,len(myfp))
    
    return save_dict
    
def save(data_file,file_dict=None):
    """Call the read function and store the resulting dict in a pickle file.
    If file_dict is given, do not call read_file and save this one instead"""
    
    if verbose: print "Now saving data..."

    try:
        output = open(data_file, 'wb')
    except (IOError, os.error):
            print "The data file %s cannot be created, stopping." % data_file
            sys.exit(1)
    else:
        # Here we save the whole dictionnary using the best available pickle format
        if file_dict is None:
            pickle.dump(read_files(), output, -1)
            print "Model information successfully stored in %s (%s)" \
            % (data_file,bytestr(os.stat(data_file).st_size))
        else:
            pickle.dump(file_dict, output, -1)
        
        output.close()
            

def load(data_file):
    """Load a file dictionnary from an archive pickle file."""
    try:
        input = open(data_file, 'rb')
    except (IOError, os.error):
            print "The data file %s cannot be read, stopping." % data_file
            sys.exit(1)
    else:
        try:
            file_dict = pickle.load(input)
        except:
            print "Error while reading %s, stopping." % data_file
            return None
        input.close()
    
    return file_dict

def restore(data_file,file_dict=None):
    """Write down a file dictionnary content into current model files.
    If file_dict is None, load the dictionnary content from the pickle
    archive."""
    
    global verbose
    
    if file_dict is None:
        ask_ok("All the current model file content will be lost. Are you sure ?")
        file_dict = load(data_file)
    
    for file, parser in file_dict.items():
        try:
            output = open(file, 'w')
            output.write(parser.to_string())
            if verbose:
                print "Content of file %s successfully written." % file
            output.close()
        except:
            print "Error while trying to write on %s, skipping." % file

def create(data_file):
    """Create the new model by loading information from the archive file and from
    current model file, comparing, updating and writing modifications."""
    
    global verbose
    global model_files
    
    ask_ok("The current model file content may be affected. Are you sure ?")
    
    ref_file_dict = load(data_file)
    new_file_dict = read_files()
    
    # Compare particles.dat parsers
    ref_pp = ref_file_dict["particles.dat"]
    new_pp = new_file_dict["particles.dat"]
    diff_pp = new_pp.compare(ref_pp)
    
    print "File particles.dat parsed and compared, found %i removed, %i new and %i modified particles." %\
        (len(diff_pp["removed"]),len(diff_pp["new"]),len(diff_pp["modified"]))

    # Compare interactions.dat parsers
    ref_ip = ref_file_dict["interactions.dat"]
    new_ip = new_file_dict["interactions.dat"]
    diff_ip = new_ip.compare(ref_ip)
    
    print "File interactions.dat parsed and compared, found %i removed, %i new and %i modified interactions." %\
        (len(diff_ip["removed"]),len(diff_ip["new"]),len(diff_ip["modified"]))
    
    # Compare param_card.dat parsers
    ref_idp = ref_file_dict["ident_card.dat"]
    new_idp = new_file_dict["ident_card.dat"]
    diff_idp = new_idp.compare(ref_idp)
    
    print "File ident_card.dat parsed and compared, found %i removed, %i new and %i modified parameters." %\
        (len(diff_idp["removed"]),len(diff_idp["new"]),len(diff_idp["modified"]))

    # Modified interactions and parameters need to be removed and added again:
    diff_ip["removed"].extend(diff_ip["modified"])
    diff_ip["new"].extend(diff_ip["modified"])
    
    diff_idp["removed"].extend(diff_idp["modified"])
    diff_idp["new"].extend(diff_idp["modified"])
    
    # The file_dict to consider is the reference one, where potentially modified files
    # are replaced
    
    file_dict = ref_file_dict
    file_dict["particles.dat"] = new_pp
    file_dict["interactions.dat"] = new_ip
    file_dict["ident_card.dat"] = new_idp
    
    # Deal with coupling modifications
    print "Updating files for coupling modifications ..."
    
    # Remove coupling variables which are not needed anymore 
    for coupl in diff_ip["removed"]:
        for var in coupl["variables"]:
            
            # First check if the corresponding variable does not appear in another coupling
            found = False
            for new_coupl in new_ip:
                if var in new_coupl["variables"]:
                    if verbose:
                        print "Variable %s also appears for another interaction and will not be removed." % var
                    found = True
                    break
            
            # otherwise, remove the variable from all model files    
            if not found:
                for file in model_files:
                    if verbose: print "Remove coupling variable %s from Fortran file %s" \
                        % (var,file)
                    rm_var(file_dict[file],var)
    
    # Add the required coupling variables which are not yet present
    for coupl in diff_ip["new"]:
        for var in coupl["variables"]:
            
            # First check if the corresponding variable does not appear in another coupling
            found = False
            for old_coupl in ref_ip:
                if var in old_coupl["variables"]:
                    if verbose:
                        print "Variable %s already appears for another interaction and will not be added." % var
                    found = True
                    break
            
            # otherwise, add the variable to the relevant model files    
            if not found:
                add_coupl_var(var,coupl,file_dict)
                    
    # Deal with parameter modifications
    print "Updating files for parameter modifications ..."

    # Added
    for param in diff_idp["new"]:
        add_param_var(param,file_dict)
    
    # Removed
    for param in diff_idp["removed"]:
        for file in model_files:
            if verbose: print "Remove external parameter %s from Fortran file %s" \
                % (param["var"],file)
            rm_var(file_dict[file],param["var"])

    if verbose: print "Writing new files ..."
    restore(data_file,file_dict)
    
def rm_var(fortran_parser,var_name):
    """Remove a variable from all statements in a given Fortran Parser. If the variable
    appears in the RHS of a definition, print a warning"""
    
    to_be_removed = []
        
    # Check all statements
    for f_statement in fortran_parser:
        
        # Check all declaration and common variables, as well as call and write variables 
        if f_statement["type"] in ("declaration","common","call","write"):
            
            for listed_var in f_statement["params"]["vars"]:
                # If var_name appears (possibly with (...)), remove it
                fnd_var = re.match("(?P<name>%s(?:\([^\)]*\))?$)" % var_name,listed_var)
                if fnd_var:
                    if f_statement["type"] in ("declaration","common"):
                        f_statement["params"]["vars"].remove(fnd_var.group("name"))
                        if len(f_statement["params"]["vars"]) == 0:
                            if f_statement not in to_be_removed:
                                to_be_removed.append(f_statement)
                    elif f_statement["type"] in ("call","write"):
                        if f_statement not in to_be_removed:
                            to_be_removed.append(f_statement)
                
                
                
        if f_statement["type"] == "definition":             
            # If var_name is in RHS, display a warning
            if re.match(".*(?:\(|\)|\s|\+|\-|\*|\/|\A)%s(?:\(|\)|\s|\+|\-|\*|\/|\Z).*" % var_name,\
                        f_statement["params"]["expr"]):
                print "Warning: variable %s appears in the RHS of variable %s definition" %\
                    (var_name,f_statement["params"]["name"])
            # Remove definition if var_name is in LHS
            if re.match("%s(\([^\)]+\))?" % var_name,f_statement["params"]["name"]):
                if f_statement not in to_be_removed:
                    to_be_removed.append(f_statement)
        
    # Remove statements no longer needed
    for f_statement in to_be_removed:
        fortran_parser.remove(f_statement)
            
            
            
        
def add_coupl_var(var_name,coupl,file_dict):
    """Add a given coupling variable associated to a particular coupling to the
    relevant model files, i.e. coupl.inc, coupl_check.inc, coupl_write.inc and add
    an empty entry to the last couplings??.f file."""
    
    global verbose
    
    # Initialize various writing variables depending on coupling type
    
    if coupl["type_var"] == "r":
        common_prefix = "real_coupl"
        decl_prefix = "double precision"
        var_decl_str = var_name
        write_statement = FortranStatement\
                          ("      write (*,1) \'%s=\', %s" % (var_name,var_name))
        check_statement = FortranStatement\
                          ("      write (1,11) \'%s \', %s" % (var_name,var_name))
        def_statements = [FortranStatement("      %s = 0" % var_name)]
    elif coupl["type_var"] == "c":
        common_prefix = "comp_coupl"
        decl_prefix = "double complex"
        var_decl_str = var_name
        write_statement = FortranStatement\
                          ("      write (*,2) \'%s=\', %s" % (var_name,var_name))
        check_statement = FortranStatement\
                          ("      write (1,12) \'%s \', %s" % (var_name,var_name))
        def_statements = [FortranStatement("      %s = 0" % var_name)]
    elif coupl["type_var"] == "dc":
        common_prefix = "chiral_coupl"
        decl_prefix = "double complex"
        var_decl_str = var_name + "(2)"
        write_statement = FortranStatement\
                          ("      write (*,4) \'%s(l)=\',%s(1),\'%s(r)=\',%s(2)"\
                           % (var_name,var_name,var_name,var_name))
        check_statement = FortranStatement\
                          ("      write (1,13) \' %s \',%s(1),%s(2)"\
                           % (var_name,var_name,var_name))
        def_statements = [FortranStatement("      %s(1) = 0" % var_name),
                         FortranStatement("      %s(2) = 0" % var_name)]
    else:
        raise ValueError, 'Unknown type for variable '+var_name
        return False
    
    # COUPL.INC
    #----------
    
    if verbose:
        print "Add coupling variable %s to Fortran file coupl.inc" % var_name
    # Read the file backward...
    for statement in file_dict["coupl.inc"][::-1]:
        
        # If the correct common block is found...
        if statement["type"] == "common":
            if re.search(common_prefix, statement["params"]["name"]):
                # Add the variable to the common block
                if var_name not in statement["params"]["vars"]:
                    statement["params"]["vars"].append(var_name)
                    # Start searching from there for a valid declaration block, then break
                    for following_statement in \
                        file_dict["coupl.inc"][file_dict["coupl.inc"].index(statement)::-1]:
                        # For the first found, add the variable and break the loop
                        if following_statement["type"] == "declaration":
                            if following_statement["params"]["type"] == decl_prefix:
                                if var_decl_str not in following_statement["params"]["vars"]:
                                    following_statement["params"]["vars"].append(var_decl_str)
                                    break
                    break
                            
    
    # COUPL_CHECK.INC
    #----------------
    
    if verbose:
        print "Add coupling variable %s to Fortran file helas_couplings.inc" % var_name
    
    # Write the coupling writing line
    if check_statement not in file_dict["helas_couplings.inc"]:
        file_dict["helas_couplings.inc"].add(check_statement)
    
    # COUPL_WRITE.INC
    #----------------
    if verbose:
        print "Add coupling variable %s to Fortran file coupl_write.inc" % var_name
        
    # If usrmod section does not exist, create it 
    head_sec = FortranStatement("      write (*,*) \' usrmod couplings: \'")
    if head_sec not in file_dict["coupl_write.inc"]:
        file_dict["coupl_write.inc"].add(head_sec)
    
    # Write the coupling writing line
    if write_statement not in file_dict["coupl_write.inc"]:
        file_dict["coupl_write.inc"].add(write_statement)

    # COUPLINGS.F
    #---------------------
    
    
    # Find the last couplings??.f file
    
    file_num = 0
    for filename in file_dict.keys():
        fnd_file = re.match(r'couplings(?P<name>[0-9]+).f',filename)
        if fnd_file and file_num < int(fnd_file.group("name")):
            file_num = int(fnd_file.group("name"))
    
    filename = "couplings%i.f" % file_num
    
    # Find the return position
    for fstatement in file_dict[filename]:
        if (fstatement["type"] == "unknown") and \
            re.match('\s*endif\s*',fstatement["init_string"],re.IGNORECASE):
            index = file_dict[filename].index(fstatement)
    
    if verbose:
        print "Add coupling variable %s to Fortran file %s" % (var_name,filename)
        
    for def_statement in def_statements:
        if def_statement not in file_dict[filename]:
            file_dict[filename].insert(index,def_statement)

def add_param_var(param,file_dict):
    """Add a given external parameter variable to the relevant model files,
    i.e. param_card.dat, input.inc, param_read.inc and param_write."""
    
    global verbose
    
    # PARAM_CARD.DAT
    #---------------
    
    if verbose:
        print "Add parameter %s to param_card.dat" % param["var"]
    
    if param["block"] == "mass":
        value = "1.0000000000000e+02"
        type = "standard"
    elif param["block"] == "decay":
        value = "1.0000000000000e+00"
        type = "decay"
    else:
        value = "0.0000000000000e+00"
        type = "standard"
    
    # Transform a ident_card member into a LHAParam object
    lhaparam = LHAParameter(param["block"],param["ids"],value,param["var"],type)
    
    # Add the LHAParam object at the right place in the LHAParser
    
    # By default put it at the end
    index = len(file_dict["param_card.dat"])
    
    for parser_member in file_dict["param_card.dat"]:
        curr_block = parser_member["block"]
        next_block = parser_member["block"]
        for next_member in file_dict["param_card.dat"]\
                [file_dict["param_card.dat"].index(parser_member)+1:]:
            if next_member["type"] is not "br":
                next_block = next_member["block"]
                break
            else: continue
        # If found the right block, and next param has a different block,
        # put it right before
        if (lhaparam["block"] == curr_block) and (next_block != curr_block):
            index = file_dict["param_card.dat"].index(next_member)
    
    file_dict["param_card.dat"].insert(index,lhaparam)
    
    # COUPL.INC & INPUT.INC
    #-----------------
    
    # Masses and widths should also appear in coupl.inc
    
    if param["block"] == "mass":
        
        if verbose:
            print "Add mass parameter %s to Fortran file coupl.inc" % param["var"]
        
        for statement in file_dict["coupl.inc"][::-1]:
            if statement["type"] == "common":
                if re.search("masses", statement["params"]["name"]):
                    if param["var"] not in statement["params"]["vars"]:
                        statement["params"]["vars"].append(param["var"])
                    for following_statement in \
                        file_dict["coupl.inc"][file_dict["coupl.inc"].index(statement)::-1]:
                        if following_statement["type"] == "declaration":
                            if following_statement["params"]["type"] == "double precision":
                                if param["var"] not in following_statement["params"]["vars"]:
                                    following_statement["params"]["vars"].append(param["var"])
                                    break
                    break
    
    
    elif param["block"] == "decay":
        
        if verbose:
            print "Add width parameter %s to Fortran file coupl.inc" % param["var"]
        
        for statement in file_dict["coupl.inc"][::-1]:
            if statement["type"] == "common":
                if re.search("widths", statement["params"]["name"]):
                    if param["var"] not in statement["params"]["vars"]:
                        statement["params"]["vars"].append(param["var"])
                    for following_statement in \
                        file_dict["coupl.inc"][file_dict["coupl.inc"].index(statement)::-1]:
                        if following_statement["type"] == "declaration":
                            if following_statement["params"]["type"] == "double precision":
                                if param["var"] not in following_statement["params"]["vars"]:
                                    following_statement["params"]["vars"].append(param["var"])
                                    break
                    break

    
    else:
        
        if verbose:
            print "Add parameter %s to Fortran file input.inc" % param["var"]
        # Read the file backward...
        for statement in file_dict["input.inc"][::-1]:
            
            # If the correct common block is found...
            if statement["type"] == "common":
                if re.search("external_params", statement["params"]["name"]):
                    # Add the variable to the common block
                    if param["var"] not in statement["params"]["vars"]:
                        statement["params"]["vars"].append(param["var"])
                    # Start searching from there for a valid declaration block, then break
                    for following_statement in \
                        file_dict["input.inc"][file_dict["input.inc"].index(statement)::-1]:
                        # For the first found, add the variable and break the loop
                        if following_statement["type"] == "declaration":
                            if following_statement["params"]["type"] == "double precision":
                                if param["var"] not in following_statement["params"]["vars"]:
                                    following_statement["params"]["vars"].append(param["var"])
                                    break
                    break
    
    # PARAM_WRITE.INC
    #----------------
    
    if verbose:
        print "Add parameter %s to Fortran file param_write.inc" % param["var"]
    
    
    
    if param["block"] == "mass":
        tag = "masses"
        write_statement = FortranStatement\
                          ("      write (*,3) \'%s=\', %s" % (param["var"],param["var"]))
    elif param["block"] == "decay":
        tag = "widths"
        write_statement = FortranStatement\
                          ("      write (*,3) \'%s=\', %s" % (param["var"],param["var"]))
    else:
        tag = "external"
        write_statement = FortranStatement\
                          ("      write (*,1) \'%s=\', %s" % (param["var"],param["var"]))
        
    # Found the right place to add an external param write statement
    index = len(file_dict["param_write.inc"])
    for statement in file_dict["param_write.inc"]:
        if statement["type"] == "write":
            if re.match(".*%s.*" % tag,statement["params"]["vars"][0], re.IGNORECASE):
                index = min(file_dict["param_write.inc"].index(statement)+2,\
                            len(file_dict["param_write.inc"]))
        
    # Write the coupling writing line
    if write_statement not in file_dict["param_write.inc"]:
        file_dict["param_write.inc"].insert(index,write_statement)
        
    
    # PARAM_READ.INC
    #----------------
    
    if verbose:
        print "Add parameter %s to Fortran file param_read.inc" % param["var"]
    
    write_statement = FortranStatement\
                          ("      call lha_get_real(npara,param,value,\"%s\",%s,\"% 16.13E\")" \
                           % (param["var"],param["var"],float(value)))
    
    if write_statement not in file_dict["param_read.inc"]:
        file_dict["param_read.inc"].append(write_statement)
            
###############################################################################
##
##  MAIN
##
###############################################################################
    
def main(argv):
    """Main function, called if script is ran from the command line.
    """
    
    global verbose
    data_file = "ref_data.pkl"
    path = "./"
    
    try:
        options, discard = getopt.getopt(argv,"hvf:p:",\
                                         ["help","verbose","file=","path="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
    for option, value in options:
        if option in ("-h", "--help"):
            usage()
            sys.exit()
        elif option in ("-v", "--verbose"):
            verbose = True        
        elif option in ("-f", "--file"):
            data_file = value         
        elif option in ("-p", "--path"):
            path = value
            try:
                os.chdir(path)
            except (IOError, os.error):
                print "Problem when trying to set the current path to: " + path
    try:
        run_mode = discard[0].lower()
        if run_mode not in ("save","create","restore"):
            raise ValueError
    except:
        usage()
        sys.exit(2)
    
    header()
    find_files()
    
    if verbose:
        print "Current path: " + os.getcwd()
        print "Selected archive file: " + data_file
        print "Running mode: " + run_mode
        print "Found %i relevant Fortran files in the current directory" % len(model_files)
    
    if run_mode == "save":
        save(data_file)
    elif run_mode == "create":
        create(data_file)
    elif run_mode == "restore":
        restore(data_file)
        
    
if __name__=="__main__":
    main(sys.argv[1:])
