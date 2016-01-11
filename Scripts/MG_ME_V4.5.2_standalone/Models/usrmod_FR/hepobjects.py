"""
This module defines the class Particle and Interaction associated to the
corresponding physical notions.

Diagnostic usage: python hepobjects.py

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
##  CLASSES
##
###############################################################################
    
class Particle(dict):
    """Particle object with attributes and associated routines"""
    
    class ParticleError(Exception):
        """Exception raised if an error occurs in the definition
        of a particle"""
        pass

    def __init__(self):
        """Initialize an empty Particle object."""
        dict.__init__(self)
        
        prop_names=["name","spin","line_type","mass","width",\
                   "color_rep","draw_name","pdg_code","kind"]
        
        for prop_name in prop_names:
            self[prop_name] = None
        
        # By default, let's assume is a particle and not an anti particle
        self["kind"] = "part" 
    
    def set(self,prop_name,value):
        """Set the value of a given particle property labeled prop_name"""
        
        # Check case by case if the input makes sense
        
        if prop_name == "name":
            
            try:
                name = value["part"]
                anti_name = value["apart"]
            except:
                raise self.ParticleError, \
                repr(value)+" is not a valid two item dictionary to define the name."
            
            if (name+'~' != anti_name) \
            and (name.replace('-','+') != anti_name) \
            and (name.replace('+','-') != anti_name) \
            and (name != anti_name):
                raise self.ParticleError, name+','+anti_name+": names do not match"
        
        elif prop_name == "spin" and value not in ('s','f','v','t'):
            raise self.ParticleError, \
            self["name"]["part"]+": "+value+" is not a valid spin tag."
            
        elif prop_name == "line_type" and value not in ('s','d','w','c'):
            raise self.ParticleError, \
            self["name"]["part"]+": "+value+" is not a valid line type tag."
            
        elif prop_name == "color_rep" and value not in ('s','t','o'):
            raise self.ParticleError, \
            self["name"]["part"]+": "+value+" is not a valid color representation tag."
        
        elif prop_name == "pdg_code":
            try:
                dummy = int(value)
            except:
                raise self.ParticleError, \
                self["name"]["part"]+": "+value+" is not a valid PDG code."
        
        # The kind tag tells if the object corresponds to a particle or its antiparticle
        elif prop_name == "kind" and value not in ("part","apart"):
            raise self.ParticleError, \
            self["name"]["part"]+": "+value+" is not a valid kind tag, only part or apart."
            
        self[prop_name]=value # If everything went well, assign the value
    
    def is_identical(self,other_part):
        """Check if two particles have strictly the same properties. Return True
        if yes, False if not"""
        
        if not isinstance(other_part,Particle):
            raise ValueError # Trying to compare apple and ...
        
        result=True
        
        for prop_name in self.keys():
            if self[prop_name] != other_part[prop_name]: result=False
        
        return result
    
    def same_name(self,other_part):
        """Check if two particles have strictly the same names and antiname (or switched).
        Return True if yes, False if not"""
        
        if not isinstance(other_part,Particle):
            raise ValueError # Trying to compare apple and ...
        
        if (self["name"]["part"] != other_part["name"]["part"] \
            and self["name"]["part"] != other_part["name"]["apart"]) or \
            (self["name"]["apart"] != other_part["name"]["apart"] \
             and self["name"]["apart"] != other_part["name"]["part"]):
            return False
        
        return True


class Interaction(dict):
    """Interaction object with attributes and associated routines"""
    
    _ref_types = (\
        {"type":"ffs","num_var":1,"type_var":"dc"},\
        {"type":"ffv","num_var":1,"type_var":"dc"},\
        {"type":"sss","num_var":1,"type_var":"c"},\
        {"type":"vvs","num_var":1,"type_var":"c"},\
        {"type":"vss","num_var":1,"type_var":"r"},\
        {"type":"vvt","num_var":1,"type_var":"r"},\
        {"type":"vvv","num_var":1,"type_var":"r"},\
        {"type":"ssss","num_var":2,"type_var":"c"},\
        {"type":"vvss","num_var":2,"type_var":"c"},\
        {"type":"vvvv","num_var":2,"type_var":"r"},\
        {"type":"vvvs","num_var":2,"type_var":"r"},\
        )
    
    class InteractionError(Exception):
        """Exception raised if an error occurs in the definition
        of an interaction"""
        pass
    
    def __init__(self):
        """Create an empty Interaction object."""
        dict.__init__(self)
        
        self["particles"] = []
        self["num_var"] = 0
        self["variables"] = []
        self["tags"] = []
        self["type"] = None
        self["type_var"] = None
        self["npart"] = 0
        
    
    def initialize(self,list_part):
        """Initialize the Interaction object using the list of particles
        entering the vertex. Return true if the vertex has been successfully
        identified, false otherwise (with an error)."""
        
        # Check input list
        try:
            for part in list_part:
                if not isinstance(part,Particle):
                    raise self.InteractionError, \
                        repr(list_part)+" is not a valid list of Particle objects."
        except:
            raise self.InteractionError, \
                        repr(list_part)+" is not a valid list of Particle objects."
                    
        self["particles"] = list_part
        self["npart"] = len(list_part)
        
        # Identify vertex type and set relevant variable once done
        my_type = ''
        for part in self["particles"]:
            my_type = my_type+part["spin"]
        permu_type = self._permutate(my_type)
        
        for type in self._ref_types:
            if type["type"] in permu_type:
                self["type"] = type["type"]
                self["type_var"] = type["type_var"]
                self["num_var"] = type["num_var"]
                return True
        raise self.InteractionError, \
                    "Type "+my_type+" of vertex "+\
                    ",".join([part["name"]["part"] for part in list_part])+" is unknown"
        return False
    
    def _permutate(self,seq):
        """permutate a sequence and return a list of all permutations"""
        if not seq:
            return [seq] # is an empty sequence
        else:
            temp = []
            for k in range(len(seq)):
                part = seq[:k] + seq[k+1:]
                for m in self._permutate(part):
                    temp.append(seq[k:k+1] + m)
            return temp
        
    def set(self,prop_name,value):
        """Set the value of a given particle property labeled prop_name"""
        
        if prop_name == "variables":
            if not isinstance(value,list):
                raise self.InteractionError, \
                    repr(value)+" is not a valid list of variables."
            elif len(value) != self["num_var"]:
                raise self.InteractionError, \
                    "Too many/too few variables for this vertex, requires "+repr(self["num_var"])
            else:
                self["variables"] = value
        
        elif prop_name == "tags":
            if not isinstance(value,list):
                raise self.InteractionError, \
                    repr(value)+" is not a valid list of tags."
            else:
                self["tags"] = value
        
        else:
            raise self.InteractionError, \
                    "Property "+repr(prop_name)+" does not exist or cannot be set manually."
            
            
    

###############################################################################
##
##  DIAGNOSTIC
##
###############################################################################

if __name__=="__main__":
    
    print __doc__
    
    mypart = Particle()
    myinter = Interaction()
    
    print "Known Particle properties:"+",".join([key for key in mypart.keys()])

    print "Known Interaction properties:"+",".join([key for key in myinter.keys()])
    
    print "Known Interaction types:"+",".join([type["type"] for type in myinter._ref_types])

    

