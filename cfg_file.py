import sys
import os

class reader:

    class parameter:
        
        def __init__(self, paramname, paramvalue):
            self.name = paramname
            self.value = paramvalue

    def __init__(self, pathtocfgfile):

        cwd = os.getcwd() + os.sep
        self.path_to_file = cwd + pathtocfgfile

        try:
            cfgfile = file(pathtocfgfile, "rU")
            cfglines = cfgfile.readlines()
        except IOError:
            exit("Error: could not find the specified config file.")

        self.parameters = list()

        for line in cfglines:

            tokens = line.split("=")
            param = tokens[0].strip()

            try:
                    val = tokens[1].strip()
            except IndexError:
                    val = None
        
            self.parameters.append(self.parameter(param,val))
