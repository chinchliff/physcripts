#!/usr/bin/env python
"""Generate simple markdown documentation for the scripts in this directory"""

from StringIO import StringIO

if __name__ == '__main__':

    import os

    observed = set([])

    print("""## Basic documentation
    
This directory contains a bunch of scripts to automate common (or not so common) tasks, mostly targeted at bioinformatics and phylogenetic systematic research. To install the scripts, clone the git repo into a local directory and add that directory to your PATH and PYTHONPATH environment variables. E.g.:

```bash
# clone the scripts repo
cd ~/ && git clone https://github.com/chinchliff/physcripts

# add it to paths
echo 'export PATH=/Users/cody/physcripts:$PATH' >> ~/.bash_profile
echo 'export PYTHONPATH=/Users/cody/physcripts:$PYTHONPATH' >> ~/.bash_profile

# load new environment variables
source ~/.bash_profile

## Running examples

Example input files (referred to in examples below) are included in the `test_files` directory, and the example script calls are written to be run from within the physcripts directory. To run the examples (using the example inputs) outside of that directory, you'll need to edit the calls to indicate the correct path to the example input files. 

```""")

    for d in os.listdir('.'):

        failed = False

        try:
            d = d.rsplit('.py')[0]
            exec('import ' + d + ' as x')
        except:
            continue
        
        if d in observed:
            continue

        observed.add(d)
        
        if (hasattr(x, '_title')):
            print('\n### ' + x._title)
    
        print('\n`' + x.__name__ + '`')

        if (x.__doc__ is not None):
            print('\n' + x.__doc__)
            
        if (hasattr(x, '_example_args')):
            print('```bash\n')
            s = StringIO()
            s.write(d + '.py ')
            for k, v in x._example_args.items():
                s.write(k + ' ')
                if v is not None:
                    s.write(v + ' ')
            s.write('\n```')                
            print s.getvalue()            