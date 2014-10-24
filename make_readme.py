#!/usr/bin/env python
"""Generate simple markdown documentation for the scripts in this directory"""

preamble = """## Basic documentation
    
This directory contains a bunch of scripts to automate common (or not so common) tasks, mostly targeted at bioinformatics and phylogenetic systematic research. To install the scripts, clone the git repo into a local directory and add that directory to your PATH and PYTHONPATH environment variables. E.g.:

```bash
# clone the scripts repo
cd ~/ && git clone https://github.com/chinchliff/physcripts

# add it to paths
echo 'export PATH=/Users/cody/physcripts:$PATH' >> ~/.bash_profile
echo 'export PYTHONPATH=/Users/cody/physcripts:$PYTHONPATH' >> ~/.bash_profile

# load new environment variables
source ~/.bash_profile
```

## Running examples

The descriptions below include a variety of example calls to the scripts. The input files used in these examples are included in the `test_files` directory. The example calls are intended to be run from within the physcripts directory, but you can run the elsewhere by if you modify the calls to indicate the correct path to the example input files.

## Scripts

Below is a (complete) list of the available Python scripts, along with (highly incomplete) descriptive documentation about each one. Links to scripts with more complete documentation are provided below, and the documentation for these scripts precedes the docs for those with less information.

More details can usually be found in the form of comments within the scripts (in addition, of course, to my optimistically somewhat self-documenting code).

#### Scripts with more complete documentation:

"""

if __name__ == '__main__':

    import os, re
    from StringIO import StringIO

    descriptions = {}
    titles = {}

    for d in os.listdir('.'):

        failed = False

        try:
            d = d.rsplit('.py')[0]
            exec('import ' + d + ' as x')
        except:
            continue
        
        if d in descriptions:
            continue

        s = StringIO()
        s.write('\n---\n')
        
        if (hasattr(x, '_title')):
            s.write('\n### ' + x._title + '\n')
            titles[d] = x._title
    
        s.write('\nscript: `' + x.__name__ + '.py`\n')

        if (x.__doc__ is not None):
            s.write('\n' + x.__doc__ + '\n')
            
        if (hasattr(x, '_example_args')):
            s.write('```bash\n')
            s.write(d + '.py ')
            for k, v in x._example_args.items():
                s.write(k + ' ')
                if v is not None:
                    s.write(v + ' ')
            s.write('\n```\n')
        
        descriptions[d] = s.getvalue()
    
    # now print the document
    print(preamble)
    
    scripts = sorted(descriptions.keys())

    # create links to descriptions with titles    
    for d in scripts:
        if d in titles:
            print('[' + titles[d] + '](#' + re.sub("\s+", "-", titles[d].lower().strip()) + ')\n')

    # write descriptions with titles
    for d in scripts:
        if d in titles:
            print descriptions[d]

    print("#### Scripts with less complete documentation:")

    # write other (incomplete) descriptions    
    for d in scripts:
        if d not in titles:
            print descriptions[d]