#!/usr/bin/env python
"""Generate simple markdown documentation for the scripts in this directory"""

from StringIO import StringIO

if __name__ == '__main__':

    import os

    observed = set([])

    print('## Basic documentation\n')

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
    
        print("\n'''" + x.__name__ + "'''")

        if (x.__doc__ is not None):
            print('\n' + x.__doc__)
            
        if (hasattr(x, '_example_args')):
            print('')
            s = StringIO()
            s.write(d + '.py ')
            for k, v in x._example_args.items():
                s.write(k + ' ')
                if v is not None:
                    s.write(v + ' ')                
            print s.getvalue()            