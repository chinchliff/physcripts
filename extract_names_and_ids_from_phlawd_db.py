#!/usr/bin/env python
"""Given a set of names, extract information about those names from the specified PHLAWD database."""
_title = 'Get information about names from a PHLAWD database'
_example_args = {
    '-t': 'test_files/tree_with_50_tips.tre',
    '-n': 'test_files/names_50.txt',
    '-b': '0.1',
    '-s': None}

if __name__ == '__main__': # if the script is being executed at the command line

    import argparse, sqlite3

    # define command line arguments using the argparse module,
    # which will simplify a lot of annoying work for us
    description = __doc__

    parser = argparse.ArgumentParser(description=description)

    # tell argparse to open the phlawd db file using the sqlite3.connect() function
    parser.add_argument('-d', '--phlawd-db', type=sqlite3.connect, nargs=1, \
        required=True, help='The phlawd database from which to extract information.')

    # tell argparse to open the names file as a python file object for reading
    parser.add_argument('-n', '--names', type=argparse.FileType('r'), nargs=1, \
        required=True, help='A file containing names on separate lines--each name will '
                            'be submitted as a separate query to the phlawd db.')

    # process the incoming command line args. this does the tedious work of validating the
    # incoming arguments and opening the corresponding files. if the incoming arguments
    # can't be processed the way we specified above, argparse will tell the user so and exit
    args = parser.parse_args()
    
    # set the phlawd_db variable to the file opened by argparse
    # (the db connection object is the first item in a list of length 1)
    phlawd_db = args.phlawd_db[0]
    c = phlawd_db.cursor()
    
    # read all the names into a list, stripping leading/trailing whitespace
    names = [l.strip() for l in args.names[0]]

    # write the header line for the output
    print 'accepted_name,input_name,name_class,ncbi_id'
    
    # initialize a dictionary to remember accepted names for ids so we only have to
    # look up the accepted name for each id once
    accepted_names_by_ncbi_id = {}
    
    for n in names:
        
        # query the database using this name. the 'execute' method will replace the '?'
        # characters in the base_query with the values from the data tuple (in order)
        data = (n,)
        r = c.execute('SELECT name_class, ncbi_id FROM taxonomy '
                      'WHERE name == ?', data)
        
        # get the data resulting from the query. this will return a list of tuples, each
        # containing the values of the fields in the select statement, in the same order
        # in which they were specified, so:
        # result[0][0] => 'name' field of first row
        # result[0][1] => 'name_class' field of first row
        # result[1][2] => 'ncbi_id' field of second row, etc.
        result = r.fetchall()

        # if there were no hits, result will be an empty list. in that case, just create
        # an empty dummy row so we get a blank line in the output for this name
        if len(result) < 1:
            result.append(('', ''))
        
        # process all the hits for this name... don't think there should just be more than
        # one hit per name, but might be edge cases...
        for row in result:
        
            # if we haven't already seen this id, then look up the accepted name for it
            id = row[1]
            if id not in accepted_names_by_ncbi_id:
                r = c.execute('SELECT name FROM taxonomy WHERE '
                              'name_class == "scientific name" and ncbi_id == ?', (id,))
                
                # set accepted name to blank string unless we found a match
                d = r.fetchone()
                accepted_name = d[0] if d is not None else ''
                
                # add the name to the dict under its id
                accepted_names_by_ncbi_id[id] = accepted_name
            
            # add the search name and the accepted name to the results
            row = (n, accepted_names_by_ncbi_id[id]) + row         
            
            # print this result row to the screen, as a comma-separated list of string
            # representations of its elements
            print ','.join([str(s) for s in row])
            
            