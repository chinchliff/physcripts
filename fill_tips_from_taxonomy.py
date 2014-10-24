#!/usr/bin/env python

if __name__ == '__main__':

    import os, sys, newick3, phylo3, sqlite3

    def get_child_taxa(parent_name, rank, db_cursor):
        
        db_cursor.execute("SELECT left_value, right_value FROM taxonomy WHERE name LIKE ?", (parent_name.strip(),))
        left_value = None
        right_value = None
        try:
            left_value, right_value = db_cursor.fetchone()
        except TypeError:
            print("\nCould not find " + parent_name + "\n")
            return []
    
        print(parent_name)

        db_cursor.execute("SELECT name FROM taxonomy WHERE left_value > ? AND " \
            "right_value < ? AND name_class == 'scientific name' AND node_rank LIKE ?", (left_value, right_value, rank.strip()))
        return [row[0] for row in db_cursor.fetchall()]
    #    for row in db_cursor.fetchall():
    #        print row
    
    #    return []

    if len(sys.argv) < 4:
        print("usage: fill_tips_from_taxonomy.py <treefile> <taxonomydb> <outfile>")
        sys.exit(0)

    tree = None
    with open(sys.argv[1],"r") as intree_file:
        tree = newick3.parse(intree_file.readline())
    
        conn = sqlite3.connect(sys.argv[2])
        c = conn.cursor()
    
        for tip in tree.leaves():
            for child_name in get_child_taxa(tip.label, "genus", c):
                child = phylo3.Node()
                child.label = child_name
                child.istip = True
    #            print child.label
                tip.add_child(child)
    #            print([t.label for t in tip.children])
    
    with open(sys.argv[3],"w") as outfile:
    #    print(newick3.to_string(tree)+";")
        outfile.write(newick3.to_string(tree)+";")
