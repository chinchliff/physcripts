#!/usr/bin/env python3

if __name__ == '__main__':

    import os, sys

    if len(sys.argv) < 3:
        print("usage: filter_list.py <inputlist> [<accept=valuesfile> | <exclude=valuesfile>] <keycolumn>")
        sys.exit(0)

    key_column_label = sys.argv[3].strip()

    filter_arg_parts = [s.strip() for s in sys.argv[2].split("=")]
    mode = filter_arg_parts[0]
    filter_fname = filter_arg_parts[1]
    if mode not in {"accept","exclude"}:
        print("mode but be set to \"accept\" or \"exclude\", not \"" + mode + "\"")
        sys.exit(0)

    filter_values = set((s.strip() for s in open(filter_arg_parts[1])))

    key_column_index = -1
    on_first_line = True
    with open(sys.argv[1], "r") as list_infile:

        for line in list_infile.readlines():
                
            parts = [s.strip() for s in line.split(",")]
        
            if on_first_line:
                on_first_line = False
                for i, val in enumerate(parts):
                    if val == key_column_label:
                        key_column_index = i
                    
                if key_column_index == -1:
                    print("the key column \"" + key_column_label + "\" must be present in the input list") 
                    sys.exit(0)

                print(line.strip())
                continue

            if (mode == "accept"):
                if parts[key_column_index].strip() in filter_values:
                    print(line.strip())
            else:
                if parts[key_column_index].strip() not in filter_values:
                    print(line.strip())