#!/usr/bin/env python3

import os, sys

if len(sys.argv) < 3:
    print("usage: combine_lists.py <input1> <input2> <keycolumn> [include=<col1,col2,...>] [exclude=<col1,col2,...>]")
    sys.exit(0)

key_column_label = sys.argv[3].strip()

# if the user has specified columns to include/exclude
colnames_to_include = set()
colnames_to_exclude = set()
user_include_colnames_set = False
user_exclude_colnames_set = False
if len(sys.argv) > 4:
    for arg in sys.argv[4:]:
        parts = arg.split("=")
        if parts[0] == "include":
            user_include_colnames_set = True
            for c in parts[1].split(","):
                colnames_to_include.add(c)
        elif parts[0] == "exclude":
            user_exclude_colnames_set = True
            for c in parts[1].split(","):
                colnames_to_exclude.add(c)

recorded_colnames = set()

# read the first list
list1 = {}
list1_colnames = []
list1_key_column_index = -1
on_first_line = True
with open(sys.argv[1], "r") as list1_infile:

    for line in list1_infile.readlines():
                
        parts = [s.strip() for s in line.split(",")]
        
        if on_first_line:
            on_first_line = False
            for i, val in enumerate(parts):
                if val == key_column_label:
                    list1_key_column_index = i
                    
            if list1_key_column_index == -1:
                print("the key column \"" + key_column_label + "\" must be present in both input lists") 
                sys.exit(0)

            list1_colnames = parts
            continue

        # make a new entry in the dict for this row
        key_value = parts[list1_key_column_index]
        list1[key_value] = {}

        # add all the other known values to the dict
        for i, val in enumerate(parts):
            if i != list1_key_column_index:
                
                cur_colname = list1_colnames[i]
                
                if user_include_colnames_set and cur_colname not in colnames_to_include:
                    continue
                
                if user_exclude_colnames_set and cur_colname in colnames_to_exclude:
                    continue
                    
                recorded_colnames.add(cur_colname) # not very efficient to attempt to add it on every line
                
                list1[key_value][cur_colname] = val

# read the second list
#list2 = {}
list2_colnames = []
list2_key_column_index = -1
on_first_line = True
with open(sys.argv[2], "r") as list2_infile:

    for line in list2_infile.readlines():
                
        parts = [s.strip() for s in line.split(",")]
        
        if on_first_line:
            on_first_line = False
            for i, val in enumerate(parts):
                if val == key_column_label:
                    list2_key_column_index = i
                    
            if list2_key_column_index == -1:
                print("the key column \"" + key_column_label + "\" must be present in both input lists") 
                sys.exit(0)

            list2_colnames = parts
            continue
            
        # make a new entry in the dict if we need one
        key_value = parts[list2_key_column_index]
        if key_value not in list1:
            list1[key_value] = {}

        # add all the other known values to the dict -- WILL OVERWRITE previously existing list values
        for i, val in enumerate(parts):
            if i != list2_key_column_index:
                
                cur_colname = list2_colnames[i]
                
                if user_include_colnames_set and cur_colname not in colnames_to_include:
                    continue
                
                if user_exclude_colnames_set and cur_colname in colnames_to_exclude:
                    continue

                recorded_colnames.add(cur_colname) # not very efficient to attempt to add it on every line

                list1[key_value][cur_colname] = val

#non_key_column_labels = list(set(list1_colnames + list2_colnames) - set([key_column_label,]))
#non_key_column_labels.sort()

recorded_colnames = list(recorded_colnames)
recorded_colnames.sort()

print(",".join([key_column_label,] + list(recorded_colnames)))
for key_value, other_columns in list1.items():
    other_values = []
#    for v in non_key_column_labels:
    for v in recorded_colnames:
        if v in other_columns:
            other_values.append(other_columns[v])
        else:
            other_values.append("")
    print(",".join([key_value,] + other_values))
