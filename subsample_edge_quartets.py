#!/usr/bin/env python

import argparse, newick3, os, phylo3, random, shutil, subprocess, sys, time
from multiprocessing import Lock, Manager, Pool, Queue

description = """Consider each node in the rooted tree to identify a bipartition, which is represented in
the tree as the outgoing edge connecting the node to its parent. In a fully bifurcating tree,
each node connected to this edge (the child and the parent) will have two other connected edges.
If we unroot the tree and consider the outgoing direction to be away from the bipartition of
interest, then leaf sets descendant from these outgoing edges form the four sets of the quartet
induced by the bipartition. We perform a number of replicated tests consisting of randomly drawing
one taxon from each of these four quartets, reconstructing the topology for the four random tips
using the sequence data with sequence data from the original alignment, and we record the topology
for each replicate. The resulting topology sets are used to calculate the ICA score for the 
bipartition."""

DEFAULT_RAXML = "raxmlHPC-AVX"
SECONDS_PER_MINUTE = 60
MINUTES_PER_HOUR = 60
SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR

def process_replicate(replicate):

    os.chdir(temp_wd)

    # just alias dictionary elements for convenience
    queue = replicate["queue"]
    lock = replicate["lock"]
    using_partitions = replicate["using_partitions"]
    node_id = replicate["node_id"]
    replicate_id = replicate["replicate_id"]
    raxml_path = replicate["raxml_path"]
    seqs = replicate["seqs"]

    result = {}

    # generate a label that will be unique within this run (but probably not among runs!)
    unique_label = node_id + "." + replicate_id
        
    # generate labels for temp files
    if using_partitions:
        # make a copy of the partitions file
        temp_part_fname = "temp_parts." + unique_label
        subprocess.call("cp temp_parts " + temp_part_fname, shell=True)
    
    temp_aln_fname = "temp_inseqs." + unique_label
    temp_aln_test_read_label = "temp_read_aln." + unique_label
    temp_ml_search_label = "temp_tree_search." + unique_label

    # write the alignment
    with open(temp_aln_fname,"w") as outfile:
        outfile.write("4 " + str(len(seqs[seqs.keys()[0]])) + "\n")
        for subtree_name, seq in seqs.iteritems():
            outfile.write(subtree_name + " " + seq + "\n")
    
    # test alignment readability by raxml, also filters entirely missing columns
    raxml_args = [raxml_path, 
        "-s", temp_aln_fname, 
        "-n", temp_aln_test_read_label, 
        "-m", "GTRCAT", 
        "-f", "c",
        "-$",  # silent alignment validation mode, currently on chinchliff branch
#        "--silent", # silent alignment validation mode, waiting for standard-raxml to work
        ">", "/dev/null"]
    if using_partitions:
        raxml_args.append("-q")
        raxml_args.append(temp_part_fname)

    subprocess.call(" ".join(raxml_args), shell=True)
    
    if os.path.exists(temp_aln_fname + ".reduced"):
        temp_aln_fname = temp_aln_fname + ".reduced"
        
    if using_partitions and os.path.exists(temp_part_fname + ".reduced"):
        temp_part_fname = temp_part_fname + ".reduced"

    # do the treesearch using the filtered data
    raxml_args = [raxml_path, \
        "-s", temp_aln_fname, \
        "-n", temp_ml_search_label, \
        "-m", "GTRCAT", \
        "-p", "123", \
        "-F",
        ">", "/dev/null"]
    if using_partitions:
        raxml_args.append("-q")
        raxml_args.append(temp_part_fname)

    result["raxml_args"] = " ".join(raxml_args)
    subprocess.call(result["raxml_args"], shell=True)

    result["label"] = unique_label
    queue.put(result)
    
    # increment counter and update user feedback
    replicate["n_completed"].value += 1
    lock.acquire()
    sys.stdout.write("\r"+str(replicate["n_completed"].value) + " / " + str(nreps))
    sys.stdout.flush()
    lock.release()

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-t", "--tree", type=file, nargs=1, required=True, help="The input tree. Must be rooted and fully bifurcating.")

    parser.add_argument("-n", "--alignment", type=file, nargs=1, required=True, help="Alignment file in \"relaxed phylip\" format, as used by RAxML.")

    parser.add_argument("-#", "--number-of-reps", type=int, nargs=1, required=True, help="The number of replicate quartet topology searches to be performed at each node.")

    parser.add_argument("-T", "--number-of-threads", type=int, nargs=1, required=True, help="The number of parallel threads to be used for quartet topology searches.")
    
    parser.add_argument("-q", "--partitions", type=file, nargs=1, help="Partitions file in RAxML format. If omitted then the entire alignment will be treated as one partition for all quartet replicate topology searches.")

    parser.add_argument("-o", "--results-dir", type=os.path.expanduser, nargs=1, help="A directory to which output files will be saved. If not supplied, the current working directory will be used.")

    parser.add_argument("-e", "--temp-dir", type=os.path.expanduser, nargs=1, help="A directory to which temporary files will be saved. If not supplied, a \"temp\" directory will be created in the current working directory.")
    
    parser.add_argument("-g", "--topology-sets-dir", type=os.path.expanduser, nargs=1, help="A directory to which topology sets will be saved. If not supplied, a directory will be created inside in the temp dir")

    parser.add_argument("-s", "--start-node-number", type=int, nargs=1, help="An integer denoting the node to which to start from. Nodes will be read from topologically identical (including isomorphism!) input trees in deterministic order, so this argument may be used to restart at an intermediate position (in case the previous run was canceled before completion, for example).")
    
    parser.add_argument("-v", "--verbose", action="store_true", help="Provide more verbose output if specified.")
    
    parser.add_argument("-X", "--raxml-executable", nargs=1, help="The name (or absolute path) of the NON-PTHREADS raxml executable to be used for inferring quartet topology replicates. If this argument is not supplied, then the name '"+ DEFAULT_RAXML + "' will be used. IMPORTANT NOTE: using a raxml version with silent alignment validation (i.e. which supports the `--silent` argument) is likely to drastically improve runtimes. The latest version from http://github.com/stamatak/standard-RAxML has this feature.")

    args = parser.parse_args()
    
    results_dir = os.path.abspath(args.results_dir[0]) if args.results_dir != None else os.path.abspath(".")
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    
    tree_result_file_path = results_dir + "/RESULT.labeled.tre"
    score_result_file_path = results_dir + "/node_scores.csv"

    temp_wd = args.temp_dir[0] if args.temp_dir != None else os.path.abspath("./temp")
    if not os.path.exists(temp_wd):
        os.mkdir(temp_wd)

    topology_dir = args.topology_sets_dir[0] if args.topology_sets_dir != None else temp_wd + "/topology_sets"
    if os.path.exists(topology_dir):
        shutil.rmtree(topology_dir)
    os.mkdir(topology_dir)

    calc_start_k = args.start_node_number if args.start_node_number != None else 1

    using_partitions = False
    if args.partitions is not None:
        using_partitions = True
        parts_file = args.partitions[0]
    
    raxml_path = args.raxml_executable if args.raxml_executable is not None else DEFAULT_RAXML
    
    nprocs = args.number_of_threads[0]
    nreps = args.number_of_reps[0]

    # shared object access for multithreading
    manager = Manager()
    lock = manager.Lock()

    # read the alignment into a dict, assumes phylip format with seqs unbroken on lines
    aln = {}
    alnfile = args.alignment[0]
    print "reading alignment from " + alnfile.name
    firstline = True
    for line in alnfile:
        if firstline:
            firstline = False
            continue
            
        parts = line.split()
        if len(parts) > 1:
            aln[parts[0]] = parts[1]
    args.alignment[0].close() 

    # get the tree to subsample
    tree = None
    treefile = args.tree[0]
    print "reading tree from " + treefile.name
    line = None
    while line != "":
        line = treefile.readline()
        try:
            tree = newick3.parse(line)
            break
        except AttributeError:
            pass
    if tree == None:
        sys.exit("Could not find a tree in the treefile: " + treefile.name)
    args.tree[0].close()
    leaves = tree.leaves()

    if args.verbose:
        print("tree has " + str(len(leaves)) + " leaves")

    os.chdir(temp_wd)

    # k is the node counter
    k = 1

    # if we are starting at the beginning, initialize the results file (otherwise assume it's already there and don't overwrite it) 
    if not calc_start_k > k:
        with open(score_result_file_path, "w") as resultsfile:
            resultsfile.write("node_label,obs_freq_of_test_bipart,ica\n")

    # process the nodes in the tree
    starttime = time.time()
    root_bipart_label = None
    for node in tree.iternodes():

        # skip tips and root
        if node.istip or node.parent == None:
            if (args.verbose):
                print("\nskipping " + ("tip " + node.label if node.label != None else "root"))
            continue

        # record the node label in the tree, these are required for user to match scores with corresponding branches
        node.label = str(k)

        if (k < calc_start_k):
            # skip nodes if we have a specified start node (i.e. not the root) and we haven't hit it yet
            k += 1
            continue
            
        else:
            # provide user feedback before incrementing            
            if k > calc_start_k:
                mean_time_secs = (time.time() - starttime) / float(k - calc_start_k)
                if mean_time_secs > 60:
                    if mean_time_secs > SECONDS_PER_HOUR: # more than one hour (yikes!) 
                        mean_time_units = "hours"
                        mean_time_scalar = SECONDS_PER_HOUR
                    else: # between 1 and 60 minutes
                        mean_time_units = "minutes"
                        mean_time_scalar = SECONDS_PER_MINUTE                    
                else: # less than 60 seconds
                    mean_time_units = "seconds"
                    mean_time_scalar = 1
                 
                # adjust for the duplicate bipart at the root (until we hit it, then stop adjusting)
                adj = -1 if root_bipart_label == None else 0
                est_remaining_time_secs = mean_time_secs * (len(leaves) - k + adj)
                 
                if est_remaining_time_secs > SECONDS_PER_MINUTE:
                    if est_remaining_time_secs > SECONDS_PER_HOUR:
                        est_remaining_time_units = "hours"
                        est_remaining_time_scalar = SECONDS_PER_HOUR
                    else: # between 1 and 60 minutes
                        est_remaining_time_units = "minutes"
                        est_remaining_time_scalar = SECONDS_PER_MINUTE                    
                else: # less than 60 seconds
                    est_remaining_time_units = "seconds"
                    est_remaining_time_scalar = 1
                
                time_string = " | average node time {:.2} {}".format(mean_time_secs / mean_time_scalar, mean_time_units) + \
                    " | est. remaining time {:.2f} {}".format(est_remaining_time_secs / est_remaining_time_scalar, est_remaining_time_units)
            else:
                time_string = ""
            print("\nprocessing node " + str(k) + time_string)
            k += 1

        # debug code
#        for i, child in enumerate(node.children):
#            print(" child " + str(i) + " [" + ", ".join([l.label for l in child.leaves()[0:10]]) + "]" + (" + " + (str(len(child.leaves())-10) + " more") if (len(child.leaves())-10) > 0 else ""))

        # require a bifurcating tree
#        assert(len(node.children) == 2)
        if len(node.children) != 2:
            print("Node %s does not have exactly 2 children. It will be skipped." % k)
            continue 

        # get leaf sets for the four connected subtrees
        leafsets = {}

        # two daughter subtrees
        leafsets["R1"] = set(node.children[0].label if node.istip else [l.label for l in node.children[0].leaves()])
        leafsets["R2"] = set(node.children[1].label if node.istip else [l.label for l in node.children[1].leaves()])

        # sibling/parent subtrees
        is_other_side_of_root = False # used when we hit the root for the second time
        skip_tip_child_of_root = False # used when one of the children of the root node is a tip
        tip_child_label = None
        for sib in node.parent.children:
            if sib != node:

                # if one of the subtrees is the root, skip over it
                if len(sib.leaves()) + len(node.leaves()) == len(leaves):

                    # if we already processed this bipart (on other side of the root), don't do it again
                    if (root_bipart_label != None):
                        is_other_side_of_root = True
                        break

                    # get the subtrees opposite the root
                    if len(sib.children) == 2:
                        leafsets["L1"] = set(sib.children[0].label if sib.children[0].istip else [l.label for l in sib.children[0].leaves()])
                        leafsets["L2"] = set(sib.children[1].label if sib.children[1].istip else [l.label for l in sib.children[1].leaves()])
                    elif len(sib.children) == 0:
                        skip_tip_child_of_root = True
                        tip_child_label = sib.label
                    else:
                        print("Node %s does not have exactly 2 children. It will be skipped." % k)
                        continue

                    # remember that we've already done the root, so we can skip it when we hit the other side
                    root_bipart_label = node.label

                # otherwise not at root, all connected subtrees have children
                else:

                    # sibling subtree
                    leafsets["L1"] = set([l.label for l in sib.leaves()])

                    # the rest of the tree
                    leafsets["L2"] = set()
                    for label in [l.label for l in leaves]:
                        if label not in leafsets["R1"] and \
                            label not in leafsets["R2"] and \
                            label not in leafsets["L1"]:
                                leafsets["L2"].add(label)
        
        if skip_tip_child_of_root:
            print("not calculating ica for tip child '" + tip_child_label + "' of the root (ica is 1.0, as for all tips).")
            continue

        # if we already processed the bipart at the root and this is the other side of that
        if is_other_side_of_root:
            print("\nskipping second instance of root-adjacent bipartition (it was already processed at node " + \
                    root_bipart_label + ").")
            node.label = root_bipart_label
            continue

        # sanity check
        t = set()
        for leafset in leafsets.itervalues():
            assert len(leafset) > 0
            t.update(leafset)
        print("t: " + ",".join(sorted(list(t))))
        print("leaves: " + ",".join(sorted([l.label for l in leaves])))
        assert len(t) == len(leaves)
        del(t)

        # randomly subsample an exemplar tip from each subtree
        replicates = []
        n_completed = manager.Value("i", 0, "lock")
        results_queue = manager.Queue()
        for j in range(nreps):

            rep = {}
            rep["queue"] = results_queue
            rep["lock"] = lock
            rep["using_partitions"] = using_partitions
            rep["node_id"] = node.label
            rep["replicate_id"] = str(j)
            rep["n_completed"] = n_completed
            rep["raxml_path"] = raxml_path
            rep["seqs"] = {}
            for subtree_name, leaf_names in leafsets.iteritems():
                while subtree_name not in rep["seqs"]: 
                    leafname = random.sample(leaf_names, 1)[0] #[0].label
                    if leafname in aln:
                        rep["seqs"][subtree_name] = aln[leafname]
                        if args.verbose:
                            print("using exemplar " + leafname + " for " + subtree_name)
                    else:
                        print("\nWARNING: name " + leafname + " not in alignment")

            replicates.append(rep)

        # clear any lingering files (e.g. from previous runs) that could interfere with raxml
        # redirecting stderr because it prints a bunch of failed calls not sure why as the
        # command seems to be working as expected...
        subprocess.call("rm *." + node.label + ".* 2> /dev/null", shell=True)

        # copy in original partitions file, should not change throughout run
        if using_partitions:
            subprocess.call("cp " + parts_file_path + " temp_parts", shell=True)

        # run the raxml calls in parallel
        # now designate multiprocessing resource pool.
        # important to do this outside the node loop as regular garbage collecting does not seem
        # to apply to the threads! also, set maxtasksperchild to release memory and files!
        pool = Pool(processes=nprocs, maxtasksperchild=1)
        pool.map(process_replicate, replicates)
        pool.close()
        pool.join()
        del(pool)
        ### map(process_replicate, replicates) # use for testing
        print("")

        # now proces the results. first open a file to hold topologies
        count_bipart_observed = 0
        topo_file_name = topology_dir + "/" + node.label + ".obs_topologies.txt"
        with open(topo_file_name, "w") as topo_file:

            while not results_queue.empty():
                result = results_queue.get()

                # attempt to open the raxml result
                raxml_result_tree_file_path = "RAxML_result.temp_tree_search." + result["label"]
                if os.path.exists(raxml_result_tree_file_path):
                    with open(raxml_result_tree_file_path, "r") as result:
                        result_tree = newick3.parse(result.readline())
                else:
                    print("WARNING: raxml did not complete successfully. The failed command was:\n\n" + result["raxml_args"] + "\n")
                    continue

                # write the result topology to the set of observed topologies for this node
                topo_file.write(newick3.to_string(result_tree)+";\n")

                # check to see if it contains the original bipartition
                for testnode in result_tree.iternodes():

                    # skip tips and root
                    if testnode.istip or testnode.parent == None:
                        continue

                    # assumes root bipart is not in pxbp output, so we can just check included names (excluded are implied)
                    testnode_leafnames = [l.label for l in testnode.leaves()]
                    if ("L1" in testnode_leafnames and "L2" in testnode_leafnames) or \
                            ("R1" in testnode_leafnames and "R2" in testnode_leafnames):
                        count_bipart_observed += 1

        # get the ICA score from phyx
        pxbp_outfile = "temp_pxbp_out." + node.label
        pxbp_args = ["pxbp", "-t", topo_file_name, ">", pxbp_outfile]
        subprocess.call(" ".join(pxbp_args),shell=True)
        ica = "-1" # if we don't find a score in the pxbp output then this bipart is never observed
        with open(pxbp_outfile,"r") as pxbp_result:
            for line in pxbp_result:
                parts = line.split()
                if len(parts) > 2:

                    # tests assume that root bipart is not in list, thus all biparts have exactly two names (the ingroup)
                    names = set(parts[0:2]) 
                    if ("L1" in names and "L2" in names) or ("R1" in names and "R2" in names):
                        ica = parts[-1] # ica score should be last item on line
                        break
    
        freq = str(count_bipart_observed / float(nreps))

        # write the scores to the file
        with open(score_result_file_path, "a") as results_file:
            results_file.write(",".join([node.label, freq, ica]) + "\n")

        # write the tree with all processed nodes labeled
        with open(tree_result_file_path,"w") as tree_file_path:
            tree_file_path.write(newick3.to_string(tree)+";")

        # clean up
        del(results_queue)
        del(n_completed)

#        del(n_completed)
#        del(count_bipart_observed)
        # redirecting stderr because it prints a bunch of failed calls. not sure why as the command seems to be working as expected...
        subprocess.call("rm *." + node.label + ".* 2> /dev/null", shell=True)
    
    print("\ndone.\nscores written to: " + score_result_file_path + \
        "\nlabeled tree written to: " + tree_result_file_path + \
        "\ntotal time {:.2f}".format((time.time() - starttime) / 60 / 60) + " hours")
