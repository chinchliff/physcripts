from copy import deepcopy as copy

#import sets
PREORDER = -99999; POSTORDER = 123456
BRANCHLENGTH = 0; INTERNODES = 1

class Node:
    def __init__(self):
        self.data = {}
        self.isroot = False
        self.istip = False
        self.label = None
        self.length = 0
        self.parent = None
        self.children = [] # should probably be using a set
        self.nchildren = 0
        self.excluded_dists = []
        self.comment = None

#    @property
#    def is_tip(self):
#        return self.istip
#    
#    @is_tip.setter
#    def is_tip(self, is_tip):
#        self.istip = is_tip        

    def order_subtrees_by_size(self, n2s=None, recurse=False, reverse=False):
        if n2s is None:
            n2s = node2size(self)
        if not self.istip:
            v = [ (n2s[c], c.label, c) for c in self.children ]
            v.sort()
            if reverse:
                v.reverse()
            self.children = [ x[-1] for x in v ]
            if recurse:
                for c in self.children:
                    c.order_subtrees_by_size(n2s, recurse=True, reverse=reverse)

    def add_child(self, child):
        assert child not in self.children
        self.children.append(child)
        child.parent = self
#        self.nchildren += 1
        self.nchildren = len(self.children)
        self.istip = False

    def remove_child(self, child):
        assert child in self.children
        self.children.remove(child)
        child.parent = None
#        self.nchildren -= 1
        self.nchildren = len(self.children)

##     def leaves(self, v=None):
##         if v is None:
##             v = []
##         if not self.children:
##             v.append(self)
##         else:
##             for child in self.children:
##                 child.leaves(v)
##         return v

    def leaves(self):
        return [ n for n in self.iternodes() if n.istip ]

    def get_node_for_name(self, inname):
        for node in self.iternodes():
            if node.label == inname:
                return node
        return None

    def iternodes(self, order=PREORDER, v=None):
        '''returns a list of nodes descendant from self - including self'''
        if order == PREORDER:
            yield self
        for child in self.children:
            for d in child.iternodes(order):
                yield d
        if order == POSTORDER:
            yield self

    def descendants(self, order=PREORDER, v=None):
        '''returns a list of nodes descendant from self - not including self!'''
        if v is None:
            v = []
        assert order in (PREORDER, POSTORDER)
        for child in self.children:
            if order == PREORDER:
                v.append(child)
            else:
                v.insert(0, child)
            if child.children:
                child.descendants(order, v)
        return v

    def find_descendant(self, label):
        if label == self.label:
            return self
        else:
            for child in self.children:
                n = child.find_descendant(label)
                if n:
                    return n
        return None

    def prune(self, logfile=None):
        
        p = self.parent
        if p == None:
            raise ValueError("attempt to remove final tip (or all tips) from tree")
        
        if (logfile != None):
            logfile.write("removing " + self.label + "\n")

        p.remove_child(self)
        p.nchildren = len(p.children)
        if p.nchildren == 0:
            p.istip = True
 
        root_joint = None
        while len(p.children) == 1:
            
            only_sib = p.children[0]
            
            if self.istip:

                if len(only_sib.children) > 0:
                    # if the only remaining sister node has children, add them to the parent
                    for gc in only_sib.children:
                        p.add_child(gc)
                    p.remove_child(only_sib)
                    p.length += only_sib.length
                    continue
                elif not only_sib.istip:
                    # the only remaining sib has no children and is not a tip, i.e. it is an empty subclade, so...
                    p.prune()
                    break
            
            pp = p.parent
            if pp == None:
                root_joint = p
                break
        
            pp.remove_child(p)
            pp.add_child(only_sib)
            only_sib.length += p.length
            p = pp
             
        if root_joint != None:
            while len(root_joint.children) < 2:
                # prune knuckles at the root of the tree if necessary
                only_child = root_joint.children[0]
                only_child.parent = None
                only_child.isroot = True
                root_joint = only_child
        
        return p
    
    def _calc_depth(self):
        '''recursively calculate the depth of this node'''
        if len(self.children) < 1:
            return self.length
        else:
            return self.length + max([c.depth for c in self.children])

    @property
    def depth(self):
        '''return the depth of this node in the tree'''
        # currently just calculating the depth for every call. should really store this
        # and just update it when necessary
        return self._calc_depth()

    def graft(self, node):
        parent = self.parent
        parent.remove_child(self)
        n = Node()
        n.add_child(self)
        n.add_child(node)
        parent.add_child(n)

    def leaf_distances(self, store=None, measure=BRANCHLENGTH):
        '''for each internal node, calculate the distance to each leaf, measured in branch length or internodes'''
        if store is None:
            store = {}
        leaf2len = {}
        if self.children:
            for child in self.children:
                if measure == BRANCHLENGTH:
                    assert child.length is not None
                    dist = child.length
                elif measure == INTERNODES:
                    dist = 1
                else:
                    raise "InvalidMeasure"
                child.leaf_distances(store, measure)
                if child.istip:
                    leaf2len[child.label] = dist
                else:
                    for k, v in store[child].items():
                        leaf2len[k] = v + dist
        else:
            leaf2len[self] = {self.label: 0}
        store[self] = leaf2len
        return store

    def rootpath(self):
        n = self
        while 1:
            yield n
            if n.parent:
                n = n.parent
            else:
                break
            
    def subtree_mapping(self, labels, clean=False):
        """
        find the set of nodes in 'labels', and create a new tree
        representing the subtree connecting them.  nodes are assumed to be
        non-nested.

        return value is a mapping of old nodes to new nodes and vice versa.
        """
        d = {}
        oldtips = [ x for x in self.leaves() if x.label in labels ]
        for tip in oldtips:
            path = list(tip.rootpath())
            for node in path:
                if node not in d:
                    newnode = Node()
                    newnode.istip = node.istip
                    newnode.length = node.length
                    newnode.label = node.label
                    d[node] = newnode
                    d[newnode] = node
                else:
                    newnode = d[node]

                for child in node.children:
                    if child in d:
                        newchild = d[child]
                        if newchild not in newnode.children:
                            newnode.add_child(newchild)
        d["oldroot"] = self
        d["newroot"] = d[self]
        if clean:
            n = d["newroot"]
            while 1:
                if n.nchildren == 1:
                    oldnode = d[n]
                    del d[oldnode]; del d[n]
                    child = n.children[0]
                    child.parent = None
                    child.isroot = True
                    d["newroot"] = child
                    d["oldroot"] = d[child]
                    n = child
                else:
                    break
                    
            for tip in oldtips:
                newnode = d[tip]
                while 1:
                    newnode = newnode.parent
                    oldnode = d[newnode]
                    if newnode.nchildren == 1:
                        child = newnode.children[0]
                        if newnode.length:
                            child.length += newnode.length
                        newnode.remove_child(child)
                        if newnode.parent:
                            parent = newnode.parent
                            parent.remove_child(newnode)
                            parent.add_child(child)
                        del d[oldnode]; del d[newnode]
                    if not newnode.parent:
                        break
            
        return d
        
    def get_path_to_root(self):
        path = []
        parent = self.parent
        while parent != None:
            path.append(parent);
#                if parent.parent != None:
            parent = parent.parent
#           else:
#                break
        return path
        
def get_tree_rooted_on(target_node):
	# we will reroot on branch subtending the target

	if target_node.parent is None:
		# the target is already the root. rerooting here would just add a knuckle
		return target_node

	if target_node.parent.parent is None and len(target_node.parent.children) < 3:
		# the target node is already a child of the root AND the root has only two
		# children. rerooting will produce an identical topology
		return target_node

	# get the branch length of the branch where the root will go
	# we will divide it evenly across both children of the new root
	new_bl = target_node.length / 2

	# get the original parent of the target node before it is lost
	old_parent = target_node.parent
	old_parent.remove_child(target_node)

	# this branch length will go on the next descendant branch when we
	# reverse all the parent rels to the old root 
	temp_bl = old_parent.length

	# make a new root node, attach the target_node as an immediate descendant
	new_root = Node()
	new_root.add_child(target_node)
	target_node.length = new_bl

	# get the target's original grandparent, this will be the starting point for reversing
	# the relationships of the target's ancestors to make them descendants of the new root
	parent = old_parent.parent
	
	# reattach the former target's parent as a sibling of the target node
	new_root.add_child(old_parent)
	old_parent.length = new_bl

	if parent is not None:
		# if the tree had a polytomy at the root and we rerooted on one of its child branches,
		# then we have basically just let the original root node become  a child of the new
		# root containing all the other children of the original polytomy. in this case, the
		# `parent` var will be None and the topology finished. If not, then we need to...

		# traverse the target's ancestors, reversing their rels to make them descendants of new root
		prev_child = old_parent
		while parent != None:

			# remember the grandparent
			gp = parent.parent

			# switch the direction of this relationship
			parent.remove_child(prev_child)
			prev_child.add_child(parent)
	
			next_bl = parent.length
			parent.length = temp_bl
			temp_bl = next_bl
	
			last_known_ancestor = prev_child # save this in case we hit the root node next
			prev_child = parent
			parent = gp
	
		# now the current previous_child variable should contain the old root node (node with no parent)
		# remove all its children and attach them to the last ancestor
		n = len(prev_child.children)
		for root_child in prev_child.children:
			last_known_ancestor.add_child(root_child)
			prev_child.remove_child(root_child)
		
			# finally, get rid of the old root node
			last_known_ancestor.remove_child(prev_child)
	
			root_child.length += prev_child.length / n
	
	new_root.isroot = True
	return new_root

def node2size(node, d=None):
    "map node and descendants to number of descendant tips"
    if d is None:
        d = {}
    size = int(node.istip)
    if not node.istip:
        for child in node.children:
            node2size(child, d)
            size += d[child]
    d[node] = size
    return d

def reroot(oldroot, newroot): # needs some work
    oldroot.isroot = False
    newroot.isroot = True
    v = []
    n = newroot
    while 1:
        v.append(n)
        if not n.parent: break
        n = n.parent
    #print [ x.label for x in v ]
    v.reverse()
    for i, cp in enumerate(v[:-1]):
        node = v[i+1]
        # node is current node; cp is current parent
        #print node.label, cp.label
        cp.remove_child(node)
        node.add_child(cp)
        cp.length = node.length
    return newroot

# DEPRECATED, use get_mrca()
def getMRCA(tree, innames):
	return get_mrca(tree, innames)

def get_mrca(tree, innames):

    # find a name in the tree to start from
    i = 0
    firstnode = None
    while firstnode == None:
        # if we don't match any names
        if i >= len(innames):
            return None
        firstnode = tree.get_node_for_name(innames[i])
        i += 1

    # if we matched a name but there aren't any others to check
    if i == len(innames):
        return firstnode
    
    guide_path = firstnode.get_path_to_root()
    mrca_index = 0
    
    # compare first path against the others
    for j in range(i, len(innames)):

        child_node = tree.get_node_for_name(innames[j])
        if child_node != None:
            compare_path = child_node.get_path_to_root()
 
#            print "looking for ancestor with " + innames[j] # testing
            cand_mrca_index = get_first_mrca_index_for_paths(guide_path, compare_path)
            if cand_mrca_index > mrca_index:
                mrca_index = cand_mrca_index
#                print "mrca index", mrca_index # testing

#    print len(guide_path) # testing
    return guide_path[mrca_index]
            

def get_first_mrca_index_for_paths(guide_path, compare_path):
    for comp_node in compare_path:
        for k in range(len(guide_path)):
#            print "testing ", guide_path[k]
#            print "against ", comp_node
            if guide_path[k] == comp_node:
#                print "they match "
                return k
    
#        candidate_mrca = getMRCATraverse

#def getMRCA(innames, tree):
#    mrca = None
#    if len(innames) == 1:
#        return None
#    else:
#        outgroup = []
#        for name in innames:
#            for i in range(len(tree.leaves())):
#                if tree.leaves()[i].label == name:
#                    outgroup.append(tree.leaves()[i])
#                    #print tree.leaves[i].label
#        cur2 = None
#        tempmrca = None
#        cur1 = outgroup.pop()
#        while len(outgroup)>0:
#            cur2 = outgroup.pop()
#            tempmrca = getMRCATraverse(cur1,cur2)
#            cur1 = tempmrca
#        mrca = cur1
#    return mrca

def getMRCATraverse(curn1, curn2):
    mrca = None
    #get path to root for first node
    path1 = []
    parent = curn1
    path1.append(parent)
    while parent != None:
        path1.append(parent);
        if parent.parent != None:
            parent = parent.parent
        else:
            break
    #find first match between this node and the first one
    parent = curn2
    x = True;
    while x == True:
        for i in range(len(path1)):
            if parent == path1[i]:
                mrca = parent
                x = False
                break
        parent = parent.parent
    return mrca
    
def getMRCATraverseFromPath(path1, curn2):
    mrca = None
    #find first match between this node and the first one
    parent = curn2
    x = True;
    while x == True:
        for i in range(len(path1)):
            if parent == path1[i]:
                mrca = parent
                x = False
                break
        parent = parent.parent
    return mrca   
    
