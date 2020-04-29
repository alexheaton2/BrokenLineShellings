def search(V, vFav, pivots, limit, misses, w):
    # INPUTS:
    #     "V" is a list of tuples. Each tuple is a vertex corresponding to a basis element.
    #         e.g. V = [[0, 0, 1, 1], [0, 1, 0, 1], [0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 1, 0], [1, 1, 0, 0]]
    #     "vFav" is your favorite vertex, whose normal cone all linear functionals will reside.
    #         e.g. vFav = [1,0,0,1]
    #     "pivots" is a list of positions where you would like to pivot, the first possible choice is "1" not "0"
    #         e.g. pivots = [1,3,4] but not pivots = [0,3,4]
    #     "limit" is the maximum number of new \ell(v_{j+1}) to collect at each pivot location
    #     "misses" is the number of consecutive failed attempts to get a new, valid \ell(v_{j+1}) before you give up
    #     "w" is the relative weight of \ell(v_j) compared to the convex combination of all v_i - v_1
    #         -higher value of "w" means a more localized search
    #         -a value of "w=0" means you destroy the dependency of \ell(v_{j+1}) on the previous \ell(v_j)
    import random
    if len(pivots) > 0:
        if pivots[0]==0:
            print "Please don't pivot at 0. The first place you can pivot is at 1."
            return 0
    VT = [tuple(v) for v in V]
    vFav = tuple(vFav)
    import time; t0 = time.time()
    normals = get_normals(VT); t1 = time.time();
    L0 = get_L_minimizing_new(vFav, VT)
    Ls = {} # better name for the linear fncls dictionary, whose keys are "sign patterns" and values are "linear fncls"
    t2 = time.time()
    sweeps, Ls = get_sweeps_new(vFav,L0,Ls,pivots,VT,normals,limit,misses,w)
    t3 = time.time(); print "time to find all the distinct, valid sweeps:   {}".format(t3-t2)
    print "we have {} sweeps total.".format(len(sweeps)); print;
    import time; t1 = time.time();
    distinct_sweeps, distinct_posets, distinct_tables = get_distinct_sweeps(sweeps, VT, Ls)
    t2 = time.time(); print "time needed to sort only the distinct posets from all sweeps:   {}".format(t2-t1);
    print "Of these {} sweeps, only {} gave rise to non-isomorphic posets".format( len(sweeps), len(distinct_posets) ); print;
    return [distinct_sweeps, distinct_posets, distinct_tables]

def update_search(result, V, vFav, pivots, limit, misses, w):
    # "result" should be a list [distinct_sweeps, distinct_posets, distinct_tables]
    import random
    if pivots[0]==0:
        print "Please don't pivot at 0. The first place you can pivot is at 1."
        return 0
    VT = [tuple(v) for v in V]
    vFav = tuple(vFav)
    import time;
    t0 = time.time()
    normals = get_normals(VT); t1 = time.time();
    L0 = get_L_minimizing_new(vFav, VT)
    Ls = {}
    t2 = time.time()
    sweeps, Ls = get_sweeps_new(vFav,L0,Ls,pivots,VT,normals,limit,misses,w)
    t3 = time.time(); print "time to find all the distinct, valid sweeps:   {}".format(t3-t2)
    print "we have {} sweeps total.".format(len(sweeps)); print;
    import time; t1 = time.time();
    distinct_sweeps, distinct_posets, distinct_tables = get_distinct_sweeps(sweeps, VT, Ls)
    distinct_sweeps, distinct_posets, distinct_tables = update_distinct_sweeps(result, [distinct_sweeps, distinct_posets, distinct_tables], VT, Ls)
    t2 = time.time(); print "time needed to sort only the distinct posets from all sweeps:   {}".format(t2-t1);
    print "Of these {} sweeps, only {} gave rise to non-isomorphic posets".format( len(sweeps), len(distinct_posets) ); print;
    return [distinct_sweeps, distinct_posets, distinct_tables]

def display_results(result, items=[]):
    # "result" should be a list [distinct_sweeps, distinct_posets, distinct_tables], output from search_new(...)
    # "items" should be a list of which posets you want to print.
    #         if there are 14 posets, you may not want to print them all, but rather just the first three...
    #         then set "items = [0,1,2]" and this prints just the first 3 out of how many there are
    if len(items)==0:
        items = range(len(result[1])) # the number of distinct_posets
    for i in items:
        if i >= len(result[1]): # user asked for too many posets, more than exist
            print "There are only {} posets. Set items=[0,1,...] until something less than the number of posets available.".format(len(result[1]))
            break
        print "Below is poset number {} and its associated data.".format(i)
        poset = result[1][i] # the ith distinct poset
        M = len(poset) / 4
        HD = poset.hasse_diagram()
        plt = plot(HD, vertex_labels=True, layout="acyclic", figsize=[M,4*M])
        show(plt)
        # now print the ith distinct table, associated to the previously printed poset
        print result[2][i]; print; print;
    return

def update_distinct_sweeps(old_result, new_result, VT, Ls):
    new_sweeps = old_result[0]
    new_posets = old_result[1]
    new_tables = old_result[2]
    sweeps = new_result[0] # new sweeps to check against old sweeps
    for i,sweep in enumerate(sweeps):
        #new_table, new_poset = get_table_and_poset_from_walk(walk, VT, Ls)
        new_table = new_result[2][i]
        new_poset = new_result[1][i]
        new = True
        for old_poset in new_posets:
            if old_poset.is_isomorphic(new_poset):
                new = False
        if new:
            # add them at the same time, so ith entry of "new_walks" corresponds to ith entry of "new_posets"
            new_posets.append(new_poset)
            new_sweeps.append(sweep)
            new_tables.append(new_table)
    return new_sweeps, new_posets, new_tables

def get_valid_Ls_new(limit, misses, w, vFav, VT, normals, Lold, Ls, step):
    neighbors_signs = []
    since_found = 0
    keep_going = True
    total_count = 0
    while keep_going:
        since_found += 1
        total_count += 1
        if since_found > misses:
            keep_going = False
        if len(neighbors_signs) < 1:
            keep_going = True
        if len(neighbors_signs) > limit:
            keep_going = False
        if total_count > misses:
            keep_going = False
        l = get_convex_comb_minimizing_check(w, Ls[Lold], vFav, VT)
        signsl = get_sign_pattern(l,normals)
        if signsl != Lold: # not the original region
            if signsl not in neighbors_signs: # not a region we already found
                if same_cut(l, Ls[Lold], step, VT):
                    since_found = 0
                    Ls[signsl] = l
                    neighbors_signs.append(signsl)
    return neighbors_signs, Ls

def get_convex_comb_minimizing_check(w, v_old, vFav, VT):
    # "w" is how much MORE of \ell(v_j) should be there versus \ell(v_{j+1})
    # "v_old" is Ls[Lold] a linear functional vector used previously
    import random
    n = len(VT)
    found = False # make sure our new linear functional minimizes "vFav"
    while not found:
        cs = get_convex_coeffs(n-1)
        vecs = [vector(v) - vector(vFav) for v in VT if v!=vFav]
        L = sum([cs[i]*vecs[i] for i in range(n-1)])
        # now we normalize and add in "v_old"
        v_old = v_old / norm(v_old) # unit vector
        L = L / norm(L) # unit vector
        Lnew = w*v_old + L
        indices = vertex_index_order(VT,Lnew)
        indFav = VT.index(vFav)
        if indices[0] == indFav:
            found = True
            return Lnew

def get_convex_coeffs(N):
    import random
    cs = [random.random() for i in range(N)]
    C = sum(cs)
    cs = [ci/C for ci in cs]
    return cs

def get_L_minimizing_new(vFav, VT):
    # OUTPUT: "L" a generic linear functional, minimizing "vFav"
    #         "L" should be perpendicular to the all ones vector, since the entire matroid polytope lives in that plane
    # INPUT: "vFav" is a tuple, the desired vertex. Should be in "VT"
    #         "VT" is a list of tuples, the vertices
    dim = len(VT[0]) # ambient dimension
    found = False
    T = SphericalDistribution(dimension=dim)
    while not found:
        L = T.get_random_element()
        # project ALONG the all ones vector
        L = vector(L) - (vector(L)*vector(RDF,[1/sqrt(dim)]*dim))*vector(RDF,[1/sqrt(dim)]*dim)
        indices = vertex_index_order(VT,L)
        indFav = VT.index(vFav)
        if indices[0] == indFav:
            found = True
            return L

def get_sweeps_new(vFav,L0,Ls,pivots,VT,normals,limit,misses,w):
    # OUTPUT: "sweeps" a list of lists, which are each a "sweep"
    #             a "sweep" is a list of sign patterns "signs", of length = len(normals) <= binomial(|VT|,2)
    #                         each "signs" is a tuple (-1,-1,1,-1,1,1,...) recording sign of dot products with each normal
    #                         These are keys to the GLOBAL dictionary "Ls" whose keys are "signs" and values are linear functionals "L"
    signsL0 = get_sign_pattern(L0,normals)
    Ls[signsL0] = L0
    limit = N
    sweepsBuild = [ [] for i in range(len(VT)+1) ] # sweepsBuild[-1] will be our final answer "sweeps"
    sweepsBuild[1] = [[signsL0]] #sweepsBuild[0] = [signsL0]
    for i in range(1, len(VT)):
        if i not in pivots:
            # just repeat the last element again
            for sweep in sweepsBuild[i]:
                Llast = sweep[-1]
                new = copy(sweep)
                new.append(Llast)
                sweepsBuild[i+1].append(new)
        if i in pivots:
            # for each old branch, create a bunch of new branches...
            for sweep in sweepsBuild[i]:
                Lold = sweep[-1]
                neighbors, Ls = get_valid_Ls_new(limit, misses, w, vFav, VT, normals, Lold, Ls, i)
                neighbors.append(Lold)
                for Lnew in neighbors:
                    new = copy(sweep)
                    new.append(Lnew)
                    sweepsBuild[i+1].append(new)
    sweeps = sweepsBuild[len(VT)] # the last list of lists is our final answer!
    return sweeps, Ls

def get_distinct_sweeps(walks, VT, signPatterns):
    # ADDED this fcn on January 22, 2020
    # sort and find only distinct walks first, before creating all the tables and data and plots
    # OUTPUT: "new_walks" is a list of walks that induce distinct non-isomorphic posets
    #         "new_posets" is a list of posets. The ith entry of "new_posets" corresponds to the ith entry in "new_walks"
    #         "new_tables" is a list of tables. The ith entry of "new_tables" corresponds to the ith entry in "new_walks" and "new_posets"
    # INPUT: "walks" is a list of tuples, each tuple is a key to the "signPatterns" dictionary
    #        "VT" is a list of tuples, each tuple is a vertex, zeros and ones
    #        "signPatterns" is a dictionary whose keys are tuples, each tuple has ones and -ones,
    #          telling signPattern of dot product against each normal vector
    new_walks = []
    new_posets = []
    new_tables = []
    for walk in walks:
        new_table, new_poset = get_table_and_poset_from_walk(walk, VT, signPatterns)
        new = True
        for old_poset in new_posets:
            if old_poset.is_isomorphic(new_poset):
                new = False
        if new:
            # add them at the same time, so ith entry of "new_walks" corresponds to ith entry of "new_posets"
            new_posets.append(new_poset)
            new_walks.append(walk)
            new_tables.append(new_table)
    return new_walks, new_posets, new_tables

def get_table_and_poset_from_walk(walk, VT, signPatterns):
    # ADDED this fcn on January 22, 2020
    # this does NOT create the plots yet. We create those last for only the distinct posets created.
    table_rows = [["vertex", "order swept", "IP set", "linear fncl"]]
    nodes = {}
    for i,step in enumerate(walk):
        l = signPatterns[step] # get the current linear functional at this step in the walk
        indices = vertex_index_order(VT,l)
        vertexInd = indices[i] # get the current vertex
        node = IP( VT[vertexInd], l, VT )
        nodes[i] = node
        fcnl = signPatterns[step].n(digits=3)
        row = [VT[vertexInd], i, list(node), fcnl]
        table_rows.append(row)
    Table2 = table(rows=table_rows, header_row=True)
    #
    poset_relation = lambda A,B : nodes[A].issubset(nodes[B])
    poset = Poset( (nodes,poset_relation) )
    return Table2, poset

def find_indices(v):
    # "v" is a list of zeros and ones
    # OUTPUT two lists, list of indices of ones, then list of indices of zeros
    ones_indices = []
    zeros_indices = []
    for i,v_i in enumerate(v):
        if v_i==1:
            ones_indices.append(i)
        elif v_i==0:
            zeros_indices.append(i)
        else:
            print "something went wrong... there was something not 0, not 1, in our v_B coordinate vector."
    return ones_indices, zeros_indices

def switch(v, one_index, zero_index):
    # switches b with b' to see if the resulting vector is actually a vertex of our polytope
    u = list(v)
    u[one_index]=0
    u[zero_index]=1
    u = tuple(u)
    return u

def IP(v,l,vertices):
    # OUTPUT: a frozenset, the IPset created at the vertex "v" with linear functional "l"
    # "v" is a tuple like (1,0,0,1,1)
    # "l" is a linear functional, a list like [3,2,1,4,5]
    # "vertices" is a list of tuples like [(1,0,0,1,1),(1,1,0,1,0),...], sometimes I call it "VT" for "vertex tuples"
    ones_indices, zeros_indices = find_indices(v)
    IPset = []
    for index_b in ones_indices:
        # decide to put it in IPset or not
        good = False
        # but if we find even one b' with l(b') < l(b) and B \setminus b \cup b' in \mathcal{B}, then change good=True
        for index_bp in zeros_indices:
            # if l(b') < l(b) becomes
            if l[index_bp] < l[index_b]:
                # and if B \setminus b \cup b' in \mathcal{B}
                potential_v = switch(v, index_b, index_bp)
                if potential_v in vertices:
                    good = True
        if good:
            IPset.append(index_b)
    IPset = frozenset( IPset )
    return IPset

def same_cut(L1,L2, step, VT):
    # "L1" and "L2" are two different linear functionals
    # "step" is where to check if they induce the same cut
    #     ...check that the elements of indices[0:step] are the same as SETS for "L1" and for "L2"
    indices1 = vertex_index_order(VT,L1)
    indices2 = vertex_index_order(VT,L2)
    okay = True
    for i in range(step):
        if indices1[i] not in indices2[0:step]:
            okay = False
    return okay

def is_minimal(v0,VT,l):
    # "v0" is the vertex we want all our linear functionals to minimize
    # "VT" is a list of tuples. Each tuple is a vertex of the matroid polytope.
    # "l" is a linear functional
    v0 = tuple(v0)
    v0_index = VT.index(v0)
    indices = vertex_index_order(VT,l)
    if indices[0] == v0_index:
        return True
    else:
        return False

def get_sign_pattern(l,normals):
    # "l" is a unit vector sampled from the sphere
    # "normals" are all the normal vectors from hyperplanes agreeing on two vertices of our matroid polytope.
    #     ...the sign patterns arising from dot products against these normals determine which subcone our "l" belongs to.
    # OUTPUT: "signs" will be a tuple of {1,-1}^(binomial(|V|,2) = len(normals))
    dotProducts = [vector(l)*vector(L) for L in normals]
    signs = tuple([sign(dot) for dot in dotProducts])
    return signs

def get_normals(VT):
    # INPUT:     let "VL" be P.vertices_list(), namely a list of lists like [[0, 1, 1], [1, 0, 1], [1, 1, 0]]
    # OUTPUT:    a list of normal vectors, ready to be dot-producted with points on the unit sphere
    VT = [tuple(v) for v in VT] # make sure they are tuples
    choose2 = [s.list() for s in Subsets(VT,2)]
    normals = [] # now create the normal vectors defining the planes orthogonal to "final - initial" for two vertices "v1 - v2"
    for verts in choose2:
        v1, v2 = verts
        normal = vector(v1) - vector(v2)
        if normal not in normals and -normal not in normals:
            normals.append(normal)
    return normals

def vertex_index_order(VT,L):
    # INPUT:     "L" is a linear functional, a vector
    #            "VT" is a list of tuples. Each tuple is a vertex of the matroid polytope.
    # OUTPUT:  "indices" is a list of indices of vertices, ordered lowest to highest by "L"
    L = vector(L) # make sure its a vector
    indices_and_vals = [(i,L*vector(v)) for (i,v) in enumerate(VT)]
    ans = sorted(indices_and_vals, key=lambda x: x[1]) # sort by the "value" of "L" acting on "v" above.
    indices = [tup[0] for tup in ans] # return a list of the indices only, discard the "values"
    return indices
