#!/home/davidec/miniconda3/envs/birc-project/bin/python
import msprime, tskit, pyslim
import numpy as np


# with open("params.txt") as f: 
#     for line in f:
#         param = line.strip().split(" ")
#         if "MUT_RATE_AMP" in param[1]:
#             Mut_rate = float(param[1].split("=")[1])
#             break

def ts_processer(ts_file_path, recapitation=False): ################################# numinds
    orig_ts = tskit.load(ts_file_path)
    if recapitation == False:
        rts = orig_ts
    else: #no recapitation needed right now 
        return "write function without recapitation"
    #Simplification: discarding less informing nodes from tree
    rng = np.random.default_rng(seed=3)
    alive_inds = pyslim.individuals_alive_at(rts, 0)
    keep_indivs = rng.choice(alive_inds, 100, replace=False) #################################
    keep_nodes = []
    for i in keep_indivs:
        keep_nodes.extend(rts.individual(i).nodes)

    sts = rts.simplify(keep_nodes, keep_input_roots=True) 

    print(f"SIMPLIFICATION:\nBefore, there were {rts.num_samples} sample nodes (and {rts.num_individuals} individuals)\n"
        f"in the tree sequence, and now there are {sts.num_samples} sample nodes\n"
        f"(and {sts.num_individuals} individuals).\n--------------------------\n")
    
    #Adding neutral mutation types after the simulation has occured
    #ACTUALLY WE HAVE ONLY NEUTRAL MUTATION TYPES IN OUR SIMULATION <====== DAVIDE
    next_id = pyslim.next_slim_mutation_id(sts)
    ts = msprime.sim_mutations(
            sts,
            rate=Mut_rate,
            model=msprime.SLiMMutationModel(type=0, next_id=next_id), #using mutation type 0 because its not used in SLiM 
            keep=True, #keeps other mut types 
    )

    print(f"ADDING NEUTRAL MUTATIONS:\nThe tree sequence now has {ts.num_mutations} mutations,\n"
        f"and mean pairwise nucleotide diversity is {ts.diversity():0.3e}.\n--------------------------\n")

    return ts

import sys # Kasper
_, slim_tree_file, processed_tree_file = sys.argv # Kasper
import re # Kasper
Mut_rate = float(re.search(r'u_([^_]+)', slim_tree_file).group(1)) # Kasper

print(Mut_rate)

processed_ts = ts_processer(slim_tree_file) # Kasper

processed_ts.dump(processed_tree_file) # Kasper
