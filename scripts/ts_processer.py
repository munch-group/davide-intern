#!/home/davidec/miniconda3/envs/birc-project/bin/python
import msprime, tskit, pyslim
import numpy as np
import pandas as pd
import sys 


Mut_rate = 1e-8 #Mut rate used for adding neutral mutations

def ts_processer(ts_file_path, recapitation=False): ################################# numinds
#Function to perform simplification and adding neutral mutations
    orig_ts = tskit.load(ts_file_path)
    if recapitation == False:
        rts = orig_ts
    else: #no recapitation needed right now 
        return "write function without recapitation"
    #Simplification: discarding less informing nodes from tree
    rng = np.random.default_rng(seed=3)
    alive_inds = pyslim.individuals_alive_at(rts, 0)
    keep_indivs = rng.choice(alive_inds, 100, replace=False) 
    keep_nodes = []
    for i in keep_indivs:
        keep_nodes.extend(rts.individual(i).nodes)

    sts = rts.simplify(keep_nodes, keep_input_roots=True) 

    print(f"SIMPLIFICATION:\nBefore, there were {rts.num_samples} sample nodes (and {rts.num_individuals} individuals)\n"
        f"in the tree sequence, and now there are {sts.num_samples} sample nodes\n"
        f"(and {sts.num_individuals} individuals).\n--------------------------\n")
    
    #Adding neutral mutation types after the simulation has occured
    next_id = pyslim.next_slim_mutation_id(sts)
    ts = msprime.sim_mutations(
            sts,
            rate=Mut_rate,
            # model=msprime.SLiMMutationModel(type=0, next_id=next_id), #using mutation type 0 because its not used in SLiM 
            model=msprime.SLiMMutationModel(type=1, next_id=next_id), #using mutation type 0 because its not used in SLiM 
            # keep=True, #keeps other mut types 
            keep=False, #keeps other mut types 
    )

    print(f"ADDING NEUTRAL MUTATIONS:\nThe tree sequence now has {ts.num_mutations} mutations,\n"
        f"and mean pairwise nucleotide diversity is {ts.diversity():0.3e}.\n--------------------------\n")

    return ts

_, slim_tree_file, processed_tree_file, windows_table_file, trees_table_file = sys.argv # Kasper

params = [c for c in slim_tree_file.split("/")[-1].split("_")[:-1] if c != ""] #getting params from file name
print(params)
params_dict = dict(zip(params[::2], map(float, params[1::2])))

def ts_to_window_stats_df(ts):
    L = int(ts.sequence_length)
    windows = np.linspace(0, L, num=L//100_000)
    df = pd.DataFrame({
        'pos': [(x+y) // 2  for x,y in zip(windows[:-1], windows[1:])],
        'start': windows[:-1],
        'end': windows[1:],
        'pi': ts.diversity(windows=windows),
        'tajimas_d': ts.Tajimas_D(windows=windows),
        **params_dict
    })   
    return df 

def ts_to_tree_stats_df(ts):
    breakpoints = list(ts.breakpoints())
    df = pd.DataFrame({
        'pos': [(x+y) // 2  for x,y in zip(breakpoints[:-1], breakpoints[1:])],
        'start': breakpoints[:-1],
        'end': breakpoints[1:],        
        'branch_length': [tree.total_branch_length for tree in ts.trees()],
        **params_dict
    })   
    return df 


#Calling function to tree sequence file
processed_ts = ts_processer(slim_tree_file) 

#Writing processed tree sequence file
processed_ts.dump(processed_tree_file) 

table = ts_to_window_stats_df(processed_ts)

table.to_hdf(windows_table_file, key="df", format="table")

table = ts_to_tree_stats_df(processed_ts)

table.to_hdf(trees_table_file, key="df", format="table")
