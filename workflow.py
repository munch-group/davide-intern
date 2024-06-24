
from gwf import Workflow

gwf = Workflow(defaults={'account': 'davide_intern'})

N = 10000 #Population size
u = 5e-6 #mutation rate
A = 3 #number of amplicons
g = 60000 #N * 10
# rec_rates = [5e-8, 5e-10]
rec_rates = [5e-8]
variances = [2, 20, 200]

slim_model_file = 'scripts/davides_drive_model.slim'
slim_output_dir = f'steps/slim_tree_seqs/{g}/'
processed_output_dir = f'steps/processed_tree_seqs/{g}/'


for i in range(20):
    for r in rec_rates:
        for v in variances:

            label = f'A_{A}__N_{N}__u_{u:.10f}__r_{r:.10f}__S_{v}_{i}'

            output_prefix = slim_output_dir + label
            slim_tree_file = output_prefix + '.trees'
            slim_plot_file = output_prefix + '.png'

            gwf.target('slim_'+label, inputs=[slim_model_file], outputs=[slim_tree_file], walltime='10:00:00', memory = '10gb') << f"""

            mkdir -p {slim_output_dir}
            slim -d 'NR_AMPLICONS={A}' -d 'MUT_RATE_AMP={u}' -d 'REC_RATE={r}' -d 'VARIANCE_SEED={v}' \
                -d 'POPULATION_SIZE={N}' -d 'OUTPUT_PREFIX="{output_prefix}"' -d 'GENERATIONS={g}' \
                 -d "SEL_COEF=0" -d "SWEEP_POS=15e6" -d "SWEEP_START={g-500}" -d 'TMPDIR="steps/tmp"' \
                scripts/davides_drive_model.slim \
            """

            processed_trees_file = processed_output_dir + label + '.trees'
            win_table_file = processed_output_dir + label + '_windowstats' + '.h5'
            tree_table_file = processed_output_dir + label + '_treestats' + '.h5'

            gwf.target('trees_'+label, inputs=[slim_tree_file], outputs=[processed_trees_file, win_table_file, tree_table_file]) << f"""
            mkdir -p {processed_output_dir} 
            python scripts/ts_processer.py {slim_tree_file} {processed_trees_file} {win_table_file} {tree_table_file}
            """

    ###########################
    # regular sweep

    # label = f'sweep__N_{N}__r_{r:.10f}_{i}'
    label = f'A_{A}__N_{N}__u_0__r_{r:.10f}__S_0.1_{i}'

    output_prefix = slim_output_dir + label
    slim_tree_file = output_prefix + '.trees'
    slim_plot_file = output_prefix + '.png'

    gwf.target('slim_'+label, inputs=[slim_model_file], outputs=[slim_tree_file], walltime='10:00:00', memory = '10gb') << f"""

    mkdir -p {slim_output_dir}
    mkdir -p steps/tmp
    slim -d 'NR_AMPLICONS={A}' -d 'MUT_RATE_AMP=0' -d 'REC_RATE={r}' -d 'VARIANCE_SEED={v}' \
        -d 'POPULATION_SIZE={N}' -d 'OUTPUT_PREFIX="{output_prefix}"' -d 'GENERATIONS={g}' \
        -d "SEL_COEF=0.1" -d "SWEEP_POS=15e6" -d "SWEEP_START={g-500}" -d 'TMPDIR="steps/tmp"' \
        scripts/davides_drive_model.slim \
    """

    processed_trees_file = processed_output_dir +  label + '.trees'
    win_table_file = processed_output_dir + label + '_windowstats' + '.h5'
    tree_table_file = processed_output_dir + label + '_treestats' + '.h5'

    gwf.target('trees_'+label, inputs=[slim_tree_file], outputs=[processed_trees_file, win_table_file, tree_table_file]) << f"""
    mkdir -p {processed_output_dir} 
    python scripts/ts_processer.py {slim_tree_file} {processed_trees_file} {win_table_file} {tree_table_file}
    """