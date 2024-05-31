
from gwf import Workflow

gwf = Workflow(defaults={'account': 'davide_intern'})

slim_model_file = 'scripts/davides_drive_model.slim'
slim_output_dir = 'steps/slim_tree_seqs/'
processed_output_dir = 'steps/processed_tree_seqs/'

N = 100
u = 5e-6
A = 3
g = 6000
rec_rates = [0, 5e-6, 5e-8, 5e-10]
vanriances = [2, 20, 200]

for i in range(10):
    for r in rec_rates:
        for v in vanriances:

            label = f'A_{A}__N_{N}__u_{u:.10f}__r_{r:.10f}__S_{v}_{i}'

            output_prefix = slim_output_dir + label
            slim_tree_file = output_prefix + '.trees'
            slim_plot_file = output_prefix + '.png'

            gwf.target('slim_'+label, inputs=[slim_model_file], outputs=[slim_tree_file]) << f"""

            mkdir -p {slim_output_dir}
            slim -d 'NR_AMPLICONS={A}' -d 'MUT_RATE_AMP={u}' -d 'REC_RATE={r}' -d 'VARIANCE_SEED={v}' \
                -d 'POPULATION_SIZE={N}' -d 'OUTPUT_PREFIX="{output_prefix}"' -d 'GENERATIONS={g}' \
                scripts/davides_drive_model.slim \
            """

            processed_trees_file = processed_output_dir + label + '.trees'
            table_file = processed_output_dir + label + '.h5'

            gwf.target('trees_'+label, inputs=[slim_tree_file], outputs=[processed_trees_file]) << f"""

            mkdir -p {processed_output_dir}
            python scripts/ts_processer.py {slim_tree_file} {processed_trees_file} {table_file}

            """
