from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SingleNodeLauncher
from parsl.providers import SlurmProvider
from parsl.addresses import address_by_hostname
import os

config = Config(
    executors=[
        HighThroughputExecutor(
            cores_per_worker=4,
            mem_per_worker=40,
            max_workers=4,
            worker_debug=True,
            address=address_by_hostname(),
            provider=SlurmProvider(
                'daenerys',
                worker_init=("conda activate /cephfs/users/jbreynier/conda/parsl_env2 ; "
                            "export PYTHONPATH='{}:{{PYTHONPATH}}'").format(os.getcwd()),
                init_blocks=1,
                max_blocks=3,
                min_blocks=0,
                nodes_per_block=1,
                walltime='99:00:00',
                scheduler_options='#SBATCH --exclude=kg15-11 --cpus-per-task=16 --mem=160gb --time=99:00:00',
            ),
        ),
    ],
    checkpoint_mode='task_exit'
)