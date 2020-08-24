from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SingleNodeLauncher
from parsl.providers import TorqueProvider
from parsl.addresses import address_by_hostname
import os

config = Config(
    executors=[
        HighThroughputExecutor(
            cores_per_worker=1,
            mem_per_worker=4,
            max_workers=12,
            worker_debug=True,
            address=address_by_hostname(),
            provider=TorqueProvider(
                launcher=SingleNodeLauncher(),
                worker_init=("module load gcc/6.2.0 perl/5.24.0 zlib/1.2.8 meme/5.0.5 bedtools/2.27.1 R/3.6.1 miniconda3/4.7.10; "
                            "source activate parsl_env; export PYTHONPATH='{}:{{PYTHONPATH}}'").format(os.getcwd()),
                init_blocks=1,
                max_blocks=10,
                min_blocks=0,
                nodes_per_block=1,
                walltime='99:00:00',
                scheduler_options='#PBS -l mem=48gb,nodes=1:ppn=12'
            ),
        ),
    ],
    checkpoint_mode='task_exit'
)