from parsl.providers import LocalProvider
from parsl.addresses import address_by_hostname
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
import os

config = Config(
    executors=[
        HighThroughputExecutor(
            address=address_by_hostname(),
            cores_per_worker=4,
            provider=LocalProvider(
                worker_init=("source activate /cephfs/users/jbreynier/conda/parsl_env2 ; "
                "export PYTHONPATH='{}:{{PYTHONPATH}}'").format("/cephfs/users/jbreynier/code/CDR3merge"),
                init_blocks=1,
                max_blocks=10,
            ),
        )
    ],
    checkpoint_mode='task_exit'
)