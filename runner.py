import sys
import os
import glob
import getopt
import classes
import exceptions
import programs
import parsl
from parsl import load
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SingleNodeLauncher
from parsl.providers import TorqueProvider
from parsl.addresses import address_by_hostname
# from parsl.dataflow.dflow import DataFlowKernel

# --genome_generate g
# --sample_names s
# --genome_fasta f
# --genome_gtf t
# --genome_dir G
# --bam_dir B
# --align a
# --fastq_dir F
# --num_threads n
# --output_dir O
# --help h
# --run_mode m
# --genome_version v
# --TRUST4_dir T
# add list algorithms?

# Not necessary: --VDJC_ref r

config = Config(
    executors=[
        HighThroughputExecutor(
            label='4core_worker',
            cores_per_worker=4,
            mem_per_worker=16,
            max_workers=10,
            worker_debug=True,
            address=address_by_hostname(),
            provider=TorqueProvider(
                launcher=SingleNodeLauncher(),
                worker_init=("module load gcc/6.2.0 zlib/1.2.11 bzip2/1.0.6 xz/5.2.2 "
                            "perl/5.24.0 STAR/2.6.1d java-jdk/1.8.0_92 "
                            "mixcr/2.1.5 trust samtools/1.3.1 miniconda3/4.7.10 ; "
                            "source activate parsl_env; export PYTHONPATH='{}:{{PYTHONPATH}}'").format(os.getcwd()),
                init_blocks=1,
                max_blocks=5,
                min_blocks=0,
                nodes_per_block=1,
                walltime='99:00:00',
                scheduler_options='#PBS -l mem=48gb,nodes=1:ppn=12'
            ),
        ),
        HighThroughputExecutor(
            label='4core_largeMEM',
            mem_per_worker=40,
            cores_per_worker=4,
            max_workers=3,
            worker_debug=True,
            address=address_by_hostname(),
            provider=TorqueProvider(
                launcher=SingleNodeLauncher(),
                worker_init=("module load gcc/6.2.0 zlib/1.2.11 bzip2/1.0.6 xz/5.2.2 "
                            "perl/5.24.0 STAR/2.6.1d java-jdk/1.8.0_92 "
                            "mixcr/2.1.5 trust samtools/1.3.1 miniconda3/4.7.10 ; "
                            "source activate parsl_env; export PYTHONPATH='{}:{{PYTHONPATH}}'").format(os.getcwd()),
                init_blocks=1,
                max_blocks=5,
                min_blocks=0,
                nodes_per_block=1,
                walltime='99:00:00',
                #change to pmem??
                scheduler_options='#PBS -l mem=120gb,nodes=1:ppn=12'
            ),
        ),
        HighThroughputExecutor(
            label='1core_worker',
            mem_per_worker=5,
            cores_per_worker=1,
            max_workers=10,
            worker_debug=True,
            address=address_by_hostname(),
            provider=TorqueProvider(
                launcher=SingleNodeLauncher(),
                worker_init=("module load gcc/6.2.0 zlib/1.2.11 bzip2/1.0.6 xz/5.2.2 "
                            "perl/5.24.0 STAR/2.6.1d java-jdk/1.8.0_92 "
                            "mixcr/2.1.5 trust samtools/1.3.1 miniconda3/4.7.10 ; "
                            "source activate parsl_env; export PYTHONPATH='{}:{{PYTHONPATH}}'").format(os.getcwd()),
                init_blocks=1,
                max_blocks=5,
                min_blocks=0,
                nodes_per_block=1,
                walltime='99:00:00',
                scheduler_options='#PBS -l mem=20gb,nodes=1:ppn=4'
            ),
        )
    ],
    checkpoint_mode='task_exit'
)

load(config)

# Same as TRUST_script except added num threads (idk why I didn't modify original)
def main():
    '''
    Extracts all the parameters from the command line and creates a Pipeline
    object and a Genome object to run all the scripts.
    '''
    try:
        opts, args = getopt.getopt(sys.argv[1:],
                                "G:B:F:T:O:gf:t:an:m:p:v:s:h",
                                ["genome_dir=",
                                    "bam_dir=",
                                    "fastq_dir=",
                                    "TRUST4_dir=",
                                    "output_dir=",
                                    "genome_generate",
                                    "genome_fasta=",
                                    "genome_gtf=",
                                    "align",
                                    "num_threads=",
                                    "run_mode=",
                                    "run_program=",
                                    "genome_version=",
                                    "sample_names=",
                                    "help"
                                    ]
                                )
    except getopt.GetoptError as e:
        print(e)
        sys.exit(2)
    
    if len(args) > 0:
        print(args)
        message = "Error: non-paired arguments are not allowed."
        raise exceptions.WrongArgumentError(message)
    
    pipeline = classes.Pipeline()
    genome = classes.Genome()
    genome_generate = False
    align = False
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            description()
            sys.exit()
        elif opt in ("-G", "--genome_dir"):
            genome.set_dir(arg)
        elif opt in ("-B", "--bam_dir"):
            pipeline.set_bam_dir(arg)
        elif opt in ("-F", "--fastq_dir"):
            pipeline.set_fastq_dir(arg)
        elif opt in ("-T", "--TRUST4_dir"):
            pipeline.set_TRUST4_dir(arg)
        elif opt in ("-V", "--VDJer_dir"):
            pipeline.set_VDJer_dir(arg)
        elif opt in ("-O", "--output_dir"):
            pipeline.set_output_dir(arg)
        elif opt in ("-g", "--genome_generate"):
            genome_generate = True
        elif opt in ("-f", "--genome_fasta"):
            genome.set_fasta(arg)
        elif opt in ("-t", "--genome_gtf"):
            genome.set_gtf(arg)
        elif opt in ("-a", "--align"):
            align = True
        elif opt in ("-n", "--num_threads"):
            pipeline.set_num_threads(arg)
        elif opt in ("-m", "--run_mode"):
            pipeline.set_run_mode(arg)
        elif opt in ("-p", "--run_program"):
            pipeline.set_run_program(arg)
        elif opt in ("-v", "--genome_version"):
            genome.set_version(arg)
        elif opt in ("-s", "--sample_names"):
            pipeline.set_sample_names(arg)
        else:
            message = "Error: {opt} is not a valid option".format(opt=opt)
            raise exceptions.WrongArgumentError(message)
    
    check_attributes(genome, pipeline, genome_generate, align)
    run_pipeline(genome, pipeline, genome_generate, align)
    
def check_attributes(genome, pipeline, genome_generate, align):
    '''Checks if all the attributes required for a specific function have been input.'''
    if genome_generate:
        for genome_attribute in ["gtf", "fasta", "dir"]:
            if not hasattr(genome, genome_attribute):
                raise exceptions.MissingArgumentError("--genome_" + genome_attribute)
    if align:
        if not hasattr(genome, "dir"):
            raise exceptions.MissingArgumentError("--genome_dir")
        for pipeline_attribute in ["bam_dir", "fastq_dir"]:
            if not hasattr(pipeline, pipeline_attribute):
                raise exceptions.MissingArgumentError("--" + pipeline_attribute)
    
    dict_requirements = {"MiXCR": ["fastq_dir"],
                        "TRUST3": ["bam_dir"],
                        "TRUST4": ["bam_dir", "TRUST4_dir"],
                        "VDJer": ["bam_dir", "VDJer_dir"]}
    for program in pipeline.run_program:
        for pipeline_attribute in dict_requirements[program]:
            if not hasattr(pipeline, pipeline_attribute):
                raise exceptions.MissingArgumentError("--" + pipeline_attribute)

# No longer need genome when calling all the functions: fix?

def run_pipeline(genome, pipeline, genome_generate, align):
    '''
    Orchestrates the entire TCR/BCR pipeline for each sample (either input from the 
    command line or inferred from the fastq/bam directories). Each program is only 
    run if the output directory for that specific program/sample/run-mode is non-empty.
    Sample names that do not have corresponding fastq/bam directories are ignored.
    '''
    # all_futures = []
    # better to run with Data Futures or Appfutures??
    if genome_generate:
        if not os.path.isfile(genome.fasta + ".fai"):
            input_align = programs.STAR_genome(genome,
                                pipeline,
                                inputs=[programs.SAM_index(pipeline, genome.fasta, "faidx")])
        else:
            input_align = [programs.STAR_genome(genome, pipeline)]
    else:
        input_align = []
    
    # all_futures += input_align

    dict_fastq_futures = {}
    dict_bam_futures = {}
    #Check if fastq dir exists for sample names
    if align:
        if hasattr(pipeline, "sample_names"):
            list_samples = pipeline.sample_names
        else:
            list_samples = [path.split("/")[-2] for path in glob.glob(pipeline.fastq_dir + "*/")]
        for sample in list_samples:
            align_sample = True
            if not os.path.isdir(pipeline.bam_dir + sample):
                if sample not in os.listdir(pipeline.fastq_dir):
                    pipeline.write_log(("Sample {sample} does not have a corresponding "
                                        "fastq dir, alignment not run.").format(sample=sample))
                    align_sample = False
                else:
                    os.mkdir(pipeline.bam_dir + sample)
            elif not os.listdir(pipeline.bam_dir + sample):
                if sample not in os.listdir(pipeline.fastq_dir):
                    pipeline.write_log(("Sample {sample} does not have a corresponding "
                                        "fastq dir, alignment not run.").format(sample=sample))
                    align_sample = False
            else:
                pipeline.write_log(("Sample {sample} already has a non-empty "
                                    "alignment directory, alignment not run.").format(sample=sample))
                align_sample = False
            if align_sample:
                list_gz = list(filter(lambda x: x.endswith(".gz"),
                                                glob.glob(pipeline.fastq_dir + sample + "/*")))
                if list_gz:
                    dict_fastq_futures[sample] = [programs.gunzip(pipeline, list_gz)]
                    dict_bam_futures[sample] = [programs.STAR_align(genome,
                                                                    pipeline,
                                                                    sample,
                                                                    inputs=input_align + dict_fastq_futures[sample])]
                else:
                    dict_bam_futures[sample] = [programs.STAR_align(genome,
                                                                    pipeline,
                                                                    sample,
                                                                    inputs=input_align)]
    
    #merge fastq samples + bam samples??
    dict_subdirectory = {"MiXCR": True, "TRUST3": True, "TRUST4": False, "VDJer": False}
    if hasattr(pipeline, "sample_names"):
        list_samples = pipeline.sample_names
    else:
        if not hasattr(pipeline, "bam_dir"):
            list_samples = [path.split("/")[-2] for path in glob.glob(pipeline.fastq_dir + "*/")]
        else:
            list_samples = [path.split("/")[-2] for path in glob.glob(pipeline.bam_dir + "*/")]
    for sample in list_samples:
        sample_output = pipeline.output_dir + sample  + "/"
        if not os.path.isdir(sample_output):
            os.mkdir(sample_output)
        for program in pipeline.run_program:
            program_input = []
            if program == "MiXCR":
                if sample not in os.listdir(pipeline.fastq_dir):
                    pipeline.write_log(("Sample {sample} does not have a corresponding "
                                        "fastq dir, {program} not run.").format(sample=sample,
                                                                                program=program))
                    continue
                if sample in dict_fastq_futures:
                    program_input += dict_fastq_futures[sample]
                else:
                    list_gz = list(filter(lambda x: x.endswith(".gz"),
                                                glob.glob(pipeline.fastq_dir + sample + "/*")))
                    if list_gz:
                        program_input += [programs.gunzip(pipeline, list_gz)]
            else:
                if sample not in os.listdir(pipeline.bam_dir):
                    pipeline.write_log(("Sample {sample} does not have a corresponding "
                                        "bam dir, {program} not run.").format(sample=sample,
                                                                            program=program))
                    continue
                if sample in dict_bam_futures:
                    program_input += dict_bam_futures[sample]
            program_function = getattr(programs, program)
            if dict_subdirectory[program]:
                if not os.path.isdir(sample_output + program):
                    os.mkdir(sample_output + program)
                run_parameters = list_run_parameters(pipeline, sample_output + "/" + program + "/")
                if program == "MiXCR":
                    #Check if main MiXCR dir is not empty
                    if [mixcr_file for mixcr_file 
                        in glob.glob(sample_output + "/" + program + "/*")
                        if os.path.isfile(mixcr_file)]:
                        specific_input = program_input
                    else:
                        specific_input = [program_function(pipeline, genome, sample, inputs=program_input)]
                    program_specific_function = getattr(programs, program + "_specific")
                    program_specific_function(pipeline, genome, sample, run_parameters, inputs=specific_input)
                if program == "TRUST3":
                    bam_file = (pipeline.bam_dir + sample + "/" + sample
                                + ".Aligned.sortedByCoord.out.bam")
                    if not os.path.isfile(bam_file + ".bai"):
                        program_input = [programs.SAM_index(pipeline, bam_file, "index", inputs=program_input)]
                    program_function(pipeline, genome, sample, run_parameters, inputs=program_input)
            else:
                if not os.path.isdir(sample_output + program):
                    os.mkdir(sample_output + program)
                    program_function(pipeline, genome, sample, inputs=program_input)
                elif not os.listdir(sample_output + program):
                    program_function(pipeline, genome, sample, inputs=program_input)
                else:
                    pipeline.write_log(("Sample {sample} already has a non-empty "
                                        "{program} directory, {program} not run.").format(sample=sample,
                                                                                        program=program))
    parsl.wait_for_current_tasks()


def list_run_parameters(pipeline, program_dir):
    '''
    Removes parameters (TCR, BCR-heavy, BCR-light) that already have
    a non-empty output directory specifically for this sample/program/parameter
    combination. Creates a directory for samples that don't have such a directory.
    '''
    list_remove = []
    for parameter in pipeline.run_mode:
        if os.path.isdir(program_dir + parameter):
            if os.listdir(program_dir + parameter):
                list_remove.append(parameter)
                pipeline.write_log(("Sample {sample} already has a non-empty "
                                    "{program} directory, "
                                    "{program} not run.").format(sample=program_dir.split("/")[-3],
                                                                program=program_dir.split("/")[-2] + " " + parameter))
        else:
            os.mkdir(program_dir + parameter)
    return [parameter for parameter in pipeline.run_mode if parameter not in list_remove]

def description():
    '''Prints command line parameters and exits the program.'''
    manual = ("\nOPTIONS:\n\t -h or --help : display the manual\n"
        "\t -d or --g_dir : specify the genome directory path\n"
        "\t -f or --g_fa: specify the hg19 genome fasta file path\n"
        "\t -g or --g_gtf : specify the hg19 genome GTF file path\n"
        "\t -r or --r1r2_dir : specify the fastq directory for the sample\n"
        "\t -a or --align_dir : specify the directory of "
        "the aligned bam file\n"
        "\t -o or --out : specify the final output directory for TRUST\n"
        "\t -p or --param : specify additional parameters for TRUST\n"
        "\t -s or --log_suffix : specify suffix for log file name\n"  
        "\t -l or --loop : input for -r/--r1r2_dir and -a/--align_dir \n"
        "\t are directories containing different samples to loop over\n"
        "\t -t or --num_threads : specify the number of CPUs\n"
        )
    print(manual)

if __name__ == "__main__":
    main()

