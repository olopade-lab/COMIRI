import sys
import os
import glob
import getopt
import classes
import exceptions
import programs
import importlib
import parsl
import argparse
import json
from shutil import which
from parsl import load, python_app, bash_app

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


# Add abspath??

def main():
    '''
    Extracts all the parameters from the command line and creates a Pipeline
    object and a Genome object to run all the scripts.
    '''
    parser = argparse.ArgumentParser(description='''CDR3merge: a parallelized consensus caller 
                                                    for immune receptor identification in bulk RNA-seq''',
                                    epilog='''If you have any questions, please submit an issue at: 
                                            https://github.com/olopade-lab/CDR3merge/issues''')
    parser.add_argument('-f', '--fastq_path', type=str, required=True, help='''Path of the input fastq files''')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='''Path of the output directory''')                                                                         
    parser.add_argument('-g', '--genome_version', type=str, required=False, choices=['hg19', 'hg38'], default='hg38', help='''Genome version to use''')
    parser.add_argument('-p', '--run_program', type=str, nargs='+', required=False, 
                        default=['TRUST3', 'TRUST4', 'MiXCR', 'VDJer', 'CATT'], choices=['TRUST3', 'TRUST4', 'MiXCR', 'VDJer', 'CATT'], help='''List of programs to run''')
    parser.add_argument('-m', '--run_mode', type=str, required=False, default='docker', 
                        choices=['docker', 'singularity', 'local'], help='''Run mode for pipeline''')
    parser.add_argument('-r', '--receptor', type=str, nargs='+', required=False, 
                        default=['TCR', 'BCR-heavy', 'BCR-light'], choices=['TCR', 'BCR-heavy', 'BCR-light'], help='''List of receptors to extract''')
    parser.add_argument('-c', '--config', type=str, default='local', help='''Config file to use for running Parsl''')
    parser.add_argument('-j', '--json', type=str, default='local', help='''JSON file to use for running Parsl''')
    parser.add_argument('-s', '--single_end', action='store_true', help='''Run in single-end mode''')
    # parser.add_argument('-G', '--genome_dir', type=str, required=True, help='''Path to the directory containing 
    #                                                                             the genome .fasta/.fa and .gtf files''')
    # parser.add_argument('-A', '--align_dir', type=str, required=True, help='''Path to the directory containing 
    #                                                                             the alignment .bam and .bai files for each sample''')
    # parser.add_argument('-g', '--genome_generate', action='store_true', help='''Run STAR genomeGenerate''')
    # parser.add_argument('-a', '--align', action='store_true', help='''Run STAR align''')
    args = parser.parse_args()
    load_config(args.config)
    pipeline = classes.Pipeline(args.run_program, args.run_mode, args.receptor, args.single_end)
    genome = classes.Genome(args.genome_version)
    pipeline.set_output_dir(args.output_dir)
    pipeline.set_fastq_dict(args.fastq_path)

    if ('TRUST3' in pipeline.run_program or 'TRUST4' in pipeline.run_program 
        or 'VDJer' in pipeline.run_program):
        if 'VDJer' in pipeline.run_program:
            genome_download_future = genome.download_genome(pipeline.output_dir, True)
        else:
            genome_download_future = genome.download_genome(pipeline.output_dir, False)
    else:
        genome_download_future = []
    
    check_attributes(genome, pipeline)
    run_pipeline2(pipeline, genome, genome_download_future)

def load_config(config_name):
    base_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1])
    config = os.path.join(base_dir, 'configs', '{}.py'.format(config_name))
    if not os.path.isfile(config):
        raise exceptions.IncorrectPathError("Cannot find the config file <{config}.py>.".format(config=config_name))
    try:
        spec = importlib.util.spec_from_file_location('', config)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        parsl.load(module.config)
    except Exception as e:
        raise exceptions.IncorrectInputFiles(("Could not load specified config from <{config}.py> :"
                                            "\n {exception}.").format(config=config_name, exception=e))

def check_attributes(genome, pipeline):
    '''Raises warning if there are incompatibility issues between different parameters chosen.'''
    dict_allowed_receptors = {"TRUST3": ['TCR', 'BCR-heavy', 'BCR-light'],
                                "TRUST4": ['TCR', 'BCR-heavy', 'BCR-light'],
                                "MiXCR": ['TCR', 'BCR-heavy', 'BCR-light'],
                                "CATT": ['TCR', 'BCR-heavy'],
                                "VDJer": ['BCR-heavy', 'BCR-light']}
    if "VDJer" in pipeline.run_program:
        if genome.version == "hg19":
            print("Warning: VDJer can only be used with samples aligned to hg38, hg38 genome will be used for VDJer specifically.")
        if pipeline.single_end:
            print("Warning: VDJer does not support single-end reads. VDJer will not be run for these samples.")
    
    for receptor in pipeline.receptor:
        for program in pipeline.run_program:
            if receptor not in dict_allowed_receptors[program]:
                print(f"Warning: {program} does not extract {receptor} CDR3 sequences.")


def run_pipeline2(pipeline, genome, download_future):
    '''
    Orchestrates the entire TCR/BCR pipeline for each sample (either input from the 
    command line or inferred from the fastq/bam directories). Each program is only 
    run if the output directory for that specific program/sample/run-mode is non-empty.
    Sample names that do not have corresponding fastq/bam directories are ignored.
    '''
    for sample in pipeline.fastq_dict.keys():
        dict_subdirectory = {"MiXCR": True,
                            "TRUST3": True,
                            "TRUST4": False,
                            "VDJer": True,
                            "CATT": True}
        sample_output = pipeline.output_dir + sample + "/"
        if not os.path.isdir(sample_output):
            os.mkdir(sample_output)
        if ('TRUST3' in pipeline.run_program or 'TRUST4' in pipeline.run_program):
            os.makedirs(sample_output + "STAR_align/", exist_ok=True)
            align_dir = sample_output + "STAR_align/" + genome.version
            if not os.path.isdir(align_dir):
                os.mkdir(align_dir)
                align_future = [programs.STAR_align(pipeline, genome, sample, download_future)]
            elif not os.listdir(align_dir):
                align_future = [programs.STAR_align(pipeline, genome, sample, download_future)]
            else:
                print(f"Non-empty {genome.version} STAR alignment directory for {sample}: alignment already run.")
                align_future = []
        if 'VDJer' in pipeline.run_program:
            if genome.version != "hg38":
                VDJer_align_dir = sample_output + "STAR_align/hg38"
                if not os.path.isdir(VDJer_align_dir):
                    os.mkdir(VDJer_align_dir)
                    VDJer_align_future = [programs.STAR_align(pipeline, genome, sample, True, download_future)]
                elif not os.listdir(VDJer_align_dir):
                    VDJer_align_future = [programs.STAR_align(pipeline, genome, sample, True, download_future)]
                else:
                    print(f"Non-empty hg38 STAR alignment directory for {sample}: alignment already run.")
                    VDJer_align_future = []
            else:
                VDJer_align_future = align_future
            
        for program in pipeline.run_program:
            if program in ["TRUST3", "TRUST4"]:
                program_input = align_future
            elif program == "VDJer":
                program_input = VDJer_align_future
            else:
                program_input = []
            if dict_subdirectory[program]:
                if not os.path.isdir(sample_output + program):
                    os.mkdir(sample_output + program)
                run_parameters = list_run_parameters(pipeline, sample_output+"/"+program+"/")
                if program == "MiXCR":
                    #Check if main MiXCR dir is not empty
                    # Make more modular ?? Check for specific files ?
                    if [mixcr_file for mixcr_file
                        in glob.glob(sample_output+"/"+program+"/*")
                        if os.path.isfile(mixcr_file)]:
                        specific_input = program_input
                    else:
                        specific_input = [programs.MiXCR_align(pipeline, genome, sample, inputs=program_input)]
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
                                        "{program} directory: {program} already run.").format(sample=sample,
                                                                                        program=program))
    parsl.wait_for_current_tasks()



            if not os.path.isdir(sample_output + sample):
                os.mkdir(pipeline.bam_dir + sample)

            elif not os.listdir(pipeline.bam_dir + sample):

            else:


            align_future = [programs.STAR_align(pipeline, genome, sample, download_future)]
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
    for parameter in pipeline.receptor:
        if os.path.isdir(program_dir + parameter):
            if os.listdir(program_dir + parameter):
                list_remove.append(parameter)
                print(("Sample {sample} already has a non-empty "
                                    "{program} directory: "
                                    "{program} already run.").format(sample=program_dir.split("/")[-3],
                                                                program=program_dir.split("/")[-2] + " " + parameter))
        else:
            os.mkdir(program_dir + parameter)
    return [parameter for parameter in pipeline.receptor if parameter not in list_remove]


if __name__ == "__main__":
    main()


