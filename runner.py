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

# Need to:
# Add more specific subparameters for running each program
# Add parameters for merging / consensus steps
# Add --noMerging and --noConsensus flags (?)


# Add defaults: https://stackoverflow.com/questions/12151306/argparse-way-to-include-default-values-in-help

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
    # Make json optional i.e. doesn't look for one if not used
    parser.add_argument('-j', '--json', type=str, default = None, help='''JSON file to use for running Parsl''')
    parser.add_argument('-s', '--single_end', action='store_true', help='''Run in single-end mode''')
    args = parser.parse_args()
    print(args)
    load_config(args.config)
    pipeline = classes.Pipeline(args.run_program, args.run_mode, args.receptor, args.single_end)
    genome = classes.Genome(args.genome_version)
    pipeline.set_output_dir(args.output_dir)
    pipeline.set_fastq_dict(args.fastq_path)
    # Fix JSON input: not required

    pipeline.extract_json(args.json)

    if ('TRUST3' in pipeline.run_program or 'TRUST4' in pipeline.run_program 
        or 'VDJer' in pipeline.run_program):
        if 'VDJer' in pipeline.run_program:
            genome_download_future = genome.download_genome(pipeline.output_dir, True)
        else:
            genome_download_future = genome.download_genome(pipeline.output_dir, False)
    else:
        genome_download_future = []
    
    check_attributes(genome, pipeline)
    run_pipeline(pipeline, genome, genome_download_future)

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

def make_bash_app(pipeline, func_name):
    program = func_name.split("_")[0]
    if program in pipeline.config_dict.keys():
        if "executor" in pipeline.config_dict[program].keys():
            executor = [pipeline.config_dict[program]["executor"]]
        else:
            executor = "all"
    else:
        executor = "all"
    program_function = getattr(programs, func_name)
    return bash_app(program_function, executors=executor)

# CHECK FOR FINAL FILE IN EACH CASE!!

# FIX Deleting files!! <- remove directories also

def run_pipeline(pipeline, genome, download_future):
    '''
    Orchestrates the entire TCR/BCR pipeline for each sample (either input from the 
    command line or inferred from the fastq/bam directories). Each program is only 
    run if the output directory for that specific program/sample/run-mode is non-empty.
    Sample names that do not have corresponding fastq/bam directories are ignored.
    '''
    dict_subdirectory = {"MiXCR": True,
                            "TRUST3": True,
                            "TRUST4": False,
                            "VDJer": True,
                            "CATT": True}
    for sample in pipeline.fastq_dict.keys():
        sample_output = pipeline.output_dir + sample + "/"
        if not os.path.isdir(sample_output):
            os.mkdir(sample_output)
        # STAR_align(pipeline, genome, sample_name, inputs=[], VDJer=False)
        if ('TRUST3' in pipeline.run_program or 'TRUST4' in pipeline.run_program):
            os.makedirs(sample_output + "STAR_align/", exist_ok=True)
            align_dir = sample_output + "STAR_align/" + genome.version
            if not os.path.isdir(align_dir):
                os.mkdir(align_dir)
                align_future = [programs.STAR_align(pipeline, genome, sample, inputs=download_future)]
            elif not os.listdir(align_dir):
                align_future = [programs.STAR_align(pipeline, genome, sample, inputs=download_future)]
            else:
                if os.path.isfile(align_dir+"/"+sample+".Aligned.sortedByCoord.out.bam"):
                    print(f"Non-empty {genome.version} STAR alignment directory for {sample}: alignment already run.")
                    align_future = []
                else:
                    align_files = glob.glob(align_dir+"/*")
                    for rm_file in align_files:
                        os.remove(rm_file)
                    align_future = [programs.STAR_align(pipeline, genome, sample, inputs=download_future)]
        if 'VDJer' in pipeline.run_program:
            if genome.version != "hg38":
                VDJer_align_dir = sample_output + "STAR_align/hg38"
                if not os.path.isdir(VDJer_align_dir):
                    os.mkdir(VDJer_align_dir)
                    VDJer_align_future = [programs.STAR_align(pipeline, genome, sample, inputs=download_future, VDJer=True)]
                elif not os.listdir(VDJer_align_dir):
                    VDJer_align_future = [programs.STAR_align(pipeline, genome, sample, inputs=download_future, VDJer=True)]
                else:
                    if os.path.isfile(VDJer_align_dir+"/"+sample+".Aligned.sortedByCoord.out.bam"):
                        print(f"Non-empty hg38 STAR alignment directory for {sample}: alignment already run.")
                        VDJer_align_future = []
                    else:
                        VDJer_align_files = glob.glob(VDJer_align_dir+"/*")
                        for rm_file in VDJer_align_files:
                            os.remove(rm_file)
                        VDJer_align_future = [programs.STAR_align(pipeline, genome, sample, inputs=download_future, VDJer=True)]
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
                run_parameters = list_run_parameters(pipeline, sample_output+"/"+program+"/", program == "CATT")
                if program == "MiXCR":
                    # Check if main MiXCR dir is not empty
                    # Make more modular ?? Check for specific files ?
                    # MiXCR_align
                    # MiXCR_assemblePartial
                    # MiXCR_extendAlignments
                    # MiXCR_assemble
                    # MiXCR_exportClones
                    rerun_exportClones = True
                    mixcr_intermediate_files = [mixcr_file for mixcr_file in glob.glob(sample_output+"/MiXCR/*") 
                                                    if os.path.isfile(mixcr_file)]
                    if mixcr_intermediate_files:
                            if os.path.isfile(sample_output+"/MiXCR/"+sample+"_clones.clns"):
                                print(f"Non-empty MiXCR directory for {sample}: align/assemblePartial/extendAlignments/assemble already run.")
                                rerun_exportClones = False
                                exportClones_input = []
                            else:
                                for rm_file in mixcr_intermediate_files:
                                    os.remove(rm_file)
                                exportClones_input = run_mixcr_prep(pipeline, genome, program_input, sample)
                    else:
                        exportClones_input = run_mixcr_prep(pipeline, genome, program_input, sample)
                    if not rerun_exportClones:
                        # rerunning because previous files have been modified
                        for parameter in run_parameters:
                            programs.MiXCR_exportClones(pipeline, genome, sample, parameter, inputs=exportClones_input)
                    else:
                        for parameter in pipeline.receptor:
                            mixcr_files = glob.glob(sample_output+"/MiXCR/"+parameter+"/*")
                            for rm_file in mixcr_files:
                                os.remove(rm_file)
                            programs.MiXCR_exportClones(pipeline, genome, sample, parameter, inputs=exportClones_input)
                if program == "TRUST3":
                    bam_file = align_dir + "/" + sample + ".Aligned.sortedByCoord.out.bam"
                    if not os.path.isfile(bam_file + ".bai"):
                        program_input = [programs.samtools_index(pipeline, bam_file, "index", inputs=program_input)]
                if program != "MiXCR":
                    program_function = getattr(programs, program)
                    for parameter in run_parameters:
                        program_function(pipeline, genome, sample, parameter, inputs=program_input)
            else:
                program_function = getattr(programs, program)
                if not os.path.isdir(sample_output + program):
                        os.mkdir(sample_output + program)
                        program_function(pipeline, genome, sample, inputs=program_input)
                elif not os.listdir(sample_output + program):
                    program_function(pipeline, genome, sample, inputs=program_input)
                else:
                    if program == "TRUST4":
                        output_file = sample+"_report.tsv"
                    if os.path.isfile(sample_output+program+"/"+output_file):
                        print(f"Non-empty {program} directory for {sample}: {program} already run.")
                    else:
                        all_output_files = glob.glob(sample_output+program+"/*")
                        for rm_file in all_output_files:
                            os.remove(rm_file)
                        program_function(pipeline, genome, sample, inputs=program_input)
    parsl.wait_for_current_tasks()

def run_mixcr_prep(pipeline, genome, program_input, sample):
    dict_assemblePartial_futures = {}
    dict_assemblePartial_futures[1] = [programs.MiXCR_align(pipeline, genome, sample, inputs=program_input)]
    for num in range(1, pipeline.rescue_num+1):
        if num == pipeline.rescue_num:
            extendAlignments_input = [programs.MiXCR_assemblePartial(pipeline, genome, sample, str(num), inputs=dict_assemblePartial_futures[num])]
            # print(extendAlignments_input[0].result())
            print("rescue_num is "+str(num)+"\n")
        else:
            dict_assemblePartial_futures[num+1] = [programs.MiXCR_assemblePartial(pipeline, genome, sample, str(num), inputs=dict_assemblePartial_futures[num])]
            # print(assemblePartial_input[0].result())
            print("rescue_num is "+str(num)+"\n")
    assemble_input = [programs.MiXCR_extendAlignments(pipeline, genome, sample, str(pipeline.rescue_num), inputs=extendAlignments_input)]
    return [programs.MiXCR_assemble(pipeline, genome, sample, inputs=assemble_input)]

# def run_pipeline(genome, pipeline, genome_generate, align):
#     '''
#     Orchestrates the entire TCR/BCR pipeline for each sample (either input from the 
#     command line or inferred from the fastq/bam directories). Each program is only 
#     run if the output directory for that specific program/sample/run-mode is non-empty.
#     Sample names that do not have corresponding fastq/bam directories are ignored.
#     '''
#     # all_futures = []
#     # better to run with Data Futures or Appfutures??
#     if genome_generate:
#         if not os.path.isfile(genome.fasta + ".fai"):
#             input_align = programs.STAR_genome(genome,
#                                 pipeline,
#                                 inputs=[programs.SAM_index(pipeline, genome.fasta, "faidx")])
#         else:
#             input_align = [programs.STAR_genome(genome, pipeline)]
#     else:
#         input_align = []
    
#     # all_futures += input_align

#     dict_fastq_futures = {}
#     dict_bam_futures = {}
#     #Check if fastq dir exists for sample names
#     if align:
#         if hasattr(pipeline, "sample_names"):
#             list_samples = pipeline.sample_names
#         else:
#             list_samples = [path.split("/")[-2] for path in glob.glob(pipeline.fastq_dir + "*/")]
#         for sample in list_samples:
#             align_sample = True
#             if not os.path.isdir(pipeline.bam_dir + sample):
#                 if sample not in os.listdir(pipeline.fastq_dir):
#                     pipeline.write_log(("Sample {sample} does not have a corresponding "
#                                         "fastq dir, alignment not run.").format(sample=sample))
#                     align_sample = False
#                 else:
#                     os.mkdir(pipeline.bam_dir + sample)
#             elif not os.listdir(pipeline.bam_dir + sample):
#                 if sample not in os.listdir(pipeline.fastq_dir):
#                     pipeline.write_log(("Sample {sample} does not have a corresponding "
#                                         "fastq dir, alignment not run.").format(sample=sample))
#                     align_sample = False
#             else:
#                 pipeline.write_log(("Sample {sample} already has a non-empty "
#                                     "alignment directory, alignment not run.").format(sample=sample))
#                 align_sample = False
#             if align_sample:
#                 list_gz = list(filter(lambda x: x.endswith(".gz"),
#                                                 glob.glob(pipeline.fastq_dir + sample + "/*")))
#                 if list_gz:
#                     dict_fastq_futures[sample] = [programs.gunzip(pipeline, list_gz)]
#                     dict_bam_futures[sample] = [programs.STAR_align(genome,
#                                                                     pipeline,
#                                                                     sample,
#                                                                     inputs=input_align + dict_fastq_futures[sample])]
#                 else:
#                     dict_bam_futures[sample] = [programs.STAR_align(genome,
#                                                                     pipeline,
#                                                                     sample,
#                                                                     inputs=input_align)]
    
#     #merge fastq samples + bam samples??
#     dict_subdirectory = {"MiXCR": True, "TRUST3": True, "TRUST4": False, "VDJer": False}
#     if hasattr(pipeline, "sample_names"):
#         list_samples = pipeline.sample_names
#     else:
#         if not hasattr(pipeline, "bam_dir"):
#             list_samples = [path.split("/")[-2] for path in glob.glob(pipeline.fastq_dir + "*/")]
#         else:
#             list_samples = [path.split("/")[-2] for path in glob.glob(pipeline.bam_dir + "*/")]
#     for sample in list_samples:
#         sample_output = pipeline.output_dir + sample  + "/"
#         if not os.path.isdir(sample_output):
#             os.mkdir(sample_output)
#         for program in pipeline.run_program:
#             program_input = []
#             if program == "MiXCR":
#                 if sample not in os.listdir(pipeline.fastq_dir):
#                     pipeline.write_log(("Sample {sample} does not have a corresponding "
#                                         "fastq dir, {program} not run.").format(sample=sample,
#                                                                                 program=program))
#                     continue
#                 if sample in dict_fastq_futures:
#                     program_input += dict_fastq_futures[sample]
#                 else:
#                     list_gz = list(filter(lambda x: x.endswith(".gz"),
#                                                 glob.glob(pipeline.fastq_dir + sample + "/*")))
#                     if list_gz:
#                         program_input += [programs.gunzip(pipeline, list_gz)]
#             else:
#                 if sample not in os.listdir(pipeline.bam_dir):
#                     pipeline.write_log(("Sample {sample} does not have a corresponding "
#                                         "bam dir, {program} not run.").format(sample=sample,
#                                                                             program=program))
#                     continue
#                 if sample in dict_bam_futures:
#                     program_input += dict_bam_futures[sample]
#             program_function = getattr(programs, program)
#             if dict_subdirectory[program]:
#                 if not os.path.isdir(sample_output + program):
#                     os.mkdir(sample_output + program)
#                 run_parameters = list_run_parameters(pipeline, sample_output + "/" + program + "/")
#                 if program == "MiXCR":
#                     #Check if main MiXCR dir is not empty
#                     if [mixcr_file for mixcr_file 
#                         in glob.glob(sample_output + "/" + program + "/*")
#                         if os.path.isfile(mixcr_file)]:
#                         specific_input = program_input
#                     else:
#                         specific_input = [program_function(pipeline, genome, sample, inputs=program_input)]
#                     program_specific_function = getattr(programs, program + "_specific")
#                     program_specific_function(pipeline, genome, sample, run_parameters, inputs=specific_input)
#                 if program == "TRUST3":
#                     bam_file = (pipeline.bam_dir + sample + "/" + sample
#                                 + ".Aligned.sortedByCoord.out.bam")
#                     if not os.path.isfile(bam_file + ".bai"):
#                         program_input = [programs.SAM_index(pipeline, bam_file, "index", inputs=program_input)]
#                     program_function(pipeline, genome, sample, run_parameters, inputs=program_input)
#             else:
#                 if not os.path.isdir(sample_output + program):
#                     os.mkdir(sample_output + program)
#                     program_function(pipeline, genome, sample, inputs=program_input)
#                 elif not os.listdir(sample_output + program):
#                     program_function(pipeline, genome, sample, inputs=program_input)
#                 else:
#                     pipeline.write_log(("Sample {sample} already has a non-empty "
#                                         "{program} directory, {program} not run.").format(sample=sample,
#                                                                                         program=program))
#     parsl.wait_for_current_tasks()


def list_run_parameters(pipeline, program_dir, CATT=False):
    '''
    Removes parameters (TCR, BCR-heavy, BCR-light) that already have
    a non-empty output directory specifically for this sample/program/parameter
    combination. Creates a directory for samples that don't have such a directory.
    '''
    list_remove = []
    list_receptors = []
    program_name = program_dir.split("/")[-2]
    sample_name = program_dir.split("/")[-3]
    dict_CATT = {
        "TCR": ["TRA", "TRB"],
        "BCR-heavy": ["IGH"],
        "BCR-light": []
    }
    dict_allowed_receptors = {"TRUST3": ['TCR', 'BCR-heavy', 'BCR-light'],
                                "TRUST4": ['TCR', 'BCR-heavy', 'BCR-light'],
                                "MiXCR": ['TCR', 'BCR-heavy', 'BCR-light'],
                                "CATT": ['TCR', 'BCR-heavy'],
                                "VDJer": ['BCR-heavy', 'BCR-light']}
    dict_file_to_check = {
        "MiXCR": "",
        "TRUST3": "",
        "TRUST4": "",
        "CATT": "",
        "VDJer": ""
    }
    # Add the check output files option

    if CATT:
        for parameter in pipeline.receptor:
            list_receptors += dict_CATT[parameter]
    else:
        list_receptors = [receptor for receptor in pipeline.receptor 
                            if receptor in dict_allowed_receptors[program_name]]

    for parameter in list_receptors:
        if os.path.isdir(program_dir + parameter):
            if os.listdir(program_dir + parameter):
                if os.path.isfile(program_dir+parameter+"/"+dict_file_to_check[program_name]):
                    list_remove.append(parameter)
                    print(("Sample {sample} already has a non-empty "
                            "{program} directory: "
                            "{program_specific} already run.").format(sample=sample_name,
                                                            program=program_name+" "+parameter,
                                                            program_specific=(program_name+" "+parameter 
                                                                if program_name != "MiXCR" else "MiXCR exportClones "+parameter)))
                else:
                    program_files = glob.glob(program_dir+parameter+"/*")
                    for rm_file in program_files:
                        os.remove(rm_file)
        else:
            os.mkdir(program_dir + parameter)
    print([parameter for parameter in list_receptors if parameter not in list_remove])
    return [parameter for parameter in list_receptors if parameter not in list_remove]

if __name__ == "__main__":
    main()


