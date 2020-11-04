import classes
from parsl import python_app, bash_app

@bash_app
def samtools_index(pipeline, bam_or_fasta, samtools_command, inputs=[]):
    import classes

    # docker container: biocontainers/samtools or mgibio/samtools
    # use --user root
     
    if pipeline.run_mode == "docker":
        samtools_index_script = ("docker run --rm -v {input_dir}:/input "
                        "{image} {samtools_command} "
                        "/input/{input_file}").format(
                            input_dir=bam_or_fasta.rsplit("/", 1)[0],
                            image=pipeline.docker_dict["samtools"],
                            samtools_command=samtools_command,
                            input_file=bam_or_fasta.rsplit("/", 1)[-1])

    elif pipeline.run_mode == "singularity":
        pass

    elif pipeline.run_mode == "local":
        samtools_index_script = ("{command} {path} {samtools_command} "
                            "{bam_or_fasta}").format(command=classes.get_command(pipeline, "samtools"),
                                                path=classes.get_path(pipeline, "samtools"),
                                                samtools_command=samtools_command,
                                                bam_or_fasta=bam_or_fasta)
    pipeline.write_log("Generating index file for {f}".format(f=bam_or_fasta.split("/")[-1]))
    return samtools_index_script

# @bash_app
# def STAR_align(pipeline, genome, sample_name, inputs=[], VDJer=False, stderr='std.err', stdout='std.out'):
#     return 'echo "Hello World!"'

@bash_app
def STAR_align(pipeline, genome, sample_name, inputs=[], VDJer=False):
    import classes
    # pipeline.write_log("STAR")
    # docker container: mgibio/star
    if VDJer:
        genome_version = "hg38"
    else:
        genome_version = genome.version
    num_threads = classes.get_threads(pipeline, "STAR")
    specific_output = pipeline.output_dir+sample_name+"/STAR_align/"+genome_version
    specific_input = pipeline.fastq_dict[sample_name][0].rsplit("/", 1)[0]
    
    # Added limitBAM
    if pipeline.run_mode == "docker":
        STAR_align_script = ("docker run --rm -v {specific_output}:/output "
                        "-v {specific_input}:/input "
                        "-v {genome_dir}:/genome {image} "
                        "STAR --genomeDir /genome "
                        "--readFilesIn {fastq} "
                        "--outSAMtype BAM SortedByCoordinate "
                        "--outFileNamePrefix /output/{sample_name}. "
                        "--alignEndsType EndToEnd "
                        "--outSAMunmapped Within "
                        "--runThreadN {num_threads}").format(specific_output = specific_output,
                            specific_input = specific_input,
                            genome_dir=genome.dir+genome_version,
                            image=pipeline.docker_dict["STAR"],
                            fastq=" ".join(["/input/"+path.rsplit("/", 1)[-1] for path in pipeline.fastq_dict[sample_name]]),
                            sample_name=sample_name,
                            num_threads=num_threads)

    elif pipeline.run_mode == "singularity":
        pass

    elif pipeline.run_mode == "local":
        STAR_align_script = ("{command} {path} STAR --genomeDir {genome_dir} "
                        "--readFilesIn {fastq} "
                        "--outSAMtype BAM SortedByCoordinate "
                        "--outFileNamePrefix {align_dir}. "
                        "--alignEndsType EndToEnd "
                        "--outSAMunmapped Within "
                        "--runThreadN {num_threads}").format(
                            command=classes.get_command(pipeline, "STAR"),
                            path=classes.get_path(pipeline, "STAR"),
                            genome_dir=genome.dir+genome_version,
                            fastq=" ".join(pipeline.fastq_dict[sample_name]),
                            align_dir=specific_output+"/"+sample_name,
                            num_threads=num_threads)

    pipeline.write_log("Alignment of {sample} with {version}".format(sample=sample_name, version=genome_version))
    return STAR_align_script


# CHANGE TO individual app for each parameter !!
@bash_app
def TRUST3(pipeline, genome, sample_name, parameter, inputs=[]):
    import classes
    # docker container: mgibio/trust

    dict_parameters = {"TCR": "",
                        "BCR-heavy": "-B",
                        "BCR-light": "-B -L"}

    num_threads = classes.get_threads(pipeline, "TRUST3")
    specific_input = pipeline.output_dir + sample_name + "/STAR_align/" + genome.version
    bam_file = sample_name + ".Aligned.sortedByCoord.out.bam"
    specific_output = pipeline.output_dir + sample_name + "/TRUST3/" + parameter

    # TRUST runs with python 2 !! : "module load gcc/6.2.0 python/2.7.13 trust && "
    
    if pipeline.run_mode == "docker":
        TRUST3_script = ("docker run --user root --rm -v {specific_output}:/output "
                        "-v {specific_input}:/input {image} "
                        "trust -f /input/{align_bam} -I 200 -g {genome} "
                        "-o /output/ {param} -E -c -n {num_threads}").format(
                            specific_output=specific_output,
                            specific_input=specific_input,
                            image=pipeline.docker_dict["TRUST3"],
                            align_bam=bam_file,
                            genome=genome.version,
                            param=dict_parameters[parameter],
                            num_threads=num_threads)

    elif pipeline.run_mode == "singularity":
        pass

    elif pipeline.run_mode == "local":
        TRUST3_script = ("{command} {path} -f {align_bam} -I 200 -g {genome} "
                        "-o {out} {param} -E -c -n {num_threads}").format(
                            command=classes.get_command(pipeline, "TRUST3"),
                            path=classes.get_path(pipeline, "TRUST3"),
                            align_bam=specific_input+"/"+bam_file,
                            genome=genome.version,
                            out=specific_output+"/",
                            param=dict_parameters[parameter],
                            num_threads=num_threads)
    
    pipeline.write_log(("Running TRUST3 {parameter} "
                        "for {sample}").format(parameter=parameter, sample=sample_name))
    return TRUST3_script

@bash_app
def MiXCR_align(pipeline, genome, sample_name, inputs=[]):
    import classes
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"
    num_threads = classes.get_threads(pipeline, "MiXCR")
    specific_input = pipeline.fastq_dict[sample_name][0].rsplit("/", 1)[0]
    if pipeline.run_mode == "docker":
        MiXCR_align_script = ("docker run --rm -v {specific_output}:/output "
                            "-v {specific_input}:/input {image} "
                            "mixcr align -s hs -p rna-seq -t {threads} "
                            "-OallowPartialAlignments=true "
                            "{fastq} /output/{vdjca}").format(specific_output=MiXCR_output,
                                                    specific_input=specific_input,
                                                    image=pipeline.docker_dict["MiXCR"],
                                                    threads=num_threads,
                                                    fastq=" ".join(["/input/"+path.rsplit("/", 1)[-1] for path in pipeline.fastq_dict[sample_name]]),
                                                    vdjca=sample_name+".vdjca")
    elif pipeline.run_mode == "singularity":
        pass

    elif pipeline.run_mode == "local":
        MiXCR_align_script = ("{command} mixcr align -s hs -p rna-seq -t {threads} "
                        "-OallowPartialAlignments=true "
                        "{fastq} {vdjca}").format(command=classes.get_command(pipeline, "MiXCR"),
                                                    path=classes.get_path(pipeline, "MiXCR"),
                                                    threads=pipeline.num_threads,
                                                    fastq=" ".join(pipeline.fastq_dict[sample_name]),
                                                    vdjca=MiXCR_output+sample_name+".vdjca")

    pipeline.write_log("Running MiXCR align for {sample}".format(sample=sample_name))
    return MiXCR_align_script

@bash_app
def MiXCR_assemblePartial(pipeline, genome, sample_name, rescue_num, inputs=[]):
    import classes
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"
    if pipeline.run_mode == "docker":
        MiXCR_assemblePartial_script = ("docker run --rm -v {specific_output}:/output "
                            "{image} mixcr assemblePartial /output/{vdjca} "
                            "/output/{rescued}").format(specific_output=MiXCR_output,
                                                    image=pipeline.docker_dict["MiXCR"],
                                                    vdjca=sample_name+(".vdjca" if rescue_num == "1" else f"_rescued{str(int(rescue_num)-1)}.vdjca"),
                                                    rescued=sample_name+"_rescued"+rescue_num+".vdjca")
    
    elif pipeline.run_mode == "singularity":
        pass

    elif pipeline.run_mode == "local":
        MiXCR_assemblePartial_script = ("{command} {path} assemblePartial {vcdja} "
                            "{rescued}").format(command=classes.get_command(pipeline, "MiXCR"),
                                                path=classes.get_path(pipeline, "MiXCR"),
                                                vcdja=MiXCR_output+sample_name+(".vdjca" if rescue_num == "1" else "_rescued1.vdjca"),
                                                rescued=MiXCR_output+sample_name+"_rescued"+rescue_num+".vdjca")

    pipeline.write_log("Running MiXCR assemblePartial {num} for {sample}".format(num=rescue_num,
                                                                            sample=sample_name))
    return MiXCR_assemblePartial_script

@bash_app
def MiXCR_extendAlignments(pipeline, genome, sample_name, rescue_num, inputs=[]):
    import classes
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"
    if pipeline.run_mode == "docker":
        MiXCR_extendAlignments_script = ("docker run --rm -v {specific_output}:/output "
                            "{image} mixcr extend /output/{rescued} "
                            "/output/{extended}").format(specific_output=MiXCR_output,
                                                    image=pipeline.docker_dict["MiXCR"],
                                                    rescued=sample_name+f"_rescued{rescue_num}.vdjca",
                                                    extended=sample_name+"_rescued_extended.vdjca")
    
    elif pipeline.run_mode == "singularity":
        pass

    elif pipeline.run_mode == "local":
        MiXCR_extendAlignments_script = ("{command} {path} extend {rescued} "
                                "{extended}").format(command=classes.get_command(pipeline, "MiXCR"),
                                                    path=classes.get_path(pipeline, "MiXCR"),
                                                    rescued=MiXCR_output+sample_name+f"_rescued{rescue_num}.vdjca",
                                                    extended=MiXCR_output+sample_name+"_rescued_extended.vdjca")

    pipeline.write_log("Running MiXCR extendAlignments for {sample}".format(sample=sample_name))
    return MiXCR_extendAlignments_script

@bash_app
def MiXCR_assemble(pipeline, genome, sample_name, inputs=[]):
    import classes
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"
    if pipeline.run_mode == "docker":
        MiXCR_assemble_script = ("docker run --rm -v {specific_output}:/output "
                            "{image} mixcr assemble /output/{extended} "
                            "/output/{clones}").format(specific_output=MiXCR_output,
                                                    image=pipeline.docker_dict["MiXCR"],
                                                    extended=sample_name+"_rescued_extended.vdjca",
                                                    clones=sample_name+"_clones.clns")
    
    elif pipeline.run_mode == "singularity":
        pass

    elif pipeline.run_mode == "local":
        MiXCR_assemble_script = ("{command} {path} assemble {extended} "
                                "{clones}").format(command=classes.get_command(pipeline, "MiXCR"),
                                                    path=classes.get_path(pipeline, "MiXCR"),
                                                    extended=MiXCR_output+sample_name+"_rescued_extended.vdjca",
                                                    clones=MiXCR_output+sample_name+"_clones.clns")

    pipeline.write_log("Running MiXCR assemble for {sample}".format(sample=sample_name))
    return MiXCR_assemble_script

@bash_app
def MiXCR_exportClones(pipeline, genome, sample_name, parameter, inputs=[]):
    import classes
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"
    dict_parameters = {"TCR": "TRA,TRB,TRG,TRD",
                        "BCR-heavy": "IGH",
                        "BCR-light": "IGL,IGK"}
    
    if pipeline.run_mode == "docker":
        MiXCR_exportClones_script = ("docker run --rm -v {specific_output}:/output "
                            "{image} mixcr exportClones -c {param} -o -t /output/{clones} "
                                "/output/{txt_clones}").format(specific_output=MiXCR_output,
                                                    image=pipeline.docker_dict["MiXCR"],
                                                    clones=sample_name+"_clones.clns",
                                                    param=dict_parameters[parameter],
                                                    txt_clones=parameter+"/"+sample_name+"_clones.txt")
    
    elif pipeline.run_mode == "singularity":
        pass

    elif pipeline.run_mode == "local":
        MiXCR_exportClones_script = ("{command} {path} exportClones -c {param} -o -t {clones} "
                                    "{txt_clones}").format(command=classes.get_command(pipeline, "MiXCR"),
                                                        path=classes.get_path(pipeline, "MiXCR"),
                                                        clones=MiXCR_output+sample_name+"_clones.clns",
                                                        param=dict_parameters[parameter],
                                                        txt_clones=MiXCR_output+parameter+"/"+sample_name+"_clones.txt")
    pipeline.write_log("Running MiXCR exportClones {param} for {sample}".format(param=parameter, sample=sample_name))
    return MiXCR_exportClones_script

@bash_app
def TRUST4(pipeline, genome, sample_name, inputs=[]):
    import classes
    # Urgent, just docker for now...
    num_threads = classes.get_threads(pipeline, "TRUST4")
    specific_input = pipeline.output_dir + sample_name + "/STAR_align/" + genome.version
    bam_file = sample_name + ".Aligned.sortedByCoord.out.bam"
    specific_output = pipeline.output_dir + sample_name + "/TRUST4"
    TRUST4_script = ("docker run --rm -v {input_dir}:/input "
                    "-v {output_dir}:/output {image} "
                    "run-trust4 -b /input/{bam_file} -f "
                    "/TRUST4/{genome_version}_bcrtcr.fa "
                    "--ref /TRUST4/human_IMGT+C.fa -t {threads} "
                    "-o /output/{output_prefix}").format(output_dir=specific_output,
                                                                input_dir=specific_input,
                                                                image=pipeline.docker_dict["TRUST4"],
                                                                threads=num_threads,
                                                                bam_file=bam_file,
                                                                genome_version=genome.version,
                                                                output_prefix=sample_name)
    pipeline.write_log("Running TRUST4 for {sample}".format(sample=sample_name))
    return TRUST4_script

# ADD Bamtools app + "sorted"

@python_app
def bamtools(pipeline, genome, sample_name, input=[]):
    import subprocess
    import exceptions
    import classes
    
    # docker container : biocontainers/bamtools
    # docker pull biocontainers/bamtools:v2.5.1dfsg-3-deb_cv1
    
    bam_file = (pipeline.bam_dir + sample_name + "/" + sample_name
                + ".Aligned.sortedByCoord.out.bam")
    bamtools_script = f"bamtools stats -insert -in {bam_file}"
    with open(pipeline.output_dir+sample_name+"/VDJer/"+sample_name+"_stats.txt", "w") as output:
        bamtools_sub = subprocess.Popen(bamtools_script.split(),
                    stdout=output,
                    stderr=subprocess.PIPE,
                    shell=True,
                    universal_newlines=True)
        _, errs = bamtools_sub.communicate()
        pipeline.write_out_err(errs)

@python_app
def VDJer(pipeline, genome, sample_name, inputs=[]):
    import subprocess
    import exceptions
    import classes

    # docker container (?): m0zack/vdjer:1.0 
    # downloads STAR, maybe too much??
    
    # Use /bin/bash -c for pipe output


    bam_file = (pipeline.bam_dir + sample_name + "/" + sample_name
                + ".Aligned.sortedByCoord.out.bam")
    output_prefix = pipeline.output_dir+sample_name+"/VDJer/"

    with open(output_prefix+sample_name+"_stats.txt", "r") as bam_data:
        line = bam_data.readline()
        while not line.startswith("Median insert size (absolute value):"):
            line = bam_data.readline()
        median_insert = str(int(float(line.split(":")[-1])))

    chain = ["replace this"]
    VDJer_script = ("vdjer --in {bam_file} --t 8 --ins {median_insert} "
                    "--chain {chain} --ref-dir {pipeline.VDJer_dir}vdjer_human_references/igh "
                    "> vdjer.sam 2> vdjer.log").format(bam_file=bam_file,
                                                    median_insert=median_insert,
                                                    chain=chain[0],
                                                    VDJer_dict=pipeline.VDJer_dir)
    pass

#For CATT: have TRA/TRB directories, or the TCR directory structure?
#I think separate is better... Might overlap file names (e.g. backup.jl)

@bash_app
def CATT(pipeline, genome, sample_name, parameter, inputs=[]):
    import classes
    import os
    specific_output = pipeline.output_dir + sample_name + "/CATT/" + parameter
    specific_input = pipeline.fastq_dict[sample_name][0].rsplit("/", 1)[0]
    UID = os.getuid()
    num_threads = classes.get_threads(pipeline, "CATT")

    if pipeline.run_mode == "docker":
        docker_str = ("docker run --rm -v {specific_output}:/output "
                        "-v {specific_input}:/input -w /output "
                        "-u {UID} {image} ").format(specific_output=specific_output,
                                                    specific_input=specific_input,
                                                    UID=UID,
                                                    image=pipeline.docker_dict["CATT"])
        if pipeline.single_end:
            CATT_script = docker_str + ("catt --chain {chain} -f /input/{fastq} "
                                        "-o /output/{out_name} -t {threads} --bowt {threads}").format(chain=parameter,
                                                            fastq=pipeline.fastq_dict[sample_name][0].rsplit("/", 1)[-1],
                                                            out_name=sample_name,
                                                            threads=num_threads)
        else:
            CATT_script = docker_str + ("catt --chain {chain} --f1 /input/{R1_fastq} --f2 /input/{R2_fastq} "
                                        "-o /output/{out_name} -t {threads} --bowt {threads}").format(chain=parameter,
                                                            R1_fastq=pipeline.fastq_dict[sample_name][0].rsplit("/", 1)[-1],
                                                            R2_fastq=pipeline.fastq_dict[sample_name][1].rsplit("/", 1)[-1],
                                                            out_name=sample_name,
                                                            threads=num_threads)
    
    elif pipeline.run_mode == "singularity":
        pass
    
    elif pipeline.run_mode == "local":
        if pipeline.single_end:
            CATT_script = ("{command} {path} --chain {chain} -f {fastq} "
                            "-o {out_name} -t {threads} --bowt {threads}").format(command=classes.get_command(pipeline, "MiXCR"),
                                                    path=classes.get_path(pipeline, "MiXCR"),
                                                    chain=parameter,
                                                    fastq=pipeline.fastq_dict[sample_name][0],
                                                    out_name=specific_output+"/"+sample_name+"_"+parameter,
                                                    threads=num_threads)
        else:
            CATT_script = ("{command} {path} --chain {chain} --f1 {R1_fastq} --f2 {R2_fastq} "
                            "-o {out_name} -t {threads} --bowt {threads}").format(command=classes.get_command(pipeline, "MiXCR"),
                                                    path=classes.get_path(pipeline, "MiXCR"),
                                                    chain=parameter,
                                                    R1_fastq=pipeline.fastq_dict[sample_name][0],
                                                    R2_fastq=pipeline.fastq_dict[sample_name][1],
                                                    out_name=specific_output+"/"+sample_name+"_"+parameter,
                                                    threads=num_threads)
    pipeline.write_log("Running CATT {chain} for {sample}".format(chain=parameter, sample=sample_name))
    return CATT_script







