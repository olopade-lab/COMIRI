import classes
from parsl import python_app, bash_app


# NEED TO REMOVE ALL DEPENDENCIES + PARAMETER CHECKING (SHOULD BE ALL IN runner.py)
# MiXCR: IGH, IGL, IGK, TRA, TRB, TRG, TRD
# TRUST: -B, -B -L

@bash_app(executors=['1core_worker'])
def SAM_index(pipeline, bam_or_fasta, command, inputs=[]):
    import classes
    SAM_index_script = ("samtools {command} "
                    "{bam_or_fasta}").format(bam_or_fasta=bam_or_fasta,
                                            command=command)
    pipeline.write_log("Generating index file for {f}".format(f=bam_or_fasta.split("/")[-1]))
    return SAM_index_script

@bash_app(executors=['4core_largeMEM'])
def STAR_align(pipeline, genome, sample_name, VDJer=False, inputs=[]):
    import classes
    if VDJer:
        genome_version = "hg38"
    else:
        genome_version = genome.version
    
    STAR_align_script = ("STAR --genomeDir {g_dir} "
                        "--readFilesIn {fastq} "
                        "--outSAMtype BAM SortedByCoordinate "
                        "--outFileNamePrefix {align_dir}. "
                        "--alignEndsType EndToEnd "
                        "--outSAMunmapped Within "
                        "--runThreadN {num_threads}").format(g_dir=genome.dir+genome_version,
                            fastq=" ".join(pipeline.fastq_dict[sample_name])
                            align_dir=pipeline.output_dir+sample_name+"/STAR_align/"+genome_version,
                            num_threads=pipeline.config_dict["STAR"]["threads"])
    pipeline.write_log("Alignment of {sample}".format(sample=sample_name))
    return STAR_align_script


# CHANGE TO individual app for each parameter !!
@bash_app(executors=['4core_worker'])
def TRUST3(pipeline, genome, sample_name, parameter, inputs=[]):
    import classes
    dict_parameters = {"TCR": "",
                        "BCR-heavy": "-B",
                        "BCR-light": "-B -L"}
    
    # list_param = param.split(",")
    # # Move this up to initial function
    # for i in range(len(list_param)):
    #     if list_param[i] == "B":
    #         list_param[i] = "-B"
    #     elif list_param[i] == "BL":
    #         list_param[i] = "-B -L"
    #     elif list_param[i] == "none":
    #         list_param[i] = ""
    #     else:
    #         message = ("The argument <{param}> "
    #             "is not allowed for TRUST").format(param=param)
    #         raise WrongArgumentError(message)

    bam_file = (pipeline.bam_dir + sample_name + "/" + sample_name
                + ".Aligned.sortedByCoord.out.bam")

    #have to wait for alignment to run this!

    # if not os.path.isfile(bam_file + ".bai"):
    #     
    #         ()
    #     inputs.append(SAM_index(pipeline, bam_file, "index"))

    specific_output = pipeline.output_dir + sample_name + "/TRUST3/" + parameter + "/"
    # Error handling for trust directory (check if directory exists way before?)
    # ext = p.replace(" ", "").replace("-", "_")
    # Need to be hg38 to match with VDJer?
    # check sample name: if there 
    # TRUST runs with python 2 !! : "module load gcc/6.2.0 python/2.7.13 trust && "
    TRUST3_script = ("trust -f {align_bam} -I 200 -g {genome} "
                    "-o {out} {param} -E -c -n {num_threads}").format(align_bam=bam_file,
                        genome=genome.version,
                        out=specific_output,
                        param=dict_parameters[parameter],
                        num_threads=pipeline.num_threads)
    pipeline.write_log(("Running TRUST3 {parameter} "
                        "for {sample}").format(parameter=parameter, sample=sample_name))
    return TRUST3_script

@bash_app(executors=['4core_worker'])
def MiXCR_align(pipeline, genome, sample_name, inputs=[]):
    import classes
    r1_cat = "replace later"
    r2_cat = "replace later"
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"
    MiXCR_align_script = ("mixcr align -s hs -p rna-seq -t {threads} "
                        "-OallowPartialAlignments=true "
                        "{r1} {r2} {vdjca}").format(threads=pipeline.num_threads,
                                                    r1=r1_cat,
                                                    r2=r2_cat,
                                                    vdjca=MiXCR_output + sample_name + ".vdjca")
    pipeline.write_log("Running MiXCR align for {sample}".format(sample=sample_name))
    return MiXCR_align_script

@bash_app(executors=['4core_worker'])
def MiXCR_assemblePartial(pipeline, genome, sample_name, rescue_num, inputs=[]):
    import classes
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"
    MiXCR_assemblePartial_script = ("mixcr assemblePartial {vcdja} "
                            "{rescued}").format(vcdja=MiXCR_output
                                                        + sample_name
                                                        + (".vdjca" if rescue_num == "1" else "_rescued1.vdjca"),
                                                rescued=MiXCR_output + sample_name + "_rescued" + rescue_num + ".vdjca")
    pipeline.write_log("Running MiXCR assemblePartial {num} for {sample}".format(num=rescue_num,
                                                                            sample=sample_name))
    return MiXCR_assemblePartial_script

@bash_app(executors=['4core_worker'])
def MiXCR_extendAlignments(pipeline, genome, sample_name, inputs=[]):
    import classes
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"
    MiXCR_extendAlignments_script = ("mixcr extendAlignments {rescued} "
                                "{extended}").format(rescued=MiXCR_output + sample_name + "_rescued2.vdjca",
                                                    extended=MiXCR_output + sample_name + "_rescued2_extended.vdjca")
    pipeline.write_log("Running MiXCR extendAlignments for {sample}".format(sample=sample_name))
    return MiXCR_extendAlignments_script

@bash_app(executors=['4core_worker'])
def MiXCR_assemble(pipeline, genome, sample_name, inputs=[]):
    import classes
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"
    MiXCR_assemble_script = ("mixcr assemble {extended} "
                                "{clones}").format(extended=MiXCR_output + sample_name + "_rescued2_extended.vdjca",
                                                    clones=MiXCR_output + sample_name + "_clones.clns")
    pipeline.write_log("Running MiXCR assemble for {sample}".format(sample=sample_name))
    return MiXCR_assemble_script

@bash_app(executors=['4core_worker'])
def MiXCR_exportClones(pipeline, genome, sample_name, parameter, inputs=[]):
    import classes
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"
    dict_parameters = {"TCR": "TRA,TRB,TRG,TRD",
                        "BCR-heavy": "IGH",
                        "BCR-light": "IGL,IGK"}
    specific_output = pipeline.output_dir + sample_name + "/MiXCR/" + parameter + "/"
    # Need to change the chain names to be user input (for B/T cells and long chain in TRUST):
    MiXCR_exportClones_script = ("mixcr exportClones -c {param} -o -t {clones} "
                                "{txt_clones}").format(clones=MiXCR_output + sample_name + "_clones.clns",
                                                    param=dict_parameters[parameter],
                                                    txt_clones=specific_output + sample_name + "_clones.txt")
    pipeline.write_log("Running MiXCR exportClones {param} for {sample}".format(param=parameter, sample=sample_name))
    return MiXCR_exportClones_script

@bash_app(executors=['1core_worker'])
def TRUST4(pipeline, genome, sample_name, inputs=[]):
    import classes
    bam_file = (pipeline.bam_dir + sample_name + "/" + sample_name
                + ".Aligned.sortedByCoord.out.bam")
    TRUST4_script = ("{TRUST4_dir}run-trust4 -b {bam_file} -f "
                    "{TRUST4_dir}{genome_version}_bcrtcr.fa "
                    "--ref {TRUST4_dir}human_IMGT+C.fa "
                    "-o {output_dir}").format(TRUST4_dir=pipeline.TRUST4_dir,
                                                                bam_file=bam_file,
                                                                genome_version=genome.version,
                                                                output_dir=pipeline.output_dir+sample_name+"/TRUST4/"+sample_name)
    pipeline.write_log("Running TRUST4 for {sample}".format(sample=sample_name))
    return TRUST4_script

# ADD Bamtools app + "sorted"

@python_app(executors=['1core_worker'])
def bamtools(pipeline, genome, sample_name, input=[]):
    import subprocess
    import exceptions
    import classes
    
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

@python_app(executors=['4core_worker'])
def VDJer(pipeline, genome, sample_name, inputs=[]):
    import subprocess
    import exceptions
    import classes

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


@bash_app(executors=['4core_worker'])
def CATT(pipeline, genome, sample_name, parameter, inputs=[]):
    import classes
    import os
    specific_output = pipeline.output_dir + sample_name + "/CATT/" + parameter
    specific_input = pipeline.fastq_dict[sample_name][0].rsplit("/", 1)[0]
    UID = os.getuid()
    if pipeline.run_mode == "docker":
        if pipeline.single_end:
            CATT_script = ("docker run --rm -v {specific_output}:/output -v {specific_input}:/input "
                            "-w /output -u {UID} guobioinfolab/catt catt -f {fastq} "
                            "-o {sample_name} -t {threads}").format(specific_output=specific_output,
                                                            specific_input=specific_input,
                                                            UID=UID,
                                                            fastq=pipeline.fastq_dict[sample_name][0],
                                                            sample_name=sample_name,
                                                            threads=pipeline.config_dict["CATT"]["threads"])
        else:
            CATT_script = ("docker run --rm -v {specific_output}:/output -v {specific_input}:/input "
                            "-w /output -u {UID} guobioinfolab/catt catt --f1 {R1_fastq} --f2 {R2_fastq} "
                            "-o {sample_name} -t {threads}").format(specific_output=specific_output,
                                                            specific_input=specific_input,
                                                            UID=UID,
                                                            R1_fastq=pipeline.fastq_dict[sample_name][0],
                                                            R2_fastq=pipeline.fastq_dict[sample_name][1],
                                                            sample_name=sample_name,
                                                            threads=pipeline.config_dict["CATT"]["threads"])
    elif pipeline.run_mode == "singularity":
        pass
    elif pipeline.run_mode == "local":
        pass
    return CATT_script







