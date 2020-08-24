import os
import subprocess
import glob
import re
import classes
import exceptions
from parsl import python_app, bash_app


# NEED TO REMOVE ALL DEPENDENCIES + PARAMETER CHECKING (SHOULD BE ALL IN runner.py)
# MiXCR: IGH, IGL, IGK, TRA, TRB, TRG, TRD
# TRUST: -B, -B -L


@python_app(executors=['4core_largeMEM'])
def STAR_genome(genome, pipeline, inputs=[]):
    import subprocess
    import exceptions
    import classes
    
    # waits for star index to complete:
    # if not os.path.isfile(genome.fasta + ".fai"):
    #     SAM_index(pipeline, genome.fasta, "faidx").result()
    STAR_genome_script = ("STAR --runMode genomeGenerate "
                        "--genomeDir {g_dir} "
                        "--genomeFastaFiles {g_fa} "
                        "--sjdbGTFfile {g_gtf} "
                        "--runThreadN {num_threads}").format(g_dir=genome.dir,
                                                        g_fa=genome.fasta,
                                                        g_gtf=genome.gtf,
                                                        num_threads=pipeline.num_threads)
    pipeline.write_log("Generating genome")
    # return bash_STAR_genome(STAR_genome_script, outputs=[genome.dir])
    STAR_genome_sub = subprocess.Popen(STAR_genome_script.split(),
                    stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE,
                    universal_newlines = True)
    outs, errs = STAR_genome_sub.communicate()
    pipeline.write_out_err(outs)
    pipeline.write_out_err(errs)
    
    if STAR_genome_sub.returncode != 0:
        pipeline.write_log("ERROR")
        raise exceptions.SubprocessError("STAR genome generate")
    return genome.dir

# @bash_app
# def bash_STAR_genome(script, stderr='std.err', stdout='std.out', outputs=[]):
#     return script

@python_app(executors=['1core_worker'])
def SAM_index(pipeline, bam_or_fasta, command, inputs=[]):
    import subprocess
    import exceptions
    import classes
    
    SAM_index_script = ("samtools {command} "
                    "{bam_or_fasta}").format(bam_or_fasta=bam_or_fasta,
                                            command=command)
    pipeline.write_log("Generating index file for {f}".format(f=bam_or_fasta.split("/")[-1]))
    # bash_SAM_index(SAM_index_script)
    SAM_index_sub = subprocess.Popen(SAM_index_script.split(),
                    stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE,
                    universal_newlines = True)
    outs, errs = SAM_index_sub.communicate()
    pipeline.write_out_err(outs)
    pipeline.write_out_err(errs)
    
    if SAM_index_sub.returncode != 0:
        pipeline.write_log("ERROR")
        raise exceptions.SubprocessError("SAM index {f}".format(f=bam_or_fasta.split("/")[-1]))

# @bash_app
# def bash_SAM_index(script, stderr='std.err', stdout='std.out', outputs=[]):
#     return script


# No need to parallelize gunzip?
@python_app(executors=['1core_worker'])
def gunzip(pipeline, list_gz):
    import subprocess
    import exceptions
    import classes
    # DO for loop sample
    for gz_file in list_gz:
        gunzip_sub = subprocess.Popen(["gunzip", gz_file],
                        stdout = subprocess.PIPE,
                        stderr = subprocess.PIPE,
                        universal_newlines = True)
        outs, errs = gunzip_sub.communicate()
        pipeline.write_out_err(outs)
        pipeline.write_out_err(errs)
        if gunzip_sub.returncode != 0:
            pipeline.write_log("ERROR")
            raise exceptions.SubprocessError("gunzip")

@python_app(executors=['4core_largeMEM'])
def STAR_align(genome, pipeline, sample_name, inputs=[]):
    import subprocess
    import exceptions
    import classes
    import glob
    import re

    # place this somewhere else?
    # For now, just ignore gunzip...
    # list_gz = list(filter(lambda x: x.endswith(".gz"),
    #                         glob.glob(pipeline.fastq_dir + sample_name + "/*")))
    # for gz_file in list_gz:
    #     gunzip(pipeline, gz_file)

    list_r1 = glob.glob(pipeline.fastq_dir + sample_name + "/*_R1_*")
    list_r2 = glob.glob(pipeline.fastq_dir + sample_name + "/*_R2_*")

    list_r1 = list(filter(lambda x: "_L" in x, list_r1))
    list_r2 = list(filter(lambda x: "_L" in x, list_r2))

    # Check this somewhere else?
    if len(list_r1) != len(list_r2):
        message = ("Error: the directory contains "
                "{lenr1} R1 file(s) and {lenr2} R2 "
                "file(s).").format(lenr1=len(list_r1), lenr2=len(list_r2))
        raise exceptions.IncorrectInputFiles(message)

    list_r1.sort(key = lambda x: re.search('_L(.+?)_R1', x).group(1))
    list_r2.sort(key = lambda x: re.search('_L(.+?)_R2', x).group(1))
    
    r1 = ','.join(list_r1)
    r2 = ','.join(list_r2)

    # if os.path.exists(align_dir):
    #     shutil.rmtree(align_dir)
    # os.makedirs(align_dir)

    # Maybe add a script to check if STAR alignment files are there already?

    STAR_align_script = ("STAR --genomeDir {g_dir} "
                        "--readFilesIn {r1} {r2} "
                        "--outSAMtype BAM SortedByCoordinate "
                        "--outFileNamePrefix {align_dir}. "
                        "--alignEndsType EndToEnd "
                        "--outSAMunmapped Within "
                        "--runThreadN {num_threads}").format(g_dir=genome.dir,
                            r1=r1,
                            r2=r2,
                            align_dir=pipeline.bam_dir + sample_name + "/" + sample_name,
                            num_threads=pipeline.num_threads)
    pipeline.write_log("Alignment of {sample}".format(sample=sample_name))
    STAR_align_sub = subprocess.Popen(STAR_align_script.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True)
    outs, errs = STAR_align_sub.communicate()
    pipeline.write_out_err(outs)
    pipeline.write_out_err(errs)
    if STAR_align_sub.returncode != 0:
        pipeline.write_log("ERROR")
        raise exceptions.SubprocessError("STAR align {sample}".format(sample=sample_name))

@python_app(executors=['4core_worker'])
def TRUST3(pipeline, genome, sample_name, list_parameters, inputs=[]):
    import subprocess
    import exceptions
    import classes
    import os
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

    TRUST_output = pipeline.output_dir + sample_name + "/TRUST3/"
    # Error handling for trust directory (check if directory exists way before?)
    for parameter in list_parameters:
        # ext = p.replace(" ", "").replace("-", "_")
        specific_output = TRUST_output + parameter + "/"
        # Need to be hg38 to match with VDJer?
        # check sample name: if there 
        TRUST_script = ("module load gcc/6.2.0 python/2.7.13 trust && "
                        "trust -f {align_bam} -I 200 -g {genome} "
                        "-o {out} {param} -E -c -n {num_threads}").format(align_bam=bam_file,
                            genome=genome.version,
                            out=specific_output,
                            param=dict_parameters[parameter],
                            num_threads=pipeline.num_threads)
        pipeline.write_log(("Running TRUST3 {parameter} "
                            "for {sample}").format(parameter=parameter, sample=sample_name))
        TRUST_sub = subprocess.Popen(TRUST_script,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    shell=True,
                    universal_newlines=True)
        outs, errs = TRUST_sub.communicate()
        pipeline.write_out_err(outs)
        pipeline.write_out_err(errs)
        if TRUST_sub.returncode != 0:
            pipeline.write_log("ERROR")
            raise exceptions.SubprocessError("TRUST3 {parameter} {sample}".format(parameter=parameter,
                                                                                    sample=sample_name))

@python_app(executors=['4core_worker'])
def MiXCR(pipeline, genome, sample_name, inputs=[]):
    import subprocess
    import exceptions
    import classes
    import os
    import glob
    
    
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"
    # Why am I checking for "_L"? Ofc cause I'm adding concatenated file in same directory!!
    # Actually maybe not since cat fastq files were sent to results file
    list_r1 = list(filter(lambda x: "_L" in x, glob.glob(pipeline.fastq_dir + sample_name + "/*_R1_*")))
    list_r2 = list(filter(lambda x: "_L" in x, glob.glob(pipeline.fastq_dir + sample_name + "/*_R2_*")))
    list_r1.sort(key = lambda x: re.search('_L(.+?)_R1', x).group(1))
    list_r2.sort(key = lambda x: re.search('_L(.+?)_R2', x).group(1))
    # Duplicate code!
    if len(list_r1) != len(list_r2):
        message = ("Error: the directory contains "
                "{lenr1} R1 file(s) and {lenr2} R2 "
                "file(s).").format(lenr1=len(list_r1), lenr2=len(list_r2))
        raise exceptions.IncorrectInputFiles(message)
    
    # Check number of R1 and R2 files, if 1 each, no need for cat
    r1_cat = pipeline.fastq_dir + sample_name + "/" + sample_name + "_R1.fastq"
    r2_cat = pipeline.fastq_dir + sample_name + "/" + sample_name + "_R2.fastq"
    if not os.path.isfile(r1_cat):
        with open(r1_cat, "w") as cat_fastq_r1:
            cat_sub_r1 = subprocess.Popen(["cat"] + list_r1,
                        stdout = cat_fastq_r1,
                        stderr = subprocess.PIPE,
                        universal_newlines = True)
            _, errs = cat_sub_r1.communicate()
            pipeline.write_out_err(errs)
        if cat_sub_r1.returncode != 0:
            pipeline.write_log("ERROR")
            raise exceptions.SubprocessError("cat fastq files {sample}".format(sample=sample_name))
    if not os.path.isfile(r2_cat):
        with open(r2_cat, "w") as cat_fastq_r2:
            cat_sub_r2 = subprocess.Popen(["cat"] + list_r2,
                        stdout = cat_fastq_r2,
                        stderr = subprocess.PIPE,
                        universal_newlines = True)
            _, errs = cat_sub_r2.communicate()
            pipeline.write_out_err(errs)
        
        if cat_sub_r2.returncode != 0:
            pipeline.write_log("ERROR")
            raise exceptions.SubprocessError("cat fastq files {sample}".format(sample=sample_name))
    
    mixcr_align_script = ("mixcr align -s hs -p rna-seq -t {threads} "
                        "-OallowPartialAlignments=true "
                        "{r1} {r2} {vdjca}").format(threads=pipeline.num_threads,
                                                    r1=r1_cat,
                                                    r2=r2_cat,
                                                    vdjca=MiXCR_output + sample_name + ".vdjca")
    pipeline.write_log("Running MiXCR align for {sample}".format(sample=sample_name))
    mixcr_align_sub = subprocess.Popen(mixcr_align_script.split(),
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                universal_newlines = True)
    outs, errs = mixcr_align_sub.communicate()
    pipeline.write_out_err(outs)
    pipeline.write_out_err(errs)
    if mixcr_align_sub.returncode != 0:
        pipeline.write_log("ERROR")
        raise exceptions.SubprocessError("MiXCR align {sample}".format(sample=sample_name))
    
    # Change "rescue_num" to be user input?
    for rescue_num in ["1", "2"]:
        mixcr_assemblePartial_script = ("mixcr assemblePartial {vcdja} "
                                "{rescued}").format(vcdja=MiXCR_output
                                                            + sample_name
                                                            + (".vdjca" if rescue_num == "1" else "_rescued1.vdjca"),
                                                    rescued=MiXCR_output + sample_name + "_rescued" + rescue_num + ".vdjca")
        pipeline.write_log("Running MiXCR assemblePartial {num} for {sample}".format(num=rescue_num,
                                                                                sample=sample_name))
        mixcr_assemblePartial_sub = subprocess.Popen(mixcr_assemblePartial_script.split(),
                    stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE,
                    universal_newlines = True)
        outs, errs = mixcr_assemblePartial_sub.communicate()
        pipeline.write_out_err(outs)
        pipeline.write_out_err(errs)
        if mixcr_assemblePartial_sub.returncode != 0:
            pipeline.write_log("ERROR")
            raise exceptions.SubprocessError("MiXCR assemblePartial {num} for {sample}".format(num=rescue_num,
                                                                                sample=sample_name))
    
    mixcr_extend_script = ("mixcr extendAlignments {rescued} "
                                "{extended}").format(rescued=MiXCR_output + sample_name + "_rescued2.vdjca",
                                                    extended=MiXCR_output + sample_name + "_rescued2_extended.vdjca")
    pipeline.write_log("Running MiXCR extendAlignments for {sample}".format(sample=sample_name))
    mixcr_extend_sub = subprocess.Popen(mixcr_extend_script.split(),
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                universal_newlines = True)
    outs, errs = mixcr_extend_sub.communicate()
    pipeline.write_out_err(outs)
    pipeline.write_out_err(errs)
    if mixcr_extend_sub.returncode != 0:
        pipeline.write_log("ERROR")
        raise exceptions.SubprocessError("MiXCR extendAlignments {sample}".format(sample=sample_name))

    mixcr_assemble_script = ("mixcr assemble {extended} "
                                "{clones}").format(extended=MiXCR_output + sample_name + "_rescued2_extended.vdjca",
                                                    clones=MiXCR_output + sample_name + "_clones.clns")
    
    pipeline.write_log("Running MiXCR assemble for {sample}".format(sample=sample_name))
    mixcr_assemble_sub = subprocess.Popen(mixcr_assemble_script.split(),
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                universal_newlines = True)
    outs, errs = mixcr_assemble_sub.communicate()
    pipeline.write_out_err(outs)
    pipeline.write_out_err(errs)
    if mixcr_assemble_sub.returncode != 0:
        pipeline.write_log("ERROR")
        raise exceptions.SubprocessError("MiXCR assemble {sample}".format(sample=sample_name))

@python_app(executors=['4core_worker'])
def MiXCR_specific(pipeline, genome, sample_name, run_parameters, inputs=[]):
    import subprocess
    import exceptions
    import classes
    
    
    dict_parameters = {"TCR": "TRA,TRB,TRG,TRD",
                        "BCR-heavy": "IGH",
                        "BCR-light": "IGL,IGK"}
    MiXCR_output = pipeline.output_dir + sample_name + "/MiXCR/"

    for parameter in run_parameters:
        specific_output = MiXCR_output + parameter + "/"
        # Need to change the chain names to be user input (for B/T cells and long chain in TRUST):
        mixcr_export_script = ("mixcr exportClones -c {param} -o -t {clones} "
                                    "{txt_clones}").format(clones=MiXCR_output + sample_name + "_clones.clns",
                                                        param=dict_parameters[parameter],
                                                        txt_clones=specific_output + sample_name + "_clones.txt")
        pipeline.write_log("Running MiXCR exportClones {param} for {sample}".format(param=parameter, sample=sample_name))
        mixcr_export_sub = subprocess.Popen(mixcr_export_script.split(),
                    stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE,
                    universal_newlines = True)
        outs, errs = mixcr_export_sub.communicate()
        pipeline.write_out_err(outs)
        pipeline.write_out_err(errs)
        if mixcr_export_sub.returncode != 0:
            pipeline.write_log("ERROR")
            raise exceptions.SubprocessError("MiXCR exportClones {param} {sample}".format(param=parameter, sample=sample_name))

@python_app(executors=['1core_worker'])
def TRUST4(pipeline, genome, sample_name, inputs=[]):
    import subprocess
    import exceptions
    import classes

    bam_file = (pipeline.bam_dir + sample_name + "/" + sample_name
                + ".Aligned.sortedByCoord.out.bam")
    TRUST4_script = ("module load gcc/6.2.0 zlib/1.2.11 bzip2/1.0.6 "
                    "xz/5.2.2 perl/5.24.0 python/3.6.0 samtools/1.9 "
                    "glibc/2.14 && {TRUST4_dir}run-trust4 -b {bam_file} -f "
                    "{TRUST4_dir}{genome_version}_bcrtcr.fa "
                    "--ref {TRUST4_dir}human_IMGT+C.fa "
                    "-o {output_dir}").format(TRUST4_dir=pipeline.TRUST4_dir,
                                                                bam_file=bam_file,
                                                                genome_version=genome.version,
                                                                output_dir=pipeline.output_dir+sample_name+"/TRUST4/"+sample_name)
    pipeline.write_log("Running TRUST4 for {sample}".format(sample=sample_name))
    TRUST4_sub = subprocess.Popen(TRUST4_script,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=True,
                universal_newlines=True)
    outs, errs = TRUST4_sub.communicate()
    pipeline.write_out_err(outs)
    pipeline.write_out_err(errs)
    if TRUST4_sub.returncode != 0:
        pipeline.write_log("ERROR")
        raise exceptions.SubprocessError("TRUST4 {sample}".format(sample=sample_name))

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
    if bamtools_sub.returncode != 0:
        pipeline.write_log("ERROR")
        raise exceptions.SubprocessError("bamtools stats {sample}".format(sample=sample_name))

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
    
    for parameter in run_parameters:

    VDJer_script = (f"vdjer --in {bam_file} --t 8 --ins {median_insert} "
                    f"--chain {chain[0]} --ref-dir {pipeline.VDJer_dir}vdjer_human_references/igh > vdjer.sam 2> vdjer.log
    pass
