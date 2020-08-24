from datetime import datetime
import exceptions
import os

# TCR,BCR-heavy,BCR-light
# TO ADD TO PIPELINE:
# r1r2_dir (fastq)
# num threads
# bam directory
# outpt

# get_state, set_state
# put pipeline path in the executor function

class Pipeline:
    def __init__(self):
        self.num_threads = 1
        self.run_mode = ["TCR", "BCR-heavy", "BCR-light"]
        self.run_program = ["MiXCR", "TRUST3", "TRUST4", "VDJer"]
        # self.align = False
    
    def time_stamp(self):
        now = datetime.now()
        return(now.strftime("%m/%d/%Y, %H:%M:%S"))

    def write_log(self, message):
        with open(self.output_dir + "TCR-BCR_log.txt", "a+") as log:
            log.write(self.time_stamp() + ":  " + message + "\n")
        
    def write_out_err(self, message):
        with open(self.output_dir + "TCR-BCR_out_err.txt", "a+") as err:
            err.write(message)
    
    def set_output_dir(self, output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except:
            raise exceptions.IncorrectPathError(output_dir)
        if output_dir[-1] != "/":
                output_dir += "/"
        self.output_dir = output_dir
    
    def set_run_mode(self, run_parameters_str):
        list_parameters = []
        for parameter in run_parameters_str.split(","):
            if parameter not in ["TCR", "BCR-heavy", "BCR-light"]:
                message = ("The following mode <{parameter}> "
                            "for --run_mode is incorrect.").format(parameter=parameter)
                raise exceptions.WrongArgumentError(message)
            else:
                list_parameters.append(parameter)
        self.run_mode = list_parameters
    
    def set_run_program(self, run_programs_str):
        list_programs = []
        for program in run_programs_str.split(","):
            if program not in ["MiXCR", "TRUST3", "TRUST4", "VDJer"]:
                message = ("The following program <{prog}> "
                            "for --run_program is incorrect.").format(prog=program)
                raise exceptions.WrongArgumentError(message)
            else:
                list_programs.append(program)
        self.run_program = list_programs

    def set_num_threads(self, num_threads_str):
        try:
            self.num_threads = int(num_threads_str)
        except:
            message = ("Error: the number of threads <{num_threads}> "
                    "is not an integer.").format(num_threads=num_threads_str)
            raise exceptions.WrongArgumentError(message)
    
    def set_fastq_dir(self, fastq_dir):
        if not os.path.isdir(fastq_dir):
            raise exceptions.IncorrectPathError(fastq_dir)
        if fastq_dir[-1] != "/":
            fastq_dir += "/"
        self.fastq_dir = fastq_dir

    def set_bam_dir(self, bam_dir):
        try:
            os.makedirs(bam_dir, exist_ok=True)
        except:
            raise exceptions.IncorrectPathError(bam_dir)
        if bam_dir[-1] != "/":
            bam_dir += "/"
        self.bam_dir = bam_dir

    def set_TRUST4_dir(self, TRUST4_dir):
        if not os.path.isdir(TRUST4_dir):
            raise exceptions.IncorrectPathError(TRUST4_dir)
        if TRUST4_dir[-1] != "/":
            TRUST4_dir += "/"
        self.TRUST4_dir = TRUST4_dir

    def set_VDJer_dir(self, VDJer_dir):
        if not os.path.isdir(VDJer_dir):
            raise exceptions.IncorrectPathError(VDJer_dir)
        if VDJer_dir[-1] != "/":
            VDJer_dir += "/"
        self.VDJer_dir = VDJer_dir

    def set_sample_names(self, sample_names_str):
        self.sample_names = sample_names_str.split(",")

    def __getstate__(self):
        return self.__dict__
    
    def __setstate__(self, d):
        self.__dict__ = d


class Genome:
    def __init__(self):
        self.version = "hg19"
        # self.generate = False
    
    def set_fasta(self, genome_fasta):
        if not os.path.isfile(genome_fasta):
            raise exceptions.IncorrectPathError(genome_fasta)
        self.fasta = genome_fasta

    def set_gtf(self, genome_gtf):
        if not os.path.isfile(genome_gtf):
            raise exceptions.IncorrectPathError(genome_gtf)
        self.gtf = genome_gtf
    
    def set_version(self, genome_version):
        if not genome_version in ["hg19", "hg38"]:
            message = "Error: the genome version must be hg19 or hg38."
            raise exceptions.WrongArgumentError(message)
        self.version = genome_version
    
    def set_dir(self, genome_dir):
        if not os.path.isdir(genome_dir):
            raise exceptions.IncorrectPathError(genome_dir)
        if genome_dir[-1] != "/":
            genome_dir += "/"
        self.dir = genome_dir
    # Cannot do genome generate here, since that will cause cyclical dependencies
    
    def __getstate__(self):
        return self.__dict__
    
    def __setstate__(self, d):
        self.__dict__ = d