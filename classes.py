from datetime import datetime
import exceptions
import os
from wcmatch import glob as wcmatchglob
from parsl import bash_app
import json
import multiprocessing


# TCR,BCR-heavy,BCR-light
# TO ADD TO PIPELINE:
# r1r2_dir (fastq)
# num threads
# bam directory
# outpt

# get_state, set_state
# put pipeline path in the executor function


@bash_app
def aws_STARGenome(url, output_dir):
    return f"aws s3 sync {url} {output_dir} --no-sign-request"

class Pipeline:

    def __init__(self, run_program, run_mode, receptor, single_end):
        self.run_program = run_program
        self.run_mode = run_mode
        self.receptor = receptor
        self.single_end = single_end
        self.rescue_num = 2
        self.docker_dict = {
            "MiXCR": "milaboratory/mixcr@sha256:26284253404ccc467ca39c60660b80eabe6f3bfa24063fc454860fa98c646c92",
            "TRUST3": "mgibio/trust@sha256:283c1102334a0c3fb51dd25b6cd053ece28bc30748045984ce190778588c7242",
            "TRUST4": "olopadelab/trust4@sha256:c05fb0b75268f558a039a4adac3a935c921b52159877ea45d9532f686703f1e4",
            "CATT": "guobioinfolab/catt:1.8.1",
            "VDJer": "",
            "STAR": "mgibio/star:2.5.2b",
            "bamtools": "",
            "samtools": "mgibio/samtools:1.9"
        }

    def time_stamp(self):
        now = datetime.now()
        return(now.strftime("%m/%d/%Y, %H:%M:%S"))

    def write_log(self, message):
        with open(self.output_dir + "CDR3merge_log.txt", "a+") as log:
            log.write(self.time_stamp() + ":  " + message + "\n")
    
    def set_output_dir(self, output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except:
            raise exceptions.IncorrectPathError(output_dir)
        if output_dir[-1] != "/":
                output_dir += "/"
        self.output_dir = output_dir
    
    def set_fastq_dict(self, fastq_path):
        try:
            base_dir = fastq_path.rsplit("/", 2)[0] + "/"
            name_format = fastq_path.rsplit("/", 1)[-1]
        except:
            raise exceptions.IncorrectPathError(fastq_path)
        list_files = wcmatchglob.glob(fastq_path, flags=wcmatchglob.BRACE)
        if not list_files:
            raise exceptions.IncorrectPathError(fastq_path)
        list_samples = list(set([file.split("/")[-2] for file in list_files]))
        self.fastq_dict = {}
        if self.single_end:
            for sample in list_samples:
                fastq_specific = wcmatchglob.glob(base_dir+sample+"/"+name_format, flags=wcmatchglob.BRACE)
                if len(fastq_specific) > 1:
                    print(f"More than 1 matching fastq file for sample {sample} in single_end mode. Sample ignored.")
                else:
                    self.fastq_dict[sample] = fastq_specific
        else:
            for sample in list_samples:
                fastq_specific = wcmatchglob.glob(base_dir+sample+"/"+name_format, flags=wcmatchglob.BRACE)
                if len(fastq_specific) > 2:
                    print(f"More than 2 matching fastq files for sample {sample} in paired_end mode. Sample ignored.")
                elif len(fastq_specific) < 2:
                    print(f"Only than 1 matching fastq file for sample {sample} in paired_end mode. Sample ignored.")
                else:
                    self.fastq_dict[sample] = fastq_specific
    
    def extract_json(self, json_name):
        if json_name is None:
            self.config_dict = {}
        else:
            base_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1])
            json_file = os.path.join(base_dir, 'configs', '{}.json'.format(json_name))
            if not os.path.isfile(json_file):
                raise exceptions.IncorrectPathError("Cannot find the json file <{json_name}.json>.".format(json_name=json_name))
            try:
                json_fileobject = open(json_file, "r")
                self.config_dict = json.load(json_fileobject)
            except Exception as e:
                raise exceptions.IncorrectInputFiles(("Could not load json from <{json_name}.json> :"
                                                    "\n {exception}.").format(json_name=json_name, exception=e))

    def __getstate__(self):
        return self.__dict__
    
    def __setstate__(self, d):
        self.__dict__ = d


class Genome:
    def __init__(self, genome_version):
        self.version = genome_version
    
    def download_genome(self, output_dir, run_VDJer):
        list_futures = []
        self.dir = output_dir + "genomes/"
        if self.version == 'hg19':
            specific_dir = self.dir + "hg19/"
            os.makedirs(specific_dir, exist_ok=True)
            if not os.path.isfile(specific_dir+"genome.fa"):
                list_futures.append(aws_STARGenome("s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/STARIndex/", specific_dir))
            else:
                print("Non-empty hg19 genome directory: hg19 already downloaded.")
        if self.version == 'hg38' or run_VDJer:
            specific_dir = self.dir + "hg38/"
            os.makedirs(specific_dir, exist_ok=True)
            if not os.path.isfile(specific_dir+"genome.fa"):
                list_futures.append(aws_STARGenome("s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh38/Sequence/STARIndex/", specific_dir))
            else:
                print("Non-empty hg38 genome directory: hg38 already downloaded.")
        return list_futures
    
    def __getstate__(self):
        return self.__dict__
    
    def __setstate__(self, d):
        self.__dict__ = d

def get_threads(pipeline, program):
    if program in pipeline.config_dict.keys():
        if "threads" in pipeline.config_dict[program].keys():
            num_threads = pipeline.config_dict[program]["threads"]
        else:
            num_threads = os.environ.get('PARSL_CORES', multiprocessing.cpu_count())
    else:
        num_threads = os.environ.get('PARSL_CORES', multiprocessing.cpu_count())
    pipeline.write_log("getting threads"+str(num_threads))
    return num_threads

def get_command(pipeline, program):
    if program in pipeline.config_dict.keys():
        if "command" in pipeline.config_dict[program].keys():
            command = pipeline.config_dict[program]["command"] + " &"
        else:
            command = ""
    else:
        command = ""
    return command

def get_path(pipeline, program):
    dict_path = {
        "MiXCR": "mixcr",
        "TRUST3": "trust",
        "TRUST4": "run-trust4",
        "CATT": "catt",
        "VDJer": "vdjer",
        "STAR": "STAR",
        "samtools": "samtools",
        "bamtools": "bamtools"
    }

    if program in pipeline.config_dict.keys():
        if "path" in pipeline.config_dict[program].keys():
            path = pipeline.config_dict[program]["path"]
        else:
            path = dict_path[program]
    else:
        path = dict_path[program]
    return path