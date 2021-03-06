import os
import glob
import parsl
from parsl import load
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SingleNodeLauncher
from parsl.providers import TorqueProvider
from parsl.addresses import address_by_hostname
from parsl import python_app

# NEED to check which chain for comparing the two cdr3 sequences

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
                worker_init=("module load gcc/6.2.0 miniconda3/4.7.10; "
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

load(config)


@python_app
def merge_results(sample_dir, attribute):
    import pandas as pd
    import merge_functions
    import os
    sample_name = sample_dir.split("/")[-1]
    dict_TRUST3 = {"TCR": ".Aligned.sortedByCoord.out.bam-TCR-ALL.fa",
                    "BCR-heavy": ".Aligned.sortedByCoord.out.bam-BCR-HeavyChain.fa",
                    "BCR-light": ".Aligned.sortedByCoord.out.bam-BCR-LightChain.fa"}
    dict_TRUST4 = {"TCR": "T",
                    "BCR-heavy": "IGH",
                    "BCR-light": ("IGL", "IGK")}
    df_merged = merge_functions.create_df(["MiXCR", "TRUST3", "TRUST4", "CATT"])
    for program in os.listdir(sample_dir):
        if program == "MiXCR":
            try:
                df_MiXCR = pd.read_csv(sample_dir+"/MiXCR/"+attribute+"/"
                                        +sample_name+"_clones.txt",
                                        sep="\t",
                                        usecols=["cloneCount",
                                                "allVHitsWithScore",
                                                "allDHitsWithScore",
                                                "allJHitsWithScore",
                                                "allCHitsWithScore",
                                                "nSeqCDR3",
                                                "aaSeqCDR3"],
                                        header=0,
                                        dtype="object")
                df_MiXCR.columns = ["count_MiXCR",
                                    "Vgene_MiXCR",
                                    "Dgene_MiXCR",
                                    "Jgene_MiXCR",
                                    "Cgene_MiXCR",
                                    "cdr3dna_MiXCR",
                                    "cdr3aa_MiXCR"]
                for gene in ["Vgene", "Dgene", "Jgene", "Cgene"]:
                    df_MiXCR[gene+"_MiXCR"] = df_MiXCR[gene+"_MiXCR"].str.split(",")
                    df_MiXCR.loc[df_MiXCR[gene+"_MiXCR"].isnull(),gene+"_MiXCR"] = df_MiXCR.loc[df_MiXCR[gene+"_MiXCR"].isnull(),gene+"_MiXCR"].apply(lambda x: [])
                    df_MiXCR[gene+"_MiXCR"] = df_MiXCR[gene+"_MiXCR"].apply(merge_functions.update_gene_name)
                df_MiXCR_merged = merge_functions.merge_within(df_MiXCR, "MiXCR")
                if df_merged.empty:
                    df_merged = df_merged.append(df_MiXCR_merged, ignore_index=True)
                else:
                    df_merged = merge_functions.merge_other(df_merged, df_MiXCR_merged, "MiXCR")
            except OSError:
                print("MiXCR {attribute} directory for {sample} does not have a clones.txt file".format(attribute=attribute,
                                                                                                        sample=sample_name))

        if program == "TRUST3":
            try:
                with open(sample_dir+"/TRUST3/"+attribute+"/"+sample_name+dict_TRUST3[attribute], "r") as fasta_file:
                    lines_fasta = fasta_file.readlines()
                    list_output = []
                    for line in lines_fasta:
                        if line.startswith(">"):
                            list_line = [element.replace(">", "").replace("|", "").replace("\n", "").split("=")[-1]
                                            for element in line.split("+")
                                            if (not element.startswith("est_clonal_exp=")
                                                and not element.startswith("seq_length=")
                                                and not element.startswith("minus_log_Eval=")
                                                and not element.startswith("est_lib_size="))]
                            list_reportgene = list(map(lambda x: x.split("*")[0].replace("/", ""), list_line[4].split("_")))
                            list_line[2] = list(set([gene.split("*")[0].replace("/", "") for gene in list_line[2].split("_") if gene != ""]
                                                + list(filter(lambda gene: gene[3] == 'V', list_reportgene))))
                            list_line[3] = list(set([gene.split("*")[0].replace("/", "") for gene in list_line[3].split("_") if gene != ""]
                                                + list(filter(lambda gene: gene[3] == 'J', list_reportgene))))
                            # if list_line[3] == [""]:
                            #     list_line[3] = float("NaN")
                            list_line.pop(4)
                            list_line = list_line[1:] + [[], []]
                            list_output.append(list_line)
                    df_TRUST3 = pd.DataFrame(list_output,
                                            columns=['count_TRUST3',
                                                    'Vgene_TRUST3',
                                                    'Jgene_TRUST3',
                                                    'cdr3aa_TRUST3',
                                                    'cdr3dna_TRUST3',
                                                    'Dgene_TRUST3',
                                                    'Cgene_TRUST3'],
                                            dtype="object")
                df_TRUST3_merged = merge_functions.merge_within(df_TRUST3, "TRUST3")
                if df_merged.empty:
                    df_merged = df_merged.append(df_TRUST3_merged, ignore_index=True)
                else:
                    df_merged = merge_functions.merge_other(df_merged, df_TRUST3_merged, "TRUST3")
            except OSError:
                print("TRUST3 {attribute} directory for {sample} does not have a {filetype} file".format(attribute=attribute,
                                                                                                        sample=sample_name,
                                                                                                        filetype=dict_TRUST3[attribute]))

        if program == "TRUST4":
            try:
                print(sample_dir+"/TRUST4/"+sample_name+"_report.tsv")
                df_TRUST4_all = pd.read_csv(sample_dir+"/TRUST4/"+sample_name+"_report.tsv",
                                            sep="\t",
                                            usecols=['#count',
                                                    'CDR3nt',
                                                    'CDR3aa',
                                                    'V',
                                                    'D',
                                                    'J',
                                                    'C'])
                df_TRUST4 = df_TRUST4_all[df_TRUST4_all[['V','J','C','D']].apply(lambda x: x.str.startswith(dict_TRUST4[attribute])).any(axis=1)]
                df_TRUST4 = df_TRUST4[(df_TRUST4.V != "*")
                                & (df_TRUST4.J != "*")
                                & (df_TRUST4.CDR3aa != "out_of_frame")
                                & (df_TRUST4.CDR3aa != "partial")]
                df_TRUST4.columns = ["count_TRUST4",
                                    "cdr3dna_TRUST4",
                                    "cdr3aa_TRUST4",
                                    "Vgene_TRUST4",
                                    "Dgene_TRUST4",
                                    "Jgene_TRUST4",
                                    "Cgene_TRUST4"]
                for gene in ["Vgene", "Dgene", "Jgene", "Cgene"]:
                    df_TRUST4[gene+"_TRUST4"].replace(to_replace="*", value="", inplace=True)
                    df_TRUST4[gene+"_TRUST4"] = df_TRUST4[gene+"_TRUST4"].str.split(",")
                    # df_TRUST4.loc[df_TRUST4[gene+"_TRUST4"]==["*"],gene+"_TRUST4"] = df_TRUST4.loc[df_TRUST4[gene+"_TRUST4"].isnull(),gene+"_TRUST4"].apply(lambda x: [])
                    # df_TRUST4[gene+"_TRUST4"].replace(to_replace=["*"], value=[], inplace=True)
                    df_TRUST4[gene+"_TRUST4"] = df_TRUST4[gene+"_TRUST4"].apply(merge_functions.update_gene_name)
                df_TRUST4_merged = merge_functions.merge_within(df_TRUST4, "TRUST4")
                if df_merged.empty:
                    df_merged = df_merged.append(df_TRUST4_merged, ignore_index=True)
                else:
                    df_merged = merge_functions.merge_other(df_merged, df_TRUST4_merged, "TRUST4")
            except OSError:
                print("TRUST4 {attribute} directory for {sample} does not have a _report.tsv file".format(attribute=attribute,
                                                                                                        sample=sample_name))

        if program == "CATT":
            if attribute == "TCR":
                try:
                    df_CATT_TRA = pd.read_csv(sample_dir+"/CATT/"+sample_name+".TRA.CDR3.CATT.csv",
                                                sep=",",
                                                usecols=['AAseq',
                                                        'NNseq',
                                                        'Vregion',
                                                        'Jregion',
                                                        'Frequency'])
                    df_CATT_TRA.columns = ["cdr3aa_CATT",
                                        "cdr3dna_CATT",
                                        "Vgene_CATT",
                                        "Jgene_CATT",
                                        "count_CATT"]
                    df_CATT_TRA["Dgene_CATT"] = ""
                    df_CATT_TRA["Cgene_CATT"] = ""
                    df_CATT_TRB = pd.read_csv(sample_dir+"/CATT/"+sample_name+".TRB.CDR3.CATT.csv",
                                                sep=",",
                                                usecols=['AAseq',
                                                        'NNseq',
                                                        'Vregion',
                                                        'Jregion',
                                                        'Dregion',
                                                        'Frequency'])
                    df_CATT_TRA.columns = ["cdr3aa_CATT",
                                        "cdr3dna_CATT",
                                        "Vgene_CATT",
                                        "Jgene_CATT",
                                        "Dgene_CATT",
                                        "count_CATT"]
                    df_CATT_TRB["Cgene_CATT"] = ""
                    df_CATT = pd.concat([df_CATT_TRA, df_CATT_TRB])
                    for gene in ["Vgene", "Dgene", "Jgene", "Cgene"]:
                        df_CATT[gene+"_CATT"] = df_CATT[gene+"_CATT"].apply(merge_functions.fix_CATT_gene_name)
                        df_CATT[gene+"_CATT"] = df_CATT[gene+"_CATT"].str.split(",")
                        df_CATT[gene+"_CATT"] = df_CATT[gene+"_CATT"].apply(merge_functions.update_gene_name)
                    df_CATT_merged = merge_functions.merge_within(df_CATT, "CATT")
                    if df_merged.empty:
                        df_merged = df_merged.append(df_CATT_merged, ignore_index=True)
                    else:
                        df_merged = merge_functions.merge_other(df_merged, df_CATT_merged, "CATT")
                except OSError:
                    print("CATT {attribute} directory for {sample} does not have a .CDR3.CATT.csv file".format(attribute=attribute,
                                                                                                            sample=sample_name))

    for column_name in df_merged.columns:
        df_merged[column_name] = df_merged[column_name].apply(lambda ls: ','.join(map(str, ls)) if isinstance(ls, list) else ls)
    df_merged["sample_name"] = sample_name
    df_merged.to_csv(sample_dir+"/"+sample_name+"_"+attribute+".csv", index=False)


for sample_dir in glob.glob("/scratch/jbreynier/NG_RNAseq_BCR-TCR/*/"):
    sample_dir = sample_dir[:-1]
    dict_TRUST3 = {"TCR": ".Aligned.sortedByCoord.out.bam-TCR-ALL.fa",
                    "BCR-heavy": ".Aligned.sortedByCoord.out.bam-BCR-HeavyChain.fa",
                    "BCR-light": ".Aligned.sortedByCoord.out.bam-BCR-LightChain.fa"}
    dict_TRUST4 = {"TCR": "T",
                    "BCR-heavy": "IGH",
                    "BCR-light": ("IGL", "IGK")}
    for attribute in ["TCR"]:
        if not os.path.isfile(sample_dir+"/"+sample_dir.split("/")[-1]+"_"+attribute+".csv"):
            merge_results(sample_dir, attribute)

parsl.wait_for_current_tasks()