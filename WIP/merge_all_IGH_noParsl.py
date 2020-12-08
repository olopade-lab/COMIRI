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


# IMPORTANT: need to convert all nulls to empty 
# Check whether I should always include the 

# NEED TO DO: final editing of output
# Start/end amino acids for BCR heavy + light

# Fix this: ask Yonglan which is the best way to merge?
# def verify_start_end(row, aa_start, aa_end, program):
#     if not row['cdr3aa_'+program].startswith(tuple(aa_start)):
#         sub_index = float("inf")
#         for subsequence in aa_start:
#             if row['cdr3aa_'+program].rfind(subsequence) != -1 and row['cdr3aa_'+program].rfind(subsequence) < sub_index:
#                 sub_index = row['cdr3aa_'+program].rfind(subsequence)
#         if sub_index != float("inf"):
#             # WHY +2 here??
#             row['cdr3aa_'+program] = row['cdr3aa_'+program][sub_index+2:]
#             row['cdr3dna_'+program] = row['cdr3dna_'+program][(sub_index+2)*3:]
#         else:
#             row['cdr3aa_'+program] = None
#             row['cdr3dna_'+program] = None
#             return row
#     if not row['cdr3aa_'+program].endswith(tuple(aa_end)):
#         sub_index = -1
#         for subsequence in aa_end:
#             if row['cdr3aa_'+program].rfind(subsequence) > sub_index:
#                 sub_index = row['cdr3aa_'+program].rfind(subsequence)
#         if sub_index != -1:
#             row['cdr3aa_'+program] = row['cdr3aa_'+program][:sub_index+2]
#             row['cdr3dna_'+program] = row['cdr3dna_'+program][:(sub_index+2)*3]
#         else:
#             row['cdr3aa_'+program] = None
#             row['cdr3dna_'+program] = None
#             return row
#     return row

# def verify_canonical(row, program):
#     match = re.search(r"C[^_]*[FW]", row['cdr3aa_'+program])
#     if not match:
#         row['cdr3aa_'+program] = None
#         row['cdr3dna_'+program] = None
#     else:
#         row['cdr3aa_'+program] = match.group()
#         row['cdr3dna_'+program] = row['cdr3dna_'+program][match.span()[0]*3:match.span()[1]*3]
#     return row


# def update_gene_name(list_genes):
#     # if pd.isnull(list_genes):
#     #     return list_genes
#     # else:
#     if list_genes == [""]:
#         return []
#     else:
#         return list(map(lambda x: x.split("*")[0], list_genes))

# def pairwise2_list(list_prev, seq):
#     align_score = -math.inf
#     best_sequence = ""
#     for prev_seq in list_prev:
#         current_score = pairwise2.align.localxs(prev_seq, seq, -1, -1, score_only=True)
#         if current_score > align_score:
#             align_score = current_score
#             best_sequence = prev_seq
#     return (align_score, best_sequence)
    
# def multiple_pairwise2(row_prev_seq, list_curr, list_columns):
#     align_score = -math.inf
#     best_prev = ""
#     best_curr = ""
#     list_curr = set(list_curr)
#     for column in list_columns:
#         # quick fix: check if instance of list
#         if isinstance(row_prev_seq[column], list):
#             list_prev = set(row_prev_seq[column])
#             for prev_seq in list_prev:
#                 for curr_seq in list_curr:
#                     current_score = pairwise2.align.localxs(prev_seq, curr_seq, -1, -1, score_only=True)
#                     if current_score > align_score:
#                         align_score = current_score
#                         best_prev = prev_seq
#                         best_curr = curr_seq
#     if align_score > -math.inf:
#         return (align_score, best_prev, best_curr)
#     else:
#         return float("NaN")


# def merge_within(input_df, program):
#     df_program = input_df.apply(verify_canonical, args=[program], axis=1)
#     df_program.dropna(subset=['cdr3aa_'+program, 'cdr3dna_'+program], inplace=True)
#     df_program = df_program[~(df_program['Vgene_'+program].isnull()) & 
#                         ~(df_program['Jgene_'+program].isnull()) & 
#                         (df_program['cdr3aa_'+program].str.len() > 6)]
#     df_program_merged = pd.DataFrame(columns=['count_'+program,
#                                     'Vgene_'+program,
#                                     'Dgene_'+program,
#                                     'Jgene_'+program,
#                                     'Cgene_'+program,
#                                     'cdr3dna_'+program,
#                                     'cdr3aa_'+program])
#     for _, program_row in df_program.iterrows():
#         mapping = df_program_merged['cdr3dna_'+program].apply(lambda cdr3dna_program: pairwise2_list(cdr3dna_program, program_row['cdr3dna_'+program]))
#         mapping_values = mapping.map(lambda result: result[0] if not pd.isna(result) else result)
#         mapping_sequences = mapping.map(lambda result: result[1] if not pd.isna(result) else result)
#         if not mapping_values.dropna().empty:
#             if mapping_values.max() >= (min(len(mapping_sequences[mapping_values.idxmax()]), len(program_row['cdr3dna_'+program])) - 2):
#                 df_program_merged.iloc[mapping_values.idxmax()]['cdr3dna_'+program].append(program_row['cdr3dna_'+program])
#                 df_program_merged.iloc[mapping_values.idxmax()]['cdr3aa_'+program].append(program_row['cdr3aa_'+program])
#                 df_program_merged.iloc[mapping_values.idxmax()]['count_'+program].append(program_row['count_'+program])
#                 for gene in ['V', 'J', 'D', 'C']:
#                     df_program_merged.iloc[mapping_values.idxmax(),
#                                             df_program_merged.columns.get_loc(gene+'gene_'+program)] = list(set(df_program_merged.iloc[mapping_values.idxmax()][gene+'gene_'+program]
#                                                                                                                 + program_row[gene+'gene_'+program]))
#             else:
#                 program_row['count_'+program] = [program_row['count_'+program]]
#                 program_row['cdr3aa_'+program] = [program_row['cdr3aa_'+program]]
#                 program_row['cdr3dna_'+program] = [program_row['cdr3dna_'+program]]
#                 df_program_merged = df_program_merged.append(program_row, ignore_index=True)
#         else:
#             program_row['count_'+program] = [program_row['count_'+program]]
#             program_row['cdr3aa_'+program] = [program_row['cdr3aa_'+program]]
#             program_row['cdr3dna_'+program] = [program_row['cdr3dna_'+program]]
#             df_program_merged = df_program_merged.append(program_row, ignore_index=True)
#     return df_program_merged

# def merge_other(df_merged, df_program, program):
#     list_columns = [column for column in df_merged.columns 
#                         if (column.startswith("cdr3dna") and column not in df_program.columns)]
#     list_columns = [column for column in list_columns if not df_merged[column].isnull().values.all()]
#     for _, program_row in df_program.iterrows():
#         mapping = df_merged.apply(lambda cdr3dna_program: multiple_pairwise2(cdr3dna_program, program_row['cdr3dna_'+program], list_columns), axis=1)
#         mapping_values = mapping.map(lambda result: result[0] if not pd.isna(result) else result)
#         mapping_prev = mapping.map(lambda result: result[1] if not pd.isna(result) else result)
#         mapping_curr = mapping.map(lambda result: result[2] if not pd.isna(result) else result)
#         max_row = mapping_values.idxmax()
#         if mapping_values[max_row] >= (min(len(mapping_prev[max_row]), len(mapping_curr[max_row])) - 2):
#             for column in df_program.columns:  
#                 df_merged.iloc[max_row][column] = program_row[column]
#         else:
#             df_merged = df_merged.append(program_row, ignore_index=True)
#     return df_merged

# def create_df(list_programs):
#     list_columns = [['count_'+program,
#                     'Vgene_'+program,
#                     'Dgene_'+program,
#                     'Jgene_'+program,
#                     'Cgene_'+program,
#                     'cdr3dna_'+program,
#                     'cdr3aa_'+program] for program in list_programs]
#     list_columns = [column for sublist in list_columns for column in sublist]
#     list_columns.sort(reverse=True)
#     return pd.DataFrame(columns=list_columns)


# list_start_aa = ["CA", "CG", "CI", "CL", "CS", "CT", "CV"]
# list_end_aa = ["AF", "FF", "HF", "IF", "IW", "LF", "MF", "NF", "QF", "RF", "SF", "TF", "VF", "YF"]

# Change to user defined inputs: sample names, data type
# Change to empty list instead of null

# config = Config(
#     executors=[
#         HighThroughputExecutor(
#             cores_per_worker=1,
#             mem_per_worker=12,
#             max_workers=4,
#             worker_debug=True,
#             address=address_by_hostname(),
#             provider=TorqueProvider(
#                 launcher=SingleNodeLauncher(),
#                 worker_init=("module load gcc/6.2.0 miniconda3/4.7.10; "
#                             "source activate parsl_env; export PYTHONPATH='{}:{{PYTHONPATH}}'").format(os.getcwd()),
#                 init_blocks=1,
#                 max_blocks=10,
#                 min_blocks=0,
#                 nodes_per_block=1,
#                 walltime='99:00:00',
#                 scheduler_options='#PBS -l mem=60gb,nodes=1:ppn=4'
#             ),
#         ),
#     ],
#     checkpoint_mode='task_exit'
# )

# load(config)

#@python_app
def merge_results(sample_dir, attribute):
    import pandas as pd
    import merge_functions
    import os
    sample_name = sample_dir.split("/")[-1]
    dict_TRUST3 = {"TCR": ".Aligned.sortedByCoord.out.bam-TCR-ALL.fa",
                    "TRA": ".Aligned.sortedByCoord.out.bam-TCR-ALL.fa",
                    "TRB": ".Aligned.sortedByCoord.out.bam-TCR-ALL.fa",
                    "BCR-heavy": ".Aligned.sortedByCoord.out.bam-BCR-HeavyChain.fa",
                    "BCR-light": ".Aligned.sortedByCoord.out.bam-BCR-LightChain.fa"}
    dict_TRUST4 = {"TCR": "T",
                    "TRA": "TRA",
                    "TRB": "TRB",
                    "BCR-heavy": "IGH",
                    "BCR-light": ("IGL", "IGK")}
    dict_CATT = {"TRA": "TRA",
                "TRB": "TRB",
                "BCR-heavy": "IGH"
    }
    df_merged = merge_functions.create_df(["MiXCR", "TRUST3", "TRUST4", "CATT"])
    for program in os.listdir(sample_dir):
        if program == "MiXCR":
            print("MiXCR")
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
            print("TRUST3")
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
                            list_genes = list_line[1] + list_line[2]
                            if [x for x in list_genes if x.startswith(attribute)]:
                                list_output.append(list_line)
                    print(list_output)
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
            print("TRUST4")
            try:
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
            print("CATT")
            if dict_CATT[attribute] in ["TRA", "IGH"]:
                try:
                    df_CATT = pd.read_csv(sample_dir+"/CATT/"+dict_CATT[attribute]+"/"+sample_name+"."+dict_CATT[attribute]+".CDR3.CATT.csv",
                                                sep=",",
                                                usecols=['AAseq',
                                                        'NNseq',
                                                        'Vregion',
                                                        'Jregion',
                                                        'Frequency'])
                    df_CATT.columns = ["cdr3aa_CATT",
                                        "cdr3dna_CATT",
                                        "Vgene_CATT",
                                        "Jgene_CATT",
                                        "count_CATT"]
                    df_CATT["Dgene_CATT"] = ""
                    df_CATT["Cgene_CATT"] = ""
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
            elif dict_CATT[attribute] == "TRB":
                try:
                    df_CATT = pd.read_csv(sample_dir+"/CATT/"+attribute+"/"+sample_name+".TRB.CDR3.CATT.csv",
                                                sep=",",
                                                usecols=['AAseq',
                                                        'NNseq',
                                                        'Vregion',
                                                        'Jregion',
                                                        'Dregion',
                                                        'Frequency'])
                    df_CATT.columns = ["cdr3aa_CATT",
                                        "cdr3dna_CATT",
                                        "Vgene_CATT",
                                        "Jgene_CATT",
                                        "Dgene_CATT",
                                        "count_CATT"]
                    df_CATT["Cgene_CATT"] = ""
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


for sample_dir in glob.glob("/cephfs/users/jbreynier/CDR3merge_outputs/FZ-83/"):
    sample_dir = sample_dir[:-1]
    dict_TRUST3 = {"TCR": ".Aligned.sortedByCoord.out.bam-TCR-ALL.fa",
                    "TRA": ".Aligned.sortedByCoord.out.bam-TCR-ALL.fa",
                    "TRB": ".Aligned.sortedByCoord.out.bam-TCR-ALL.fa",
                    "BCR-heavy": ".Aligned.sortedByCoord.out.bam-BCR-HeavyChain.fa",
                    "BCR-light": ".Aligned.sortedByCoord.out.bam-BCR-LightChain.fa"}
    dict_TRUST4 = {"TCR": "T",
                    "BCR-heavy": "IGH",
                    "BCR-light": ("IGL", "IGK")}
    for attribute in ["BCR-heavy"]:
        if not os.path.isfile(sample_dir+"/"+sample_dir.split("/")[-1]+"_"+attribute+".csv"):
            merge_results(sample_dir, attribute)
