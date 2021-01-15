import os
import glob
import pandas as pd
import re

# NEED to check which chain for comparing the two cdr3 sequences

def merge_results(sample_dir, attribute):
    sample_name = sample_dir.split("/")[-1]
    dict_TRUST3 = {"TCR": ".Aligned.sortedByCoord.out.bam-TCR-ALL.fa",
                    "BCR-heavy": ".Aligned.sortedByCoord.out.bam-BCR-HeavyChain.fa",
                    "BCR-light": ".Aligned.sortedByCoord.out.bam-BCR-LightChain.fa"}
    dict_TRUST4 = {"TCR": "T",
                    "BCR-heavy": "IGH",
                    "BCR-light": ("IGL", "IGK")}
    df_merged = pd.DataFrame()
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
                    df_MiXCR[gene+"_MiXCR"] = df_MiXCR[gene+"_MiXCR"].apply(update_gene_name)
                df_MiXCR_merged = merge_within(df_MiXCR, "MiXCR")
                if df_merged.empty:
                    df_merged = df_merged.append(df_MiXCR_merged, ignore_index=True)
                else:
                    df_merged = merge_other(df_merged, df_MiXCR_merged)
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
                df_TRUST3_merged = merge_within(df_TRUST3, "TRUST3")
                if df_merged.empty:
                    df_merged = df_merged.append(df_TRUST3_merged, ignore_index=True)
                else:
                    df_merged = merge_other(df_merged, df_TRUST3_merged)
            except OSError:
                print("TRUST3 {attribute} directory for {sample} does not have a {filetype} file".format(attribute=attribute,
                                                                                                        sample=sample_name,
                                                                                                        filetype=dict_TRUST3[attribute]))

        if program == "TRUST4":
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
                    df_TRUST4[gene+"_TRUST4"] = df_TRUST4[gene+"_TRUST4"].apply(update_gene_name)
                df_TRUST4_merged = merge_within(df_TRUST4, "TRUST4")
                if df_merged.empty:
                    df_merged = df_merged.append(df_TRUST4_merged, ignore_index=True)
                else:
                    df_merged = merge_other(df_merged, df_TRUST4_merged)
            except OSError:
                print("TRUST4 {attribute} directory for {sample} does not have a _report.tsv file".format(attribute=attribute,
                                                                                                        sample=sample_name))

        if program == "CATT":
            if attribute in ["TCR", "BCR-heavy"]:
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
                        df_CATT_TRB.columns = ["cdr3aa_CATT",
                                            "cdr3dna_CATT",
                                            "Vgene_CATT",
                                            "Jgene_CATT",
                                            "Dgene_CATT",
                                            "count_CATT"]
                        df_CATT_TRB["Cgene_CATT"] = ""
                        df_CATT = pd.concat([df_CATT_TRA, df_CATT_TRB])
                    except OSError:
                        print("CATT {attribute} directory for {sample} does not have a proper TRA and/or TRB .CDR3.CATT.csv file".format(attribute=attribute,
                                                                                                                sample=sample_name))
                else:
                    try:
                        df_CATT = pd.read_csv(sample_dir+"/CATT/IGH/"+sample_name+".IGH.CDR3.CATT.csv",
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
                    except OSError:
                        print("CATT {attribute} directory for {sample} does not have a proper .CDR3.CATT.csv file".format(attribute=attribute,
                                                                                                                sample=sample_name))
                # FIX THIS: df_CATT referenced even if file doesn't exist

                for gene in ["Vgene", "Dgene", "Jgene", "Cgene"]:
                    df_CATT[gene+"_CATT"] = df_CATT[gene+"_CATT"].apply(fix_CATT_gene_name)
                    df_CATT[gene+"_CATT"] = df_CATT[gene+"_CATT"].str.split(",")
                    df_CATT[gene+"_CATT"] = df_CATT[gene+"_CATT"].apply(update_gene_name)
                df_CATT_merged = merge_within(df_CATT, "CATT")
                if df_merged.empty:
                    df_merged = df_merged.append(df_CATT_merged, ignore_index=True)
                else:
                    df_merged = merge_other(df_merged, df_CATT_merged)

    for column_name in df_merged.columns:
        df_merged[column_name] = df_merged[column_name].apply(lambda ls: ','.join(map(str, ls)) if isinstance(ls, list) else ls)
    df_merged["sample_name"] = sample_name
    df_merged.drop(labels="cdr3dna", axis=1, inplace=True)
    df_merged.to_csv(sample_dir+"/"+sample_name+"_"+attribute+"_simple.csv", index=False)

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

def update_gene_name(list_genes):
    if list_genes == [""]:
        return []
    else:
        return list(map(lambda x: x.split("*")[0], list_genes))

def fix_CATT_gene_name(gene_str):
    if gene_str == "None":
        gene_str = ""
    if "|" in gene_str:
        return gene_str.split("|")[1]
    else:
        return gene_str

def verify_canonical(row, program):
    match = re.search(r"C[^_]*[FW]", row['cdr3aa_'+program])
    if not match:
        row['cdr3aa_'+program] = None
        row['cdr3dna_'+program] = None
    else:
        row['cdr3aa_'+program] = match.group()
        row['cdr3dna_'+program] = row['cdr3dna_'+program][match.span()[0]*3:match.span()[1]*3]
    return row

def merge_within(input_df, program):
    df_program = input_df.apply(verify_canonical, args=[program], axis=1)
    df_program.dropna(subset=['cdr3aa_'+program, 'cdr3dna_'+program], inplace=True)
    if program == "CATT":
        df_program = df_program[(~(df_program['Vgene_'+program].isnull()) | 
                        ~(df_program['Jgene_'+program].isnull()) ) & 
                        (df_program['cdr3aa_'+program].str.len() > 6)]
    else:
        df_program = df_program[~(df_program['Vgene_'+program].isnull()) & 
                            ~(df_program['Jgene_'+program].isnull()) & 
                            (df_program['cdr3aa_'+program].str.len() > 6)]
    # Required, otherwise sum() skips over "nuisance columns"
    dict_apply = {'count_'+program: "sum",
                    'Vgene_'+program: "sum",
                    'Dgene_'+program: "sum",
                    'Jgene_'+program: "sum",
                    'Cgene_'+program: "sum"
    }
    df_program = df_program.groupby(['cdr3dna_'+program, 'cdr3aa_'+program], as_index=False).agg(dict_apply)
    for gene in ['Vgene_', 'Dgene_', 'Jgene_', 'Cgene_']:
        df_program[gene+program] = df_program[gene+program].apply(lambda row: list(set(row)))
    df_program['cdr3dna'] = df_program['cdr3dna_'+program]
    return df_program

def merge_other(df_merged, df_program):
    return df_merged.merge(df_program, on="cdr3dna", how="outer")

# More complex version, if want to keep dna and aa separate for merging:

# def merge_rows(groupby_df, program):
#     groupby
#     groupby_df = groupby_df.sum()
#     for col in ['Vgene_', 'Dgene_', 'Jgene_', 'Cgene_']:
#         groupby_df[col+program] = set(groupby_df[col+program])
#     return groupby_df

for sample_dir in glob.glob("/cephfs/users/jbreynier/CDR3merge_outputs/FZ-116/"):
    sample_dir = sample_dir[:-1]
    for attribute in ["BCR-heavy"]:
        if not os.path.isfile(sample_dir+"/"+sample_dir.split("/")[-1]+"_"+attribute+"_simple.csv"):
            merge_results(sample_dir, attribute)