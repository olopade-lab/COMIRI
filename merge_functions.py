import pandas as pd
import re
from Bio import pairwise2
import math

def verify_start_end(row, aa_start, aa_end, program):
    if not row['cdr3aa_'+program].startswith(tuple(aa_start)):
        sub_index = float("inf")
        for subsequence in aa_start:
            if row['cdr3aa_'+program].rfind(subsequence) != -1 and row['cdr3aa_'+program].rfind(subsequence) < sub_index:
                sub_index = row['cdr3aa_'+program].rfind(subsequence)
        if sub_index != float("inf"):
            # WHY +2 here??
            row['cdr3aa_'+program] = row['cdr3aa_'+program][sub_index+2:]
            row['cdr3dna_'+program] = row['cdr3dna_'+program][(sub_index+2)*3:]
        else:
            row['cdr3aa_'+program] = None
            row['cdr3dna_'+program] = None
            return row
    if not row['cdr3aa_'+program].endswith(tuple(aa_end)):
        sub_index = -1
        for subsequence in aa_end:
            if row['cdr3aa_'+program].rfind(subsequence) > sub_index:
                sub_index = row['cdr3aa_'+program].rfind(subsequence)
        if sub_index != -1:
            row['cdr3aa_'+program] = row['cdr3aa_'+program][:sub_index+2]
            row['cdr3dna_'+program] = row['cdr3dna_'+program][:(sub_index+2)*3]
        else:
            row['cdr3aa_'+program] = None
            row['cdr3dna_'+program] = None
            return row
    return row

def verify_canonical(row, program):
    match = re.search(r"C[^_]*[FW]", row['cdr3aa_'+program])
    if not match:
        row['cdr3aa_'+program] = None
        row['cdr3dna_'+program] = None
    else:
        row['cdr3aa_'+program] = match.group()
        row['cdr3dna_'+program] = row['cdr3dna_'+program][match.span()[0]*3:match.span()[1]*3]
    return row


def update_gene_name(list_genes):
    # if pd.isnull(list_genes):
    #     return list_genes
    # else:
    if list_genes == [""]:
        return []
    else:
        return list(map(lambda x: x.split("*")[0], list_genes))

def fix_CATT_gene_name(gene_str):
    if "|" in gene_str:
        return gene_str.split("|")[1]
    else:
        return gene_str

def pairwise2_list(list_prev, seq):
    align_score = -math.inf
    best_sequence = ""
    for prev_seq in list_prev:
        current_score = pairwise2.align.localxs(prev_seq, seq, -1, -1, score_only=True)
        if current_score > align_score:
            align_score = current_score
            best_sequence = prev_seq
    return (align_score, best_sequence)
    
def multiple_pairwise2(row_prev_seq, list_curr, list_columns):
    align_score = -math.inf
    best_prev = ""
    best_curr = ""
    list_curr = set(list_curr)
    for column in list_columns:
        # quick fix: check if instance of list
        if isinstance(row_prev_seq[column], list):
            list_prev = set(row_prev_seq[column])
            for prev_seq in list_prev:
                for curr_seq in list_curr:
                    current_score = pairwise2.align.localxs(prev_seq, curr_seq, -1, -1, score_only=True)
                    if current_score > align_score:
                        align_score = current_score
                        best_prev = prev_seq
                        best_curr = curr_seq
    if align_score > -math.inf:
        return (align_score, best_prev, best_curr)
    else:
        return float("NaN")


def merge_within(input_df, program):
    df_program = input_df.apply(verify_canonical, args=[program], axis=1)
    df_program.dropna(subset=['cdr3aa_'+program, 'cdr3dna_'+program], inplace=True)
    # CATT often has "None" for gene names
    if program == "CATT":
        df_program = df_program[(~(df_program['Vgene_'+program].isnull()) | 
                        ~(df_program['Jgene_'+program].isnull()) ) & 
                        (df_program['cdr3aa_'+program].str.len() > 6)]
    else:
        df_program = df_program[~(df_program['Vgene_'+program].isnull()) & 
                            ~(df_program['Jgene_'+program].isnull()) & 
                            (df_program['cdr3aa_'+program].str.len() > 6)]
    df_program_merged = pd.DataFrame(columns=['count_'+program,
                                    'Vgene_'+program,
                                    'Dgene_'+program,
                                    'Jgene_'+program,
                                    'Cgene_'+program,
                                    'cdr3dna_'+program,
                                    'cdr3aa_'+program])
    for _, program_row in df_program.iterrows():
        mapping = df_program_merged['cdr3dna_'+program].apply(lambda cdr3dna_program: pairwise2_list(cdr3dna_program, program_row['cdr3dna_'+program]))
        mapping_values = mapping.map(lambda result: result[0] if not pd.isna(result) else result)
        mapping_sequences = mapping.map(lambda result: result[1] if not pd.isna(result) else result)
        if not mapping_values.dropna().empty:
            if mapping_values.max() >= (min(len(mapping_sequences[mapping_values.idxmax()]), len(program_row['cdr3dna_'+program])) - 2):
                df_program_merged.iloc[mapping_values.idxmax()]['cdr3dna_'+program].append(program_row['cdr3dna_'+program])
                df_program_merged.iloc[mapping_values.idxmax()]['cdr3aa_'+program].append(program_row['cdr3aa_'+program])
                df_program_merged.iloc[mapping_values.idxmax()]['count_'+program].append(program_row['count_'+program])
                for gene in ['V', 'J', 'D', 'C']:
                    df_program_merged.iloc[mapping_values.idxmax(),
                                            df_program_merged.columns.get_loc(gene+'gene_'+program)] = list(set(df_program_merged.iloc[mapping_values.idxmax()][gene+'gene_'+program]
                                                                                                                + program_row[gene+'gene_'+program]))
            else:
                program_row['count_'+program] = [program_row['count_'+program]]
                program_row['cdr3aa_'+program] = [program_row['cdr3aa_'+program]]
                program_row['cdr3dna_'+program] = [program_row['cdr3dna_'+program]]
                df_program_merged = df_program_merged.append(program_row, ignore_index=True)
        else:
            program_row['count_'+program] = [program_row['count_'+program]]
            program_row['cdr3aa_'+program] = [program_row['cdr3aa_'+program]]
            program_row['cdr3dna_'+program] = [program_row['cdr3dna_'+program]]
            df_program_merged = df_program_merged.append(program_row, ignore_index=True)
    return df_program_merged

def merge_other(df_merged, df_program, program):
    list_columns = [column for column in df_merged.columns 
                        if (column.startswith("cdr3dna") and column not in df_program.columns)]
    list_columns = [column for column in list_columns if not df_merged[column].isnull().values.all()]
    for _, program_row in df_program.iterrows():
        mapping = df_merged.apply(lambda cdr3dna_program: multiple_pairwise2(cdr3dna_program, program_row['cdr3dna_'+program], list_columns), axis=1)
        mapping_values = mapping.map(lambda result: result[0] if not pd.isna(result) else result)
        mapping_prev = mapping.map(lambda result: result[1] if not pd.isna(result) else result)
        mapping_curr = mapping.map(lambda result: result[2] if not pd.isna(result) else result)
        max_row = mapping_values.idxmax()
        if mapping_values[max_row] >= (min(len(mapping_prev[max_row]), len(mapping_curr[max_row])) - 2):
            for column in df_program.columns:  
                df_merged.iloc[max_row][column] = program_row[column]
        else:
            df_merged = df_merged.append(program_row, ignore_index=True)
    return df_merged

def create_df(list_programs):
    list_columns = [['count_'+program,
                    'Vgene_'+program,
                    'Dgene_'+program,
                    'Jgene_'+program,
                    'Cgene_'+program,
                    'cdr3dna_'+program,
                    'cdr3aa_'+program] for program in list_programs]
    list_columns = [column for sublist in list_columns for column in sublist]
    list_columns.sort(reverse=True)
    return pd.DataFrame(columns=list_columns)