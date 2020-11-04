import pandas as pd
import statistics
from Bio.Seq import Seq
from collections import Counter
import random
# choose between sequences by looking at count if there is a tie


def consensus(row):
    new_row = {"sample_name": row["sample_name"]}
    for col in ["cdr3dna", "Vgene", "Jgene", "Cgene", "Dgene"]:
        all_results = []
        for program in ["MiXCR", "TRUST3", "TRUST4", "CATT"]:
            if isinstance(row[col+"_"+program], list):
                all_results += row[col+"_"+program]
        if not len(all_results) == 0:
            counter = Counter(all_results)
            _,val = counter.most_common(1)[0]
            modes = [x for x,y in counter.items() if y == val]
            # freqs = groupby(Counter(all_results).most_common(), lambda x:x[1])
            # modes = [val for val,count in next(freqs)[1]]
            if col.endswith("gene"):
                # simplified_gene_list = list(set([gene.split("-")[0] for gene in modes]))
                simplified_gene_list = modes
                if len(simplified_gene_list) == 1:
                    new_row[col] = simplified_gene_list
                else:
                    simplified_gene_list_filtered = list(filter(lambda x: x[-1].isdigit(), simplified_gene_list))
                    if len(simplified_gene_list) == 0:
                        random.shuffle(simplified_gene_list)
                        new_row[col] = simplified_gene_list
                    else:
                        random.shuffle(simplified_gene_list_filtered)
                        new_row[col] = simplified_gene_list_filtered
            else:
                if len(modes) > 1:
                    list_len = [len(seq) for seq in modes]
                    if list_len.count(max(list_len)) == 1:
                        new_row[col] = [modes[list_len.index(max(list_len))]]
                        print(new_row[col])
                    else:
                        all_count = []
                        total_count_mode = []
                        for program in ["MiXCR", "TRUST3", "TRUST4", "CATT"]:
                            if isinstance(row["count_"+program], list):
                                all_count += [int(count) for count in row["count_"+program]]
                        for mode in modes:
                            mode_indices = [i for i, seq in enumerate(all_results) if seq == mode]
                            mode_count = 0
                            for index in mode_indices:
                                mode_count += all_count[index]
                            total_count_mode.append(mode_count)
                        if total_count_mode.count(max(total_count_mode)) == 1:
                            new_row[col] = [modes[list_len.index(max(list_len))]]
                        else:
                            final_list_modes = []
                            count_indices = [i for i, count in enumerate(total_count_mode) if count == max(total_count_mode)]
                            for index in count_indices:
                                final_list_modes.append(modes[index])
                            print(f"Ambiguous choice of sequence for sample {new_row['sample_name']}, options are: {str(final_list_modes)}")
                            new_row[col] = [random.choice(final_list_modes)]
                else:
                    new_row[col] = modes
                new_row['cdr3aa'] = str(Seq(new_row[col][0]).translate())
        else:
            new_row[col] = float("NaN")
    all_count = []
    for program in ["MiXCR", "TRUST3", "TRUST4", "CATT"]:
        if isinstance(row["count_"+program], list):
            all_count += [sum([int(count) for count in row["count_"+program]])]
    new_row["count"] = int(round(sum(all_count) / len(all_count), 0))
    return pd.Series(new_row)

def quorum(row, list_algorithms, num_votes):
    return ["Y" if not row["cdr3dna_"+algorithm].isna() else "N" for algorithm in list_algorithms].count("Y") >= num_votes


df_merged = pd.read_csv("/Users/jbreynier/results_benchmarking/RNASeq_SPX6730-2_TRA.csv")
df_merged.loc[:, df_merged.columns != 'sample_name'] = df_merged.loc[:, df_merged.columns != 'sample_name'].apply(lambda x: x.str.split(","), axis=1)

df_merged = df_merged[(~df_merged["cdr3dna_MiXCR"].isna() & ~df_merged["cdr3dna_TRUST3"].isna()) | 
                        (~df_merged["cdr3dna_MiXCR"].isna() & ~df_merged["cdr3dna_TRUST4"].isna()) |
                        (~df_merged["cdr3dna_TRUST3"].isna() & ~df_merged["cdr3dna_TRUST4"].isna())]
df_consensus = df_merged.apply(consensus, axis=1)
for column_name in df_consensus.columns:
    df_consensus[column_name] = df_consensus[column_name].apply(lambda ls: ','.join(map(str, ls)) if isinstance(ls, list) else ls)
df_consensus.to_csv("/Users/jbreynier/results_benchmarking/RNASeq_SPX6730-2_TRA_consensus.csv", index=False)

