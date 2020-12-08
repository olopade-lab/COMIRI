import pandas as pd
import re

def verify_canonical(row):
    match = re.search(r"C[^_]*[FW]", row['cdr3aa'])
    if not match:
        row['cdr3aa'] = None
        # row['cdr3dna'] = None
    else:
        row['cdr3aa'] = match.group()
        # row['cdr3dna'] = row['cdr3dna'][match.span()[0]*3:match.span()[1]*3]
    return row

# Do we care about amino acid or DNA?

list_rows = []

for sample in ["SPX6730", "SPX8151"]:
    for receptor in ["TRB"]:
        df_TRUST4_all = pd.read_csv(f"/Users/jbreynier/results_benchmarking/TRUST4/RNASeq_{sample}-2_report.tsv",
                                    sep="\t",
                                    usecols=['CDR3aa',
                                            'V',
                                            'D',
                                            'J',
                                            'C'])
        df_TRUST4 = df_TRUST4_all[df_TRUST4_all[['V','J','C','D']].apply(lambda x: x.str.startswith(receptor)).any(axis=1)]
        df_TRUST4 = df_TRUST4[(df_TRUST4.CDR3aa != "out_of_frame") & (df_TRUST4.CDR3aa != "partial")]
        print(len(df_TRUST4))