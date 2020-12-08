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
        row = [sample, receptor]
        df_consensus = pd.read_csv(f"/Users/jbreynier/results_benchmarking/RNASeq_{sample}-2_{receptor}_consensus_q4.csv")
        df_TCR_seq = pd.read_csv(f"/Users/jbreynier/results_benchmarking/{sample}_{receptor}.clonotypes.{receptor}.txt", sep="\t", usecols=["aaSeqCDR3"])
        df_TCR_seq.columns = ["cdr3aa"]
        df_TCR_seq.drop_duplicates(inplace=True)
        df_TCR_seq = df_TCR_seq.apply(verify_canonical, axis=1)
        df_TCR_seq.dropna(inplace=True)
        # print(f"{sample} {receptor}:")
        row.append(len(df_consensus))
        row.append(len(pd.merge(df_consensus, df_TCR_seq, on=['cdr3aa'], how='inner')))


        #TRUST3
        with open("/Users/jbreynier/results_benchmarking/TRUST3/RNASeq_SPX8151-2.Aligned.sortedByCoord.out.bam-TCR-ALL.fa") as fasta_file:
            list_
            for line in lines_fasta:
                if line.startswith(">"):
                    print(line.split("+")[-3])



        #TRUST4
        df_TRUST4_all = pd.read_csv(f"/Users/jbreynier/results_benchmarking/TRUST4/RNASeq_{sample}-2_report.tsv",
                                    sep="\t",
                                    usecols=['CDR3aa',
                                            'V',
                                            'D',
                                            'J',
                                            'C'])
        df_TRUST4 = df_TRUST4_all[df_TRUST4_all[['V','J','C','D']].apply(lambda x: x.str.startswith(receptor)).any(axis=1)]
        df_TRUST4 = df_TRUST4[(df_TRUST4.CDR3aa != "out_of_frame") & (df_TRUST4.CDR3aa != "partial")]
        df_TRUST4.drop(columns=["V", "D", "J", "C"], inplace=True)
        df_TRUST4.columns = ["cdr3aa"]
        df_TRUST4.drop_duplicates(inplace=True)
        row.append(len(df_TRUST4))
        row.append(len(pd.merge(df_TRUST4, df_TCR_seq, on=['cdr3aa'], how='inner')))

        #MiXCR
        df_MiXCR = pd.read_csv(f"/Users/jbreynier/results_benchmarking/MiXCR/RNASeq_{sample}-2_{receptor}_clones.txt",
                                sep="\t",
                                usecols=["aaSeqCDR3"],
                                header=0,
                                dtype="object")
        df_MiXCR.columns = ["cdr3aa"]
        df_MiXCR.drop_duplicates(inplace=True)
        row.append(len(df_MiXCR))
        row.append(len(pd.merge(df_MiXCR, df_TCR_seq, on=['cdr3aa'], how='inner')))
        
        #CATT
        df_CATT = pd.read_csv(f"/Users/jbreynier/results_benchmarking/CATT/RNASeq_{sample}-2.{receptor}.CDR3.CATT.csv",
                                sep=",", usecols=['AAseq'])
        df_CATT.columns = ["cdr3aa"]
        df_CATT.drop_duplicates(inplace=True)
        row.append(len(df_CATT))
        row.append(len(pd.merge(df_CATT, df_TCR_seq, on=['cdr3aa'], how='inner')))

        print(row)
        list_rows.append(row)


df_output = pd.DataFrame(list_rows, columns=["sample", "receptor", 
                                        "consensus_total", "consensus_correct", 
                                        "TRUST4_total", "TRUST4_correct",
                                        "MiXCR_total", "MiXCR_correct",
                                        "CATT_total", "CATT_correct"])
print(df_output)

df_output.to_csv("/Users/jbreynier/results_benchmarking/results_total_q4.csv", index=False)