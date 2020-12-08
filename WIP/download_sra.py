import subprocess
import pandas as pd
import sys
import os

path_runinfo = sys.argv[1]
output_path = sys.argv[2]
df_runinfo = pd.read_csv(path_runinfo)
for index, row in df_runinfo.iterrows():
    if row["SampleName"] in ["GSM4455935", "GSM4455933"]:
        specific_output = os.path.join(output_path, row["SampleName"])
        if not os.path.isdir(specific_output):
            os.mkdir(specific_output)
        wget_process = subprocess.Popen(["wget", "-P", specific_output, row["download_path"]])
        wget_process.wait()

