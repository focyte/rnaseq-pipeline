import pandas as pd
import sys
import os

def merge_files(file_paths, output_path):
    dfs = []
    column_names = []


    for file_path in file_paths:
        file_name = os.path.basename(file_path)
        df = pd.read_csv(file_path, sep='\t', skiprows=1, usecols=[0, 6])
        dfs.append(df)
        column_names.append(file_name)

    merged_df = pd.concat(dfs, axis=1, keys=column_names)
    merged_df.columns = merged_df.columns.droplevel(1)  # Drop the second level of the multi-level columns
    merged_df.to_csv(output_path, index=False)
    print(merged_df)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py file1_path file2_path ... output_path")
        sys.exit(1)

    file_paths = sys.argv[1:-1]
    output_path = sys.argv[-1]

    merge_files(file_paths, output_path)
