import subprocess
import pandas as pd
import os
import glob

def run_iupred2a(fasta_path, mode='long', use_anchor=False):
    if use_anchor and mode != 'long':
        raise ValueError("Anchor mode must use 'long' as the prediction type.")

    if use_anchor:
        cmd = [
            'python3',
            '/home/mpradhan007/Academic/Research_Projects/Intern_Research/iupred2a/iupred2a.py',
            '-a',
            fasta_path,
            'long'
        ]
    else:
        cmd = [
            'python3',
            '/home/mpradhan007/Academic/Research_Projects/Intern_Research/iupred2a/iupred2a.py',
            fasta_path,
            mode
        ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout

    if not output or "Usage:" in output or "not found" in output:
        raise RuntimeError(f"IUPred2A failed for {os.path.basename(fasta_path)} in mode '{mode}':\n{result.stderr or output}")

    lines = output.strip().split('\n')
    data = []

    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        fields = line.strip().split()
        if len(fields) < 3:
            continue
        pos, aa, score = fields[0], fields[1], fields[2]
        label = 'anchor' if use_anchor else mode
        data.append({
            'position': int(pos),
            'amino_acid': aa,
            f'iupred2a_{label}_score': float(score)
        })

    return pd.DataFrame(data)


def extract_ids_from_filename(filename):
    basename = os.path.basename(filename)
    pdb_chain = basename.replace('.fasta', '')
    if '_' not in pdb_chain:
        raise ValueError(f"Filename {basename} does not follow the expected 'pdb_chain.fasta' format.")
    pdb_id, chain_id = pdb_chain.split('_')
    if chain_id == 'SX':
        print("incorrect chain_id, skipping")
    return pdb_id.lower(), chain_id.upper()


if __name__ == "__main__":
    input_dir = "/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/raw/iupred_fasta"
    output_csv = "/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/iupred2a_scores_all.csv"

    all_files = glob.glob(os.path.join(input_dir, "*.fasta"))
    all_dfs = []

    for fasta_file in all_files:
        pdb_id, chain_id = extract_ids_from_filename(fasta_file)

        try:
            df_long = run_iupred2a(fasta_file, mode='long')
            df_short = run_iupred2a(fasta_file, mode='short')
            df_anchor = run_iupred2a(fasta_file, mode='long', use_anchor=True)

            df = df_long.merge(df_short, on=['position', 'amino_acid']) \
                        .merge(df_anchor, on=['position', 'amino_acid'])

            df['pdb_id'] = pdb_id
            df['chain_id'] = chain_id

            all_dfs.append(df)
        except Exception as e:
            print(f"Failed on {fasta_file}: {e}")

    if all_dfs:
        final_df = pd.concat(all_dfs, ignore_index=True)
        final_df.to_csv(output_csv, index=False)
        print(f"Saved combined IUPred2A scores to {output_csv}")
    else:
        print("No data processed.")
