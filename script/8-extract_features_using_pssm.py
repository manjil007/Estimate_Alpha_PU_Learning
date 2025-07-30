import pandas as pd
import os
import glob

def parse_pssm_ascii(pssm_path):
    with open(pssm_path, "r") as f:
        lines = f.readlines()

    start = None
    for i, line in enumerate(lines):
        if line.startswith("Last position-specific scoring matrix"):
            start = i + 2
            break

    if start is None:
        raise ValueError(f"Could not find scoring matrix in PSSM file: {pssm_path}")

    aa_order = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    pssm_data = []
    for line in lines[start:]:
        if line.strip() == "":
            break
        fields = line.strip().split()
        if len(fields) < 22:
            continue

        pos = int(fields[0])
        aa = fields[1]
        scores = list(map(int, fields[2:22]))
        pssm_data.append({
            "position": pos,
            "amino_acid": aa,
            **{f"PSSM_{aa_col}": score for aa_col, score in zip(aa_order, scores)}
        })

    return pd.DataFrame(pssm_data)

if __name__ == "__main__":
    pssm_dir = "/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/non_biolip_pssm_1"
    output_csv = "/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/non_biolip_pssm_features.csv"

    all_pssm_dfs = []
    pssm_files = glob.glob(os.path.join(pssm_dir, "*.pssm"))


    for pssm_path in pssm_files:
        try:
            
            pdb_chain = pssm_path
            if "_" not in pdb_chain:
                print(f"Skipping malformed filename: {pdb_chain}")
                continue
            _, _, _, _, _, pdb_id, chain = pdb_chain.split("_")
            _, pdb_id = pdb_id.split("/")
            chain = chain.split(".")[0]  # Remove file extension if present
            pdb_id = pdb_id.lower()
            chain = chain.upper()

            df_pssm = parse_pssm_ascii(pssm_path)
            df_pssm["pdb_id"] = pdb_id
            df_pssm["chain"] = chain

            all_pssm_dfs.append(df_pssm)

        except Exception as e:
            print(f"Failed to parse {pssm_path}: {e}")

    if all_pssm_dfs:
        final_df = pd.concat(all_pssm_dfs, ignore_index=True)
        final_df.to_csv(output_csv, index=False)
        print(f"\n✅ All PSSM features saved to: {output_csv}")
    else:
        print("⚠️ No PSSM data parsed.")
