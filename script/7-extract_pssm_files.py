import subprocess
import os
import glob

def run_psiblast(fasta_path, db_path, output_dir, num_iterations=3):
    pdb_id = os.path.splitext(os.path.basename(fasta_path))[0]
    pssm_ascii = os.path.join(output_dir, f"{pdb_id}.pssm")

    if os.path.exists(pssm_ascii):
        print(f"Skipping {pdb_id}, PSSM already exists.")
        return

    cmd = [
        "psiblast",
        "-query", fasta_path,
        "-db", db_path,
        "-num_iterations", str(num_iterations),
        "-out_ascii_pssm", pssm_ascii,
        "-evalue", "0.001"
    ]
    print(f"Running PSI-BLAST for {pdb_id}...")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"PSI-BLAST failed for {pdb_id}: {e}")

if __name__ == "__main__":
    fasta_dir = "/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/raw/iupred_fasta"
    blast_db = "/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/raw/blast_db/uniprot_sprot_db"
    output_dir = "/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/pssm"
    os.makedirs(output_dir, exist_ok=True)

    fasta_files = glob.glob(os.path.join(fasta_dir, "*.fasta"))[38000:40662]
    # fasta_files = glob.glob(os.path.join(fasta_dir, "*.fasta"))[13000:17000]
    for fasta_path in fasta_files:

        run_psiblast(fasta_path, blast_db, output_dir)


