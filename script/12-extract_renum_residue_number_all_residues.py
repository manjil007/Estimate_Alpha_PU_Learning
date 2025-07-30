from pathlib import Path
import csv
import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.PDB import MMCIFParser, PPBuilder
from Bio.Align import substitution_matrices
import csv, sys
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.PDB import MMCIFParser, Polypeptide
from pathlib import Path
from Bio import SeqIO, pairwise2
from Bio.PDB import MMCIFParser, PPBuilder, PDBList
from Bio.Align import substitution_matrices

# ---------------------------------------------------------------------------
# 1.  I/O locations
# ---------------------------------------------------------------------------
BIO_LIP_FILE = Path("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/raw/BioLiP_nr.txt")
CIF_DIR      = Path("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/raw/structures_cif")
FASTA_DIR    = Path("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/raw/iupred_fasta")
OUT_FILE     = Path("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/processed_mapping/map_pdb_biolip_renum_residue_number.csv")
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 2.  Load BioLiP and list unique (pdb_id, chain) pairs
# ---------------------------------------------------------------------------
cols = [
    "pdb_id","receptor_chain","resolution","binding_site_id","ligand_id",
    "ligand_chain","ligand_serial_number","binding_residues_pdb",
    "binding_residues_renum","catalytic_residues_pdb","catalytic_residues_renum",
    "ec_number","go_terms","binding_affinity_literature","binding_affinity_moad",
    "binding_affinity_pdbbind","binding_affinity_bindingdb","uniprot_id",
    "pubmed_id","ligand_residue_seq_number","receptor_sequence"
]

biolip_df = pd.read_csv(BIO_LIP_FILE, sep="\t", header=None, names=cols)
chain_pairs = (
    biolip_df[["pdb_id", "receptor_chain"]]
    .drop_duplicates()
    .itertuples(index=False, name=None)
)

chain_pairs = list(chain_pairs)  

chain_pairs = chain_pairs # Limit to first 1000 pairs for testing

def fetch_structure(pdb_id: str, cif_dir="cifs") -> Path:
    """Download the mmCIF if it is not already on disk and return its path."""
    Path(cif_dir).mkdir(exist_ok=True)
    cif_path = Path(cif_dir) / f"{pdb_id.lower()}.cif"
    if not cif_path.exists():
        PDBList().retrieve_pdb_file(
            pdb_code=pdb_id, file_format="mmCif",
            pdir=cif_dir, overwrite=False
        )
    return cif_path

def chain_sequence(cif_path, chain_id: str) -> str:
    """
    Accepts either a str or Path; returns the one-letter sequence
    for `chain_id` in the mmCIF file.
    """
    path_obj = Path(cif_path)           # ← normalize to Path
    parser   = MMCIFParser(QUIET=True)
    structure = parser.get_structure(path_obj.stem, str(path_obj))

    ppb = PPBuilder()
    for model in structure:
        try:
            chain = model[chain_id]
        except KeyError:
            raise ValueError(f"Chain {chain_id} not found in {path_obj.name}")
        # Concatenate peptides that belong to this chain
        seq = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(chain))
        return seq

    raise RuntimeError("No model found in the structure")


def align_sequences(seq_pdb: str, seq_fasta: str):
    """Global alignment with BLOSUM62; returns alignment objects."""
    blosum62 = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.globalds(seq_pdb, seq_fasta, blosum62, -11, -1)
    # choose the top-scoring alignment
    return alignments[0]


def build_mapping(aln_pdb, aln_fasta, residues):
    """
    Given aligned strings and a list of (res_id, full_residue_obj),
    return a list of rows with mapping between PDB residue id and FASTA index.
    """
    mapping = []
    pdb_idx = 0   # index in original PDB sequence
    fas_idx = 0   # index in FASTA sequence
    for a_pdb, a_fas, res in zip(aln_pdb, aln_fasta, residues):
        if a_pdb != "-":   # residue exists in PDB chain
            res_id = residues[pdb_idx][0]      # (resseq, icode) pair or residue.id
            pdb_idx += 1
        if a_fas != "-":
            fas_idx += 1
        if a_pdb != "-" and a_fas != "-":
            mapping.append((res_id, fas_idx))   # 1-based indices
    return mapping

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# ───────────────────────────────────────────────────────────────────────────
# 0-bis.  Build expected sequence lengths once  (outside the loop)
# ───────────────────────────────────────────────────────────────────────────
seq_len = {
    (row.pdb_id.lower(), row.receptor_chain): len(row.receptor_sequence)
    for row in biolip_df.itertuples(index=False)
}

# Decide whether we are appending or creating the CSV
mode   = "a" if OUT_FILE.exists() else "w"
header = mode == "w"

# ───────────────────────────────────────────────────────────────────────────
# 3.  Write / append mapping rows
# ───────────────────────────────────────────────────────────────────────────
with OUT_FILE.open(mode, newline="") as fh:
    writer = csv.writer(fh)
    if header:
        writer.writerow([
            "pdb_id", "chain_id", "pdb_residue_number",
            "insertion_code", "aa", "renum_residue_number"
        ])

    for pdb_id, chain_id in chain_pairs:
        pdb_id_lc = pdb_id.lower()

        # ───── heavy work starts here ─────
        cif_path   = CIF_DIR / f"{pdb_id_lc}.cif"
        fasta_path = FASTA_DIR / f"{pdb_id}_{chain_id}.fasta"
        if not (cif_path.exists() and fasta_path.exists()):
            continue

        try:
            pdb_seq = chain_sequence(cif_path, chain_id)
        except Exception as err:
            print(f"[warn] {pdb_id} {chain_id}: {err}", file=sys.stderr)
            continue

        fasta_seq = str(next(SeqIO.parse(fasta_path, "fasta")).seq)
        try:
            aln_pdb, aln_fas, *_ = align_sequences(pdb_seq, fasta_seq)
        except (IndexError, ValueError, SystemError) as e:
            print(f"[warn] alignment failed ({e}): {pdb_id} {chain_id}", file=sys.stderr)
            continue

        structure = MMCIFParser(QUIET=True).get_structure(pdb_id_lc, cif_path)
        residues  = [
            (res.id, res) for res in structure[0][chain_id]
            if Polypeptide.is_aa(res, standard=True)
        ]

        id_to_aa = {res.id: d.get(res.get_resname().upper(), "X")
                    for (_, res) in residues}

        for res_id, renum in build_mapping(aln_pdb, aln_fas, residues):
            het, auth_num, icode = res_id
            icode_str  = icode.strip()
            pdb_number = f"{auth_num}{icode_str}"
            aa = id_to_aa.get(res_id, "X")
            writer.writerow([
                pdb_id_lc, chain_id, pdb_number,
                icode_str, aa, renum
            ])
