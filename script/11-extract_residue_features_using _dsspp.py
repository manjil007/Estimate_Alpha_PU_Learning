# extract_features_part1.ipynb
import pandas as pd
from pathlib import Path
from Bio.PDB import MMCIFParser, DSSP
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# ----------------------- settings -----------------------
split_id_file  = "/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/raw/pdb_ids_part3.csv"
structure_dir  = Path("/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/raw/structures_cif")
dssp_exe       = "mkdssp"
output_path    = "/home/mpradhan007/Academic/Research_Projects/Intern_Research/data/processed/residue_features_part3.csv"
# --------------------------------------------------------

warnings.simplefilter("ignore", PDBConstructionWarning)

pdb_ids = pd.read_csv(split_id_file)["pdb_id"].tolist()
parser  = MMCIFParser(QUIET=True)
rows    = []

for pdb_id in pdb_ids:
    cif_file = structure_dir / f"{pdb_id}.cif"
    if not cif_file.exists():
        print(f"[!] Missing structure: {pdb_id}")
        continue

    try:
        structure = parser.get_structure(pdb_id, cif_file)
        dssp      = DSSP(structure[0], cif_file, dssp=dssp_exe)
    except Exception as e:
        print(f"[!] Skipping {pdb_id}: DSSP failed – {e}")
        continue

    for (chain_id, res_id), data in dssp.property_dict.items():
        hetflag, resnum, icode = res_id
        residue_id = f"{resnum}{icode.strip()}" if icode.strip() else str(resnum)
        aa         = data[1]

        rows.append({
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "residue_id": residue_id,
            "residue_name": aa,
            "secondary_structure": data[2],
            "absolute_sasa": data[3],
            "relative_asa": data[4],
            "phi": data[5],
            "psi": data[6],
            "hbond_NH_O1_energy": data[7],
            "hbond_NH_O2_energy": data[8],
            "hbond_O_NH1_energy": data[9],
            "hbond_O_NH2_energy": data[10],
            "is_aromatic": int(aa in {"F","Y","W","H"}),
            "is_polar":     int(aa in {"S","T","N","Q","Y","C"}),
            "is_charged":   int(aa in {"R","K","D","E","H"}),
            "is_hydrophobic": int(aa in {"A","V","I","L","M","F","W","Y"}),
        })

df = pd.DataFrame(rows)
df.to_csv(output_path, index=False)
print(f"[✓] Saved {len(df)} residue rows → {output_path}")
