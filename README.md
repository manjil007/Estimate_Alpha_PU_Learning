# PU-Learning for Ligand Binding Site Prediction

This project applies **Positive and Unlabeled Learning (PU Learning)** techniques to identify ligand-binding residues in proteins. The pipeline integrates sequence, structural, and evolutionary features extracted from BioLiP and related resources to build machine learning models (including PULSNAR) for residue-level binding site prediction.

---

## 🔧 Setup Instructions

### 1. Clone the Repository
```bash
git clone https://github.com/manjil007/Estimate_Alpha_PU_Learning.git
cd Estimate_Alpha_PU_Learning
```

### 2. Install Python Dependencies
Make sure you have Python 3.8 or higher installed. Then run:
```bash
pip install -r requirements.txt
```

### 3. Install PULSNAR
Follow the official installation guide provided at the [PULSNAR GitHub repository](https://github.com/unmtransinfo/PULSNAR). This may include:
```bash
git clone https://github.com/unmtransinfo/PULSNAR.git
cd PULSNAR
pip install .
```

### 4. Install IUPred2A
Download and install IUPred2A from its [official GitHub repository](https://github.com/cbalbin-bio/iupred-parser). Follow the installation steps to ensure the CLI tools like `iupred2a.py` are executable and available in your `$PATH`.

---

## 📂 Data Setup

Before running any scripts, please perform the following:

1. **Download BioLiP Data**  
   Download `BioLiP_nr.txt` from the [BioLiP website](https://zhanggroup.org/BioLiP/) and save it in the following folder:
   ```
   data/raw/BioLiP_nr.txt
   ```

2. **Create Processed Directory**  
   Create a directory to store intermediate processed files:
   ```bash
   mkdir -p data/processed
   ```

---

## 🚀 Execution Order

All scripts and notebooks are numbered (e.g., `1_extract_features.py`, `2_run_pulsnar.py`, etc.) and must be executed in **chronological order of their prefix numbers** to ensure correct data flow and model training.

---

## 📌 Notes

- Make sure all downloaded datasets and installed tools are in place **before** running the pipeline.
- Intermediate files are stored in `data/processed` and will be reused by subsequent scripts.
- The final outputs include balanced datasets, calibrated probabilities, and PU-learning performance summaries.

---

## 📞 Contact

For questions or contributions, please contact Manjil Pradhan at [mpradhan@unm.edu or pradhanmanjil292@gmail.com].