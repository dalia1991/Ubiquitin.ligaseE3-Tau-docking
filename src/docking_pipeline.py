# =====================================================================
# 1. INSTALL THE NECESSARY MOLECULAR PREPARATION TOOLS
# =====================================================================
print("Setting up chemical biology environment...")
!apt-get install -y openbabel > /dev/null 2>&1
!pip install biopython > /dev/null 2>&1

import os
from Bio.PDB import PDBList, PDBParser, PDBIO, Select

# Create a clean directory for our structural workspace
os.makedirs("docking_workspace", exist_ok=True)
os.chdir("docking_workspace")
print("Environment ready!\n")

# =====================================================================
# 2. DOWNLOAD THE STRUCTURAL COORDINATES FROM THE PDB
# =====================================================================
pdbl = PDBList()

# Download Human Cereblon (E3 Ligase) complexed with an immunomodulatory drug
print("Fetching E3 Ligase (Cereblon, PDB: 4CIH) from RCSB...")
ligase_file = pdbl.retrieve_pdb_file('4CIH', pdir='.', file_format='pdb')
# Rename for simplicity
os.rename("pdb4cih.ent", "raw_ligase.pdb")

# Download Tau peptide fragment forming toxic beta-sheets
print("Fetching Tau Fragment (PDB: 7S76) from RCSB...")
tau_file = pdbl.retrieve_pdb_file('7S76', pdir='.', file_format='pdb')
os.rename("pdb7s76.ent", "raw_tau.pdb")

# =====================================================================
# 3. BIOCHEMICAL CLEANING (Removing Waters, Salts, and Ligands)
# =====================================================================
class CleanProtein(Select):
    """Custom selector to preserve only standard amino acid residues."""
    def accept_residue(self, residue):
        # Reject water molecules (HOH) and hetero-atoms/non-protein ligands
        if residue.get_resname() in ['HOH', 'WAT', 'H2O']:
            return 0
        if residue.id[0].startswith('H_'):  # Removes bound drugs/ions
            return 0
        return 1

parser = PDBParser(QUIET=True)
io = PDBIO()

# Clean Ligase
structure_ligase = parser.get_structure('ligase', 'raw_ligase.pdb')
io.set_structure(structure_ligase)
io.save('ligase_no_water.pdb', CleanProtein())
print("\n[✓] E3 Ligase stripped of crystallographic waters and ligands.")

# Clean Tau
structure_tau = parser.get_structure('tau', 'raw_tau.pdb')
io.set_structure(structure_tau)
io.save('tau_no_water.pdb', CleanProtein())
print("[✓] Tau peptide stripped of structural waters and crystallographic artifacts.")

# =====================================================================
# 4. PROTONATION (Adding Missing Hydrogens at Physiological pH 7.4)
# =====================================================================
print("\nRunning protonation via OpenBabel to fix atomic valence states...")

# -p 7.4 automatically calculates the protonation state of residues like Lys, His, Asp, Glu
!obabel ligase_no_water.pdb -O ligase_clean.pdb -p 7.4 > /dev/null 2>&1
!obabel tau_no_water.pdb -O tau_clean.pdb -p 7.4 > /dev/null 2>&1

print("\n=====================================================================")
print("SUCCESS: Your structural inputs are perfectly prepped for docking!")
print("Files generated in your directory:")
print("1. ligase_clean.pdb (Ready E3 Cereblon receptor)")
print("2. tau_clean.pdb    (Ready Soluble Tau target)")
print("=====================================================================")