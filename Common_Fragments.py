from __future__ import print_function
from rdkit import Chem
import os
import matplotlib.pyplot as plt
import numpy as np
from rdkit import RDConfig
from rdkit.Chem import FragmentCatalog
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Draw

### Creating a directory for output
work_dir = os.getcwd()
out_dir = work_dir + "/output/"
try:
    os.mkdir(work_dir + "/output")
except:
    FileExistsError
### Original Smiles string used to search PubChem
orig_smiles = 'CC1=CC(=O)OC2=C1C=C3C=C(OC3=C2C)C'

### Making a list of the smiles
### strings and the IDs
### from the downloaded file

smiles_list = []
id_list = []
for line in open("./search_result.smi"):
	smiles_list.append(line.split()[0])
	id_list.append(line.split()[1].rstrip())

### Creating a list of functional groups
fName=os.path.join(RDConfig.RDDataDir,'FunctionalGroups.txt')
fparams = FragmentCatalog.FragCatParams(1,6,fName)
print("We are considering " + str(fparams.GetNumFuncGroups()) + " functional groups")

### Identifying fragments in original smiles string
fcat=FragmentCatalog.FragCatalog(fparams)
fcgen=FragmentCatalog.FragCatGenerator()
m = Chem.MolFromSmiles(orig_smiles)
total_orig_frags = fcgen.AddFragsFromMol(m,fcat)
print("There are " + str(total_orig_frags) + " fragments in the original molecule")

### Creating the list of fragments in the original molecule
unique_fragments = []
for i in range(0, total_orig_frags):
    if fcat.GetEntryDescription(i) not in unique_fragments:
        unique_fragments.append(fcat.GetEntryDescription(i))

### Creating a matrix of 0s that has dimensions:
### Rows = 100; number of molecules
### Columns = total_orig_frags; number of fragments found in original smiles (236)
frag_matrix = []
for i in range(0, len(smiles_list)):
    to_add = []
    for j in range(0, total_orig_frags):
        to_add.append(0)
    frag_matrix.append(to_add)

### Identifying common fragment between downloaded and original molecule
for i in range(0, len(smiles_list)):
    m = Chem.MolFromSmiles(smiles_list[i])
    fcat=FragmentCatalog.FragCatalog(fparams)
    fcgen=FragmentCatalog.FragCatGenerator()
    total_frags = fcgen.AddFragsFromMol(m,fcat)
    for j in range(0, total_frags):
        if fcat.GetEntryDescription(j) in unique_fragments:
            frag_index = unique_fragments.index(fcat.GetEntryDescription(j))
            frag_matrix[i][frag_index] = 1


### Creating a list for the number of common fragments
common_frags = []
for i in range(0, len(smiles_list)):
    common_frags.append(sum(frag_matrix[i]))

### Preparing a bar plot to see the common fragments	
plt.bar(range(1, len(smiles_list) + 1), common_frags)
plt.xticks(np.arange(1, len(smiles_list) + 1, step=5))
plt.yticks(np.arange(0, len(unique_fragments), step=10))
plt.ylabel('Number of Common Fragments')
plt.xlabel('Downloaded Smiles ID')
plt.savefig(out_dir + 'Common_Fragments.png')
plt.show()

### Creating a list of Tanimoto Fingerprint Similarities
### for each of the compounds
tanimoto_sim = []
m1 = Chem.MolFromSmiles(orig_smiles)
fp1 = FingerprintMols.FingerprintMol(m1)
for i in smiles_list:
    m2 = Chem.MolFromSmiles(i)
    fp2 = FingerprintMols.FingerprintMol(m2)
    tanimoto_sim.append(round(DataStructs.FingerprintSimilarity(fp1,fp2), 2))

### Preparing a bar plot to see the Tanimoto Similarites
plt.bar(range(1, len(smiles_list) + 1), tanimoto_sim)
plt.xticks(np.arange(1, len(smiles_list) + 1, step=5))
plt.yticks(np.arange(0, 1.2, step=0.2))
plt.ylabel('Tanimoto_Similarity')
plt.xlabel('Downloaded_Smiles_ID')
plt.savefig(out_dir + 'Tanimoto_Similarites.png')
plt.show()

### Creating a list of Morgan Fingerprint Similarities
### for each of the compounds
morgan_sim = []
m1 = Chem.MolFromSmiles(orig_smiles)
fp1 = AllChem.GetMorganFingerprint(m1,2)
for i in smiles_list:
    m2 = Chem.MolFromSmiles(i)
    fp2 = AllChem.GetMorganFingerprint(m2,2)
    morgan_sim.append(round(DataStructs.DiceSimilarity(fp1,fp2), 2))

### Preparing a bar plot to see the Morgan Similarites
plt.bar(range(1, len(smiles_list) + 1), morgan_sim)
plt.xticks(np.arange(1, len(smiles_list) + 1, step=5))
plt.yticks(np.arange(0, 1.2, step=0.2))
plt.ylabel('Morgan_Similarity')
plt.xlabel('Downloaded_Smiles_ID')
plt.savefig(out_dir + 'Morgan_Similarites.png')
plt.show()

### Plotting Morgan vs Tanimoto Fingerprints
plt.scatter(morgan_sim, tanimoto_sim)
plt.xticks(np.arange(0, 1.2, step=0.2))
plt.yticks(np.arange(0, 1.2, step=0.2))
plt.ylabel('Tanimoto_Similarity')
plt.xlabel('Morgan_Similarity')
plt.savefig(out_dir + 'Tanimoto_vs_Morgan.png')
plt.show()

### Using MurckoScaffold to identify common substructures
### cores variable will be used to develop unique substructures/clusters

cores = []
cores_list = []
core_id_list = []
for i in range(0, len(id_list)):
	m1 = Chem.MolFromSmiles(smiles_list[i])
	core = MurckoScaffold.GetScaffoldForMol(m1)
	cores_list.append(str(Chem.MolToSmiles(core)))
	if str(Chem.MolToSmiles(core)) not in cores:
		cores.append(str(Chem.MolToSmiles(core)))
	core_id_list.append(cores.index(str(Chem.MolToSmiles(core))))

### Renumber Core-ID values
renumber_cores = []
for i in core_id_list:
    renumber_cores.append(i+1)
### Preparing a CSV file to view smiles with
### fingerprint and substructure cluster information
f1 = open(out_dir + "Tab_Delimited_File.txt", 'a')
f1.write("Smiles_CMPD" + "\t" + "Common_Fragments" + "\t" + "ID" + "\t" + "ID_Substructure" + "\t" + "Tanimoto_Similarity" + "\t" + "Morgan_Similarity" + "\n")
f1.close()
for row in zip(smiles_list, common_frags, id_list, renumber_cores, tanimoto_sim, morgan_sim):
    row_list = (str(list(row)))
    row_list = row_list.replace("[", "")
    row_list = row_list.replace("]", "")
    row_list = row_list.replace("'", "")
    row_list = row_list.replace(", ", "\t")
    f1 = open(out_dir + "Tab_Delimited_File.txt", 'a')
    f1.write(row_list + "\n")
    f1.close()

### Preparing a bar plot to see the unique cores
core_counts = []
core_id_for_barplot = []
for i in range(0, len(cores)):
	core_counts.append(core_id_list.count(i))
	core_id_for_barplot.append(i + 1)	
plt.bar(core_id_for_barplot, core_counts)
plt.xticks(np.arange(1, len(cores) + 1, step=2))
plt.yticks(np.arange(0, 20, step=2))
plt.ylabel('Frequency per Cluster')
plt.xlabel('Substructure_Cluster_ID')
plt.savefig(out_dir + 'Bar-Plot-Substructures.png')
plt.show()

### Showing substructure for clusters with more than 2 members
cores_for_image = []
core_id_for_image = []
for i in core_id_for_barplot:
    j = i - 1
    if core_counts[j] > 2:
        cores_for_image.append(cores[j])
        core_id_for_image.append('Core-' + str(i))
core_set = []
for i in range(0, len(cores_for_image)):
    m = Chem.MolFromSmiles(cores_for_image[i])
    core_set.append(m)
    
### Adding in original smiles string
core_set.append(Chem.MolFromSmiles(orig_smiles))
core_id_for_image.append('Original')

### Creating grid to view Substructures
img = Draw.MolsToGridImage(core_set,molsPerRow=4,subImgSize=(200,200),legends=[x for x in core_id_for_image])
img.save(out_dir + 'Popular-SubStructures.png')   
img.show()