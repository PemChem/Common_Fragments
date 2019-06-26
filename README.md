The python script developed is called Common_Fragments.py

There are two dependancies:
	1. A file with a list of smiles strings called search_result.smi 
		a. must placed in the same folder
	2. An original smiles string used to make all comparisons
		a. There is a variable on line 22 of Common_Fragments.py called orig_smiles
		b. This can be editted

The script can be run in a python/terminal/unix window with the command:

	python Common_Fragments.py

Alternatively, the script can be pasted into an interactive python interface like spyder and run from there.

Goals of the script:
	1. Identify the common fragments between a reference molecule and a list of molecules
	2. Plot the total number of common fragments between each molecule and the reference
	3. Compute Tanimoto and Morgan similarities between each molecule and the reference
	4. Plot those similarites
	5. Use the MurckoScaffold substructure routine to identify common substructures among the list of molecules
		a. This can provide an additional element to describe diversity in the list of molecules
		b. Assign each molecule in the list to a common substructure

Ultimate Goal:
	To make an educated guess of which molecules to select for synthesis based on similarity and properties of the original molecule.

Environmental characteristics are described in the file called "requirements.txt"
