from pdbfixer import PDBFixer
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from openmm.app import PDBFile
def fix_protein(filename='',addHs_pH=7.4,output='',try_renumberResidues=False):

    fix = PDBFixer(filename=filename)
    fix.findMissingResidues()
    fix.findNonstandardResidues()
    fix.replaceNonstandardResidues()
    fix.removeHeterogens(True)
    fix.findMissingAtoms()
    fix.addMissingAtoms()
    fix.addMissingHydrogens(addHs_pH)
    PDBFile.writeFile(fix.topology, fix.positions, open(output, 'w'))

    if try_renumberResidues == True:
        try:
            original=mda.Universe(filename)
            from_fix=mda.Universe(output)

            resNum=[res.resid for res in original.residues]
            for idx,res in enumerate(from_fix.residues):
                res.resid = resNum[idx]

            save=PDB.PDBWriter(filename=output)
            save.write(from_fix)
            save.close()
        except Exception:
            print('Not possible to renumber residues, check excepton for extra details')
