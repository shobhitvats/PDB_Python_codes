import os
import sys
import urllib.request
import numpy as np
import Bio
import Bio.PDB
from Bio.PDB import PDBParser, PICIO, PDBIO
import Bio.SeqRecord


# rsync -rlpt -v -z --delete --port=33444 \
# rsync.rcsb.org::ftp_data/structures/divided/pdb/ ./pdb


# def download_read_pdb(pdbcode, datadir, keepfile=True):
#     """
#     Downloads a PDB file from the Internet and saves it in a data directory.
#     Then it reads and returns the structure inside.
#     :param pdbcode: The standard PDB ID e.g. '3ICB'
#     :param datadir: The directory where the downloaded file will be saved
#     :param keepfile: if False, then the downloaded file will be deleted (default: keep the downloaded file)
#     :return: a Bio.PDB Structure object or None if something went wrong
#     """
#     pdbfilenm = download_pdb(pdbcode, datadir)
#     if pdbfilenm is None:
#         return None
#     struct = read_pdb(pdbcode, pdbfilenm)
#     if not keepfile:
#         os.remove(pdbfilenm)
#     return struct

def download_pdb(pdbcode, datadir, downloadurl="http://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
        Note that the unencrypted HTTP protocol is used by default
        to avoid spurious OpenSSL errors...
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        # all sorts of things could have gone wrong...
        print(str(err), file=sys.stderr)
        return None
    
# def read_pdb(pdbcode, pdbfilenm):
#     """
#     Read a PDB structure from a file.
#     :param pdbcode: A PDB ID string
#     :param pdbfilenm: The PDB file
#     :return: a Bio.PDB.Structure object or None if something went wrong
#     """
#     try:
#         pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
#         struct = pdbparser.get_structure(pdbcode, pdbfilenm)
#         return struct
#     except Exception as err:
#         print(str(err), file=sys.stderr)
#         return None 

# def extract_seqrecords(pdbcode, struct):
#     """
#     Extracts the sequence records from a Bio.PDB structure.
#     :param pdbcode: the PDB ID of the structure, needed to add a sequence ID to the result
#     :param struct: a Bio.PDB.Structure object
#     :return: a list of Bio.SeqRecord objects
#     """
#     ppb = Bio.PDB.PPBuilder()
#     seqrecords = []
#     for i, chain in enumerate(struct.get_chains()):
#         # extract and store sequences as list of SeqRecord objects
#         pps = ppb.build_peptides(chain)    # polypeptides
#         seq = pps[0].get_sequence() # just take the first, hope there's no chain break
#         seqid = pdbcode + chain.id
#         seqrec = Bio.SeqRecord.SeqRecord(seq, id=seqid, 
#             description="Sequence #{}, {}".format(i+1, seqid))
#         seqrecords.append(seqrec)
#     return seqrecords

# def get_calphas(struct):
#     """
#     Extracts the C-alpha atoms from a PDB structure.
#     :param struct: A Bio.PDB.Structure object.
#     :return: A list of Bio.PDB.Atom objects representing the C-alpha atoms in `struct`.
#     """
#     calphas = [ atom for atom in struct.get_atoms() if atom.get_fullname() == " CA " ]
#     return calphas

# pdbcode = "1fat.pdb"
# datadir = "./PDB"

download_pdb("1fat","./PDB","http://files.rcsb.org/download/")
    
# parser = PDBParser()
# structure = parser.get_structure("./PDB", "1fou.pdb")

# structure.atom_to_internal_coordinates()
# chain = list(structure.get_chains())[0]
# ic_chain = chain.internal_coord
# dihedrals = ic_chain.dihedra

# newdihedrals = np.random.uniform(-np.pi, np.pi, size = (56, 2))

# modified = PICIO.read_PIC('internal')
# modified.internal_to_atom_coordinates()

# io = set_structure()
# io.set_structure(modified)
# io.save('1fou_modified.pdb', preserve_atom_numbering = True) 

