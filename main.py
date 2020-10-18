import urllib.request
import os
import csv
from datetime import date
from pathlib import Path
import difflib
from Bio.PDB import PDBList
from Bio.PDB.mmcifio import MMCIFIO
import glob
from Bio.PDB import *

def save_file(text, filename):
    file = open(filename, 'w')
    file.write(text)
    file.close()

def download_non_redundant_set():
    resolution = '3.0A'
    base_url = 'http://rna.bgsu.edu/rna3dhub/nrlist/download'
    release = 'current'
    url = '/'.join([base_url, release, resolution])
    response = urllib.request.urlopen(url)
    try:
        info = response.info()['Content-Disposition']
        filename = info.split('filename=',1)[1]
    except:
        today = date.today()
        timestamp = today.strftime("%d-%m-%Y")
        filename = 'output_d' + timestamp + '.csv'
    data = response.read()
    text = data.decode('utf-8')

    path_to_location = Path('./RNA_SETS')
    save_file(text, path_to_location / filename)


def pick_file(n):
    base_path = Path('./RNA_SETS')
    list_of_files = list(base_path.glob('*.csv'))
    list_of_files_sorted = sorted(list_of_files, key=os.path.getmtime)
    picked_file = list_of_files_sorted[-n]
    return picked_file

def parse_output_file(filename):
    list_of_structures = []
    path_to_file = Path(filename)
    with open(path_to_file, newline='') as csvfile:
        csv_file_content = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in csv_file_content:
            list_of_structures.append('+' + row[1][1:7])
    return list_of_structures

def find_difference(list_of_structures_new, list_of_structures_old):
    f = open('files_to_update.txt', 'w')
    for line in difflib.unified_diff(list_of_structures_old, list_of_structures_new, fromfile='file1', tofile='file2', lineterm='', n=0):
        for prefix in ('---', '+++', '@@'):
            if line.startswith(prefix):
                break
        else:
            f.write(line +'\n')
    f.close()

def parse_output_file1():
    output_data_folder = Path("./")
    if (is_non_zero_file(output_data_folder / 'files_to_update.txt')):
        file = open('files_to_update.txt', mode = 'r')
    elif (is_non_zero_file(output_data_folder / 'init_set.txt')):
        file = open('init_set.txt', mode = 'r')
    lines = file.readlines()
    file.close()
    for line in lines:
        if line[0] =='+':
             download_PDB_structures(line.replace("+","")[:4])



def download_PDB_structures(pdb_ID):
    pdb_data_folder = Path("PDB_files_raw/")
    if not os.path.exists(pdb_data_folder):
        os.makedirs(pdb_data_folder)
    pdb_ID = pdb_ID.lower()
    if (not is_non_zero_file(pdb_data_folder / (pdb_ID + ".pdb"))) and (not is_non_zero_file(pdb_data_folder / (pdb_ID + ".cif"))):
        PDBList().retrieve_pdb_file(pdb_ID,pdir=pdb_data_folder, file_format='mmCif')

def is_non_zero_file(file_path):
    return os.path.isfile(file_path) and os.path.getsize(file_path) > 0

def merge_dictionary(dict1, dict2):
    res = {**dict1, **dict2}
    return res

def standardize_models():
    non_standard_residues_A = ['A23', 'A2L', 'A2M', 'A39', 'A3P', 'A44', 'A5O', 'A6A', 'A7E', 'A9Z',
                               'ADI', 'ADP', 'AET', 'AMD', 'AMO', 'AP7', 'AVC', 'MA6', 'MAD', 'MGQ',
                               'MIA', 'MTU', 'M7A', '26A', '2MA', '6IA', '6MA', '6MC', '6MP', '6MT',
                               '6MZ', '6NW', 'F3N', 'N79', 'RIA', 'V3L', 'ZAD', '31H', '31M', '7AT',
                               'O2Z', 'SRA', '00A', '45A', '8AN', 'LCA', 'P5P', 'PPU', 'PR5', 'PU',
                               'T6A', 'TBN', 'TXD', 'TXP', '12A', '1MA', '5FA']

    non_standard_residues_G = ['A6G', 'E6G', 'E7G', 'EQ4', 'IG', 'IMP', 'M2G', 'MGT', 'MGV', 'MHG',
                               'QUO', 'YG', 'YYG', '23G', '2EG', '2MG', '2SG', 'B8K', 'B8W', 'B9B',
                               'BGH', 'N6G', 'RFJ', 'ZGU', '7MG', 'CG1', 'G1G', 'G25', 'G2L', 'G46',
                               'G48', 'G7M', 'GAO', 'GDO', 'GDP', 'GH3', 'GNG', 'GOM', 'GRB', 'GTP',
                               'KAG', 'KAK', 'O2G', 'OMG', '8AA', '8OS', 'LG', 'PGP', 'P7G', 'TPG',
                               'TG', 'XTS', '102', '18M', '1MG']

    non_standard_residues_C = ['A5M', 'A6C', 'E3C', 'IC', 'M4C', 'M5M', '6OO', 'B8Q', 'B8T', 'B9H',
                               'JMH', 'N5M', 'RPC', 'RSP', 'RSQ', 'ZBC', 'ZCY', '73W', 'C25', 'C2L',
                               'C31', 'C43', 'C5L', 'CBV', 'CCC', 'CH', 'CSF', 'OMC', 'S4C', '4OC',
                               'LC', 'LHH', 'LV2', 'PMT', 'TC', '10C', '1SC', '5HM', '5IC', '5MC']

    non_standard_residues_U =  ['A6U', 'IU', 'I4U', 'MEP', 'MNU', 'U25', 'U2L', 'U2P', 'U31', 'U34',
                                'U36', 'U37', 'U8U', 'UAR', 'UBB', 'UBD', 'UD5', 'UPV', 'UR3', 'URD',
                                'US5', 'UZR', 'UMO', 'U23', '2AU', '2MU', '2OM', 'B8H', 'FHU', 'FNU',
                                'F2T', 'RUS', 'ZBU', '3AU', '3ME', '3MU', '3TD', '70U', '75B', 'CNU',
                                'OMU', 'ONE', 'S4U', 'SSU', 'SUR', '4SU', '85Y', 'DHU', 'H2U', 'LHU',
                                'PSU', 'PYO', 'P4U', 'T31', '125', '126', '127', '1RN', '5BU', '5FU',
                                '5MU', '9QV']

    list_of_nucleotides = ['A', 'C', 'G', 'U']
    list_of_nucleotides = [*list_of_nucleotides, *non_standard_residues_A, *non_standard_residues_G,
                           *non_standard_residues_C, *non_standard_residues_U]

    dict_non_standard_residues_A = {non_standard_residues_A[i]: 'A' for i in range(0, len(non_standard_residues_A))}
    dict_non_standard_residues_G = {non_standard_residues_G[i]: 'G' for i in range(0, len(non_standard_residues_G))}
    dict_non_standard_residues_C = {non_standard_residues_C[i]: 'C' for i in range(0, len(non_standard_residues_C))}
    dict_non_standard_residues_U = {non_standard_residues_U[i]: 'U' for i in range(0, len(non_standard_residues_U))}
    dict_list = [dict_non_standard_residues_A, dict_non_standard_residues_G, dict_non_standard_residues_C, dict_non_standard_residues_U]
    dict_non_standard_residues = {}
    for dictionary in dict_list:
        dict_non_standard_residues = merge_dictionary(dict_non_standard_residues, dictionary)
    acceptable_atoms = ['C2', 'C4', 'C6', 'C8', 'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9', 'O2', 'O4', 'O6',
                        'C1\'', 'C2\'', 'C3\'', 'C4\'', 'C5\'', 'O2\'', 'O3\'', 'O4\'', 'O5\'', 'OP1', 'OP2', 'P' ]

    for files in glob.glob('PDB_files_raw/*.cif'):
        residue_to_remove= []
        chain_to_remove = []
        atom_to_remove = []
        pdb_data_folder = Path("PDB_files_raw/")
        structure_name = files[len(files)-8:-4]
        if (not is_non_zero_file(pdb_data_folder / (structure_name + ".pdb"))) and (not is_non_zero_file(pdb_data_folder / (structure_name + ".cif"))):
            try:
                structure = get_structure(pdb_data_folder, structure_name)
            except:
                continue
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.resname not in list_of_nucleotides:

                            residue_to_remove.append((chain.id, residue.id))
                        if residue.resname in dict_non_standard_residues:

                            residue.resname = dict_non_standard_residues[residue.resname]
                            residue_id = list(residue.id)
                            residue_id[0] = ' '
                            residue.id = tuple(residue_id)
                            for atom in residue:
                                if atom.id not in acceptable_atoms:

                                    atom_to_remove.append((chain.id, residue.id, atom.id))
                                    atom_id = list(atom.full_id[3])
                                    atom_id[0] = ' '
                                    atom_tuple = tuple(atom_id)
                                    atom.full_id = [atom.full_id[0], atom.full_id[1], atom.full_id[2], atom_tuple, atom.full_id[4]]


                    if len(chain) == 0:
                        chain_to_remove.append(chain.id)
            for atom in atom_to_remove:
                try:
                    model[atom[0]][atom[1]].detach_child(atom[2])
                except:
                    continue

            for residue in residue_to_remove:
                try:
                    model[residue[0]].detach_child(residue[1])
                except:
                    continue

            for chain in chain_to_remove:
                model.detach_child(chain)
            save_structure(structure, structure_name)



def get_structure(pdb_data_folder, structure_PDB_ID):
    """
    Function to retrieve information about structure of specified molecule
    :param pdb_data_folder: path to folder contaiting all 3D structures
    :param structure_PDB_ID: structures PDB ID
    :return: the structure object
    """

    parser_cif = MMCIFParser(QUIET=True)
    structure = parser_cif.get_structure(structure_PDB_ID, pdb_data_folder / (structure_PDB_ID + ".cif"))

    return structure




def save_structure(structure, structure_PDB_ID):
    """
    Function creates the file with junction 3D model (*.cif format) and returns its name
    :param structure: the whole PDB structure
    :param structure_PDB_ID: PDB ID of whole structure
    :return: name of the created file
    """

    pdbout_data_folder = Path("PDB_files/")
    if not os.path.exists(pdbout_data_folder):
        os.makedirs(pdbout_data_folder)
    name_of_file = structure_PDB_ID + ".cif"

    io=MMCIFIO()
    io.set_structure(structure)
    io.save('./PDB_files/' + name_of_file)

    return name_of_file


if __name__ == "__main__":
    if not os.path.exists(Path('./RNA_SETS')):
        os.makedirs(Path('RNA_SETS'))
    download_non_redundant_set()
    base_path = Path('./RNA_SETS')
    if len(list(base_path.glob('*.csv'))) > 1:
        list_of_structures_new = parse_output_file(pick_file(1))
        list_of_structures_old = parse_output_file(pick_file(2))
        find_difference(list_of_structures_new, list_of_structures_old)
    else:
        save_file('\n'.join(parse_output_file(pick_file(1))), 'init_set.txt')
    parse_output_file1()
    standardize_models()


