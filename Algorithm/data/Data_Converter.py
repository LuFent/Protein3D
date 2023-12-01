import pandas as pd
from pandas import DataFrame as df
import os
from Bio import PDB as pdb

class Data_Converter():

    def __init__(self, hyperparameters):
        self._directory = hyperparameters.proj_dir + hyperparameters.data_dir
        self._dataName = hyperparameters.data_name
        self._dataExtention = hyperparameters.data_ext

        self._structure = None

        self._path_data = self._directory  + "/" + self._dataName + self._dataExtention
        self._pdb_id = self._directory  + "/" + self._dataName


    def get_structure(self):
        if(self._dataExtention == ".pdb"):
            parcer = pdb.PDBParser()
        elif(self._dataExtention == ".cif"):
            parcer = pdb.MMCIFParser()

        try:
            structure = parcer.get_structure(self._pdb_id, self._path_data)
        except FileNotFoundError:
            print(FileNotFoundError)

        parcer_resolution = parcer.header["resolution"]
        parcer_keywords = parcer.header["keywords"]

        print("Sample name: " + str(parcer))
        print("parcer_resolution = ", parcer_resolution)
        print("parcer_keywords = ", parcer_keywords)
        print("structure = ", structure)

        self._structure = structure
        return structure

    def get_alpha_carbon_chain(self):
        alphaSkeleton = {}

        for model in self._structure.get_list():
            for chain in model.get_list():
                ch_alpha_C = []

                for residue in chain.get_list():
                    #gprint("residute = ", residue.__repr__())
                    if residue.has_id("CA"):
                        ca = residue["CA"]
                        ch_alpha_C.append(ca)

                alphaSkeleton[chain.get_id()] = ch_alpha_C
        print("alphaSkeleton = ", alphaSkeleton)
        return alphaSkeleton





"""
        #self._targetExtention = hyperparameters.target_ext
        #self._path_target = self._directory  + "/" + self._dataName + self._targetExtention

    def get_curr_dir(self):
        return os.getcwd()

    def convert_pdb_to_cdv(self):

        with open(self._path_data) as file:
            file_csv = df([line.strip().split() for line in file])
            file_csv.to_csv(self._path_data.split(".")[0] + self._targetExtention)
            file.close()

    def get_titles(self, data):
        return list(data.drop_duplicates())

    def parce_by_titles(self):

        file = pd.DataFrame(pd.read_csv(self._path_target))
        file.drop(file.columns[[0]], axis=1, inplace=True)      #удаляем первый столбец с нумировкой

        self._titles = self.get_titles(file[file.columns[0]])         #список заголовков


        for title in self._titles:
            mask = file[file.columns[0]] == title
            file_atom = pd.DataFrame(file[mask])                #выделяем строки с нужным заголовком
            file_atom.to_csv(self._directory  + "/" + self._dataName + "_{}".format(title) + self._targetExtention)

    def convert_data(self):
        self.convert_pdb_to_cdv()
        self.parce_by_titles()

    def read_atoms(self):
        tmp_path = self._directory  + "/" + self._dataName + "_ATOM" + self._targetExtention

        file = pd.DataFrame(pd.read_csv(tmp_path))
        file.drop(file.columns[[0]], axis=1, inplace=True) #удаляем нумерацию лишнюю
        data = file.dropna(axis=1, how='all')
        print(data)

        compiles_atoms = []
        compiles_molecules = []
        for index, row in data.iterrows():
            if(compiles_atoms == [] or
                    compiles_atoms[len(compiles_atoms)-1].get_identifier() == row.iloc[4]): # смотрим на смену идентификатора
                compiles_atoms.append(Atom.Atom(row.iloc))
            else:
                print("molecule ended ", len(compiles_atoms))
                compiles_molecules.append(Molecule.Molecule(compiles_atoms, compiles_atoms[0].get_identifier()))
                compiles_atoms = []

        print("molecule ended ", len(compiles_atoms))
        compiles_molecules.append(Molecule.Molecule(compiles_atoms, compiles_atoms[0].get_identifier()))
        del compiles_atoms

        print(compiles_molecules)
        return compiles_molecules
"""
