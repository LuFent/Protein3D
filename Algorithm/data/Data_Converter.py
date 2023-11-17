import pandas as pd
from pandas import DataFrame as df
import os

from Arcitecture import Molecule, Atom

class Data_Converter():

    def __init__(self, hyperparameters):
        self._directory = hyperparameters.proj_dir + hyperparameters.data_dir
        self._dataName = hyperparameters.data_name
        self._dataExtention = hyperparameters.data_ext
        self._targetExtention = hyperparameters.target_ext
        self._path_data = self._directory  + "/" + self._dataName + self._dataExtention
        self._path_target = self._directory  + "/" + self._dataName + self._targetExtention

    def Get_Curr_Dir(self):
        return os.getcwd()

    def Convert_PDB_to_CSV(self):

        with open(self._path_data) as file:
            file_csv = df([line.strip().split() for line in file])
            file_csv.to_csv(self._path_data.split(".")[0] + self._targetExtention)
            file.close()

    def Get_Titles(self, data):
        return list(data.drop_duplicates())

    def Parce_by_Titles(self):

        file = pd.DataFrame(pd.read_csv(self._path_target))
        file.drop(file.columns[[0]], axis=1, inplace=True)      #удаляем первый столбец с нумировкой

        self._titles = self.Get_Titles(file[file.columns[0]])         #список заголовков


        for title in self._titles:
            mask = file[file.columns[0]] == title
            file_atom = pd.DataFrame(file[mask])                #выделяем строки с нужным заголовком
            file_atom.to_csv(self._directory  + "/" + self._dataName + "_{}".format(title) + self._targetExtention)

    def Convert_data(self):
        self.Convert_PDB_to_CSV()
        self.Parce_by_Titles()

    def Read_Atoms(self):
        tmp_path = self._directory  + "/" + self._dataName + "_ATOM" + self._targetExtention

        file = pd.DataFrame(pd.read_csv(tmp_path))
        file.drop(file.columns[[0, 2]], axis=1, inplace=True) #удаляем нумерацию 2 раза и заголовок
        data = file.dropna(axis=1, how='all')
        #print(data)

        compiles_atoms = []
        compiles_molecules = []
        for index, row in data.iterrows():

            if(compiles_atoms == [] or
                    compiles_atoms[len(compiles_atoms)-1].Get_identifier() == row.iloc[2]): # смотрим на смену идентификатора
                compiles_atoms.append(Atom.Atom(row.iloc))
            else:
                print("molecule ended ", len(compiles_atoms))
                compiles_molecules.append(Molecule.Molecule(compiles_atoms, compiles_atoms[0].Get_identifier()))
                compiles_atoms = []

        print("molecule ended ", len(compiles_atoms))
        compiles_molecules.append(Molecule.Molecule(compiles_atoms, compiles_atoms[0].Get_identifier()))
        del compiles_atoms

        print(compiles_molecules)
        return compiles_molecules

