from Bio import PDB as pdb

class Data_Structure():

    def __init__(self, hyperparameters):
        self._directory = hyperparameters.proj_dir + hyperparameters.data_dir
        self._dataName = hyperparameters.data_name
        self._dataExtention = hyperparameters.data_ext

        self._path_data = self._directory  + "/" + self._dataName + self._dataExtention
        self._pdb_id = self._directory  + "/" + self._dataName


    def get_structure(self):
        if(self._dataExtention == ".pdb"):
            parcer = pdb.PDBParser(QUIET=True)
        elif(self._dataExtention == ".cif"):
            parcer = pdb.MMCIFParser(QUIET=True)

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

        return structure

    def formJSON(self):
        pass