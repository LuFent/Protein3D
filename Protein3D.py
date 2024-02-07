import os

from data.Data_Parcer import Data_Structure
from data.Hyperparameters import HyperParameters
from Arcitecture import CIHBS

from Bio.PDB import Selection
import json

class Protein3D:
    def __init__(self, structure, path):
        self._structure = structure
        self._pathJSON = path


    def get_alpha_carbon_chain(self):
        return [ atom.get_serial_number() for atom in Selection.unfold_entities(self._structure, 'A') if atom.get_name() == "CA" ]

    def getCHIBS(self):

        cihbsObj = CIHBS.CIHBS()
        innerCIHBS = cihbsObj.checkInnerGroups(self._structure)
        # cihbsObj.checkMultiGroups()
        from itertools import chain
        return list(sorted(set(chain.from_iterable(innerCIHBS))))

    def sendJSON(self, show_atoms,  extra_bonds):
        """
            {
                "show_atoms" : [id атомов которые надо оставить!],
                "extra_bonds": [Bond(..), Bond(..), Bond(..), ...]
            }
        """

        data = {
                "show_atoms" : show_atoms,
                "extra_bonds" : extra_bonds
                }
        with open(self._pathJSON, "w") as outfile:
            json.dump(data, outfile)

if __name__ == '__main__':
    """
    где то тут должно реагировать на событие
    """
    #logging.basicConfig(level=logging.INFO)

    root_dir = os.path.dirname(os.path.realpath(__file__)).replace(os.sep, '/')
    print("root_dir = ", root_dir)

    #создаем параметры по умолчанию
    parser = HyperParameters().add_parameters(root_dir)
    hyperparams = parser.parse_args()
    print(hyperparams)

    #получаем структуру
    protein = Protein3D(Data_Structure(hyperparams).get_structure(), "data/output.json")

    #альфаскелет
    alphaSkeleton = protein.get_alpha_carbon_chain()

    #водородные связи
    cihbs = protein.getCHIBS()

    #формируем и сохраняем JSON
    protein.sendJSON(cihbs, [])