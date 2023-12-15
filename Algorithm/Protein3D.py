import logging
import os

from data.GraphStructure import Data_Structure
from data.Data_Parcer import HyperParameters
from Bio import PDB as pdb


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
    data_converter = Data_Structure(hyperparams)
    structure = data_converter.get_structure()
    alphaSkeleton = data_converter.get_alpha_carbon_chain()