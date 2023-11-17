import logging
import os

from data.Data_Converter import Data_Converter
from data.Data_Parcer import HyperParameters


if __name__ == "__main__":
    """
    где то тут должно реагировать на событие
    """
    # logging.basicConfig(level=logging.INFO)

    root_dir = os.path.dirname(os.path.realpath(__file__)).replace(os.sep, "/")
    print("root_dir = ", root_dir)

    # создаем параметры по умолчанию
    parser = HyperParameters().Add_Parameters(root_dir)
    hyperparams = parser.parse_args()
    print(hyperparams)

    # парсим
    data_converter = Data_Converter(hyperparams)
    data_converter.Convert_data()

    # читаем данные и составляем молекулы
    molecules = data_converter.Read_Atoms()
