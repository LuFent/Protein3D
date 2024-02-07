from test_tube import HyperOptArgumentParser

class HyperParameters():
    def __init__(self):
        self.parser = HyperOptArgumentParser(strategy='grid_search')

    #@staticmethod
    def add_parameters(self, root_dir):

        #параметры файлов
        self.parser.opt_list('--data_ext', default=".pdb", options=[".pdb1"], type=str, tunable=True)
        self.parser.opt_list('--target_ext', default=".pdb", type=str, tunable=True)
        self.parser.opt_list('--data_dir', default='/data', type=str, tunable=True)
        self.parser.opt_list('--proj_dir', default=root_dir, type=str, tunable=True)
        self.parser.opt_list('--data_name', default='8fvs', type=str, tunable=True) #8fvs  7n1l

        #параметры командной строки
        #self.parser.add_argument('--data_name', type=str)
        #self.parser.add_argument('--data_dir', default='./data', type=str)

        return self.parser