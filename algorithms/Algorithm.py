from .Mask import Mask


class Algorithm:
    def __init__(self):
        pass

    def _execute(self, s) -> Mask:
        return Mask()

    def execute(self, s) -> Mask:
        return self._execute(s)


class AlgorithmsStorage:
    def __init__(self):
        self.labels = set()
        self.code_names = set()
        self.algs = dict()

    def add_algorithm(self, alg: Algorithm):
        if not alg.code_name:
            raise ValueError("Algorithm must have an code_name")

        if not alg.label:
            raise ValueError("Algorithm must have an label")

        if alg.code_name in self.code_names:
            raise ValueError("Algorithm with such code_names already exists")

        if alg.label in self.labels:
            raise ValueError("Algorithm with such label already exists")

        self.labels.add(alg.label)
        self.code_names.add(alg.code_name)
        self.algs[alg.code_name] = alg


    def get_all_algorithms(self):
        algs = dict()
        for code_name, alg in self.algs.items():
            algs[code_name] = {"label": alg.label,
                               "icon": alg.icon}
        return algs

    def copy(self):
        algorithm_storage = AlgorithmsStorage()
        algorithm_storage.labels = self.labels
        algorithm_storage.code_names = self.code_names
        algorithm_storage.algs = self.algs
        return algorithm_storage


    def __getitem__(self, key):
        return self.algs.get(key, None)

    def __repr__(self):
        return str(self.algs)

    def __str__(self):
        return str(self.algs)

