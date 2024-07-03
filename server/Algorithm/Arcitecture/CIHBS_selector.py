from .CIHBS import BaseCIHBS
from .BaseEnums import CHIBSBond


class SelectorCIHBS(BaseCIHBS):
    def __init__(self):
        super().__init__()

    def setParameters(self, new_cihbs):
        self.new_cihbs = new_cihbs

    # получить сложные ССИВС между аминокислотами
    def get_selector_map(self):
        if self.new_cihbs:
            dictinary = {'Hydrogen': {**self.getR_O(), **self.getR_R(), **self.getOthers()},
                         'Physics': self.getNH_O(),
                         'Bridges': self.getBridges()}
            return dictinary
        raise ValueError("Возвращается пустой объект.")

    # получение физических операторов между i и i-n элементами пентофрагментов
    def getNH_O(self):
        def mapping(distance):
            return [[i[0].get_serial_number(), i[1].get_serial_number(), i[3]] for i in self.new_cihbs
                    if
                    abs(int(i[0].get_parent().get_id()[1]) - int(i[1].get_parent().get_id()[1])) == distance
                    and i[0].get_name() == 'N'
                    and i[1].get_name() == 'O']

        return {'NH_Oi_3': mapping(3),
                'NH_Oi_4': mapping(4),
                'NH_Oi_5': mapping(5),
                'NH_Oi_6': mapping(6)}

    # получение ССИВС между iм остатком и i-n кислородом пентофрагментов
    def getR_O(self):
        def mapping(distance):
            return [[i[0].get_serial_number(), i[1].get_serial_number(), i[3]] for i in self.new_cihbs
                    if
                    abs(int(i[0].get_parent().get_id()[1])) - int(i[1].get_parent().get_id()[1]) == 3
                    and len(i[0].get_name()) >= 2
                    and i[1].get_name() == 'O']

        return {'Ri_Oi_3': mapping(3),
                'Ri_Oi_4': mapping(4),
                'Ri_Oi_5': mapping(5),
                'Ri_Oi_6': mapping(6)}

    # получение ССИВС между iм остатком и i-n остатком
    def getR_R(self):
        def mapping(distance):
            return [[i[0].get_serial_number(), i[1].get_serial_number(), i[3]] for i in self.new_cihbs
                    if
                    abs(int(i[0].get_parent().get_id()[1])) - int(i[1].get_parent().get_id()[1]) == distance
                    and len(i[0].get_name()) >= 2
                    and len(i[1].get_name()) >= 2]
        return {'Ri_Ri_3': mapping(3),
                'Ri_Ri_4': mapping(4),
                'Ri_Ri_5': mapping(5),
                'Ri_Ri_6': mapping(6)}

    # получение остальных ССИВС
    def getOthers(self):
        return {'Others': [[i[0].get_serial_number(), i[1].get_serial_number(), i[3]] for i in self.new_cihbs
                           if
                           i[2] != CHIBSBond.physicalOperator.value.get("type")
                           and ((len(i[0].get_name()) >= 2 and i[1].get_name() != 'O')
                                or (len(i[0].get_name()) == 1 and len(i[1].get_name()) >= 2))]}

    def getBridges(self):
        return [[i[0].get_serial_number(), i[1].get_serial_number(), i[3]] for i in self.new_cihbs
                if i[2] == CHIBSBond.residueBridge.value.get("type")]
