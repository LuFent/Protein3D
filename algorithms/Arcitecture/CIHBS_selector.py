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
                         'Physics': self.getNH_O()}
            return dictinary
        raise ValueError("Возвращается пустой объект.")

    # получение физических операторов между i и i-n элементами пентофрагментов
    def getNH_O(self):
        return {'NH_Oi_3': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1]) - int(i[1].get_parent().get_id()[1])) == 3
                            and i[0].get_name() == 'N'
                            and i[1].get_name() == 'O'],
                'NH_Oi_4': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1]) - int(i[1].get_parent().get_id()[1])) == 4
                            and i[0].get_name() == 'N'
                            and i[1].get_name() == 'O'],
                'NH_Oi_5': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1]) - int(i[1].get_parent().get_id()[1])) == 5
                            and i[0].get_name() == 'N'
                            and i[1].get_name() == 'O'],
                'NH_Oi_6': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1]) - int(i[1].get_parent().get_id()[1])) == 6
                            and i[0].get_name() == 'N'
                            and i[1].get_name() == 'O']}

    # получение ССИВС между iм остатком и i-n кислородом пентофрагментов
    def getR_O(self):
        return {'Ri_Oi_3': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1])) - int(i[1].get_parent().get_id()[1]) == 3
                            and len(i[0].get_name()) >= 2
                            and i[1].get_name() == 'O'],
                'Ri_Oi_4': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1]) - int(i[1].get_parent().get_id()[1])) == 4
                            and len(i[0].get_name()) >= 2
                            and i[1].get_name() == 'O'],
                'Ri_Oi_5': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1]) - int(i[1].get_parent().get_id()[1])) == 5
                            and len(i[0].get_name()) >= 2
                            and i[1].get_name() == 'O'],
                'Ri_Oi_6': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1]) - int(i[1].get_parent().get_id()[1])) == 6
                            and len(i[0].get_name()) >= 2
                            and i[1].get_name() == 'O']}

    # получение ССИВС между iм остатком и i-n остатком
    def getR_R(self):
        return {'Ri_Ri_3': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1])) - int(i[1].get_parent().get_id()[1]) == 3
                            and len(i[0].get_name()) >= 2
                            and len(i[1].get_name()) >= 2],
                'Ri_Ri_4': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1]) - int(i[1].get_parent().get_id()[1])) == 4
                            and len(i[0].get_name()) >= 2
                            and len(i[1].get_name()) >= 2],
                'Ri_Ri_5': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1]) - int(i[1].get_parent().get_id()[1])) == 5
                            and len(i[0].get_name()) >= 2
                            and len(i[1].get_name()) >= 2],
                'Ri_Ri_6': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                            if
                            abs(int(i[0].get_parent().get_id()[1]) - int(i[1].get_parent().get_id()[1])) == 6
                            and len(i[0].get_name()) >= 2
                            and len(i[1].get_name()) >= 2]}

    # получение остальных ССИВС
    def getOthers(self):
        return {'Others': [[i[0].get_serial_number(), i[1].get_serial_number(), i[2].value] for i in self.new_cihbs
                           if
                           i[2] != CHIBSBond.physicalOperator
                           and ((len(i[0].get_name()) >= 2 and i[1].get_name() != 'O')
                                or (len(i[0].get_name()) == 1 and len(i[1].get_name()) >= 2))]}
