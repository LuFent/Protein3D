import math

syze_vector=2
max_vector=2
epsilon=0.00001
Vector = []             #=array[0..syze_vector] of extended;
Vectors = []            #=array[0..max_vector] of Vector;
Index_Miller = []       #=array[0..syze_vector] of integer;
Index_Miller_Brave = [] #=array[0..syze_vector+1] of integer;

class V_algebra:

    def Radian(self, alpha):
        return (math.pi/180) * alpha

    def PowerDunction(self, x, a):
        return math.exp(x*math.log(a))

    def Addc(self, i, k, max):
        i = i + k
        if (i >= 0) and (i <= max):
            return i
        if i > max:
            return i - (max + 1)
        return max + i + 1

    def Add_c(self, i, k, min, max, flag):
        i = i + k
        if (i >= min) and (i <= max):
            return 0, i     #flag = 0
        if i > max:
            return 1, (min + i) - (max + 1)     #flag = 1
        return -1, (max + 1) - (min - i)

    def Compare(self, x, y, MinX, MinY, MaxX, MaxY): #!!!!!!!!!
        return -1

    def Vect_Modul(self, array):
        return math.sqrt(sum(i^2 for i in array))

    def Scal_Mult(self, array_1, array_2):
        return sum(map(lambda x, y: x * y, array_1, array_2))

    def CheckForLess(self, array):
        return (all(x < epsilon for x in array))

    def I_Miller_Brave(self, array_1, array_2):
        array_2[0] = 2 * array_1[0] - array_1[1]
        array_2[1] = 2 * array_1[1] - array_1[0]
        array_2[2] = -(array_1[0] + array_1[1])
        array_2[3] = 3 * array_1[2]

        if(all(x < epsilon for x in array_2)):
            array_2 = [int(i / 2) for i in array_2 if i < epsilon]
        if (all(x < epsilon for x in array_2)):
            array_2 = [int(i / 3) for i in array_2 if i < epsilon]

        return array_2

    def I_Brave_Miller(self, array_1, array_2):
        array_1[0] = array_2[0] - array_2[2]
        array_1[1] = array_2[1] - array_2[2]
        array_1[2] = array_2[3]

        if (all(x < epsilon for x in array_1)):
            array_1 = [int(i / 2) for i in array_2 if i < epsilon]
        if (all(x < epsilon for x in array_1)):
            array_1 = [int(i / 3) for i in array_2 if i < epsilon]

        return array_1

    def Numb_Mult(self, k, array_1):
        return [k*i for i in array_1]
    #Numb_Mult_Int не нужен

    def Vect_Add(self, array_1, array_2):
        return list(map(lambda x, y: x + y, array_1, array_2))
    #Vect_Add_Int не нужен

    def Vector_Mult(self, array_1, array_2):
        c = []
        c[0] = array_1[1] * array_2[2] - array_1[2] * array_2[1]
        c[1] = array_1[2] * array_2[0] - array_1[0] * array_2[2]
        c[2] = array_1[0] * array_2[1] - array_1[1] * array_2[0]
        return c

    def Smesh_Mult(self, a, b, c):
        return self.Scal_Mult(a, self.Vector_Mult(b, c))

    def Kr_Vectors(self, beta, alfa, gamma, a, b, c):
        Vcskr = [[0] * 2 for i in range(2)]

        Vcskr[0][0] = a * math.sqrt(1 - math.cos(beta)**2 - (math.sin(alfa) * math.cos(gamma))**2)
        Vcskr[0][1] = a * math.sin(alfa) * math.cos(gamma)
        Vcskr[0][2] = a * math.cos(beta)
        Vcskr[1][0] = 0
        Vcskr[1][1] = b * math.sin(alfa)
        Vcskr[1][2] = b * math.cos(alfa)
        Vcskr[2][0] = 0
        Vcskr[2][1] = 0
        Vcskr[2][2] = c
        return Vcskr

    def Ob_Vectors(self, Vcskr, Vcsob):
        V = self.Smesh_Mult(Vcskr[0], Vcskr[1], Vcskr[2])
        for i in range(3):
            Vcsob[i] = self.Vector_Mult(Vcskr[(i+1)%3], Vcskr[(i+2)%3])
        return Vcsob/V

    def Kr_uvw(self, array_1, array_2, a): # аналогичен Ob_hkl
        a = self.Numb_Mult(array_1[0],array_2[0])
        a1 = self.Numb_Mult(array_1[1],array_2[1])
        a2 = self.Numb_Mult(array_1[2],array_2[2])
        a = self.Vect_Add(a1, a)
        a = self.Vect_Add(a2, a)
        return a #что вернуть?



