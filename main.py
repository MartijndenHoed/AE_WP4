# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

#bruhhh

import math

print("bruhh")


##all SI units please
g = 9.80665


RTGmass = 52.8
RTGlength = 1.07
RTGdiameter = 0.47

GStaticlong = 4.55
GStaticLat = 0.25

GDynamicLong = 1.0
GDynamicLat = 0.8

SafetyFactorLaunch = 1.25

GLongMax = (GStaticlong + GDynamicLong) * SafetyFactorLaunch
GlatMax = (GStaticLat + GDynamicLat) * SafetyFactorLaunch

PLong = GLongMax * RTGmass * g
Plat = GlatMax * RTGmass * g

class material:
    def __init__(self,E,sigmaYield,density):
        self.density = density
        self.E = E
        self.sigmaYield = sigmaYield

def generateLugConfiguration(axleMat,lugMat,W):
    #this should be enough variables to generate a thickness and calculate a weight

    return config

def calcRtr(Rax):
    # this calculates the Rtr  ratio from a given Rax ratio (with a safety margin of 1)
    x = 1-(math.pow(Rax,1.6))
    Rtr = math.pow(x,1/1.6)
    return Rtr


def calcAxleRadius(material):
    #calculate the minium radius of the (circular) axle


    stressAllow = material.sigmaYield
    axleArea = Plat / (2*stressAllow)
    r = math.sqrt(axleArea/math.pi)

    #pain
    return r

def calcMinimumFlangeThickness(material,Diameter,width):
    #loop through all possible Rtr and Rtx ratios and calculate the absolute minimum flange thickness
    stepSize = 0.01

    for Rax in range(0,1,stepSize):
        Rtr = calcRtr(Rax)

        tAxial = calcFlangeThicknessTensionAxial(material,Diameter,width,Rax)




    return t


def calcFlangeThicknessTensionAxial(material,Diameter,height,Rax):
    # calculate the flange thickness for axial loads
    k = 0.7 ##arbitrary, function of height , Diameter


    load = (0.5*Fy_max)/Rax
    t = load/(height*k*material.sigmaYield)

    return t

def calcFlangeThicknessTensionTransverse(material,Diameter,width,Rtr):
    #calculate the flange thickness for transverse loads
    return t


def calcFlangeThicknessBearingShearOut(material,Diameter,width):
    # calculate the flange thickness for bearing shear out

    return t

##calculate the thickness for a certain width, take the minimum value, compare the mass for all materials.

def linearSolve(eq1, eq2, lower, upper, cycles, steps):
    #this is a simple numerical solver
    upperAnswer = upper
    lowerAnswer = lower

    for y in range(0, cycles):

        diffSign = 0
        cycleDone = False

        for x in range(0, steps):
            if (not cycleDone):
                i = lowerAnswer + (upperAnswer - lowerAnswer) / steps * x
                # print(x)
                val1 = eq1(i)
                val2 = eq2(i)
                diff = val1 - val2
                # print(diff)
                if (diff == 0):
                    print("exact: " + str(i))
                    return i
                sign = diff / abs(diff)
                if (diffSign == 0):
                    diffSign = sign

                if (diffSign != sign):
                    cycleDone = True
                    lowerAnswer = lowerAnswer + (upperAnswer - lowerAnswer) / steps * (x - 1)
                    upperAnswer = i
                    # print("lower " + str(lowerAnswer))
                    # print("upper " +  str(upperAnswer))

    diff1 = eq2(lowerAnswer) - eq1(lowerAnswer)
    diff2 = eq2(upperAnswer) - eq1(upperAnswer)

    diff3 = upperAnswer - lowerAnswer

    return lowerAnswer + diff1 / ((diff1 - diff2) / diff3)


aluminium = Material(70000000000,200000000,2.6) ##just some test values
testR = calcAxleRadius(aluminium)
print(testR)
testAxialT = calcFlangeThicknessTensionAxial(aluminium,testR,0.4,1) #(material,Diameter,width,Rax)
print(testR)



def Ktu_val(material, Ratio1):
    if material == Material_curve1:
        Ktu = 0.9180553870465704 * (Ratio1 ** 1) + 8.18122273185643 * (Ratio1 ** 2) + -57.7056377483359 * (Ratio1 ** 3) + 231.21388530032073 * (Ratio1 ** 4) + -578.4496388173648 * (Ratio1 ** 5) + 921.2687606248473 * (Ratio1 ** 6) + -928.1664126745236 * (Ratio1 ** 7) + 570.5859381665941 * (Ratio1 ** 8) + -194.87683183789522 * (Ratio1 ** 9) + 28.300734922855554 * (Ratio1 ** 10)
    elif material == Material_curve2:
        Ktu = 1.1540799943441162 * (Ratio1 ** 1) + 3.5658936323271035 * (Ratio1 ** 2) + -13.955444513318298 * (Ratio1 ** 3) + 16.23858184466644 * (Ratio1 ** 4) + 19.999288298125066 * (Ratio1 ** 5) + -75.33757395030716 * (Ratio1 ** 6) + 82.87193625333839 * (Ratio1 ** 7) + -41.04677495808528 * (Ratio1 ** 8) + 7.779874724243314 * (Ratio1 ** 9)
    elif material == Material_curve4:
        Ktu = 1.0968401953214486 * (Ratio1 ** 1) + 5.911055078229481 * (Ratio1 ** 2) + -41.276480666988846 * (Ratio1 ** 3) + 170.3631152105354 * (Ratio1 ** 4) + -485.8393789099573 * (Ratio1 ** 5) + 946.0509482390116 * (Ratio1 ** 6) + -1188.580514624009 * (Ratio1 ** 7) + 901.8312818112099 * (Ratio1 ** 8) + -372.68816788151344 * (Ratio1 ** 9) + 64.20378959440313 * (Ratio1 ** 10)
    elif material == Material_curve5:
        Ktu = 0.8935998289127944 * (Ratio1 ** 1) + 5.273317497891158 * (Ratio1 ** 2) + -28.55089032590277 * (Ratio1 ** 3) + 52.05627723859402 * (Ratio1 ** 4) + 28.107639294058117 * (Ratio1 ** 5) + -256.71084383513744 * (Ratio1 ** 6) + 425.1831418346137 * (Ratio1 ** 7) + -339.90178976325706 * (Ratio1 ** 8) + 136.91477026647703 * (Ratio1 ** 9) + -22.29608887194499 * (Ratio1 ** 10)
    elif material == Material_curve6:
        Ktu = 0.9564655594429534 * (Ratio1 ** 1) + 5.166345227885259 * (Ratio1 ** 2) + -34.01201870580809 * (Ratio1 ** 3) + 104.20957326598887 * (Ratio1 ** 4) + -180.16817428435706 * (Ratio1 ** 5) + 184.97898791984153 * (Ratio1 ** 6) + -112.11763123206947 * (Ratio1 ** 7) + 37.11718579286003 * (Ratio1 ** 8) + -5.182081049544365 * (Ratio1 ** 9)

