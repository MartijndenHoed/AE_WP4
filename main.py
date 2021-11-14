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

TotalFz_max = PLong
TotalFx_max = Plat
TotalFy_max = Plat
TotalMy_max = 0 #0.5 * RTGlength * TotalFx_max



LugFzMax = TotalFz_max/4 + TotalMy_max*lugSpacingX*0.*0.5 #The loads on the individual lugs
LugFxMax = TotalFx_max/4 + TotalMy_max*lugSpacingZ*0.5*0.5
LugFyMax = TotalFy_max/4


#global Variables:
globalW = 0
globalFy = 0
globalD = 0
globalAav = 0
globalSigma = 0

class Material:
    def __init__(self,E,sigmaYield,density):
        self.density = density
        self.E = E
        self.sigmaYield = sigmaYield

class Config:
    def __init__(self,t,W,D,mass):
        self.t = t
        self.W = W
        self.D = D
        self.mass = mass

def iterateConfigs():
    bestMass = 1000
    config = -1
    Diameter = 2 * calcAxleRadius(aluminium)
    for x in range(0,300):
        W = x/10000
        #print("checking W: " + str(W))
        tempConfig = generateLugConfiguration(aluminium,W,0.2,Diameter)
        if(tempConfig != -1):
            tempMass = tempConfig.mass
            if(bestMass > tempMass):
                bestMass = tempMass
                config = tempConfig

    return config

def generateLugConfiguration(lugMat,W,l,Diameter):
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
    totalAxleForce = math.sqrt(LugFzMax**2 + LugFyMax**2)
    axleArea = totalAxleForce / (2*stressAllow)
    r = math.sqrt(axleArea/math.pi)

    #pain
    return r

def calcMinimumFlangeThickness(material,Diameter,width,l):
    #loop through all possible Rtr and Rtx ratios and calculate the absolute minimum flange thickness
    steps = 1000

    smallestLugT = 10000 #LARGEness
    smallestTTension = 100
    for x in range(1,steps):
        Rax = x/steps
        Rtr = calcRtr(Rax)

        tAxial = calcFlangeThicknessTensionAxial(material,Diameter,width,Rax)
        tTransverse = calcFlangeThicknessTensionTransverse(material,Diameter,width,Rtr)
        minTTension = max(tAxial,tTransverse)
        smallestTTension = min(minTTension,smallestTTension)
        #print("Rax: " + str(Rax) + "Rtr: " + str(Rtr) + "minT: " + str(minTTension))

    tBendingXZ = calcFlangeThicknessBendingXZ(material,Diameter,width,l)
    #print( tBendingXZ)

    minT = max(smallestTTension,tBendingXZ)


    return minT


def calcFlangeThicknessTensionAxial(material,D,W,Rax):
    # calculate the flange thickness for axial loads
    Ratio = W/D

    k = Kt_val_axial(Ratio)
    #print(k)

    load = (0.5*LugFyMax)/Rax
    t = load/((W-D)*k*material.sigmaYield)

    return t

def calcFlangeThicknessTensionTransverse(material,D,W,Rtr): #should be finished
    #calculate the flange thickness for transverse loads

    A1 = (1-0.5*math.sqrt(2))*(D/2) + (W-D)*0.5
    A2 = (W-D)*0.5
    A3 = (W - D) * 0.5
    A4 = (1 - 0.5 * math.sqrt(2)) * (D / 2) + (W - D) * 0.5
    Aav = 6/((3/A1) + (1/A2) + (1/A3) + (1/A4))

    load = (LugFyMax*0.5)/Rtr
    global globalW
    globalW = W
    global globalD
    globalD = D
    global globalFy
    globalFy = load
    global globalSigma
    globalSigma = material.sigmaYield
    global globalAav
    globalAav = Aav

    t = linearSolve(transverseEqLeft,transverseEqRight,0.0001,1,4,10)

    return t

def transverseEqLeft(t):
    Abr = globalD * t
    ratio = Abr / globalAav
    kt = Kt_val_transverse(ratio)
    value = globalFy/(kt)

    return value

def transverseEqRight(t):


    return t * globalD * globalSigma

def Kbrycalc(D, t, W):
    Param_eD = W/(2*D) - 0.5
    if Param_eD > 3.5:
        print("Param_eD is out of boundaries")
        return
    if t/D >= 0.6:
        K_bry = -0.023 + 3.0145753647164595 * (Param_eD ** 1) + -6.131229073233705 * (Param_eD ** 2) + 13.81165851383463 * (Param_eD ** 3) + -19.791798609852435 * (Param_eD ** 4) + 17.31017221490646 * (Param_eD ** 5) + -9.537564894148687 * (Param_eD ** 6) + 3.334368780153838 * (Param_eD ** 7) + -0.717854396732335 * (Param_eD ** 8) + 0.0868080698849476 * (Param_eD ** 9) + -0.00451187785287549 * (Param_eD ** 10)
    elif t/D > 0.4:
        K_bry = -0.035 + 3.0564119946981476 * (Param_eD ** 1) + -5.5035073918290305 * (Param_eD ** 2) + 10.327115088375328 * (Param_eD ** 3) + -12.27065902515286 * (Param_eD ** 4) + 8.624285881608284 * (Param_eD ** 5) + -3.6427509060765177 * (Param_eD ** 6) + 0.9117687857917013 * (Param_eD ** 7) + -0.12478081617805703 * (Param_eD ** 8) + 0.00719910832024418 * (Param_eD ** 9)
    elif t/D > 0.3:
        K_bry = -0.028 + 2.8922672056827325 * (Param_eD ** 1) + -5.125369774158555 * (Param_eD ** 2) + 11.274747882284034 * (Param_eD ** 3) + -16.28695512496961 * (Param_eD ** 4) + 13.988407150257245 * (Param_eD ** 5) + -7.366631725247057 * (Param_eD ** 6) + 2.409683108093695 * (Param_eD ** 7) + -0.4772937520465902 * (Param_eD ** 8) + 0.052345031985765234 * (Param_eD ** 9) + -0.002434048328197378 * (Param_eD ** 10)
    elif t/D > 0.2:
        K_bry = -0.018 + 2.924623014945673 * (Param_eD ** 1) + -5.933006018725241 * (Param_eD ** 2) + 15.286866566089987 * (Param_eD ** 3) + -25.560447367408067 * (Param_eD ** 4) + 25.26597829363366 * (Param_eD ** 5) + -15.302041990531661 * (Param_eD ** 6) + 5.76702212228471 * (Param_eD ** 7) + -1.3209525478559079 * (Param_eD ** 8) + 0.16843991622532656 * (Param_eD ** 9) + -0.009174555773228737 * (Param_eD ** 10)
    elif t/D > 0.15:
        K_bry = -0.015 + 2.6746380373328176 * (Param_eD ** 1) + -3.8095514487546716 * (Param_eD ** 2) + 8.623219742992589 * (Param_eD ** 3) + -15.750885241422282 * (Param_eD ** 4) + 17.115175048197273 * (Param_eD ** 5) + -11.18034715883286 * (Param_eD ** 6) + 4.4747996668137695 * (Param_eD ** 7) + -1.0765463102425024 * (Param_eD ** 8) + 0.14305215852155762 * (Param_eD ** 9) + -0.008072914212872394 * (Param_eD ** 10)
    elif t/D > 0.12:
        K_bry = -0.015 + 2.26199460291689 * (Param_eD ** 1) + 0.2998846940719729 * (Param_eD ** 2) + -5.896253454827716 * (Param_eD ** 3) + 8.866390732323316 * (Param_eD ** 4) + -6.659250538483327 * (Param_eD ** 5) + 2.8922806130518577 * (Param_eD ** 6) + -0.7356872271946134 * (Param_eD ** 7) + 0.10187317354853766 * (Param_eD ** 8) + -0.005936465946266028 * (Param_eD ** 9)
    elif t/D > 0.1:
        K_bry = -0.012 + 2.1958873981812173 * (Param_eD ** 1) + 1.4121860944569744 * (Param_eD ** 2) + -11.094326947388348 * (Param_eD ** 3) + 19.05242541029777 * (Param_eD ** 4) + -17.492127624958844 * (Param_eD ** 5) + 9.777977619561582 * (Param_eD ** 6) + -3.430942783487547 * (Param_eD ** 7) + 0.739575340442149 * (Param_eD ** 8) + -0.08966973195725639 * (Param_eD ** 9) + 0.004687489286406277 * (Param_eD ** 10)
    elif t/D > 0.08:
        K_bry = -0.000 + 2.047168809200642 * (Param_eD ** 1) + 2.8439054359210623 * (Param_eD ** 2) + -18.38619717436508 * (Param_eD ** 3) + 34.28307332550166 * (Param_eD ** 4) + -34.299009336084836 * (Param_eD ** 5) + 20.676428605730308 * (Param_eD ** 6) + -7.734726022194671 * (Param_eD ** 7) + 1.7588160610321202 * (Param_eD ** 8) + -0.22283882512474207 * (Param_eD ** 9) + 0.012071823552433879 * (Param_eD ** 10)
    elif t/D > 0.06:
        K_bry = -0.019 + 2.815705548826579 * (Param_eD ** 1) + -3.8386040984312917 * (Param_eD ** 2) + 0.5570512486788379 * (Param_eD ** 3) + 6.174683983615872 * (Param_eD ** 4) + -9.49021569548029 * (Param_eD ** 5) + 6.967071323446899 * (Param_eD ** 6) + -2.941962100212095 * (Param_eD ** 7) + 0.7284623056343401 * (Param_eD ** 8) + -0.09845184660891078 * (Param_eD ** 9) + 0.005614099708215215 * (Param_eD ** 10)
    else:
        print("t/D too small")
        return
    if K_bry < 0:
        K_bry = 0
    return K_bry


def calcFlangeThicknessBearingShearOut(material,Diameter,width):
    # calculate the flange thickness for bearing shear out

    return t

def calcFlangeThicknessBendingX(material,Diameter,width,l):
    # calculate the flange thickness for bending
    t = 0
    dt = 0.001
    M = LugFzMax * l
    Bending_stressX = 0
    while Bending_stressX < material.sigmaYield:
        t = t + dt
        I = t * width ** 3 / 12
        Bending_stressX = M * (width / 2) / I
        print(t)

    return t

def calcFlangeThicknessBendingZ(material,Diameter,width,l):
    # calculate the flange thickness for bending
    t = 0
    dt = 0.001
    M = LugFxMax * l
    Bending_stressZ = 0
    while Bending_stressZ < material.sigmaYield:
        t = t + dt
        I = width * t ** 3 / 12
        Bending_stressZ = M * (t / 2) / I
        print(t)

    return t


def calcFlangeThicknessBendingXZ(material,Diameter,w,l):
    M1 = LugFzMax * l
    M2 = LugFxMax * l
    sigma = material.sigmaYield
    return  (-((6*M1)/(w**2)) - math.sqrt( ((6*M1)/(w**2))**2 + 4 * sigma * ((6*M2)/(w)) ))/(-2*sigma)

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
testAxialT = calcFlangeThicknessTensionAxial(aluminium,testR,0.05,1) #(material,Diameter,width,Rax)
print(testAxialT)
#config = generateLugConfiguration(aluminium,0.03,0.1,0.0017)
config = iterateConfigs()
print("mass " + str(config.mass))
print("diameter " + str(config.D))
print("W " + str(config.W))
print("T " + str(config.t))




