import math
import matplotlib.pyplot as plt

g = 9.80665 #m/(s**2)


RTGmass = 52.8 #kg
RTGlength = 1.07 #m
RTGdiameter = 0.47 #m
RTGCGy = 0.5 * RTGdiameter

lugSpacingX = 0.4 #m
lugSpacingZ = RTGlength #m

GStaticlong = 4.55 #g
GStaticLat = 0.25 #g

GDynamicLong = 1.0 #g
GDynamicLat = 0.8 #g

SafetyFactorLaunch = 1.5

GLongMax = (GStaticlong + GDynamicLong) * SafetyFactorLaunch
GlatMax = (GStaticLat + GDynamicLat) * SafetyFactorLaunch

PLong = GLongMax * RTGmass * g
Plat = GlatMax * RTGmass * g

TotalFz_max = PLong
TotalFx_max = Plat
TotalFy_max = Plat
TotalMy_max = 0



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
    def __init__(self,E,sigmaYield,density,sigmaYieldShear):
        self.density = density
        self.E = E
        self.sigmaYield = sigmaYield
        self.sigmaYieldShear = sigmaYieldShear

class Config:
    def __init__(self,material,t,W,D,l):
        self.t = t
        self.W = W
        self.D = D
        self.mass = (W * l + 0.5 * math.pi * (0.5*W)**2 + math.pi * (D/2) ** 2) * t * material.density * 2
        self.l = l


def iterateW(lugMat,l,Diameter):
    RTGMomentArm = l - RTGCGy
    moment = TotalFz_max * RTGMomentArm
    global LugFyMax
    LugFyMax = TotalFy_max / 4 + abs(moment / (RTGlength*0.5))/4

    bestMass = math.inf
    steps = 100
    lowerW = 2 * Diameter
    upperW = 5 * Diameter

    wArr = []
    tArr = []

    for x in range(0,steps):
        W = (upperW - lowerW)/steps  * x + lowerW
        #print(W)
        #print(W/Diameter)
        t = iterateT(lugMat,W,l,Diameter)
        config = Config(lugMat,t,W,Diameter,l)
        if(config.mass < bestMass):
            bestConfig = config
            bestMass = config.mass
        wArr.append(W)
        tArr.append(config.mass)

    plt.plot(wArr,tArr)
    plt.show()

    return bestConfig

def iterateT(lugMat,W,l,Diameter):
    t=0
    stepSize = 0.0001



    stressMax = math.inf
    #while (stressMax > lugMat.sigmaYield):
    #    t += stepSize
    #    stressMax = calcTotalMaxStress(t,[Diameter,l,W,lugMat,sigmaYield])
        #print(str(flangeStress) + " " + str(bendingStress) + " " + str(stressMax))
        #print(bearingOutStress)

    t = linearSolve(calcStressDifference,zeroFunc,0.001,1,4,8,[Diameter,l,W,lugMat.sigmaYield])
    #print(str([Diameter,l,W,lugMat.sigmaYield]) + " " + str(t))
    #global wArr
    #global tArr
    #wArr.append(W)
    #tArr.append(t)



    return t

def calcTotalMaxStress(t,parameters):
    Diameter = parameters[0]
    l = parameters[1]
    W = parameters[2]
    flangeStress = calcFlangeStress(t, Diameter, W, l)
    bendingStress = calcBendingStress(t, W, l)
    bearingOutStress = calcBearingOutStress(t, Diameter, W)
    # print(flangeStress)
    return max(flangeStress, bendingStress, bearingOutStress)

def calcStressDifference(t,parameters):
    return (calcTotalMaxStress(t, parameters) - parameters[3])

def zeroFunc(x):
    return 0

def calcBearingOutStress(t,Diameter,W):
    K_bry = Kbrycalc(Diameter, t, W)
    A_br = Diameter * t
    load = TotalFy_max * 0.5
    stress = load / (K_bry * A_br)
    if (K_bry == -1):
        stress = math.inf
    return stress

def calcBendingStress(t,W,l):
    M1 = LugFzMax * l * 0.5
    M2 = LugFxMax * l * 0.5

    return (M1 * 0.5 * W)/((1/12) * t * W**3) + (M2 * 0.5 * t)/((1/12) * t**3 * W)

def calcFlangeStress(t,Diameter,width,l):
    steps = 1000

    smallestTTension = 100

    minStressAxial = math.inf
    minStressTransverse = math.inf
    minStressTotal = math.inf
    #print(str(width) + " " + str(width/Diameter) + " " + str(Kt_val_axial(width/Diameter)))

    for x in range(1, steps):
        Rax = x / steps
        Rtr = calcRtr(Rax)

        stressAxial = calcFlangeStressAxial(t,Diameter,width,Rax)
        stressTransverse = calcFlangeStressTransverse(t,Diameter,width,Rtr)

        minStressAxial = min(stressAxial,minStressAxial)
        minStressTransverse = min(stressTransverse,minStressTransverse)

        maxStressTotal = max(stressAxial,stressTransverse)
        minStressTotal = min(minStressTotal,maxStressTotal)
    return minStressTotal

def calcFlangeStressTransverse(t,D,W,Rtr): #should be finished
    #calculate the flange thickness for transverse loads

    A1 = (1-0.5*math.sqrt(2))*(D/2) + (W-D)*0.5
    A2 = (W-D)*0.5
    A3 = (W - D) * 0.5
    A4 = (1 - 0.5 * math.sqrt(2)) * (D / 2) + (W - D) * 0.5
    Aav = 6/((3/A1) + (1/A2) + (1/A3) + (1/A4))
    Abr = D*t

    #print(Aav)
    #print(Abr)
    #print(Aav/Abr)

    load = (LugFyMax*0.5)/Rtr
    k = Kt_val_transverse(Aav/Abr)
    return load/(Abr * k)


def calcFlangeStressAxial(t,D,W,Rax):
    Ratio = W/D
    k = Kt_val_axial(Ratio)
    #print(k)
    A = t*(W-D)
    load = (0.5*LugFyMax)/Rax
    return load/(A*k)

def calcRtr(Rax):
    # this calculates the Rtr  ratio from a given Rax ratio (with a safety margin of 1)
    x = 1-(math.pow(Rax,1.6))
    Rtr = math.pow(x,1/1.6)
    return Rtr


def calcAxleRadius(material):
    #calculate the minium radius of the (circular) axle


    stressAllow = material.sigmaYieldShear
    totalAxleForce = math.sqrt(LugFzMax**2 + LugFyMax**2)
    axleArea = totalAxleForce / (2*stressAllow)
    r = math.sqrt(axleArea/math.pi)

    #pain
    return r

def Kt_val_axial(R):
    Kt = 1.009 + -0.22902304135292462 * (R ** 1) + 0.35905460079484186 * (R ** 2) + -1.0729145294629987 * (R ** 3) + 1.2672573367948203 * (R ** 4) + -0.8003541023438316 * (R ** 5) + 0.3032426494074585 * (R ** 6) + -0.07136553672277213 * (R ** 7) + 0.010232930693775955 * (R ** 8) + -0.000818435700391848 * (R ** 9) + 0.000027954852209106408 * (R ** 10)
    if(R>5):
        return 0.00001
    return Kt

def Kt_val_transverse(r):
    Kt = 0.003 + 1.1711822557878742 * (r ** 1) + 2.058639823919016 * (r ** 2) + -17.838056169717532 * (r ** 3) + 62.31283878518559 * (r ** 4) + -118.14607854029032 * (r ** 5) + 129.78187628214891 * (r ** 6) + -82.82263597513247 * (r ** 7) + 28.55280116408801 * (r ** 8) + -4.118577606465635 * (r ** 9)
    if(r>1.4):
        return 1.29

    return Kt

def Kbrycalc(D, t, W):
    Param_eD = W/(2*D) - 0.5
    if t/D >= 0.6:
        K_bry = -0.023 + 3.0145753647164595 * (Param_eD ** 1) + -6.131229073233705 * (Param_eD ** 2) + 13.81165851383463 * (Param_eD ** 3) + -19.791798609852435 * (Param_eD ** 4) + 17.31017221490646 * (Param_eD ** 5) + -9.537564894148687 * (Param_eD ** 6) + 3.334368780153838 * (Param_eD ** 7) + -0.717854396732335 * (Param_eD ** 8) + 0.0868080698849476 * (Param_eD ** 9) + -0.00451187785287549 * (Param_eD ** 10)
        if Param_eD > 3.5:
            K_bry = 1.75
    elif t/D > 0.4:
        K_bry = -0.035 + 3.0564119946981476 * (Param_eD ** 1) + -5.5035073918290305 * (Param_eD ** 2) + 10.327115088375328 * (Param_eD ** 3) + -12.27065902515286 * (Param_eD ** 4) + 8.624285881608284 * (Param_eD ** 5) + -3.6427509060765177 * (Param_eD ** 6) + 0.9117687857917013 * (Param_eD ** 7) + -0.12478081617805703 * (Param_eD ** 8) + 0.00719910832024418 * (Param_eD ** 9)
        if Param_eD > 3.5:
            K_bry = 1.7
    elif t/D > 0.3:
        K_bry = -0.028 + 2.8922672056827325 * (Param_eD ** 1) + -5.125369774158555 * (Param_eD ** 2) + 11.274747882284034 * (Param_eD ** 3) + -16.28695512496961 * (Param_eD ** 4) + 13.988407150257245 * (Param_eD ** 5) + -7.366631725247057 * (Param_eD ** 6) + 2.409683108093695 * (Param_eD ** 7) + -0.4772937520465902 * (Param_eD ** 8) + 0.052345031985765234 * (Param_eD ** 9) + -0.002434048328197378 * (Param_eD ** 10)
        if Param_eD > 3.5:
            K_bry = 1.65
    elif t/D > 0.2:
        K_bry = -0.018 + 2.924623014945673 * (Param_eD ** 1) + -5.933006018725241 * (Param_eD ** 2) + 15.286866566089987 * (Param_eD ** 3) + -25.560447367408067 * (Param_eD ** 4) + 25.26597829363366 * (Param_eD ** 5) + -15.302041990531661 * (Param_eD ** 6) + 5.76702212228471 * (Param_eD ** 7) + -1.3209525478559079 * (Param_eD ** 8) + 0.16843991622532656 * (Param_eD ** 9) + -0.009174555773228737 * (Param_eD ** 10)
        if Param_eD > 3.5:
            K_bry = 1.5
    elif t/D > 0.15:
        K_bry = -0.015 + 2.6746380373328176 * (Param_eD ** 1) + -3.8095514487546716 * (Param_eD ** 2) + 8.623219742992589 * (Param_eD ** 3) + -15.750885241422282 * (Param_eD ** 4) + 17.115175048197273 * (Param_eD ** 5) + -11.18034715883286 * (Param_eD ** 6) + 4.4747996668137695 * (Param_eD ** 7) + -1.0765463102425024 * (Param_eD ** 8) + 0.14305215852155762 * (Param_eD ** 9) + -0.008072914212872394 * (Param_eD ** 10)
        if Param_eD > 3.5:
            K_bry = 1.4
    elif t/D > 0.12:
        K_bry = -0.015 + 2.26199460291689 * (Param_eD ** 1) + 0.2998846940719729 * (Param_eD ** 2) + -5.896253454827716 * (Param_eD ** 3) + 8.866390732323316 * (Param_eD ** 4) + -6.659250538483327 * (Param_eD ** 5) + 2.8922806130518577 * (Param_eD ** 6) + -0.7356872271946134 * (Param_eD ** 7) + 0.10187317354853766 * (Param_eD ** 8) + -0.005936465946266028 * (Param_eD ** 9)
        if Param_eD > 3.5:
            K_bry = 1.32
    elif t/D > 0.1:
        K_bry = -0.012 + 2.1958873981812173 * (Param_eD ** 1) + 1.4121860944569744 * (Param_eD ** 2) + -11.094326947388348 * (Param_eD ** 3) + 19.05242541029777 * (Param_eD ** 4) + -17.492127624958844 * (Param_eD ** 5) + 9.777977619561582 * (Param_eD ** 6) + -3.430942783487547 * (Param_eD ** 7) + 0.739575340442149 * (Param_eD ** 8) + -0.08966973195725639 * (Param_eD ** 9) + 0.004687489286406277 * (Param_eD ** 10)
        if Param_eD > 3.5:
            K_bry = 1.25
    elif t/D > 0.08:
        K_bry = -0.000 + 2.047168809200642 * (Param_eD ** 1) + 2.8439054359210623 * (Param_eD ** 2) + -18.38619717436508 * (Param_eD ** 3) + 34.28307332550166 * (Param_eD ** 4) + -34.299009336084836 * (Param_eD ** 5) + 20.676428605730308 * (Param_eD ** 6) + -7.734726022194671 * (Param_eD ** 7) + 1.7588160610321202 * (Param_eD ** 8) + -0.22283882512474207 * (Param_eD ** 9) + 0.012071823552433879 * (Param_eD ** 10)
        if Param_eD > 3.5:
            K_bry = 1.12
    elif t/D > 0.06:
        K_bry = -0.019 + 2.815705548826579 * (Param_eD ** 1) + -3.8386040984312917 * (Param_eD ** 2) + 0.5570512486788379 * (Param_eD ** 3) + 6.174683983615872 * (Param_eD ** 4) + -9.49021569548029 * (Param_eD ** 5) + 6.967071323446899 * (Param_eD ** 6) + -2.941962100212095 * (Param_eD ** 7) + 0.7284623056343401 * (Param_eD ** 8) + -0.09845184660891078 * (Param_eD ** 9) + 0.005614099708215215 * (Param_eD ** 10)
        if Param_eD > 3.5:
            K_bry = 1.0
    else:
        #print("t/D too small")
        #print(t / D)
        return -1
    if K_bry < 0:
        K_bry = 0

    return K_bry


def linearSolve(eq1, eq2, lower, upper, cycles, steps, parameters = []):
    upperAnswer = upper
    lowerAnswer = lower

    for y in range(0, cycles):

        diffSign = 0
        cycleDone = False

        for x in range(0, steps+1):
            if (not cycleDone):
                i = lowerAnswer + (upperAnswer - lowerAnswer) / steps * x
                # print(x)
                val1 = eq1(i, parameters)
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
                    #print("lower " + str(lowerAnswer))
                    #print("upper " +  str(upperAnswer))

    diff1 = eq2(lowerAnswer) - eq1(lowerAnswer, parameters)
    diff2 = eq2(upperAnswer) - eq1(upperAnswer, parameters)

    diff3 = upperAnswer - lowerAnswer

    return lowerAnswer + diff1 / ((diff1 - diff2) / diff3)

Al7075 = Material(71700000000, 503000000, 2810,331000000)

R = calcAxleRadius(Al7075)
#print(R)

#print(calcBendingStress(0.01,0.02,0.05))
config = iterateW(Al7075,0.02,2*R)
print("mass " + str(config.mass))
print("diameter " + str(config.D))
print("W " + str(config.W))
print("T " + str(config.t))
("L " + str(config.l))

#wArr = []
#tArr = []
#for x in range(1,100):
 #   t = x/10000
 #   diff = calcStressDifference(t,[0.0028,0.03,0.012,200000000])
#
 #   print(str(t) + " " + str(diff))
 #   wArr.append(t)
#    tArr.append(diff)
#plt.plot(wArr,tArr)
#plt.show()

print(linearSolve(calcStressDifference, zeroFunc, 0.002, 1, 10, 3, [0.00295,0.03,0.012,200000000]))

#print(calcFlangeStressAxial(0.001,2*testR,0.014,0.7))