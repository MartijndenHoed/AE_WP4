# Stuff that needs to be done:
# bearing shear out needs to be tested and integrated
# different materials need to be integrated - CARMELO



import math




##all SI units please
g = 9.80665 #m/(s**2)


RTGmass = 52.8 #kg
RTGlength = 1.07 #m
RTGdiameter = 0.47 #m

lugSpacingX = 0.4 #m
lugSpacingZ = RTGlength #m

GStaticlong = 4.55 #g
GStaticLat = 0.25 #g

GDynamicLong = 1.0 #g
GDynamicLat = 0.8 #g

SafetyFactorLaunch = 1.25*1.5

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
    def __init__(self,t,W,D,mass,l):
        self.t = t
        self.W = W
        self.D = D
        self.mass = mass
        self.l = l

def iterateConfigs(material):
    bestMass = 1000
    config = -1
    Diameter = 2 * calcAxleRadius(material)



    for y in range(0,18):
        global l
        l = (y/100) + 0.03
        print("checking L: " + str(l))
        RTGCGy = 0.5 * RTGdiameter
        RTGMomentArm = l - RTGCGy
        moment = LugFzMax * RTGMomentArm * 0.5
        global LugFyMax
        LugFyMax = TotalFy_max / 4 + abs((moment / RTGlength) * 0.5 * 0.5)
        #print(LugFyMax)
        for x in range(0,30):
            W = x/1000 + Diameter * 2
            #print("checking W: " + str(W))
            tempConfig = generateLugConfiguration(material,W,l,Diameter)
            if(tempConfig != -1):
                tempMass = tempConfig.mass
                if(bestMass > tempMass):
                    bestMass = tempMass
                    config = tempConfig

    return config

def generateLugConfiguration(lugMat,W,l,Diameter):
    #this should be enough variables to generate a thickness and calculate a weight

    if(Diameter >= W):
        return -1

    t = calcMinimumFlangeThickness(lugMat, Diameter, W,l)
    A = (W * l + 0.5 * math.pi * (0.5*W)**2 + math.pi * (Diameter/2) ** 2)
    V = A * t
    #print(str(t) + " " + str(W) + " " + str(l) + " " + str(W*t*l))
    singleMass = V * lugMat.density
    totalMass = 2 * singleMass

    config = Config(t,W,Diameter,totalMass,l)

    return config

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
    tBearing = calcFlangeThicknessBearingShearOut(material,Diameter,width)
    #print("t " + str(tBearing))


    minT = max(smallestTTension,tBendingXZ,tBearing)


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
    if(ratio > 1.4):
        kt = 1.29
    else:
        kt = Kt_val_transverse(ratio)
    value = globalFy/(kt)

    return value

def transverseEqRight(t):


    return t * globalD * globalSigma

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


def calcFlangeThicknessShearXY(material,width):
    t = max(LugFyMax/(2*material.sigmaYieldShear*width), LugFxMax/(2*material.sigmaYieldShear*width))
    return t


def calcFlangeThicknessBearingShearOut(material,Diameter,width):
    # calculate the flange thickness for bearing shear out
    P_bry = 0.
    t_step = 0.001
    t=0
    while P_bry < LugFyMax:
        t = t + t_step
        #print("t: " + str(t))
        #print("Ed: " + str( width/(2*Diameter)))
        K_bry = Kbrycalc(Diameter, t, width)
        A_br = Diameter * t
        F_ty = material.sigmaYield
        P_bry = F_ty * A_br * K_bry * 0.5
        if(K_bry == -1):
            P_bry = 0
        #print("k: " + str(K_bry))
        #print("load: " + str(P_bry))
        #print("load2: " + str(LugFyMax))
    #print("Ed: " + str( width/(2*Diameter)))
    #print(t)
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
    #a bit more efficient and combines the above two functions into one
    M1 = LugFzMax * l / 2
    M2 = LugFxMax * l / 2
    sigma = material.sigmaYield
    return  (-((6*M1)/(w**2)) - math.sqrt( ((6*M1)/(w**2))**2 + 4 * sigma * ((6*M2)/(w)) ))/(-2*sigma)


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


def numericalExtremeValue(eq1, lower, upper, cycles, steps, h=0.00001, parameters = []):
    upperAnswer = upper
    lowerAnswer = lower
    rootFound = False

    for y in range(0, cycles):

        derSign = 0
        cycleDone = False

        for x in range(0, steps + 1):
            if (not cycleDone):
                i = lowerAnswer + (upperAnswer - lowerAnswer) / steps * x
                # print(x)
                val1 = eq1(i + h,parameters)
                val2 = eq1(i,parameters)
                der = (val1 - val2) / h
                # print(der)
                if (der == 0):
                    # print("exact: " + str(i))
                    return i
                sign = der / abs(der)
                if (derSign == 0):
                    derSign = sign
                # print(str(derSign) + " " + str(sign) + " " + str(i))
                if (derSign != sign):
                    cycleDone = True
                    rootFound = True
                    lowerAnswer = lowerAnswer + (upperAnswer - lowerAnswer) / steps * (x - 1)
                    upperAnswer = i
                    # print("lower " + str(lowerAnswer))
                    # print("upper " +  str(upperAnswer))

    if (not rootFound):
        print("extreme value not found")
        return "null"

    der1 = -(eq1(lowerAnswer + h,parameters) - eq1(lowerAnswer,parameters))
    der2 = -(eq1(upperAnswer + h,parameters) - eq1(upperAnswer,parameters))

    der3 = upperAnswer - lowerAnswer

    return lowerAnswer + der1 / ((der1 - der2) / der3)

def Kty_val_transverse(Ratio1):
    # if material == Material_curve1:
    #     Ktu = 0.9180553870465704 * (Ratio1 ** 1) + 8.18122273185643 * (Ratio1 ** 2) + -57.7056377483359 * (Ratio1 ** 3) + 231.21388530032073 * (Ratio1 ** 4) + -578.4496388173648 * (Ratio1 ** 5) + 921.2687606248473 * (Ratio1 ** 6) + -928.1664126745236 * (Ratio1 ** 7) + 570.5859381665941 * (Ratio1 ** 8) + -194.87683183789522 * (Ratio1 ** 9) + 28.300734922855554 * (Ratio1 ** 10)
    # elif material == Material_curve2:
    #     Ktu = 1.1540799943441162 * (Ratio1 ** 1) + 3.5658936323271035 * (Ratio1 ** 2) + -13.955444513318298 * (Ratio1 ** 3) + 16.23858184466644 * (Ratio1 ** 4) + 19.999288298125066 * (Ratio1 ** 5) + -75.33757395030716 * (Ratio1 ** 6) + 82.87193625333839 * (Ratio1 ** 7) + -41.04677495808528 * (Ratio1 ** 8) + 7.779874724243314 * (Ratio1 ** 9)
    # elif material == Material_curve4:
    #     Ktu = 1.0968401953214486 * (Ratio1 ** 1) + 5.911055078229481 * (Ratio1 ** 2) + -41.276480666988846 * (Ratio1 ** 3) + 170.3631152105354 * (Ratio1 ** 4) + -485.8393789099573 * (Ratio1 ** 5) + 946.0509482390116 * (Ratio1 ** 6) + -1188.580514624009 * (Ratio1 ** 7) + 901.8312818112099 * (Ratio1 ** 8) + -372.68816788151344 * (Ratio1 ** 9) + 64.20378959440313 * (Ratio1 ** 10)
    # elif material == Material_curve5:
    #     Ktu = 0.8935998289127944 * (Ratio1 ** 1) + 5.273317497891158 * (Ratio1 ** 2) + -28.55089032590277 * (Ratio1 ** 3) + 52.05627723859402 * (Ratio1 ** 4) + 28.107639294058117 * (Ratio1 ** 5) + -256.71084383513744 * (Ratio1 ** 6) + 425.1831418346137 * (Ratio1 ** 7) + -339.90178976325706 * (Ratio1 ** 8) + 136.91477026647703 * (Ratio1 ** 9) + -22.29608887194499 * (Ratio1 ** 10)
    # elif material == Material_curve6:
    #     Ktu = 0.9564655594429534 * (Ratio1 ** 1) + 5.166345227885259 * (Ratio1 ** 2) + -34.01201870580809 * (Ratio1 ** 3) + 104.20957326598887 * (Ratio1 ** 4) + -180.16817428435706 * (Ratio1 ** 5) + 184.97898791984153 * (Ratio1 ** 6) + -112.11763123206947 * (Ratio1 ** 7) + 37.11718579286003 * (Ratio1 ** 8) + -5.182081049544365 * (Ratio1 ** 9)
    # elif material == Material_curve7:
    #     Ktu = 2.039426957816968 * (Ratio1 ** 1) + -13.863718880066074 * (Ratio1 ** 2) + 96.71439084070298 * (Ratio1 ** 3) + -386.6676379116315 * (Ratio1 ** 4) + 925.0396689887041 * (Ratio1 ** 5) + -1370.8650598029553 * (Ratio1 ** 6) + 1268.331094686256 * (Ratio1 ** 7) + -712.9116263332044 * (Ratio1 ** 8) + 222.7186655351037 * (Ratio1 ** 9) + -29.667334586692217 * (Ratio1 ** 10)
    # elif material == Material_curve8:

    Kty = 1.2545232944973213 * (Ratio1 ** 1) + 2.2613250246083663 * (Ratio1 ** 2) + -20.705248912030186 * (
                Ratio1 ** 3) + 83.2739255326694 * (Ratio1 ** 4) + -182.25919008412453 * (
                      Ratio1 ** 5) + 228.18661744033676 * (Ratio1 ** 6) + -163.82182854204245 * (
                      Ratio1 ** 7) + 62.861464155581686 * (Ratio1 ** 8) + -10.003844555636892 * (Ratio1 ** 9)

    return Kty

def Kt_val_axial(R):
    Kt = 1.009 + -0.22902304135292462 * (R ** 1) + 0.35905460079484186 * (R ** 2) + -1.0729145294629987 * (R ** 3) + 1.2672573367948203 * (R ** 4) + -0.8003541023438316 * (R ** 5) + 0.3032426494074585 * (R ** 6) + -0.07136553672277213 * (R ** 7) + 0.010232930693775955 * (R ** 8) + -0.000818435700391848 * (R ** 9) + 0.000027954852209106408 * (R ** 10)
    if(R>5):
        return 0.001
    return Kt

def Kt_val_transverse(r):
    Kt = 0.003 + 1.1711822557878742 * (r ** 1) + 2.058639823919016 * (r ** 2) + -17.838056169717532 * (r ** 3) + 62.31283878518559 * (r ** 4) + -118.14607854029032 * (r ** 5) + 129.78187628214891 * (r ** 6) + -82.82263597513247 * (r ** 7) + 28.55280116408801 * (r ** 8) + -4.118577606465635 * (r ** 9)
    return Kt

##Materials: (material(self, E, sigma yield, density))
#materialsList = []
Al7075 = Material(71700000000, 503000000, 2810,331000000)
#Al7050 = Material(71700000000, 490000000, 2830)
#Al2024 = Material(72400000000, 345000000, 2780)
#Al2219 = Material(73100000000, 393000000, 2840)
#Al2014 = Material(72400000000, 414000000, 2800)



#aluminium = Material(70000000000,200000000,2700) ##just some test values
#testR = calcAxleRadius(aluminium)
#print(testR)
#testAxialT = calcFlangeThicknessTensionAxial(aluminium,testR,0.05,1) #(material,Diameter,width,Rax)
#print(testAxialT)
#config = generateLugConfiguration(aluminium,0.03,0.1,0.0017)
config = iterateConfigs(Al7075)
print("mass " + str(config.mass))
print("diameter " + str(config.D))
print("W " + str(config.W))
print("T " + str(config.t))
print("L " + str(config.l))




