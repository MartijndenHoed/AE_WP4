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
        self.sigmaYield

def calcAxleRadius(material):
    #calcualte the minium radius of the (circular) axle


    stressAllow = material.sigmaYield
    axleArea = Plat / stressAllow
    r = math.sqrt(axleArea/math.pi)

    #pain
    return r

def calcFlangeThicknessTension(material,Diameter,width):
    k = 0.7 ##arbitrary, find function for this (pain)

    t = diameter/(k*material.sigmaYield,width)

    return t

def calcFlangeThicknessBearing(material,Diameter,width):



##calculate the thickness for a certain width, take the minimum value, compare the mass for all materials.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
