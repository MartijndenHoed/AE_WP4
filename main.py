# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

#bruhhh

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

PLong = GLongMax * RTGmass
Plat = GlatMax * RTGmass

class material:
    def __init__(self,E,sigmaYield):



def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
