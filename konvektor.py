# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 18:57:46 2022

@author: Jindřich Novák
"""
#print("Script started.   Importing libraries.")
import math
import pandas
import matplotlib.pyplot as plot
import matplotlib.ticker as ticker
import numpy as np
print("Libraries imported.   Importing data.")


filename = "./data/RAW_11520_2023-07-06_1200.csv" #Název datového souboru
ktstoms = 0.5144444 #Převodní poměr mezi uzly a metry za sekundu
g = 9.80665 #Tíhové zrychlení zemské
k = 0.052 #Součinitel vnitřního tření konvektivní částice
n = 0.190284 #Koeficient pro výpočet QNH
Rv = 461.5 #Měrná plynová konstanta pro vodní páru
Rd = 287.04 #Měrná plynová konstanta pro suchý vzduch
R = 8.31446 #Molární plynová konstanta
Mv = 0.018016 #Molární hmotnost vodní páry
Md = 0.0289652 #Molární hmotnost suchého vzduchu
eps = Rd/Rv #Poměr měrných plynových konstant
L = 2501000 #Měrné skupenské teplo vypařování vody za normálních podmínek
cp = 1.006 #Izobarická měrná tepelná kapacita vzduchu
cv = 0.718 #Izochorická měrná tepelná kapacita vzduchu
opt = 3.5 #Optimalizační konstanta pro výpočet rychlosti stoupání bržděné částice
kappa = cp/cv
gammaD = g/cp #Sucho-adiabatický vertikální teplotní gradient
stdRho = 1.225 #Hustota vzduchu na hladině moře dle mezinárodní standardní atmosféry
stabDev = -0.4 #Při kladných hodnotách bude zvrstvení hodnoceno jako více stabilní, při záporných jako méně stabilní
stabDif = 5 #Výškový krok mezi srovnávanými hladinami při určování stability teplotního zvrstvení
showNEL = False #Zobrazovat NEL?
simUntilCCL = True #Při hodnotě True se simulace automaticky zastaví při dostoupání částice do kondenzační hladiny
CCLlimit = 400 #Nejnižší povolená výska konvektivní kondenzační hladiny

plotTmin = -25   #Nejnižší teplota emagramu
plotTmax = "auto"   #Nejvyšší teplota emagramu. Nastavením na "auto" se automaticky nastaví nejvyšší teplota datové řady zvýšená o 1°C.
plotAmin = "auto"   #Nejnižší výška emagramu. Nastavením na "auto" se automaticky nastaví na hladinu zemského povrchu.
plotAmax = 5000   #Nejvyšší výška emagramu
dpi = 60 #Rozlišení grafů v DPI
Nbarbs = 30   #Vykreslovat šipku větru pro každý Nbarbs-tý datový bod


dt = 0.01   #Časový krok pro simulace
Simulated_seconds = 1200   #Simulovaná doba v sekundách
v0 = 0   #Počáteční vertikální rychlost konvektivní částice
virtv0 = 0   #Počáteční vertikální rychlost konvektivní částice při výpočtech s virtuální teplotou
cycles = int(Simulated_seconds/dt)   #Počet iterací simulace
speedInterval = 0.2   #Rychlostní krok, o který se zvyšuje počáteční vertikální rychlost konvektivní částice
itInt = 50000   #Interval, po kolika iteracích se printuje oznámení

if input("Use default preferences? (Y/N): ") == 'N':
    showNEL = bool(input("Draw New Equilibrium Level? (True/False): "))
    simUntilCCL = bool(input("Stop simulation when CCL is reached? (True/False): "))
    CCLlimit = float(input("Minimum CCL altitude AMSL: (float)"))
    plotTmin = float(input("Minimum emagram temperature: (float)"))
    plotTmax = float(input("Maximum emagram temperature: (float)"))
    plotAmin = float(input("Minimum emagram altitude: (float)"))
    plotAmax = float(input("Maximum emagram altitude: (float)"))
    dpi = int(input("Resolution: (int)"))
    dt = int(input("Simulation time step: (float)"))
    Simulated_seconds = int(input("Simulated time period: (int)"))

data = pandas.read_csv(filename, header=0) #Přečtení datového souboru

#Převod sloupců tabulky do jednotlivých seznamů:
Plist = list(data.Pressure) #Seznam tlaků
Alist = list(data.Altitude) #Seznam výšek
Tlist = list(data.Temperature) #Seznam teplot
Hlist = list(data.Dew_point) #Seznam rosných bodů
Dlist = list(data.Wind_direction) #Seznam směrů větru
Slist = list(data.Wind_speed) #Seznam rychlostí větru

"""
Plist = [993.0,990.0,925.0,924.0,898.0,855.0,850.0,843.0,828.0,827.0,775.0,770.0,700.0,655.0,620.0,613.0,602.0,559.0,550.0,512.0,500.0,497.0,492.0,489.0,480.0,475.0,435.0,412.0,400.0,391.0,388.0,383.0,349.0,348.0,339.0,323.0,311.0,300.0,290.0,281.0,279.0,250.0,217.0,200.0,191.0,189.0,167.0,150.0,148.0,136.0,122.0,118.0,111.0,107.0,105.0,100.0,98.0,92.0,86.0,79.0,77.0,75.0,70.8,70.0,67.0,63.0,60.0,59.0,57.0,52.0,50.0,48.3,47.0,45.0,43.0,40.0,39.0,38.0,37.0,35.0,33.0,31.0,30.0,26.3,26.0,25.0,23.0,22.0,21.0,20.0,19.6,19.0,18.1,18.0,17.0,16.0,15.6,14.0,13.0,12.2,11.0,10.0,9.0,8.0,7.0,6.3,6.2] #Seznam tlaků
Alist = [304,329,894,903,1138,1538,1586,1653,1799,1809,2338,2390,3160,3686,4121,4211,4352,4929,5053,5599,5780,5825,5901,5946,6084,6162,6808,7198,7410,7572,7626,7717,8367,8386,8565,8895,9149,9390,9615,9824,9871,10590,11493,12010,12304,12371,13158,13840,13924,14452,15131,15338,15718,15948,16065,16370,16495,16887,17305,17831,17990,18153,18510,18580,18852,19234,19537,19641,19856,20426,20670,20884,21052,21321,21603,22050,22207,22368,22532,22876,23240,23627,23830,24648,24720,24967,25491,25771,26063,26370,26498,26695,27003,27038,27404,27792,27954,28649,29125,29533,30203,30820,31518,32299,33195,33902,34010] #Seznam výšek
Tlist = [14.0,13.7,7.2,7.0,5.6,3.6,3.4,3.2,4.2,4.2,5.2,4.9,0.2,-3.0,-5.6,-6.1,-6.7,-9.3,10.3,-15.0,-16.5,-16.9,-17.0,-17.1,-18.0,-18.5,-24.3,-27.5,-29.3,-30.9,-31.4,-32.2,-38.1,-38.5,-39.8,-42.3,-44.4,-46.3,-47.6,-48.8,-49.1,-52.1,56.9,55.1,-54.8,-54.7,-56.1,-57.3,-57.6,-59.5,-61.9,-60.9,-59.1,-59.5,-59.7,-60.3,-60.5,-60.9,-61.4,-62.1,-62.3,-62.5,-62.9,-62.7,-61.9,-60.9,-60.0,-59.7,-60.2,-61.5,-62.1,-62.9,-62.8,-62.7,-62.5,-62.3,-62.2,-62.1,-62.0,-61.8,-61.6,-61.4,-61.3,-61.5,-61.3,-60.6,-59.0,-58.2,-57.4,-56.5,-56.3,-56.6,-57.1,-56.9,-55.2,-53.3,-52.5,-53.6,-54.4,-55.1,-52.1,-49.3,-46.9,-44.3,-43.8,-43.4,-43.3] #Seznam teplot
Hlist = list(data.Dew_point) #Seznam rosných bodů
Dlist = list(data.Wind_direction) #Seznam směrů větru
Slist = list(data.Wind_speed) #Seznam rychlostí větru
"""

print("Data imported.   Calculating.")

######################################################################
#Výpočty##############################################################
######################################################################

#Převod seznamu rychlostí větru z uzlů na metry za sekundu
for i in range(len(Slist)):
    Slist[int(i)] = Slist[int(i)]*ktstoms

#Přepočet úhlu směru větru na směr vektoru větru v radiánech
radDlist = [] #Seznam směrů rychlostního vektoru větru v radiánech
for j5 in range(len(Dlist)):
    rads = -1*math.radians(180-Dlist[j5])
    radDlist.append(rads)

#Rozdělení seznamů rychlostí a směrů dle výšky
S1list = [] #Seznam rychlostí větru mezi zemským povrchem a hladinou 1 km
S2list = [] #Seznam rychlostí větru mezi hladinami 1 km a 2 km
S3list = [] #Seznam rychlostí větru mezi hladinami 2 km a 3 km
S6list = [] #Seznam rychlostí větru mezi hladinami 3 km a 6 km
S9list = [] #Seznam rychlostí větru mezi hladinami 6 km a 9 km
S12list = [] #Seznam rychlostí větru mezi hladinami 9 km a 12 km
S20list = [] #Seznam rychlostí větru nad hladinou 12 km
D1list = [] #Seznam směrů větru mezi zemským povrchem a hladinou 1 km
D2list = [] #Seznam směrů větru mezi hladinami 1 km a 2 km
D3list = [] #Seznam směrů větru mezi hladinami 2 km a 3 km
D6list = [] #Seznam směrů větru mezi hladinami 3 km a 6 km
D9list = [] #Seznam směrů větru mezi hladinami 6 km a 9 km
D12list = [] #Seznam směrů větru mezi hladinami 9 km a 12 km
D20list = [] #Seznam směrů větru nad hladinou 12 km
for j6 in range(len(Slist)):
    if Alist[j6]<=1000:
        S1list.append(Slist[j6])
        D1list.append(radDlist[j6])
    elif Alist[j6]<=2000:
        S2list.append(Slist[j6])
        D2list.append(radDlist[j6])
    elif Alist[j6]<=3000:
        S3list.append(Slist[j6])
        D3list.append(radDlist[j6])
    elif Alist[j6]<=6000:
        S6list.append(Slist[j6])
        D6list.append(radDlist[j6])
    elif Alist[j6]<=9000:
        S9list.append(Slist[j6])
        D9list.append(radDlist[j6])
    elif Alist[j6]<=12000:
        S12list.append(Slist[j6])
        D12list.append(radDlist[j6])
    else:
        S20list.append(Slist[j6])
        D20list.append(radDlist[j6])

#Výpočet relativní vlhkosti, parciálních tlaků, směšovacího poměru a virtuální teploty

virtTlistK = [] #Seznam virtuálních teplot v Kelvinech
virtTlistC = [] #Seznam virtuálních teplot ve stupních Celsia
VTdiff = [] #Seznam rozdílů suché a virtuální teploty
Rlist = [] #Seznam relativních vlhkostí
elist = [] #Seznam parciálních tlaků nasycení vzduchu vodní parou
Elist = [] #Seznam maximálních parciálních tlaků nasycení vzduchu vodní parou
Wlist = [] #Seznam směšovacích poměrů

for j in range(len(Tlist)):
    e = 0.61121*math.exp((18.678-(Hlist[j]/234.5))*(Hlist[j]/(257.14+Hlist[j])))
    E = 0.61121*math.exp((18.678-(Tlist[j]/234.5))*(Tlist[j]/(257.14+Tlist[j])))
    r = e/E
    Rlist.append(r)
    elist.append(e)
    Elist.append(E)
    w = (r*eps*E)/100*(Plist[int(j)]-E)/1000
    Wlist.append(w)
    s = w/(w+1)
    virtTlistK.append((Tlist[int(j)]+273.15)*(1+((Rv/Rd)-1)*s))
    VTdiff.append(virtTlistK[int(j)]-Tlist[int(j)]-273.15)
    if VTdiff[int(j)]<0: #Virtuální teplota nesmí být nižší než suchá teplota
        print("Calculation ERROR: Virtual temperature less than dry bulb temperature!")
    if Rlist[int(j)]>1:
        print("Calculation ERROR: Relative humidity greater than 100%!")
for j2 in range(len(virtTlistK)):
    virtTlistC.append(virtTlistK[j2] - 273.15)


#Výpočet redukovaných tlaků

QNH = Plist[0]*(1+((1013.25**n)*0.0065*0.003472)*Alist[0]/(Plist[0]**n))**(1/n)
QFF = Plist[0]*math.exp(g*Alist[0]/(287.04*virtTlistK[0]))


#Vytvoření seznamu relativních vlhkostí v procentech

R100list = []
for j3 in range(len(Rlist)):
    R100list.append(Rlist[j3]*100)


#Výpočet suché adiabaty a izogramy

WTlist = [] #Seznam hodnot teploty izogramy konvektivní částice
TGDlist = [] #Seznam hodnot teploty suché adiabaty pro částici stoupající od zemského povrchu
for i2 in range(len(Alist)):
    Tgamma = Tlist[0]+(Alist[int(i2)] - Alist[0])/1000*(-gammaD)
    e = 100*Plist[int(i2)]/((eps/Wlist[0]*1000)+1)
    wTd = ((math.log(e/0.61094)/-(L/Rv))+(1/273.15))**(-1)
    WTlist.append(wTd-273.15)
    TGDlist.append(Tgamma)


#Oprava teploty izogramy

diff = Hlist[0] - WTlist[0]
for i4 in range(len(WTlist)):
    WTlist[i4] += diff


#Výpočet VKH (výstupné kondenzační hladiny)

for i3 in range(len(WTlist)):
    if TGDlist[i3]<=WTlist[i3]:
        LCLtemp = TGDlist[i3]
        LCLalt = Alist[i3]
        LCLi = i3
        break


#Výpočet KKH (konvektivní kondenzační hladiny)

for i6 in range(len(WTlist)):
    if Tlist[i6]<WTlist[i6] and Alist[i6]>CCLlimit:
        CCLtemp = Tlist[i6]
        CCLalt = Alist[i6]
        CCLi = i6
        break
    

#Výpočet konvektivní teploty a deficitu konvektivní teploty

konvTGDlist = []
for i7 in range(len(Alist)):
    Konvgamma = Tlist[CCLi]+(Alist[int(i7)] - Alist[CCLi])/1000*(-gammaD)
    konvTGDlist.append(Konvgamma)
Tkonv = konvTGDlist[0]
dTkonv = Tlist[0] - konvTGDlist[0]


#Výpočet nasyceně-adiabatických vertikálních teplotních gradientů a nasycených adiabat

#Z konvektivní kondenzační hladiny
konvTGSlistD = [CCLtemp]
for a2 in range(CCLi):
    Egs = 0.61094*math.exp((L/Rv)*(1/273.15-(1/(konvTGSlistD[a2]+273.15))))
    Wgs = (eps*Egs)/(Plist[CCLi-a2]-Egs)
    cit = 1 + ((L*10*Wgs)/(Rd*(konvTGSlistD[a2]+273.15)))
    jme = 1000*cp + ( (10*Wgs/Rv)  *  (((L)**2)/((konvTGSlistD[a2]+273.15)**2)   )    )
    salr = g*(cit/jme)
    tgs = konvTGSlistD[a2] + salr*(Alist[CCLi-a2+1]-Alist[CCLi-a2])
    konvTGSlistD.append(tgs)

konvTGSlist = []
for a3 in range(len(konvTGSlistD)):
    konvTGSlist.append(konvTGSlistD[-a3])


konvTGSlistU = [CCLtemp]
for a1 in range(CCLi, len(Tlist)):
    if konvTGSlist[a1]>-270 and a1+1 != len(Tlist):
        Egs = 0.61094*math.exp((L/Rv)*(1/273.15-(1/(konvTGSlist[a1]+273.15))))
        Wgs = (eps*Egs)/(Plist[a1]-Egs)
        cit = 1 + ((L*10*Wgs)/(Rd*(konvTGSlist[a1]+273.15)))
        jme = 1000*cp + ( (10*Wgs/Rv)  *  (((L)**2)/((konvTGSlist[a1]+273.15)**2)   )    )
        salr = g*(cit/jme)
        tgs = konvTGSlist[a1] - salr*(Alist[a1+1]-Alist[a1])
        konvTGSlist.append(tgs)
    else:
        konvTGSlist.append(-273.15)
for a4 in range(len(konvTGSlistU)):
    konvTGSlist.append(konvTGSlistU[a4])
konvTGSlist.pop(-1)
konvTGSlist.pop(-2)
konvTGSlist[0] = konvTGSlist[1]


#Z výstupné kondenzační hladiny
TGSlistD = [LCLtemp]
for a1 in range(LCLi):
    Egs = 0.61094*math.exp((L/Rv)*(1/273.15-(1/(TGSlistD[a1]+273.15))))
    Wgs = (eps*Egs)/(Plist[LCLi-a1]-Egs)
    cit = 1 + ((L*10*Wgs)/(Rd*(TGSlistD[a1]+273.15)))
    jme = 1000*cp + ( (10*Wgs/Rv)  *  (((L)**2)/((TGSlistD[a1]+273.15)**2)   )    )
    salr = g*(cit/jme)
    tgs = TGSlistD[a1] + salr*(Alist[LCLi-a1+1]-Alist[LCLi-a1])
    TGSlistD.append(tgs)

TGSlist = []
for a3 in range(len(TGSlistD)):
    TGSlist.append(TGSlistD[-a3])
if len(TGSlist)>1:
    TGSlist[0] = TGSlist[1]


TGSlistU = [LCLtemp]
for a2 in range(LCLi, len(Tlist)):
    if TGSlist[a2]>-270 and a2+1 != len(Tlist):
        Egs = 0.61094*math.exp((L/Rv)*(1/273.15-(1/(TGSlist[a2]+273.15))))
        Wgs = (eps*Egs)/(Plist[a2]-Egs)
        cit = 1 + ((L*10*Wgs)/(Rd*(TGSlist[a2]+273.15)))
        jme = 1000*cp + ( (10*Wgs/Rv)  *  (((L)**2)/((TGSlist[a2]+273.15)**2)   )    )
        salr = g*(cit/jme)
        tgs = TGSlist[a2] - salr*(Alist[a2+1]-Alist[a2])
        TGSlist.append(tgs)
    else:
        TGSlist.append(-273.15)
for a4 in range(len(TGSlistU)):
    TGSlist.append(TGSlistU[a4])
TGSlist.pop(-1)
TGSlist.pop(-2)


#Výpočet izobarické ekvivalentní teploty

ekvTlist = []
for i9 in range(len(Tlist)):
    Tekv = Tlist[i9] + (L*Wlist[i9]/(1000*cp))
    ekvTlist.append(Tekv)


#Výpočet potenciální teploty

thetaTlist = []
for i8 in range(len(Plist)):
    theta = (Tlist[i8] + 273.15)*(1000/Plist[i8])**(Rd/(cp*1000))
    thetaTlist.append(theta-273.15)


#Výpočet izobarické ekvivalentní potenciální teploty

thetaElist = []
for i10 in range(len(Plist)):
    theta = (ekvTlist[i10] + 273.15)*(1000/Plist[i10])**(Rd/(cp*1000))
    thetaElist.append(theta-273.15)


#Výpočet zrychlení (vztlaku působícího na částici)

Blist = []
for j7 in range(len(Tlist)):
    if j7 <= CCLi:
        b = g*( (konvTGDlist[j7] - Tlist[j7]) /  virtTlistK[j7]  )
    else:
        b = g*( (konvTGSlist[j7] - Tlist[j7]) /  virtTlistK[j7]  )
    Blist.append(b)


#Výpočet zrychlení (vztlaku působícího na částici) podle virtuální teploty

virtBlist = []
for j7 in range(len(Tlist)):
    if j7 <= CCLi:
        virtb = g*( (konvTGDlist[j7] - virtTlistC[j7]) /  virtTlistK[j7]  )
    else:
        virtb = g*( (konvTGSlist[j7] - virtTlistC[j7]) /  virtTlistK[j7]  )
    virtBlist.append(virtb)


#Výpočet HNV (hladiny nulového vztlaku)

HNVi = len(Blist)-1
HNValt = Alist[-1]
for b1 in range(1, len(Blist)):
    if Blist[len(Blist)-b1]>0:
        HNValt = Alist[len(Blist)-b1+1]
        HNVi = len(Blist)-b1+1
        break


#Výpočet HVK (hladiny volné konvekce)

HVKi = 0
for b2 in range(CCLi, len(Blist)):
    if Blist[b2]>0:
        HVKalt = Alist[b2]
        HVKi = b2
        break


#Výpočet CAPE

cape = 0
for b3 in range(HVKi, HNVi):
    cape += virtBlist[b3]*(Alist[b3+1]-Alist[b3])
    cape += (virtBlist[b3+1]-virtBlist[b3])*(Alist[b3+1]-Alist[b3])/2


#Výpočet NCAPE

ncape = cape/(HNValt - HVKalt)


#Výpočet CIN

cin = 0
for b4 in range(0, HVKi):
    cin += (-1)*virtBlist[b4]*(Alist[b4+1]-Alist[b4])
    cin += (-1)*(Alist[b4+1]-Alist[b4])*(virtBlist[b4+1]-virtBlist[b4])


#Výpočet TOP (hladiny vrcholu oblaku)

tp = 0
TOPi = len(Blist)-1
TOPalt = Alist[-1]
for b5 in range(HNVi, len(Alist)):
    if tp<cape:
        tp += -1*virtBlist[b5]*(Alist[b5]-Alist[b5-1])
        tp += -1*(Alist[b5]-Alist[b5-1])*(virtBlist[b5+1]-virtBlist[b5])
    else:
        TOPalt = Alist[b5]
        TOPi = b5
        break


#Výpočet nové rovnovážné hladiny po klesání z hladiny top

nelT = konvTGSlist[TOPi]
for b6 in range(0, len(Tlist)-TOPi):
    if Tlist[TOPi-b6]<=nelT:
        NELi = TOPi-b6
        NELt = Tlist[NELi]
        NELalt = Alist[NELi]
        break
    else:
        nelT = konvTGSlist[TOPi] + gammaD*(Alist[TOPi]-Alist[TOPi-b6])/1000


#Výpočet DCAPE

dcape = 0
for b3 in range(HNVi, TOPi):
    dcape += Blist[b3]*(Alist[b3+1]-Alist[b3])


#Výpočet Bruntovy-Vaisalovy frekvence

BruntVaisalaFreqlist = []
BruntVaisalaFreqSqlist = []
for f1 in range(len(Tlist)):
    try:
        dtheta = (thetaTlist[f1]-thetaTlist[f1-1])/(Alist[f1]-Alist[f1-1])
        bvf = (g/thetaTlist[f1])*dtheta
        if bvf >= 0:
            BruntVaisalaFreqlist.append(bvf**0.5)
            BruntVaisalaFreqSqlist.append(bvf)
        else:
            BruntVaisalaFreqlist.append(-1)
            BruntVaisalaFreqSqlist.append(bvf)
    except:
        BruntVaisalaFreqlist.append(-1)
        BruntVaisalaFreqSqlist.append(0)
        print("ERROR: Division by zero when calculating Brunt-Vaisala frequency.")


#Výpočet Bruntovy-Vaisalovy periody

BruntVaisalaPlist = []
for i2 in range(len(BruntVaisalaFreqlist)):
    if BruntVaisalaFreqlist[i2] != 0:
        BruntVaisalaPlist.append(1/BruntVaisalaFreqlist[i2])
    else:
        BruntVaisalaPlist.append(-1)


#Výpočet hustoty vzduchu

Rholist = []
for i3 in range(len(Tlist)):
    rho = ((Plist[i3]*100-elist[i3]*100)*Md + 100*elist[i3]*Mv)/(R*(Tlist[i3]+273.15))
    Rholist.append(rho)


#Výpočet výšky první konvenční tropopauzy

tropo1i = 0
tropo1alt = 0
tropo1temp = 0
tropo = True
for i4 in range(1, len(Tlist)):
    grad = (Tlist[i4-10]-Tlist[i4])/((Alist[i4]-Alist[i4-10])/1000)
    if grad < 2 and Plist[i4]<500:
        for i5 in range(i4, len(Alist)):
            if (Alist[i5]-Alist[i4])>2000:
                A1 = i5
                break
        Agrad = (Tlist[i4]-Tlist[A1])/((Alist[A1]-Alist[i4])/1000)
        if Agrad < 2:
            tropo1i = i4
            tropo1alt = Alist[i4]
            tropo1temp = Tlist[i4]
            break


#Výpočet výšky druhé konvenční tropopauzy

#Zjišťování, zda li druhá konvenční tropopauza vůbec existuje
for i1 in range(tropo1i, len(Tlist)):
    try:
        grad = (Tlist[i1]-Tlist[i1-1])/((Alist[i1]-Alist[i1-1])/1000)
    except:
        print("ERROR: Division by zero when calculating second tropopause.")
    if grad > 3:
        for i2 in range(i1, len(Alist)):
            if (Alist[i2]-Alist[i1])>1000:
                A2 = i2
                break
        Agrad = (Tlist[i1]-Tlist[A2])/((Alist[A2]-Alist[i1])/1000)
        if Agrad > 3:
            mintropo2i = A2
            tropo2 = True
            break
        
#Samotný výpočet, pokud bylo zjištěno, že existuje     
if tropo2 is True:   
    for i3 in range(mintropo2i, len(Tlist)):
        grad = (Tlist[i3]-Tlist[i3-1])/((Alist[i3]-Alist[i3-1])/1000)
        if grad < 2:
            for i5 in range(i3, len(Alist)):
                if (Alist[i5]-Alist[i3])>2000:
                    A1 = i5
                    break
            Agrad = (Tlist[i3]-Tlist[A1])/((Alist[A1]-Alist[i3])/1000)
            if Agrad < 2:
                tropo2i = i3
                tropo2alt = Alist[i3]
                tropo2temp = Tlist[i3]
                break


#Výpočet hořejší hodnoty tlakové škály

Pmax = -1
for i5 in range(len(Alist)):
    if Alist[i5] >= plotAmax:
        Pmax = i5
        break
if Pmax == -1:
    Pmax = len(Plist)-1


#Automatické nastavení krajních hodnot os diagramu

if plotAmin == "auto":
    plotAmin = Alist[0]
if plotTmax == "auto":
    plotTmax = max([max(konvTGDlist), max(virtTlistC), max(Tlist)])+1


print("Calculations completed.   Simulating parcel motion.")

######################################################################
#Simulace#############################################################
######################################################################

t = 0
Ai = 0
simv = v0
sima = Alist[0]
simVlist = []
simAlist = [-1]
simtlist = []
simtop = False
simmaxA = 0
simCCLi = 0
simx = 0
simy = 0
simXlist = []
simYlist = []
simrhox = 0
simrhoy = 0
simrhoXlist = []
simrhoYlist = []
simTiltlist = []
simrhoTiltlist = []
pV = Plist[0]*100*(4*math.pi/3)**kappa

print("Simulating an undamped parcel in a real environment.")


if Blist[0]>0:
    simAlist = []
    for i1 in range(cycles):
        if simv < 0 and simtop is False:
            simmaxi = Ai
            simtop = True

        if simtop is False:
            simv += Blist[Ai]*dt
        else:
            simTGD = konvTGSlist[simmaxi] + (Alist[simmaxi]-Alist[Ai])*gammaD/1000
            simB = (simTGD-Tlist[Ai])/virtTlistK[Ai]
            simv += dt*simB
        simx += dt*Slist[Ai]*math.sin(radDlist[Ai])
        simy += dt*Slist[Ai]*math.cos(radDlist[Ai])
        simrhox += dt*Slist[Ai]*math.sin(radDlist[Ai])*Rholist[Ai]/stdRho
        simrhoy += dt*Slist[Ai]*math.cos(radDlist[Ai])*Rholist[Ai]/stdRho
        sima += simv*dt
        t += dt
        if sima < Alist[0]:
            simv = 0
        if sima > Alist[Ai]:
            Ai += 1
        if Ai != 0:
            if sima < Alist[Ai-1]:
                Ai -= 1
        simVlist.append(simv)
        simAlist.append(sima)
        simXlist.append(simx)
        simYlist.append(simy)
        simrhoXlist.append(simrhox)
        simrhoYlist.append(simrhoy)
        vel = (simVlist[i1]**2 + Slist[Ai]**2)**0.5
        tilt = math.acos(abs(simVlist[i1]/vel))
        tilt = tilt*180/math.pi
        rhovel = (simVlist[i1]**2 + (Slist[Ai]*Rholist[Ai]/stdRho)**2)**0.5
        rhotilt = math.acos(abs(simVlist[i1]/rhovel))
        rhotilt = rhotilt*180/math.pi
        simTiltlist.append(tilt)
        simrhoTiltlist.append(rhotilt)
        simtlist.append(t)
        if Alist[Ai]<=CCLalt and simCCLi+10>i1:
            simCCLi = i1
        if i1 % itInt == 0:
            print(str(i1) + " simulation iterations complete.")
        if simUntilCCL is True and simAlist[i1]>CCLalt:
            break
    print(str(i1) + " simulation iterations complete.")
else:
    while max(simAlist)<CCLalt:
        v0 += speedInterval
        print("Simulating for initial speed of " + str(round(v0, 3)) + " m/s")
        t = 0
        simVlist = []
        simAlist = []
        simtlist = []
        Ai = 0
        simv = v0
        sima = Alist[0]
        simtop = False
        simmmaxA = 0
        simCCLi = 0
        simx = 0
        simy = 0
        simXlist = []
        simYlist = []
        simrhox = 0
        simrhoy = 0
        simrhoXlist = []
        simrhoYlist = []
        simTiltlist = []
        simrhoTiltlist = []
        for i1 in range(cycles):
            if simv < 0 and simtop is False:
                simmaxi = Ai
                simtop = True

            if simtop is False:
                simv += Blist[Ai]*dt
            else:
                simTGD = konvTGSlist[simmaxi] + (Alist[simmaxi]-Alist[Ai])*gammaD/1000
                simB = (simTGD-Tlist[Ai])/virtTlistK[Ai]
                simv += dt*simB
            simx += dt*Slist[Ai]*math.sin(radDlist[Ai])
            simy += dt*Slist[Ai]*math.cos(radDlist[Ai])
            simrhox += dt*Slist[Ai]*math.sin(radDlist[Ai])*Rholist[Ai]/stdRho
            simrhoy += dt*Slist[Ai]*math.cos(radDlist[Ai])*Rholist[Ai]/stdRho
            sima += simv*dt
            t += dt
            if sima < Alist[0]:
                simv = 0
            if sima > Alist[Ai]:
                Ai += 1
            if Ai != 0:
                if sima < Alist[Ai-1]:
                    Ai -= 1
            simVlist.append(simv)
            simAlist.append(sima)
            simXlist.append(simx)
            simYlist.append(simy)
            simrhoXlist.append(simrhox)
            simrhoYlist.append(simrhoy)
            vel = (simVlist[i1]**2 + Slist[Ai]**2)**0.5
            tilt = math.acos(abs(simVlist[i1]/vel))
            tilt = tilt*180/math.pi
            rhovel = (simVlist[i1]**2 + (Slist[Ai]*Rholist[Ai]/stdRho)**2)**0.5
            rhotilt = math.acos(abs(simVlist[i1]/rhovel))
            rhotilt = rhotilt*180/math.pi
            simTiltlist.append(tilt)
            simrhoTiltlist.append(rhotilt)
            simtlist.append(t)
            if Alist[Ai]<=CCLalt and simCCLi+10>i1:
                simCCLi = i1
            if i1 % itInt == 0:
                print(str(i1) + " simulation iterations complete.")
            if simUntilCCL is True and simAlist[i1]>CCLalt:
                break
        print(str(i1) + " simulation iterations complete.")


print("Simulating an undamped parcel in a virtual environment.")

if virtBlist[0]<0:
    virtv0 = v0

vt = 0
simvirttlist = []
simvirtv = virtv0
simvirta = Alist[0]
vAi = 0
simvirtVlist = []
simvirtAlist = []
simvCCLi = 0
simvx = 0
simvy = 0
simvXlist = []
simvYlist = []
simvrhox = 0
simvrhoy = 0
simvrhoXlist = []
simvrhoYlist = []
simvTiltlist = []
simvrhoTiltlist = []
simvirttop = False

for i1 in range(cycles):
    if simvirtv < 0 and simvirttop is False:
        simvirtmaxi = vAi
        simvirttop = True

    if simvirttop is False:
        simvirtv += virtBlist[vAi]*dt
    else:
        simvirtTGD = konvTGSlist[simvirtmaxi] + (Alist[simvirtmaxi]-Alist[vAi])*gammaD/1000
        simvirtB = (simvirtTGD-virtTlistC[vAi])/virtTlistK[vAi]
        simvirtv += dt*simvirtB
    simvirta += simvirtv*dt
    simvx += dt*Slist[vAi]*math.sin(radDlist[vAi])
    simvy += dt*Slist[vAi]*math.cos(radDlist[vAi])
    simvrhox += dt*Slist[vAi]*math.sin(radDlist[vAi])*Rholist[vAi]/stdRho
    simvrhoy += dt*Slist[vAi]*math.cos(radDlist[vAi])*Rholist[vAi]/stdRho
    if simvirta < Alist[0]:
        simvirtv = 0
    vt += dt
    if simvirta > Alist[vAi]:
        vAi += 1
    if vAi != 0:
        if simvirta < Alist[vAi-1]:
            vAi -= 1
    simvirtVlist.append(simvirtv)
    simvirtAlist.append(simvirta)
    simvXlist.append(simvx)
    simvYlist.append(simvy)
    simvrhoXlist.append(simvrhox)
    simvrhoYlist.append(simvrhoy)
    vvel = (simvirtVlist[i1]**2 + Slist[vAi]**2)**0.5
    vtilt = math.acos(abs(simvirtVlist[i1]/vvel))
    vtilt = vtilt*180/math.pi
    vrhovel = (simvirtVlist[i1]**2 + (Slist[vAi]*Rholist[vAi]/stdRho)**2)**0.5
    vrhotilt = math.acos(abs(simvirtVlist[i1]/vrhovel))
    vrhotilt = vrhotilt*180/math.pi
    simvTiltlist.append(vtilt)
    simvrhoTiltlist.append(vrhotilt)
    simvirttlist.append(vt)
    if Alist[vAi]<=CCLalt and simvCCLi+10>i1:
        simvCCLi = i1
    if i1 % itInt == 0:
        print(str(i1) + " simulation iterations complete.")
    if simUntilCCL is True and simvirtAlist[i1]>CCLalt:
        break
if min(simvirtAlist) < Alist[0]-50:
    vt = 0
    simvirttlist = [0]
    simvirtv = virtv0
    simvirta = Alist[0]
    vAi = 0
    simvirtVlist = [0]
    simvirtAlist = [0]
    simvCCLi = 0
    simvx = 0
    simvy = 0
    simvXlist = [0]
    simvYlist = [0]
    simvrhox = 0
    simvrhoy = 0
    simvrhoXlist = [0]
    simvrhoYlist = [0]
    simvTiltlist = [0]
    simvrhoTiltlist = [0]
    simvirttop = False
print(str(i1) + " simulation iterations complete.")

print("Simulating a real parcel in a damped environment.")

rt = 0
simrtlist = []
simrv = 0
simra = Alist[0]
rAi = 0
simrVlist = []
simrAlist = []
simrCCLi = 0
simrx = 0
simry = 0
simrXlist = []
simrYlist = []
simrrhox = 0
simrrhoy = 0
simrrhoXlist = []
simrrhoYlist = []
simrTiltlist = []
simrrhoTiltlist = []
simrirttop = False

for i1 in range(cycles):
    if simra <= CCLalt:
        simrv = g*(Tlist[0]-Tlist[rAi]-((Alist[rAi]-Alist[0])*gammaD/1000)+opt)/(k*(Tlist[rAi+25]+273.15))
    simra += simrv*dt
    if Alist[rAi]<=CCLalt and simrCCLi+10>i1:
        simrCCLi = i1
    if simra > Alist[rAi+1]:
        rAi += 1
    if simra < Alist[rAi]:
        rAi -= 1
    simrx += dt*Slist[rAi]*math.sin(radDlist[rAi])
    simry += dt*Slist[rAi]*math.cos(radDlist[rAi])
    simrrhox += dt*Slist[rAi]*math.sin(radDlist[rAi])*Rholist[rAi]/stdRho
    simrrhoy += dt*Slist[rAi]*math.cos(radDlist[rAi])*Rholist[rAi]/stdRho
    rt += dt
    simrtlist.append(rt)
    simrAlist.append(simra)
    simrVlist.append(simrv)
    simrXlist.append(simrx)
    simrYlist.append(simry)
    simrrhoXlist.append(simrrhox)
    simrrhoYlist.append(simrrhoy)
    vel = (simrVlist[i1]**2 + Slist[rAi]**2)**0.5
    rtilt = math.acos(abs(simrVlist[i1]/vel))
    rtilt = rtilt*180/math.pi
    rrhovel = (simrVlist[i1]**2 + (Slist[rAi]*Rholist[rAi]/stdRho)**2)**0.5
    rrhotilt = math.acos(abs(simrVlist[i1]/rrhovel))
    rrhotilt = rrhotilt*180/math.pi
    simrTiltlist.append(rtilt)
    simrrhoTiltlist.append(rrhotilt)
    if simra > CCLalt:
        break
    if i1 % itInt == 0:
        print(str(i1) + " simulation iterations complete.")

print(str(i1) + " simulation iterations complete.")

print("Simulations complete.   Generating plots.")

#Výpočet průměrné hodnoty stoupání
meanclimb = sum(simrVlist)/len(simrVlist)

#Přepočet polohy konvektivní částice do polárních souřadnic

simDlist = []
simAnglelist = []
for i1 in range(len(simXlist)):
    dist = (simXlist[i1]**2 + simYlist[i1]**2)**0.5
    angle = math.atan2(simXlist[i1], simYlist[i1])
    simDlist.append(dist)
    simAnglelist.append(angle)

simvDlist = []
simvAnglelist = []
for i2 in range(len(simvXlist)):
    vdist = (simvXlist[i2]**2 + simvYlist[i2]**2)**0.5
    vangle = math.atan2(simvXlist[i2], simvYlist[i2])
    simvDlist.append(vdist)
    simvAnglelist.append(vangle)

simrhoDlist = []
simrhoAnglelist = []
for i3 in range(len(simrhoXlist)):
    rhodist = (simrhoXlist[i3]**2 + simrhoYlist[i3]**2)**0.5
    rhoangle = math.atan2(simrhoXlist[i3], simrhoYlist[i3])
    simrhoDlist.append(rhodist)
    simrhoAnglelist.append(rhoangle)

simvrhoDlist = []
simvrhoAnglelist = []
for i4 in range(len(simvrhoXlist)):
    vrhodist = (simvrhoXlist[i4]**2 + simvrhoYlist[i4]**2)**0.5
    vrhoangle = math.atan2(simvrhoXlist[i4], simvrhoYlist[i4])
    simvrhoDlist.append(vrhodist)
    simvrhoAnglelist.append(vrhoangle)
    
simrDlist = []
simrAnglelist = []
for i1 in range(len(simrXlist)):
    dist = (simrXlist[i1]**2 + simrYlist[i1]**2)**0.5
    angle = math.atan2(simrXlist[i1], simrYlist[i1])
    simrDlist.append(dist)
    simrAnglelist.append(angle)

simrrhoDlist = []
simrrhoAnglelist = []
for i3 in range(len(simrrhoXlist)):
    rhodist = (simrrhoXlist[i3]**2 + simrrhoYlist[i3]**2)**0.5
    rhoangle = math.atan2(simrrhoXlist[i3], simrrhoYlist[i3])
    simrrhoDlist.append(rhodist)
    simrrhoAnglelist.append(rhoangle)

######################################################################
#Printění dat#########################################################
######################################################################

#Printění dat z povrchu zemského

print("Surface ALTITUDE:                " + str(Alist[0]) + " m AMSL")
print("Surface TEMP:                    " + str(Tlist[0]) + " °C")
print("Surface DEWPOINT:                " + str(Hlist[0]) + " °C")
print("Surface RELATIVE HUMIDITY:       " + str(round(Rlist[0] * 100, 1)) + " %")
print("Surface MIXING RATIO:            " + str(round(Wlist[0]*1000, 2)) + " g/kg")
print("Surface WIND SPEED:              " + str(round(Slist[0], 1)) + " m/s")
print("Surface WIND DIRECTION:          " + str(Dlist[0]) + "°")
print("Station QNH:                     " + str(round(QNH, 2)) + " hPa")
print("Station QFF:                     " + str(round(QFF, 2)) + " hPa")
print('\n' + "Lifting condensation level:      " + str(LCLalt) + " m AMSL = " + str(LCLalt - Alist[0]) + " m AGL, " + str(round(LCLtemp, 1)) + " °C")
print("Convective condensation level:   " + str(CCLalt) + " m AMSL = " + str(CCLalt-Alist[0]) + " m AGL, " + str(CCLtemp) + " °C")
print("Convective temperature:          " + str(round(Tkonv, 1)) + " °C")
print("T-Tkon:                          " + str(round(dTkonv, 1)) + " °C")
print("Level of free convection:        " + str(HVKalt) + " m AMSL")
print("Equilibrium level:               " + str(HNValt) + " m AMSL")
print("Top of clouds:                   " + str(TOPalt) + " m AMSL")
print("New equilibrium level:           " + str(NELalt) + " m AMSL")
print("First tropopause altitude:       " + str(tropo1alt) + " m AMSL")
if tropo2 is True:
    print("Second tropopause altitude:      " + str(tropo2alt) + " m AMSL")
print("CAPE:                            " + str(round(cape)) + " J/kg")
print("CIN:                             " + str(round(cin)) + " J/kg")
print("DCAPE:                           " + str(round(dcape)) + " J/kg")
print("NCAPE:                           " + str(round(ncape, 3)) + " N/kg")
print("Average parcel climb rate:       " + str(round(meanclimb, 3)) + " m/s")


######################################################################
#Kreslení diagramů####################################################
######################################################################

#Vytvoření pomocných seznamů

s = np.linspace(Alist[0], Alist[0], len(Tlist))
x = np.linspace(Tlist[0]+1000, Tlist[-1]-1000, len(Tlist))

#Vytvoření diagramu

z = max([Tlist[0]-plotAmax/100, min(Tlist)-50])

fig = plot.figure(figsize = (20,20), dpi = dpi)
ax1 = plot.subplot(441)
ax2 = plot.subplot(442)
ax3 = plot.subplot(443)
ax01 = ax1.twinx()
ax4 = plot.subplot(444, projection = 'polar')
ax5 = plot.subplot(445)
ax6 = plot.subplot(446)
ax05 = ax5.twiny()
ax10 = plot.subplot(447)
ax11 = plot.subplot(448)
ax12 = plot.subplot(449)
ax7 = plot.subplot(4,4,10)
ax8 = plot.subplot(4,4,11)
ax9 = plot.subplot(4,4,12)
ax13 = plot.subplot(4,4,13, projection = 'polar')
ax14 = plot.subplot(4,4,14, projection = 'polar')
ax15 = plot.subplot(4,4,15)
ax16 = plot.subplot(4,4,16)
ax1.axis([plotTmin, plotTmax, plotAmin, plotAmax])
ax2.axis([0, 100, plotAmin, plotAmax])
ax3.axis([z, thetaElist[0]+plotAmax/125, plotAmin, plotAmax])
ax5.axis([0,40,plotAmin,plotAmax])
if plotAmax>TOPalt:
    ax6.axis([Blist[TOPi]-0.05, max(Blist)+0.05, plotAmin, plotAmax])
else:
    ax6.axis([-max(Blist)-0.05, max(Blist)+0.05, plotAmin, plotAmax])
ax7.axis([0, max([max(simtlist), max(simvirttlist), max(simrtlist)]), min([min(simVlist), min(simrVlist)])*1.025, max([max(simVlist), max(simvirtVlist), max(simrVlist)])])
ax8.axis([0, max([max(simtlist), max(simvirttlist), max(simrtlist)]), Alist[0], max([max(simAlist), max(simvirtAlist), max(simrAlist)])*1.025])
ax9.set_ylim(Alist[0], max([max(simAlist), max(simvirtAlist), max(simrAlist)])*1.025)
ax10.axis([0, max(BruntVaisalaFreqlist), plotAmin, plotAmax])
ax11.axis([-max(BruntVaisalaFreqlist)**2, max(BruntVaisalaFreqlist)**2, plotAmin, plotAmax])
ax12.axis([0,max(Rholist) + 0.025, plotAmin, plotAmax])
ax05.set_xlim(0, 360)
ax05.xaxis.set_major_locator(ticker.MultipleLocator(45))
ax4.set_theta_zero_location('N')
ax4.set_theta_direction(-1)
ax13.set_theta_zero_location('N')
ax13.set_theta_direction(-1)
ax14.set_theta_zero_location('N')
ax14.set_theta_direction(-1)
ax15.axis([0, 90, Alist[0], max([max(simAlist), max(simvirtAlist)])*1.025])
ax16.axis([0, max([max(simDlist), max(simrhoDlist), max(simrDlist), max(simrrhoDlist)])*1.025, Alist[0], max(simAlist)*1.025])
print("Plots initialized.")

#Kreslení křivek

#Emagram
if cape > 0:
    ax1.fill_betweenx(Alist[HVKi:HNVi], konvTGSlist[HVKi:HNVi], Tlist[HVKi:HNVi], facecolor = 'xkcd:creme' )
ax1.plot(WTlist, Alist, color = 'xkcd:gold')
ax1.plot(TGDlist, Alist, color = 'xkcd:orange')
ax1.plot(x, s, color = 'xkcd:black')
ax1.plot(konvTGDlist, Alist, color = 'xkcd:orange')
ax1.plot(virtTlistC, Alist, 'xkcd:maroon', linestyle = 'dotted')
ax1.plot(Hlist, Alist, color = 'b')
ax1.plot(konvTGSlist, Alist, color = 'g')
ax1.plot(TGSlist, Alist, color = 'g')
ax1.plot(Tlist, Alist, color = 'r')
print("Emagram plotted.")

#Diagram relativní vlhkosti
ax2.plot(R100list, Alist, color = 'b')

#Thetagram
for j5 in range(len(thetaTlist)):
    if thetaTlist[j5]+stabDev<thetaTlist[j5-stabDif]:
        ax3.plot(z, Alist[j5], color = 'xkcd:lime', linewidth = 0, marker = '_', markersize = 30)
for j4 in range(len(thetaTlist)):
    if thetaTlist[j4]+stabDev>thetaTlist[j4-stabDif]:
        ax3.plot(z, Alist[j4], color = 'r', linewidth = 0, marker = '_', markersize = 30)
ax3.plot(Tlist, Alist, color = 'r')
ax3.plot(thetaTlist, Alist, color = 'xkcd:pink')
ax3.plot(thetaElist, Alist, color = 'xkcd:magenta')
ax3.plot(ekvTlist, Alist, color = 'r', linestyle = 'dashed')
print("Thetagram plotted.")

#Hodograf
ax4.plot(D12list, S12list, color = 'xkcd:violet')
ax4.plot(D9list, S9list, color = 'b')
ax4.plot(D6list, S6list, color = 'g')
ax4.plot(D3list, S3list, color = 'y')
ax4.plot(D2list, S2list, color = 'xkcd:orange')
ax4.plot(D1list, S1list, color = 'r')
    #ax4.plot(D20list, S20list, color = 'xkcd:violet')

#Profil větru
ax5.plot(Slist, Alist, color = 'r')
ax05.scatter(Dlist, Alist, 10, color = 'b')
for bi in range(len(Alist)):
    if bi%(int(Nbarbs*plotAmax/5000)) == 0:
        ax5.barbs(37.5, Alist[bi], Slist[bi]*math.sin(radDlist[bi])/ktstoms, Slist[bi]*math.cos(radDlist[bi])/ktstoms, barbcolor='k', flagcolor='k', length=6, linewidth=0.5)
    if Alist[bi]>plotAmax:
        break


print("Wind plotted.")

ax6.plot(Blist, Alist, color = 'xkcd:brown')
ax6.plot(virtBlist, Alist, color = 'xkcd:brown', linestyle = 'dotted')

if Blist[0]>0 or v0 != 0:
    if v0 == 0:
        ax7.plot(simtlist, simVlist, color = 'xkcd:violet')
    else:
        ax7.plot(simtlist, simVlist, color = 'xkcd:violet', linestyle = 'dashed')
        ax7.text(0, max(simVlist)*1.025, "Initial vertical speed of: " + str(round(v0, 3)) + ' $m/s$')
        #ax7.plot(simtlist, simvirtfVlist, color = 'xkcd:violet', linestyle = 'dotted')
if virtBlist[0]>0:
    ax7.plot(simvirttlist, simvirtVlist, color = 'xkcd:violet', linestyle = 'dotted')
    ax7.plot(simvirttlist[simvCCLi], simvirtVlist[simvCCLi], color = 'k', marker = '+', markersize = 12)
    ax7.text(simvirttlist[simvCCLi], simvirtVlist[simvCCLi], "   CCL")
ax7.plot(simrtlist, simrVlist, color = 'r')
ax7.plot(simrtlist[simrCCLi], simrVlist[simrCCLi], color = 'k', marker = '+', markersize = 12)
ax7.text(simrtlist[simrCCLi], simrVlist[simrCCLi], "   CCL")

if Blist[0]>0 or v0 != 0:
    if v0 == 0:
        ax8.plot(simtlist, simAlist, color = 'xkcd:violet')
    else:
        ax8.plot(simtlist, simAlist, color = 'xkcd:violet', linestyle = 'dashed')
        #ax8.plot(simvirttlist, simvirtfAlist, color = 'xkcd:violet', linestyle = 'dotted')
if virtBlist[0]>0:
    ax8.plot(simvirttlist, simvirtAlist, color = 'xkcd:violet', linestyle = 'dotted')
    ax8.plot(simvirttlist[simvCCLi], simvirtAlist[simvCCLi], color = 'k', marker = '+', markersize = 12)
    ax8.text(simvirttlist[simvCCLi], simvirtAlist[simvCCLi], "   CCL")
ax8.plot(simrtlist, simrAlist, color = 'r')
ax8.plot(simrtlist[simrCCLi], simrAlist[simrCCLi], color = 'k', marker = '+', markersize = 12)
ax8.text(simrtlist[simrCCLi], simrAlist[simrCCLi], "   CCL")

if Blist[0]>0 or v0 != 0:
    if v0 == 0:
        ax9.plot(simVlist, simAlist, color = 'xkcd:violet')
    else:
        ax9.plot(simVlist, simAlist, color = 'xkcd:violet', linestyle = 'dashed')
        #ax9.plot(simvirtfVlist, simvirtfAlist, color = 'xkcd:violet', linestyle = 'dotted')
if virtBlist[0]>0:
    ax9.plot(simvirtVlist, simvirtAlist, color = 'xkcd:violet', linestyle = 'dotted')
ax9.plot(simrVlist, simrAlist, color = 'r')

ax10.scatter(BruntVaisalaFreqlist, Alist, color = 'xkcd:pink')
ax11.scatter(BruntVaisalaFreqSqlist, Alist, color = 'xkcd:pink')
ax12.plot(Rholist, Alist, color = 'xkcd:brown')
if v0 != 0:
    ax13.plot(simAnglelist, simDlist, color = 'g', linestyle = 'dashed')
    ax14.plot(simrhoAnglelist, simrhoDlist, color = 'xkcd:orange', linestyle = 'dashed')
else:
    ax13.plot(simAnglelist, simDlist, color = 'g')
    ax14.plot(simrhoAnglelist, simrhoDlist, color = 'xkcd:orange')
if virtBlist[0]>0:
    ax13.plot(simvAnglelist, simvDlist, color = 'g', linestyle = 'dotted')
    ax14.plot(simvrhoAnglelist, simvrhoDlist, color = 'xkcd:orange', linestyle = 'dotted')
ax13.plot(simrAnglelist, simrDlist, color = 'r')
ax14.plot(simrrhoAnglelist, simrrhoDlist, color = 'xkcd:pink')


if Blist[0]>0:
    ax15.plot(simTiltlist, simAlist, color = 'g')
    ax15.plot(simrhoTiltlist, simAlist, color = 'xkcd:orange')
else:
    ax15.plot(simTiltlist, simAlist, color = 'g', linestyle = 'dashed')
    ax15.plot(simrhoTiltlist, simAlist, color = 'xkcd:orange', linestyle = 'dashed')
if virtBlist[0]>0:
    ax15.plot(simvTiltlist, simvirtAlist, color = 'g', linestyle = 'dotted')
    ax15.plot(simvrhoTiltlist, simvirtAlist, color = 'xkcd:orange', linestyle = 'dotted')
    
if Blist[0]>0:
    ax16.plot(simDlist, simAlist, color = 'g')
    ax16.plot(simrhoDlist, simAlist, color = 'xkcd:orange')
else:
    ax16.plot(simDlist, simAlist, color = 'g', linestyle = 'dashed')
    ax16.plot(simrhoDlist, simAlist, color = 'xkcd:orange', linestyle = 'dashed')
if virtBlist[0]>0:
    ax16.plot(simvDlist, simvirtAlist, color = 'g', linestyle = 'dotted')
    ax16.plot(simvrhoDlist, simvirtAlist, color = 'xkcd:orange', linestyle = 'dotted')
ax15.plot(simrTiltlist, simrAlist, color = 'r')
ax15.plot(simrrhoTiltlist, simrAlist, color = 'xkcd:pink')
ax16.plot(simrDlist, simrAlist, color = 'r')
ax16.plot(simrrhoDlist, simrAlist, color = 'xkcd:pink')

#Markery

ax1.plot(LCLtemp, LCLalt, color = 'k', marker = '_', markersize = 20)
ax1.plot(CCLtemp, CCLalt, color = 'k', marker = '_', markersize = 20)
ax1.text(LCLtemp, LCLalt, "     LCL     " + str(LCLalt) + 'm, ' + str(round(LCLtemp, 1)) + " °C", horizontalalignment='right')
ax1.text(CCLtemp, CCLalt, "     CCL     " + str(CCLalt) + 'm, ' + str(round(CCLtemp, 1)) + " °C", horizontalalignment='right')
if tropo1alt<plotAmax:
    ax1.plot(tropo1temp, tropo1alt, color = 'k', marker = '_', markersize = 20)
    ax1.text(tropo1temp, tropo1alt, "     TROPO", horizontalalignment='left')
if tropo2 is True and tropo2alt<plotAmax:
    ax1.plot(tropo2temp, tropo2alt, color = 'k', marker = '_', markersize = 20)
    ax1.text(tropo2temp, tropo2alt, "     TROPO", horizontalalignment='left')
ax4.plot(D1list[-1], S1list[-1], color = 'k', marker = 'x', markersize = 12)
ax4.plot(D2list[-1], S2list[-1], color = 'k', marker = 'x', markersize = 12)
ax4.plot(D3list[-1], S3list[-1], color = 'k', marker = 'x', markersize = 12)
ax4.plot(D6list[-1], S6list[-1], color = 'k', marker = 'x', markersize = 12)
ax4.plot(D9list[-1], S9list[-1], color = 'k', marker = 'x', markersize = 12)
ax4.plot(D12list[-1], S12list[-1], color = 'k', marker = 'x', markersize = 12)
ax4.text(D1list[-1], S1list[-1], "  1")
ax4.text(D2list[-1], S2list[-1], "  2")
ax4.text(D3list[-1], S3list[-1], "  3")
ax4.text(D6list[-1], S6list[-1], "  6")
ax4.text(D9list[-1], S9list[-1], "  9")
ax4.text(D12list[-1], S12list[-1], "  12")

if round(cape)>0:
    if plotAmax>HNValt:
        ax6.plot(Blist[HNVi], HNValt, color = 'k', marker = '_', markersize = 20)
        ax6.text(Blist[HNVi], HNValt, "EL     " + str(HNValt) + 'm    ', horizontalalignment='right')
        ax1.plot(Tlist[HNVi], HNValt, color = 'k', marker = '_', markersize = 20)
        ax1.text(Tlist[HNVi], HNValt, "EL     " + str(HNValt) + 'm    ', horizontalalignment='right')
    if plotAmax>HVKalt:
        ax6.plot(Blist[HVKi], HVKalt, color = 'k', marker = '_', markersize = 20)
        ax6.text(Blist[HVKi], HVKalt, "     LFC     " + str(HVKalt) + 'm')
        ax1.plot(Tlist[HVKi], HVKalt, color = 'k', marker = '_', markersize = 20)
        ax1.text(Tlist[HVKi], HVKalt, "     LFC     " + str(HVKalt) + 'm    ', horizontalalignment='left')
    if plotAmax>TOPalt:
        ax6.plot(Blist[TOPi], TOPalt, color = 'k', marker = '_', markersize = 20)
        ax6.text(Blist[TOPi], TOPalt, "     TOP     " + str(TOPalt) + 'm')
        ax1.plot(konvTGSlist[TOPi], TOPalt, color = 'k', marker = '_', markersize = 20)
        ax1.text(konvTGSlist[TOPi], TOPalt, "     TOP     " + str(TOPalt) + 'm', horizontalalignment='left')
    if plotAmax>NELalt and showNEL is True:
        ax6.plot(Blist[NELi], NELalt, color = 'k', marker = '_', markersize = 20)
        ax6.text(Blist[NELi], NELalt, "     NEL     " + str(NELalt) + 'm')
        ax1.plot(NELt, NELalt, color = 'k', marker = '_', markersize = 20)
        ax1.text(NELt, NELalt, "     NEL     " + str(NELalt) + 'm', horizontalalignment='left')
#ax4.plot(D20list[-1], S20list[-1], color = 'k', marker = 'x', markersize = 12)

ax7.plot(simtlist[simCCLi], simVlist[simCCLi], color = 'k', marker = '+', markersize = 12)
ax7.text(simtlist[simCCLi], simVlist[simCCLi], "   CCL")

ax8.plot(simtlist[simCCLi], simAlist[simCCLi], color = 'k', marker = '+', markersize = 12)
ax8.text(simtlist[simCCLi], simAlist[simCCLi], "   CCL")



ax9.plot(0, LCLalt, color = 'k', marker = '_', markersize = 20)
ax9.plot(0, CCLalt, color = 'k', marker = '_', markersize = 20)
ax9.text(0, LCLalt, "      LCL")
ax9.text(0, CCLalt, "      CCL")

if virtBlist[0]>0:
    ax13.plot(simvAnglelist[simvCCLi], simvDlist[simvCCLi], color = 'k', marker = 'x', markersize = 12)
    ax13.text(simvAnglelist[simvCCLi], simvDlist[simvCCLi], "   CCL")
    ax14.plot(simvrhoAnglelist[simvCCLi], simvrhoDlist[simvCCLi], color = 'k', marker = 'x', markersize = 12)
    ax14.text(simvrhoAnglelist[simvCCLi], simvrhoDlist[simvCCLi], "   CCL")
ax13.plot(simAnglelist[simCCLi], simDlist[simCCLi], color = 'k', marker = 'x', markersize = 12)
ax13.text(simAnglelist[simCCLi], simDlist[simCCLi], "   CCL")
ax14.plot(simrhoAnglelist[simCCLi], simrhoDlist[simCCLi], color = 'k', marker = 'x', markersize = 12)
ax14.text(simrhoAnglelist[simCCLi], simrhoDlist[simCCLi], "   CCL")

ax13.plot(simrAnglelist[-1], simrDlist[-1], color = 'k', marker = 'x', markersize = 12)
ax13.text(simrAnglelist[-1], simrDlist[-1], "   CCL")
ax14.plot(simrrhoAnglelist[-1], simrrhoDlist[-1], color = 'k', marker = 'x', markersize = 12)
ax14.text(simrrhoAnglelist[-1], simrrhoDlist[-1], "   CCL")

ax15.plot(simTiltlist[simCCLi], simAlist[simCCLi], color = 'k', marker = '_', markersize = 20)
ax15.text(simTiltlist[simCCLi], simAlist[simCCLi], "CCL     ", horizontalalignment = 'right')
if virtBlist[0]>0:
    ax15.text(simvTiltlist[simvCCLi], simvirtAlist[simvCCLi], "CCL     ", horizontalalignment = 'right')
    ax15.plot(simvTiltlist[simvCCLi], simvirtAlist[simvCCLi], color = 'k', marker = '_', markersize = 20)

ax16.plot(simDlist[simCCLi], simAlist[simCCLi], color = 'k', marker = '_', markersize = 20)


#Mřížky

ax1.grid()
ax2.grid()
ax3.grid()
ax5.grid()
ax6.grid()
ax7.grid()
ax8.grid()
ax9.grid()
ax10.grid()
ax11.grid()
ax12.grid()
ax15.grid()
ax16.grid()

#Tlaková osa Y

ax01.set_ylim(Plist[0], Plist[int(Pmax)]) #Nastavení krajních limitů tlakové osy
ax01.plot(0,0) #Placeholder
ax01.set_yscale('log') #Nastavení tlakové škály na logaritmickou
ax01.yaxis.set_major_locator(ticker.MultipleLocator(50)) #Interval tlakových značek

#Popisky os

ax01.set_ylabel("Pressure [$hPa$]")
ax1.set_ylabel("Altitude AMSL [$m$]")
ax1.set_xlabel("Temperature [$°C$]")
ax3.set_xlabel("Temperature [$°C$]")
ax2.set_xlabel("Relative humidity [%]")
ax5.set_xlabel("Wind speed [$m/s$]")
ax5.set_ylabel("Altitude AMSL [$m$]")
ax05.set_xlabel("Wind direction [$°$]")
ax6.set_xlabel("Parcel buoyancy [$N/kg$] | Vertical acceleration [$m/s^2$]")
ax6.set_ylabel("Altitude AMSL [$m$]")
ax7.set_xlabel("Simulated time [$s$]")
ax7.set_ylabel("Simulated parcel vertical speed [$m/s$]")
ax8.set_xlabel("Simulated time [$s$]")
ax8.set_ylabel("Simulated parcel altitude AMSL [$m$]")
ax9.set_xlabel("Simulated parcel vertical speed [$m/s$]")
ax9.set_ylabel("Simulated parcel altitude AMSL [$m$]")
ax10.set_xlabel("Brunt-Vaisala frequency [$Hz$]")
ax10.set_ylabel("Altitude AMSL [$m$]")
ax11.set_xlabel("Brunt-Vaisala frequency squared [$Hz^2$]")
ax11.set_ylabel("Altitude AMSL [$m$]")
ax12.set_xlabel("Air density [$kg/m^3$]")
ax12.set_ylabel("Altitude AMSL [$m$]")
ax15.set_xlabel("Simulated parcel tilt from vertical [$°$]")
ax15.set_ylabel("Simulated parcel altitude AMSL [$m$]")
ax16.set_xlabel("Simulated parcel distance from origin [$m$]")
ax16.set_ylabel("Simulated parcel altitude AMSL [$m$]")
print("Plots completed.")
#Vykreslení celého grafu
fig.tight_layout()
plot.show()

print('\n' + "Script finished.")
