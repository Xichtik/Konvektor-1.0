# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 18:57:46 2022

@author: Jindřich Novák
"""
print("Script started.")
import pandas
import math
import matplotlib.pyplot as plot
import matplotlib.ticker as ticker
import numpy as np
print("Libraries imported.")

class konvektor:
    def __init__(self):
        
        
        filename = "RAW_11520_2021-12-01_1200.csv" #Název datového souboru
        ktstoms = 0.5144444 #Převodní poměr mezi uzly a metry za sekundu
        g = 9.80665 #Tíhové zrychlení zemské
        n = 0.190284 #Koeficient pro výpočet QNH
        Rv = 461.5 #Měrná plynová konstanta pro vodní páru
        Rd = 287.04 #Měrná plynová konstanta pro suchý vzduch
        eps = Rd/Rv #Poměr měrných plynových konstant
        L = 2264000 #Měrné skupenské teplo vypařování vody za normálních podmínek
        cp = 1.006
        gammaD = g/cp
        stabDev = -0.4 #Při kladných hodnotách bude zvrstvení hodnoceno jako více stabilní, při záporných jako méně stabilní
        stabDif = 5 #
        
        
        plotTmin = -15
        plotTmax = "auto"   #nastavením na "auto" se automaticky nastaví nejvyšší teplota datové řady zvýšená o 1°C
        plotAmin = "auto"   #nastavením na "auto" se automaticky nastaví na hladinu zemského povrchu
        plotAmax = 3000
        dpi = 300 #Rozlišení grafů v DPI
        
        data = pandas.read_csv(filename, header=0) #Přečtení datového souboru
        #Převod sloupců tabulky do jednotlivých seznamů:
        Plist = list(data.Pressure)
        Alist = list(data.Altitude)
        Tlist = list(data.Temperature)
        Hlist = list(data.Dew_point)
        Dlist = list(data.Wind_direction)
        Slist = list(data.Wind_speed)
        print("Data imported.")
        
        ######################################################################
        #Výpočty##############################################################
        ######################################################################
        
        #Převod seznamu rychlostí větru z uzlů na metry za sekundu
        for i in range(len(Slist)):
            Slist[int(i)] = Slist[int(i)]*ktstoms
        
        
        #Výpočet relativní vlhkosti, parciálních tlaků, směšovacího poměru a virtuální teploty
        
        virtTlistK = [] #Seznam virtuálních teplot v Kelvinech
        virtTlistC = [] #Seznam virtuálních teplot ve stupních Celsia
        VTdiff = [] #Seznam rozdílů suché a virtuální teploty
        Rlist = [] #Seznam relativních vlhkostí
        elist = [] #Seznam parciálních tlaků nasycení vzduchu vodní parou
        Elist = [] #Seznam maximálních parciálních tlaků nasycení vzduchu vodní parou
        Wlist = [] #Seznam směšovacích poměrů
        
        for j in range(len(Tlist)):
            E = 0.61094*math.exp((L/Rv)*(1/273.15-(1/(Tlist[int(j)]+273.15))))
            e = 0.61094*math.exp((L/Rv)*(1/273.15-(1/(Hlist[int(j)]+273.15))))
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
            e = 100*Plist[int(i2)]/((eps/Wlist[int(0)]*1000)+1)
            wTd = ((math.log(e/0.61094)/-(L/Rv))+(1/273.15))**(-1)
            WTlist.append(wTd-273.15)
            TGDlist.append(Tgamma)


        #Oprava teploty izogramy
        
        diff = Hlist[0] - WTlist[0]
        for i4 in range(len(WTlist)):
            WTlist[i4] += diff
            
            
        #Výpočet VKH    
        
        for i3 in range(len(WTlist)):
            if TGDlist[i3]<=WTlist[i3]:
                LCLtemp = TGDlist[i3]
                LCLalt = Alist[i3]
                break
            else:
                pass
         
          
        #Výpočet KKH
            
        for i6 in range(len(WTlist)):
            if Tlist[i6]<=WTlist[i6]:
                CCLtemp = Tlist[i6]
                CCLalt = Alist[i6]
                CCLi = i6
                break
            else:
                pass
            
        
        #Výpočet konvektivní teploty a deficitu konvektivní teploty
        
        konvTGDlist = []
        for i7 in range(len(Alist)):
            Konvgamma = Tlist[CCLi]+(Alist[int(i7)] - Alist[CCLi])/1000*(-gammaD)
            konvTGDlist.append(Konvgamma)
        Tkonv = konvTGDlist[0]
        dTkonv = Tlist[0] - konvTGDlist[0]
        
    
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
        
        
        #Výpočet hořejší hodnoty tlakové škály    
            
        Pmax = -1
        for i5 in range(len(Alist)):
            if (Alist[i5] >= plotAmax):
                Pmax = i5
                break
            else:
                pass
        if Pmax == -1:    
            Pmax = len(Plist)-1


        #Automatické nastavení krajních hodnot os diagramu       

        if plotAmin == "auto": 
            plotAmin = Alist[0]
        if plotTmax == "auto":
            plotTmax = max(Tlist)+1
            
            
        print("Calculations completed.")    
            
        ######################################################################
        #Kreslení diagramu####################################################
        ######################################################################
            
        #Přepis seznamů do seznamů knihovny
        y = np.array(Alist)
        s = np.linspace(Alist[0], Alist[0], len(Tlist))
        z = np.linspace(5+min(thetaTlist)+plotAmax/100, 5+min(thetaTlist)+plotAmax/100, len(thetaTlist))
        x = np.linspace(Tlist[0]+1000, Tlist[-1]-1000, len(Tlist))
        t = np.array(Tlist)
        d = np.array(Hlist)
        w = np.array(WTlist)
        a = np.array(TGDlist)
        k = np.array(konvTGDlist)
        v = np.array(virtTlistC)
        th = np.array(thetaTlist)
        thE = np.array(thetaElist)
        eT = np.array(ekvTlist)
        r = np.array(R100list)
        ws = np.array(Slist)
        wd = np.array(Dlist)
        
        #Vytvoření diagramu
        fig, diagram = plot.subplots(2, 3, figsize=(20,16), dpi = dpi, sharey = 'row')
        diagram[0,0].axis([plotTmin, plotTmax, plotAmin, plotAmax])
        diagram[0,1].set_xlim(0, 100)
        diagram[0,2].set_xlim(max(ekvTlist)-plotAmax/100, min(thetaTlist)+plotAmax/100)
        
        print("Generating plots.")
        print("Drawing lines.")
        z = max(ekvTlist)-plotAmax/100
        
        diagram[0,0].plot(w, y, color = 'xkcd:gold')
        diagram[0,0].plot(a, y, color = 'xkcd:orange')
        diagram[0,0].plot(x, s, color = 'xkcd:black')
        diagram[0,0].plot(k, y, color = 'xkcd:orange')
        diagram[0,0].plot(v, y, 'xkcd:maroon', linestyle = 'dotted')
        diagram[0,0].plot(t, y, color = 'r')
        diagram[0,0].plot(d, y, color = 'b')
        for j5 in range(len(thetaTlist)):
            if thetaTlist[j5]+stabDev<thetaTlist[j5-stabDif]:
                diagram[0,2].plot(z, Alist[j5], color = 'xkcd:lime', linewidth = 0, marker = '_', markersize = 30)
        for j4 in range(len(thetaTlist)):
            if thetaTlist[j4]+stabDev>thetaTlist[j4-stabDif]:
                diagram[0,2].plot(z, Alist[j4], color = 'xkcd:orange', linewidth = 0, marker = '_', markersize = 30)
        diagram[0,2].plot(t, y, color = 'r')
        diagram[0,2].plot(th, y, color = 'xkcd:pink')
        diagram[0,2].plot(thE, y, color = 'xkcd:magenta')
        diagram[0,2].plot(eT, y, color = 'r', linestyle = 'dashed')
        diagram[0,1].plot(r, y, color = 'b')
        diagram[0,0].plot(LCLtemp, LCLalt, color = 'xkcd:black', marker = '_', markersize = 20)
        diagram[0,0].plot(CCLtemp, CCLalt, color = 'xkcd:black', marker = '_', markersize = 20)
        diagram[0,0].text(LCLtemp, LCLalt, "     LCL     " + str(LCLalt) + 'm, ' + str(round(LCLtemp, 1)) + " °C")
        diagram[0,0].text(CCLtemp, CCLalt, "     CCL     " + str(CCLalt) + 'm, ' + str(round(CCLtemp, 1)) + " °C")
        
        diagram[1,0].plot.subplots(subplot_kw = {'projection': 'polar'})
        diagram[1,0].plot(wd, ws, color = 'r')
        
        
        #Mřížky
        diagram[0,0].grid("True")
        diagram[0,1].grid("True")
        diagram[0,2].grid("True")
        diagram[1,0].grid("True")
        
        #Popisky os
        diagram[0,0].set_ylabel("Altitude AMSL [m]")
        diagram[0,0].set_xlabel("Temperature [°C]")
        diagram[0,2].set_xlabel("Temperature [°C]")
        diagram[0,1].set_xlabel("Relative humidity [%]")
        
        #Tlakový osa Y
        Paxis = diagram[0,0].twinx() #Vytvoření druhé osy Y
        Paxis.set_ylim(Plist[0], Plist[int(Pmax)]) #Nastavení krajních limitů tlakové osy
        Paxis.plot(0,0) #Placeholder
        Paxis.set_yscale('log') #Nastavení tlakové škály na logaritmickou
        Paxis.yaxis.set_major_locator(ticker.MultipleLocator(50)) #Interval tlakových značek
        Paxis.set_ylabel("Pressure [hPa]")
        
        fig.tight_layout()
        plot.show() #Vykreslení celého grafu
        
        ######################################################################
        #Printění dat#########################################################
        ######################################################################
        
        #Printění dat z povrchu zemského
        
        print("Surface ALTITUDE:                " + str(Alist[0]) + " m AMSL")
        print("Surface TEMP:                    " + str(Tlist[0]) + " °C")
        print("Surface DEWPOINT:                " + str(Hlist[0]) + " °C")
        print("Surface RELATIVE HUMIDITY:       " + str(round(Rlist[0] * 100, 1)) + " %")
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
        print("        QNH:                     " + str(round(QNH, 2)) + " hPa")
        print("        QFF:                     " + str(round(QFF, 2)) + " hPa")
        print('\n' + "Lifting condensation level:      " + str(LCLalt) + " m AMSL = " + str(LCLalt - Alist[0]) + " m AGL, " + str(round(LCLtemp, 1)) + " °C")
        print("Convective condensation level:   " + str(CCLalt) + " m AMSL = " + str(CCLalt-Alist[0]) + " m AGL, " + str(CCLtemp) + " °C")
        print("Convective temperature:          " + str(round(Tkonv, 1)) + " °C")
        print("T-Tkon:                          " + str(round(dTkonv, 1)) + " °C")
        
        
        
        
        
k = konvektor()