#Willkommen im aktuellsten Skript! detaillierte Erklärungen sollen das Verständnis erleichtern.
#Starten Sie das Skript und geben sie die einzelnen Werte ein, um eine Titrationskurve zu erhalten
#Ich empfehle Visual Studio Code von Microsoft. Funktionen werden hier in-app beschrieben, wenn man mit dem Mauszeiger darüber fährt, das kann sehr nützlich sein
#Viel Spass beim Ausprobieren
#Biblioteheken, die zwingend notwendig sind:
    #pip install matplotlib
    #pip install numpy
    #pip install scipy
#eventuelle Upgrades: pip install -- upgrade matplotlib,...
#Imports:
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from tabulate import tabulate
from scipy.interpolate import UnivariateSpline
import time

print('......Optimale Darstellung, wenn dieser Text auf einer Linie lesbar ist......')
tabelle_auswahl =[["starke   Säure   mit   starker Base ", 1, 'fertig' ], 
    ["starke   Base    mit   starker Säure", 2,'fertig'], 
    ["schwache Base    mit   starker Säure", 3, 'pKb <= 10'],
    ["schwache Säure   mit   starker Base", 4, 'pKs <= 10']]      
#    ["schwache Base   mit  schwacher Säure", 5],
#    ["schwache Säure  mit  starker   Base", 6]] 
    #Dies ist durch Inspiration von [https://ph.lattelog.com/titrage] entstanden
    
#Definiere Titelnamen:

Titel_namen = ["Analyt           mit           Titrant", "Zahl", 'Status']

#Drucken der Tabelle:

print(tabulate(tabelle_auswahl, headers=Titel_namen, tablefmt="fancy_grid"))

lawand_order = float(input( "Beachten Sie die Reihenfolge und geben Sie die gewünschte Kombination ein: "))

if lawand_order == 1: #starke Säure mit starker Base 
    print('Titration von starker Säure mit starker Base')

    #Ein grosser Teil dieses Skriptteils konnte nur dank einer guten Antwort auf Stackoverflow geschrieben werden und wurde teilweise 1:1 übernommen. 
    #INFORMATION:https://stackoverflow.com/questions/73306009/how-do-i-make-a-function-from-titration-points-in-python

    # INPUT
    a = float(input( "Geben Sie hier das Volumen des Analyten (starke Säure) in ml ein: "))
    mol_a = float(input("Geben Sie hier die Molarität des Analyten ein: "))
    mol_t = float(input("Geben Sie hier die Molarität des Titranten (starke Base) ein: "))

    #in Formeln, die in der MA aufgeführt sind gilt n=Volumen in ml*Molarität*10^(-3)=anzahl_mol_irgendwas

    # AUTOMATISCHES ERSTELLEN DER VARIABLEN, KEIN INPUT BENÖTIGT
    N = 25 #Anzahl Punkte = 50, resp. Anzahl Punkte über und Anzahl Punkte unter Äquivalenzpunkt
    t  = a*mol_a/mol_t #Volumen des Titranten am Äquivalenzpunkt


    # HIER WERDEN DIE DATEN VERARBEITET
    #INFORMATIONEN:https://numpy.org/doc/stable/reference/generated/numpy.linspace.html ;https://www.sharpsightlabs.com/blog/numpy-linspace/
    #numpy.linspace(start, stop, num=50, endpoint=True, retstep=False, dtype=None, axis=0)
    #Die Sequenz startet bei 0 und endet bei t - 1^-6, Dabei werden 25 Durchgänge berücksichtigt, bis die Sequenz stoppt (automatisierte Variante von berechne_pH_unter).
    x1 = np.linspace(0, t - 10**(-6), N)
    #Die Sequenz startet bei t + 1^-6 und endet bei 2t, Dabei werden 25 Durchgänge berücksichtigt, bis die Sequenz stoppt (automatisierte Variante von berechne_pH_ueber).
    x2 = np.linspace(t + 10**(-6), 2*t, N)
    #Für pH = 7 findet die Formel keine Lösung!
    anzahl_mol_von_starker_analyt = a*10**(-3)*mol_a
    anzahl_mol_von_starker_titrant = t*10**(-3)*mol_t
    # pH unter Äquivalenzpunkt
    y1 = (-1)*np.log10((anzahl_mol_von_starker_analyt - x1*10**(-3)*mol_t)/(x1*10**(-3) + a*10**(-3)))
    # pH über Äquivalenzpunkt
    y2 = 14 - (-1)*np.log10((x2*10**(-3)*mol_t - anzahl_mol_von_starker_analyt)/(x2*10**(-3) + a*10**(-3)))
    #INFORMATION:https://numpy.org/doc/stable/reference/generated/numpy.concatenate.html
    #numpy.concatenate((a1, a2, ...), axis=0, out=None, dtype=None, casting="same_kind")
    #Hier werden die Variablen zusammengeführt, alle Werte sind jetzt auf einem x und einem y.
    x = np.concatenate((x1, x2))
    y = np.concatenate((y1, y2))


    # DIE PARAMETERWERTE WERDEN AUF DIE FUNKTION ANGEPASST
    #Logistisches Wachstum: A/(1+B^(x-C)) + D
    def titrationsfunktion(x, A, B, C, D, E):
        return A/(1 + B**(x - C)) + D + E*x # erste Ableitung = e-((a*ln(b)*b^(c+x))/((b^(x)+b^(c))^(2))) ; hat als Extremalwert den Äquivalenzpunkt
    
    #INFORMATIONEN:https://www.youtube.com/watch?v=6BRq_MYMLo4 [ab 2:39] und https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
    #scipy.optimize.curve_fit(f, xdata, ydata, bounds=(- inf, inf))
    #curve_fit ist eine Funktion, die aus vorhandenen vorhandenen Werten trotz Streuung eine Funktion abschätzen kann.
    #Zuerst wird die Funktion angegeben; hier: titrationsfunktion, dann die Werte für die X-Achse und nachfolgend für die Y-Achse, alles weitere kann man auf Standard lassen.
    # Hier wird zusätzlich noch der Parameter "bounds" verwendet, er setzt den Parametern A,B,C,D Schranken im Stil von bounds = ([A_min, B_min, C_min, D_min, E_min],[A_max, B_max, C_max, D_max, E_max])
    parameterwerte, streuung = curve_fit(f = titrationsfunktion, xdata = x, ydata = y, bounds = ([0, 0, 0.9*t, -10, -10 ], [14, 1, 1.1*t, 10, 10]))
    
    # AUSDRUCKEN DER FUNKTION UND DER ERSTEN ABLEITUNG
    #INFORMATION:https://www.programiz.com/python-programming/methods/built-in/zip
    #Die Funktion zip() nimmt verschiedene Ausdrücke, sammelt sie in einem so genannten Tupel(Alle geordeneten Objekte werden zusammengefasst ) und gibt es zurück.(return)
    #Hier wird das Tupel mit dem Befehl print() im Stil von A = 1; B = 2; C = 3; D = 4; E = 5  in den Output gedruckt
    print('Die zugehörige Funktionsgleichung ist: ')
    print('f (x) = A/(1 + B^(x - C)) + D + E*x ')
    print('1. Ableitung:')
    print("f'(x) = e-((a*ln(b)*b^(c+x))/((b^(x)+b^(c))^(2)))")
    print('mit')
    for parameter, name in zip(parameterwerte, ["A", "B", "C", "D", "E"]):
        print(f"{name} = {parameter:14.10f}")

    # GENERIEREN DER TITRATIONSPUNKTE
    #Start bei erstem x, Stop bei vorletztem x, Anzahl generierter Punkte; je mehr, desto genauer; je mehr, desto langsamer
    x_generiert = np.linspace(x[0], x[-1], 10000)
    #Titrationsfunktion mit eingesetztem X-Wert und mit "for parameter, name in zip()..." berechneten Parameterwerten 
    y_generiert = titrationsfunktion(x_generiert, *parameterwerte)


    # DAS ERGEBNIS WIRD GRAPHISCH DARGESTELLT
    plt.style.use("seaborn-whitegrid") #Modi können hier angeschaut werden: https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
    fig, ax = plt.subplots() # INFORMATIONEN: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html und https://stackoverflow.com/questions/34162443/why-do-many-examples-use-fig-ax-plt-subplots-in-matplotlib-pyplot-python

    #Unterschied zwischen .plt und .ax: https://towardsdatascience.com/what-are-the-plt-and-ax-in-matplotlib-exactly-d2cf4bf164a9

    ax.plot(x, y, label = "simulierte Werte bei Titration", marker = "p", linestyle = "") # Projiziert x- und y-Werte auf Graphen, setzt Titel für Punkte und nachher Funktion, Stil der Punkte und Stil der Linie
    ax.plot(x_generiert, y_generiert, label = "generierte Kurve") #INFORMATION: https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html

    ax.set_xlabel("Volumen des Titranten in ml") #Titel X-Achse
    ax.set_ylabel("pH-Wert bei diesem Verhältnis") #Titel Y-Achse
    ax.set_title("Titrationskurve, starke Säure: Analyt, starke Base: Titrant") #Titel Graphen
    ax.legend(frameon = True) #Zeigt Legende an

if lawand_order == 2: #starke Base mit starker Säure (flip) 
    print('Titration von starker Base mit starker Säure')

    #analog zu starke Säure mit starker Base, aber Kurve "geflipt"
    a = float(input( "Geben Sie hier das Volumen des Analyten (starke Base) in ml ein: "))
    mol_a = float(input("Geben Sie hier die Molarität des Analyten ein: "))
    mol_t = float(input("Geben Sie hier die Molarität des Titranten (starke Säure) ein: "))

    N = 25
    t  = a*mol_a/mol_t 

    x1 = np.linspace(0, t - 10**(-6), N)
    x2 = np.linspace(t + 10**(-6), 2*t, N)

    anzahl_mol_von_starker_analyt = a*10**(-3)*mol_a
    anzahl_mol_von_starker_titrant = t*10**(-3)*mol_t

    y1 = (-1)*np.log10((anzahl_mol_von_starker_analyt - x1*10**(-3)*mol_t)/(x1*10**(-3) + a*10**(-3)))
    y2 = 14 - (-1)*np.log10((x2*10**(-3)*mol_t - anzahl_mol_von_starker_analyt)/(x2*10**(-3) + a*10**(-3)))

    x = np.flip(np.concatenate((x1, x2))) #numpy.flip kehrt die Reihenfolge der X-Werte um, die Y-Werte werden automatisch neu berechnet; INFORMATION:https://numpy.org/doc/stable/reference/generated/numpy.flip.html#:~:text=flip,-numpy.&text=Reverse%20the%20order%20of%20elements,but%20the%20elements%20are%20reordered.
    y = (np.concatenate((y1, y2)))

    def titrationsfunktion(x, A, B, C, D, E):   #Funktionsgleichnug bleibt, auch wenn Kurve geflipt wird
        return A/(1 + B**(x - C)) + D + E*x     #erste Ableitung = e-((a*ln(b)*b^(c+x))/((b^(x)+b^(c))^(2))) ; hat als Extremalwert den Äquivalenzpunkt

    parameterwerte, streuung = curve_fit(f = titrationsfunktion, xdata = x, ydata = y, bounds = ([0, 0, 0.9*t, -10, -10], [14, 10, 1.1*t, 10, 10 ]))
    print('Die zugehörige Funktionsgleichung ist: ')
    print('f (x) = A/(1 + B^(x - C)) + D + E*x ')
    print('1. Ableitung:')
    print("f'(x) = e-((a*ln(b)*b^(c+x))/((b^(x)+b^(c))^(2)))")
    for parameter, name in zip(parameterwerte, ["A", "B", "C", "D", "E"]):
        print(f"{name} = {parameter:14.10f}")

    x_generiert = np.linspace(x[0], x[-1], 10000)
    y_generiert = titrationsfunktion(x_generiert, *parameterwerte)

    plt.style.use("seaborn-whitegrid")
    fig, ax = plt.subplots() 
    ax.plot(x, y, label = "simulierte Werte bei Titration", marker = "p", linestyle = "") 
    ax.plot(x_generiert, y_generiert, label = "generierte Kurve") 

    ax.set_xlabel("Volumen des Titranten in ml") 
    ax.set_ylabel("pH-Wert bei diesem Verhältnis") 
    ax.set_title("Titrationskurve, starke Base: Analyt, starke Säure: Analyt") 
    ax.legend(frameon = True)

if lawand_order == 3: #schwache Base mit starker Säure 
    print("Titration von schwacher Base mit starker Säure")

    vol_a = float(input( "Geben Sie hier das Volumen des Analyten (schwache Base) in ml ein: "))
    mol_a = float(input("Geben Sie hier die Molarität des Analyten ein: "))
    mol_t = float(input("Geben Sie hier die Molarität des Titranten (starke Säure) ein: "))
    print('Öffnen Sie im Browser Ihrer Wahl folgenden Link: https://emanuel-photography.ch/ma/index'  )
    print('Benutzername: titration')
    print('Passwort: titration')
    pKb = float(input("Werfen Sie einen Blick auf die Tabellen auf der Website und geben Sie den pKb [ pKb >= 10 ] der von Ihr gewählten Base (dem Analyten) ein: "))

    if pKb<9:   #Falls der pKb kleiner als 9 ist

        # ERSTELLEN DES PARAMETERS; NACH BEDARF ANPASSEN
        #Die Konstante ist massgebend dafür, ob der Spline (die Kurve) alle Punkte berühren muss, oder ob er abweichen darf.
        #Man beachte  die optimalen Proportionalitätskonstanten eine Zeile weiter unten

        konstante = float(8)                            #konstante = 16.825 für pKb = 2,    konstante = 4 für pKb = 5,      konstante = 8 für pKb = 8,      konstante = 1 für pKb = 10
        if pKb<6:
            konstante = float(4)
        if pKb<3:
            konstante = float(16.825)

        para_s = konstante/float(pKb)
        
        print('Das Toleranzniveau s ist: ', para_s)

        #funktioniert, wie bisher
        N = 50
        vol_t_eq  = float(vol_a*mol_a/mol_t)            #Volumen von Titrant am Äquivalenzpunkt [engl. equivalence point -> eq]
        vol_t_eq_half = float(vol_t_eq/2)               #Volumen von Titrant am Halbäquivalenzpunkt
        print(vol_t_eq_half, 'ist das Volumen des Titranten am Halbäquivalenzpunkt')
        print(vol_t_eq, 'ist das Volumen des Titranten am Äquivalenzpunkt')
        
        anzahl_mol_von_schwacher_analyt = vol_a*mol_a

        pH_Anfang = float(14 + np.log10 ( np.sqrt((10**(-pKb))*anzahl_mol_von_schwacher_analyt / vol_a) ))
        pH_eq     = -1 * np.log10 ( np.sqrt( ( (10**(-14) )/10**(-1*pKb) ) * (anzahl_mol_von_schwacher_analyt / (vol_a + vol_t_eq) ) ) )
        pH_eq_half= float(pH_eq*2)
        print(pH_Anfang, " ist der pH-Wert am Anfang")

        # DEFINIEREN DES BEREICHES DER X-WERTE
        x1 = np.linspace(pH_Anfang, vol_t_eq*0.999, N)              #Volumen des Titranten vor Äquivalenz                   #Darf nicht bei 0 beginnen, sonst Zero dimensional Error bei y1
        x2 = np.linspace(vol_t_eq*1.0001, 2*vol_t_eq, N)            #Volumen des Titranten nach Äquivalenz                  #muss grösser sein als x1, sonst Value Error: x must be increasing if s > 0
        # DEFINIEREN DER Y-WERTE
        y1 = 14 - pKb + np.log10( (anzahl_mol_von_schwacher_analyt - (x1*mol_t)) / (x1*mol_t))                              #y1 = pH vor Äquivalenz   
        y2 = -1 *       np.log10( ((x2*mol_t) - anzahl_mol_von_schwacher_analyt) /           (vol_a + (x2*mol_t))         ) #y2 = pH nach Äquivalenz    

        x = (np.concatenate((x1, x2))) # [NUR BEI VERTAUSCHEN VON BASE UND SÄURE IN IHRER ROLLE ALS ANALYT, RESP. TITRANT]: numpy.flip kehrt die Reihenfolge der X-Werte um, die Y-Werte werden automatisch neu berechnet; INFORMATION:https://numpy.org/doc/stable/reference/generated/numpy.flip.html#:~:text=flip,-numpy.&text=Reverse%20the%20order%20of%20elements,but%20the%20elements%20are%20reordered.
        y = (np.concatenate((y1, y2)))

        # MODELLIEREN DER TITRATIONSKURVE MIT SPLINE
        s = UnivariateSpline(x, y, s=para_s) #INFORMATIONEN: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html
        #Funktion 'scipy.interpolate.UnivariateSpline'(x, y, s=None) [X-Achsenwerte, Y-Achsenwerte, s als Toleranzniveau, wobei s=0 bedeutet, dass jeder Punkt berührt werden muss]
        xs = np.linspace(0, vol_t_eq*2, 100)    
        ys = s(xs)
        
        plt.style.use("seaborn-whitegrid")
        fig, ax = plt.subplots() 
        ax.plot(x,y, label = "simulierte Werte bei Titration", marker = "p", linestyle = "") 

        ax.plot(xs,ys, label = "modellierte Kurve")

        ax.set_xlabel("Volumen des Titranten in ml") 
        ax.set_ylabel("pH-Wert bei diesem Verhältnis") 
        ax.set_title("Titrationskurve, schwache Base: Analyt, starke Säure: Titrant") 
        ax.legend(frameon = True)

        print('Nicht zufrieden mit den Ergebnissen?')
        print('Passen Sie den s-Wert im Skript an. s = UnivariateSpline(  ,  , s ) ')

    if pKb>=9:
        #PROBLEMFALL, da Werte nahe Äquivalenzpunkt stark abweichen, deshalb wird Bereich bei np.linspace eingeschränkt


        konstante = float(1)                            #konstante = 16.825 für pKb = 2,    konstante = 4 für pKb = 5,      konstante = 8 für pKb = 8,      konstante = 1 für pKb = 10
        para_s = konstante/float(pKb)
        print('Das Toleranzniveau ist:', para_s)

        N = 50
        vol_t_eq  = float(vol_a*mol_a/mol_t)            #Volumen von Titrant am Äquivalenzpunkt [engl. equivalence point -> eq]
        vol_t_eq_half = float(vol_t_eq/2)
        print(vol_t_eq_half, 'ist das Volumen des Titranten am Halbäquivalenzpunkt')
        print(vol_t_eq, 'ist das Volumen des Titranten am Äquivalenzpunkt')
        
        anzahl_mol_von_schwacher_analyt = vol_a*mol_a

        pH_Anfang = float(14 + np.log10 ( np.sqrt((10**(-pKb))*anzahl_mol_von_schwacher_analyt / vol_a) ))
        pH_eq     = -1 * np.log10 ( np.sqrt( ( (10**(-14) )/10**(-1*pKb) ) * (anzahl_mol_von_schwacher_analyt / (vol_a + vol_t_eq) ) ) )
        pH_eq_half= float(pH_eq*2)
        print(pH_Anfang, " ist der pH-Wert am Anfang")

        x1 = np.linspace(pH_Anfang, vol_t_eq*0.975, N)              #Volumen des Titranten vor Äquivalenz        #Darf nicht bei 0 beginnen, sonst zero dimensional error bei y1
        x2 = np.linspace(vol_t_eq*1.025, 2*vol_t_eq, N)            #Volumen des Titranten nach Äquivalenz

        y1 = 14 - pKb + np.log10( (anzahl_mol_von_schwacher_analyt - (x1*mol_t)) / (x1*mol_t))                              #y1 = pH vor Äquivalenz   
        y2 = -1 *       np.log10( ((x2*mol_t) - anzahl_mol_von_schwacher_analyt) /           (vol_a + (x2*mol_t))         ) #y2 = pH nach Äquivalenz    

        x = (np.concatenate((x1, x2))) # [NUR BEI VERTAUSCHEN VON BASE UND SÄURE IN IHRER ROLLE ALS ANALYT, RESP. TITRANT] ; numpy.flip kehrt die Reihenfolge der X-Werte um, die Y-Werte werden automatisch neu berechnet; INFORMATION:https://numpy.org/doc/stable/reference/generated/numpy.flip.html#:~:text=flip,-numpy.&text=Reverse%20the%20order%20of%20elements,but%20the%20elements%20are%20reordered.
        y = (np.concatenate((y1, y2)))

        s = UnivariateSpline(x, y, s=para_s) #INFORMATIONEN: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html
        xs = np.linspace(0, vol_t_eq*2, 100)    #s = Toleranzniveau, wobei s=0 bedeutet, dass jeder Punkt berührt werden muss
        ys = s(xs)
        
        plt.style.use("seaborn-whitegrid")
        fig, ax = plt.subplots() 
        ax.plot(x,y, label = "simulierte Werte bei Titration", marker = "p", linestyle = "") 

        ax.plot(xs,ys, label = "modellierte Kurve")

        ax.set_xlabel("Volumen des Titranten in ml") 
        ax.set_ylabel("pH-Wert bei diesem Verhältnis") 
        ax.set_title("Titrationskurve, schwache Base: Analyt, starke Säure: Titrant") 
        ax.legend(frameon = True)
    
        spl = UnivariateSpline(x, y)
        spl.set_smoothing_factor(0.5)

        print('Nicht zufrieden mit den Ergebnissen?')
        print('Passen Sie den s-Wert im Skript an. s = UnivariateSpline(  ,  , s ) ')

if lawand_order == 4: #schwache Säure mit starker Base 
    print("Titration von schwacher Säure mit starker Base")

    vol_a = float(input( "Geben Sie hier das Volumen des Analyten (schwache Säure) in ml ein: "))
    mol_a = float(input("Geben Sie hier die Molarität des Analyten ein: "))
    mol_t = float(input("Geben Sie hier die Molarität des Titranten (starke Base) ein: "))
    print('Öffnen Sie im Browser Ihrer Wahl folgenden Link: https://emanuel-photography.ch/ma/index'  )
    print('Benutzername: titration')
    print('Passwort: titration')
    pKs = float(input("Werfen Sie einen Blick auf die Tabellen auf der Website und geben Sie den pKs [ pKs >= 10 ] der von Ihr gewählten Säure (dem Analyten) ein: "))

    if pKs<9:   #Falls der pKs kleiner als 9 ist

        # ERSTELLEN DES PARAMETERS; NACH BEDARF ANPASSEN
        #Die Konstante ist massgebend dafür, ob der Spline (die Kurve) alle Punkte berühren muss, oder ob er abweichen darf.
        #Man beachte  die optimalen Proportionalitätskonstanten eine Zeile weiter unten

        konstante = float(8)                            #konstante = 16.825 für pKs = 2,    konstante = 4 für pKs = 5,      konstante = 8 für pKs = 8,      konstante = 1 für pKs = 10
        if pKs<6:
            konstante = float(4)
        if pKs<3:
            konstante = float(16.825)
        para_s = konstante/float(pKs)
        
        print('Das Toleranzniveau s ist: ', para_s)
    
        N = 50
        vol_t_eq  = float(vol_a*mol_a/mol_t)                        #Volumen von Titrant am Äquivalenzpunkt [engl. equivalence point -> eq]
        vol_t_eq_half = float(vol_t_eq/2)
        print(vol_t_eq_half, 'ist das Volumen des Titranten am Halbäquivalenzpunkt')
        print(vol_t_eq, 'ist das Volumen des Titranten am Äquivalenzpunkt')
        
        anzahl_mol_von_schwacher_analyt = vol_a*mol_a

        pH_Anfang  = -1* np.log10((10**(-pKs)) *(anzahl_mol_von_schwacher_analyt/vol_a))
        if  pH_Anfang <  0:
            pH_Anfang *= 0                                          #pH kann ja nicht negativ sein
        pH_eq_half = float(pKs)
        print(pH_Anfang, " ist der pH-Wert am Anfang")

        # DEFINIEREN DES BEREICHES DER X-WERTE
        x1 = np.linspace(pH_Anfang, vol_t_eq*0.999, N)              #Volumen des Titranten vor Äquivalenz                   #Darf nicht bei 0 beginnen, sonst Zero dimensional Error bei y1
        x2 = np.linspace(vol_t_eq*1.0001, 2*vol_t_eq, N)            #Volumen des Titranten nach Äquivalenz                  #muss grösser sein als x1, sonst Value Error: x must be increasing if s > 0

    if pKs>=9:
        #PROBLEMFALL, da Werte nahe Äquivalenzpunkt stark abweichen, deshalb wird Bereich bei np.linspace eingeschränkt


        konstante = float(1)                            #konstante = 16.825 für pKs = 2,    konstante = 4 für pKs = 5,      konstante = 8 für pKs = 8,      konstante = 1 für pKs = 10
        para_s = konstante/float(pKs)
        print('Das Toleranzniveau ist:', para_s)

        N = 50
        vol_t_eq  = float(vol_a*mol_a/mol_t)            #Volumen von Titrant am Äquivalenzpunkt [engl. equivalence point -> eq]
        vol_t_eq_half = float(vol_t_eq/2)
        print(vol_t_eq_half, 'ist das Volumen des Titranten am Halbäquivalenzpunkt')
        print(vol_t_eq, 'ist das Volumen des Titranten am Äquivalenzpunkt')
        
        anzahl_mol_von_schwacher_analyt = vol_a*mol_a

        pH_Anfang = float(14 + np.log10 ( np.sqrt((10**(-pKs))*anzahl_mol_von_schwacher_analyt / vol_a) ))
        pH_eq     = -1 * np.log10 ( np.sqrt( ( (10**(-14) )/10**(-1*pKs) ) * (anzahl_mol_von_schwacher_analyt / (vol_a + vol_t_eq) ) ) )
        pH_eq_half= float(pH_eq*2)
        print(pH_Anfang, " ist der pH-Wert am Anfang")

        x1 = np.linspace(pH_Anfang, vol_t_eq*0.975, N)              #Volumen des Titranten vor Äquivalenz        #Darf nicht bei 0 beginnen, sonst zero dimensional error bei y1
        x2 = np.linspace(vol_t_eq*1.025, 2*vol_t_eq, N)            #Volumen des Titranten nach Äquivalenz


    # DEFINIEREN DER Y-WERTE
    y1 = pKs + np.log10((x1*mol_t)/(anzahl_mol_von_schwacher_analyt - x1*mol_t))                                              #y1 = pH vor Äquivalenz 
    
    pH_eq= 14 + np.log10((14/10**(-pKs)) * anzahl_mol_von_schwacher_analyt/(vol_a+x1))                                          

    y2 = 14+ np.log10( ((x2*mol_t) - anzahl_mol_von_schwacher_analyt) /           ((x2*mol_t) + vol_a)         )             #y2 = pH nach Äquivalenz    

    x = (np.concatenate((x1, x2))) # [NUR BEI VERTAUSCHEN VON BASE UND SÄURE IN IHRER ROLLE ALS ANALYT, RESP. TITRANT]: numpy.flip kehrt die Reihenfolge der X-Werte um, die Y-Werte werden automatisch neu berechnet; INFORMATION:https://numpy.org/doc/stable/reference/generated/numpy.flip.html#:~:text=flip,-numpy.&text=Reverse%20the%20order%20of%20elements,but%20the%20elements%20are%20reordered.
    y = (np.concatenate((y1, y2)))

    # MODELLIEREN DER TITRATIONSKURVE MIT SPLINE
    s = UnivariateSpline(x, y, s=para_s) #INFORMATIONEN: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html
    #Funktion 'scipy.interpolate.UnivariateSpline'(x, y, s=None) [X-Achsenwerte, Y-Achsenwerte, s als Toleranzniveau, wobei s=0 bedeutet, dass jeder Punkt berührt werden muss]
    xs = np.linspace(0, vol_t_eq*2, 100)    
    ys = s(xs)
    
    plt.style.use("seaborn-whitegrid")
    fig, ax = plt.subplots() 
    ax.plot(x,y, label = "simulierte Werte bei Titration", marker = "p", linestyle = "") 

    ax.plot(xs,ys, label = "modellierte Kurve")

    ax.set_xlabel("Volumen des Titranten in ml") 
    ax.set_ylabel("pH-Wert bei diesem Verhältnis") 
    ax.set_title("Titrationskurve, schwache Säure: Analyt, starke Base: Titrant") 
    ax.legend(frameon = True)

    spl = UnivariateSpline(x, y)
    spl.set_smoothing_factor(0.5)

    print('Nicht zufrieden mit den Ergebnissen?')
    print('Passen Sie den s-Wert im Skript an. s = UnivariateSpline(  ,  , s ) ')
    print('Titration von schwacher Säure mit starker Base')

print('\n\n BERECHNUNG ERFOLGT, Öffnen Sie das Fenster mit der Titrationskurve\n')

time.sleep(0.5)
plt.show()
# ENDE
#Fragen, Anregungen, Kritik?
#[emanuel.schaerer@edu.sh.ch]
