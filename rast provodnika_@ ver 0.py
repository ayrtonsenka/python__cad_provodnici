import numpy as np
import sympy as sy
import scipy as scp
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from tkinter import *
from tkinter.ttk import *
from sympy.solvers.solveset import solveset, solveset_real, S
from scipy import optimize
from colorama import Fore, Style, init
import builtins
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
import pandas as pd
import math
import sys


from pandas import DataFrame

def ln(broj):
    return(np.log(broj))

def zaokr(broj1):
    return (np.round(broj1,3))

def ipis_cmplx(broj):
    broj=zaokr(broj)
    return (str(broj)+"="+str( zaokr(np.absolute(broj)))+" (< "+str(np.round(np.angle(broj, deg=True),1))+" deg)"+"\n")

def paralel_impedance(b1,b2):
    return(b1*b2)/(b1+b2)

def cmplx_broj_std_ispis (broj):
    redeo=zaokr(np.real(broj))
    imdeo=zaokr(np.imag(broj))
    if imdeo<0:
        znak=" -"
    else:
        znak=" +"
    return (str(redeo)+znak+"j "+str(np.absolute(imdeo)) +" --->  "+str( zaokr(np.absolute(broj)))+" (< "+str(np.round(np.angle(broj, deg=True),1))+" deg)")

def y_lanc (x,p):   # jednacina lancanice
    return p*np.cosh(x/p)

def a_dod(p,a,h): # dodatni raspon pri kosom rasponu
    return 2*p*np.arcsinh(h/2/p/np.sinh(a/2/p))

def a_tot(p,a,h):
    return a+abs(a_dod(p,a,h))

def  f_x(x,p): # ugib u totalnom rasponu
    return y_lanc(x,p)-p

def f_max (p,a,h): # max ugib u kosom rasponu
    return f_x(a_tot(p,a,h)/2,p)

def x_max (p,a,h): # koordinata u kojoj je ugib max
    return p*np.arcsinh(h/a)

def ugao_otklona(diametar, pritisak_vetra,poduzna_tezina):
    return np.rad2deg(np.arctan(diametar/1000*pritisak_vetra*1/poduzna_tezina))   

def D_potrebno (koef_k, L_lanca,f,sr_sklopni): # f i L u cm moraju biti
    return (koef_k*np.sqrt(f+L_lanca)+sr_sklopni)

def koef_provodnika (distanca_TRANS,distanca_VERT,ugao_vetra, sr_sklopni):
    k_niz=[]
    slovo_niz=[]
    for elementT in distanca_TRANS:      
        for elementV in distanca_VERT:  
            if (elementV==0): #horizontalni
                k=builtins.max(6,4+ugao_vetra/25)
                slovo="HORZ"
            else:   
                if np.abs(elementT*100)<sr_sklopni: #vertikalni
                    k=builtins.max(14,4+ugao_vetra/5)
                    slovo="VERT"
                if np.abs(elementT*100)>sr_sklopni: #kosi
                    k=builtins.max(7,2+ugao_vetra/10)
                    slovo="KOSI"
        k_niz=np.append(k_niz,k) 
        slovo_niz=np.append(slovo_niz,slovo) 
    return(k_niz,slovo_niz) 

pi=np.pi 
mi_nula=4*pi*1e-7
f=50 # Hz
w=2*pi*f #omega
epsilon_nula=8.854*1e-12
imag_j=0+1j

provodnici = {
    "Al/Če 120/20****26/7": (0.4846, 15.5),
    "Al/Če 150/25****26/7": (0.5935, 17.1),
    "Al/Če 240/40****26/7": (0.9682, 21.9),
    "Al/Če 360/57****26/19": (1.4175, 26.6),
    "Al/Če 360/57****26/7": (1.4146, 26.4),
    "Al/Če 490/65****54/7": (1.8305, 30.6),
    "Al/Če 505/585****78/91": (6.091, 43.1),
    "Al/Če 75/80****18/19": (0.8181, 16.1),
    "Al/Če 95/15****26/7": (0.3757, 13.6),
    "awg 126.1": (0.826, 14.5),
    "Če 150**III****19": (1.176, 15.8),
    "Če 35**III****19": (0.27958, 7.8),
    "Če 35**III****70": (0.2668, 7.5),
    "Če 50**III****7": (0.3835, 9.0),
    "Če 70**III****19": (0.513, 10.5),
    "Draka tip D": (0.33046, 10.0),
    "Draka tip G": (0.226611, 8.3),
    "OPGW tip E": (0.38945, 12.6),
    "OPGW tp B": (0.572904, 15.0)
}


def unos_podataka():
    def submit():
        naziv_provere = entry_naziv.get()

        prov1 = {
            "x1": float(entry_x1.get()),
            "z1": float(entry_z1.get()),
            "L_lanca_t1": float(entry_L_lanca1_t1.get()),
            "x2": float(entry_x2.get()),
            "z2": float(entry_z2.get()),
            "L_lanca_t2": float(entry_L_lanca1_t2.get()),
            "L": float(entry_L1.get()),
            "p": float(entry_p1.get()),
            "vetar": float(entry_vetar1.get()),
            "tip": combo_prov1.get(),
            "w": float(entry_w1.get()),
            "d": float(entry_d1.get())
        }
        prov2 = {
            "x1": float(entry_x1_2.get()),
            "z1": float(entry_z1_2.get()),
            "L_lanca_t1": float(entry_L_lanca2_t1.get()),
            "x2": float(entry_x2_2.get()),
            "z2": float(entry_z2_2.get()),
            "L_lanca_t2": float(entry_L_lanca2_t2.get()),
            "L": float(entry_L2.get()),
            "p": float(entry_p2.get()),
            "vetar": float(entry_vetar2.get()),
            "tip": combo_prov2.get(),
            "w": float(entry_w2.get()),
            "d": float(entry_d2.get())
        }
        root.prov1, root.prov2, root.naziv_provere = prov1, prov2, naziv_provere
        root.destroy()

    def on_select(combo, entry_w, entry_d):
        name = combo.get()
        w, d = provodnici[name]
        entry_w.delete(0, "end")
        entry_w.insert(0, str(w))
        entry_d.delete(0, "end")
        entry_d.insert(0, str(d))

    root = tk.Tk()
    root.title("rast provodnika - unos podataka")
    root.geometry("1000x700")   

    canvas = tk.Canvas(root)
    h_scroll = tk.Scrollbar(root, orient="horizontal", command=canvas.xview)
    v_scroll = tk.Scrollbar(root, orient="vertical", command=canvas.yview)

    canvas.configure(xscrollcommand=h_scroll.set, yscrollcommand=v_scroll.set)

    canvas.pack(side="left", fill="both", expand=True)
    h_scroll.pack(side="bottom", fill="x")
    v_scroll.pack(side="right", fill="y")

    inner_frame = ttk.Frame(canvas)
    canvas.create_window((0, 0), window=inner_frame, anchor="nw")
    inner_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))


    style = ttk.Style()
    style.configure("TLabel", font=("Segoe UI", 9))
    style.configure("TEntry", font=("Segoe UI", 9))
    style.configure("TCombobox", font=("Segoe UI", 9))

    lbl_title = ttk.Label(inner_frame, text="Dvant", font=("Segoe UI", 14, "bold"), foreground="#5DADE2")
    lbl_title.grid(column=0, row=0, sticky="w", pady=(0, 10))

    ttk.Label(inner_frame, text="Provodnik 1", font=("Segoe UI", 10, "bold")).grid(column=0, row=1, sticky="w", pady=(10, 5))

    # Red 1
    frame1 = tk.LabelFrame(inner_frame, text="Stub 1", bg="#D6EAF8", font=("Segoe UI", 9, "bold"))
    frame1.grid(column=0, row=2, columnspan=12, pady=5, padx=5, sticky="ew")
    ttk.Label(frame1, text="redukovane konzole [m]").grid(column=0, row=0, padx=8)
    ttk.Label(frame1, text="visina [m]").grid(column=1, row=0, padx=8)
    ttk.Label(frame1, text="dužina izol. lanca [cm]").grid(column=2, row=0, padx=8)
    entry_x1 = ttk.Entry(frame1, width=14); entry_x1.grid(column=0, row=1, padx=8, pady=6)
    entry_z1 = ttk.Entry(frame1, width=14); entry_z1.grid(column=1, row=1, padx=8, pady=6)
    entry_L_lanca1_t1 = ttk.Entry(frame1, width=14); entry_L_lanca1_t1.grid(column=2, row=1, padx=8, pady=6)

    # Red 2
    frame2 = tk.LabelFrame(inner_frame, text="Stub 2", bg="#D6EAF8", font=("Segoe UI", 9, "bold"))
    frame2.grid(column=0, row=3, columnspan=12, pady=5, padx=5, sticky="ew")
    ttk.Label(frame2, text="redukovane konzole [m]").grid(column=0, row=0, padx=8)
    ttk.Label(frame2, text="visina [m]").grid(column=1, row=0, padx=8)
    ttk.Label(frame2, text="dužina izol. lanca [cm]").grid(column=2, row=0, padx=8)
    entry_x2 = ttk.Entry(frame2, width=14); entry_x2.grid(column=0, row=1, padx=8, pady=6)
    entry_z2 = ttk.Entry(frame2, width=14); entry_z2.grid(column=1, row=1, padx=8, pady=6)
    entry_L_lanca1_t2 = ttk.Entry(frame2, width=14); entry_L_lanca1_t2.grid(column=2, row=1, padx=8, pady=6)

    # Red 3
    frame3 = tk.Frame(inner_frame, bg="#D6EAF8")
    frame3.grid(column=0, row=4, columnspan=12, pady=5, padx=5, sticky="ew")
    ttk.Label(frame3, text="dužina raspona [m]").grid(column=0, row=0, padx=12)
    ttk.Label(frame3, text="parametar lančanice [m]").grid(column=1, row=0, padx=12)
    entry_L1 = ttk.Entry(frame3, width=14); entry_L1.grid(column=0, row=1, padx=12, pady=6)
    entry_p1 = ttk.Entry(frame3, width=14); entry_p1.grid(column=1, row=1, padx=12, pady=6)

    # Red 4
    frame4 = tk.Frame(inner_frame, bg="#D6EAF8")
    frame4.grid(column=0, row=5, columnspan=12, pady=5, padx=5, sticky="ew")
    ttk.Label(frame4, text="pritisak vetra [daN/m²]").grid(column=0, row=0, padx=8)
    ttk.Label(frame4, text="provodnik").grid(column=1, row=0, padx=8)
    ttk.Label(frame4, text="w [daN/m]").grid(column=2, row=0, padx=8)
    ttk.Label(frame4, text="d [mm]").grid(column=3, row=0, padx=8)
    entry_vetar1 = ttk.Entry(frame4, width=14); entry_vetar1.grid(column=0, row=1, padx=8, pady=6)
    combo_prov1 = ttk.Combobox(frame4, values=list(provodnici.keys()), width=22); combo_prov1.grid(column=1, row=1, padx=8, pady=6)
    entry_w1 = ttk.Entry(frame4, width=14); entry_w1.grid(column=2, row=1, padx=8, pady=6)
    entry_d1 = ttk.Entry(frame4, width=14); entry_d1.grid(column=3, row=1, padx=8, pady=6)
    combo_prov1.bind("<<ComboboxSelected>>", lambda e: on_select(combo_prov1, entry_w1, entry_d1))

    ttk.Label(inner_frame, text="Provodnik 2", font=("Segoe UI", 10, "bold")).grid(column=0, row=6, sticky="w", pady=(20, 5))

    # Red 1
    frame1b = tk.LabelFrame(inner_frame, text="Stub 1", bg="#D6EAF8", font=("Segoe UI", 9, "bold"))
    frame1b.grid(column=0, row=7, columnspan=12, pady=5, padx=5, sticky="ew")
    ttk.Label(frame1b, text="redukovane konzole [m]").grid(column=0, row=0, padx=8)
    ttk.Label(frame1b, text="visina [m]").grid(column=1, row=0, padx=8)
    ttk.Label(frame1b, text="dužina izol. lanca [cm]").grid(column=2, row=0, padx=8)
    entry_x1_2 = ttk.Entry(frame1b, width=14); entry_x1_2.grid(column=0, row=1, padx=8, pady=6)
    entry_z1_2 = ttk.Entry(frame1b, width=14); entry_z1_2.grid(column=1, row=1, padx=8, pady=6)
    entry_L_lanca2_t1 = ttk.Entry(frame1b, width=14); entry_L_lanca2_t1.grid(column=2, row=1, padx=8, pady=6)

    # Red 2
    frame2b = tk.LabelFrame(inner_frame, text="Stub 2", bg="#D6EAF8", font=("Segoe UI", 9, "bold"))
    frame2b.grid(column=0, row=8, columnspan=12, pady=5, padx=5, sticky="ew")
    ttk.Label(frame2b, text="redukovane konzole [m]").grid(column=0, row=0, padx=8)
    ttk.Label(frame2b, text="visina [m]").grid(column=1, row=0, padx=8)
    ttk.Label(frame2b, text="dužina izol. lanca [cm]").grid(column=2, row=0, padx=8)
    entry_x2_2 = ttk.Entry(frame2b, width=14); entry_x2_2.grid(column=0, row=1, padx=8, pady=6)
    entry_z2_2 = ttk.Entry(frame2b, width=14); entry_z2_2.grid(column=1, row=1, padx=8, pady=6)
    entry_L_lanca2_t2 = ttk.Entry(frame2b, width=14); entry_L_lanca2_t2.grid(column=2, row=1, padx=8, pady=6)

    # Red 3
    frame3b = tk.Frame(inner_frame, bg="#D6EAF8")
    frame3b.grid(column=0, row=9, columnspan=12, pady=5, padx=5, sticky="ew")
    ttk.Label(frame3b, text="dužina raspona [m]").grid(column=0, row=0, padx=12)
    ttk.Label(frame3b, text="parametar lančanice [m]").grid(column=1, row=0, padx=12)
    entry_L2 = ttk.Entry(frame3b, width=14); entry_L2.grid(column=0, row=1, padx=12, pady=6)
    entry_p2 = ttk.Entry(frame3b, width=14); entry_p2.grid(column=1, row=1, padx=12, pady=6)

    # Red 4
    frame4b = tk.Frame(inner_frame, bg="#D6EAF8")
    frame4b.grid(column=0, row=10, columnspan=12, pady=5, padx=5, sticky="ew")
    ttk.Label(frame4b, text="pritisak vetra [daN/m²]").grid(column=0, row=0, padx=8)
    ttk.Label(frame4b, text="provodnik").grid(column=1, row=0, padx=8)
    ttk.Label(frame4b, text="w [daN/m]").grid(column=2, row=0, padx=8)
    ttk.Label(frame4b, text="d [mm]").grid(column=3, row=0, padx=8)
    entry_vetar2 = ttk.Entry(frame4b, width=14); entry_vetar2.grid(column=0, row=1, padx=8, pady=6)
    combo_prov2 = ttk.Combobox(frame4b, values=list(provodnici.keys()), width=22); combo_prov2.grid(column=1, row=1, padx=8, pady=6)
    entry_w2 = ttk.Entry(frame4b, width=14); entry_w2.grid(column=2, row=1, padx=8, pady=6)
    entry_d2 = ttk.Entry(frame4b, width=14); entry_d2.grid(column=3, row=1, padx=8, pady=6)
    combo_prov2.bind("<<ComboboxSelected>>", lambda e: on_select(combo_prov2, entry_w2, entry_d2))

    lbl_naziv = ttk.Label(inner_frame, text="Naziv provere:")
    lbl_naziv.grid(column=0, row=11, pady=(20, 5), sticky="w")
    entry_naziv = ttk.Entry(inner_frame, width=50)
    entry_naziv.grid(column=1, row=11, columnspan=6, pady=(20, 5), sticky="w")

    ttk.Button(inner_frame, text="Potvrdi", command=submit).grid(column=0, row=12, columnspan=12, pady=15)

    root.mainloop()
    return root.prov1, root.prov2, root.naziv_provere


if __name__ == "__main__":
    # --- POZIV GUI ---
    prov1, prov2, naziv_provere = unos_podataka()
    print(f"Naziv provere: {naziv_provere}")
    Calc_name = " Rastojanje provodnika VP .xlsx"
    Rezolocija_tacaka = 1500
    SR_za_sklopni_prenapon = 80  # cm
    broj_stubova = 2

    diam1 = prov1["d"]
    w1 = provodnici[prov1["tip"]][0]        
    vetar1 = prov1["vetar"]
    param1 = prov1["p"]
    point1 = [prov1["x1"], 0, prov1["z1"]]
    point2 = [prov1["x2"], prov1["L"], prov1["z2"]]
    L_lanca_prov_1_tacka_1 = prov1["L_lanca_t1"]
    L_lanca_prov_1_tacka_2 = prov1["L_lanca_t2"]

    diam2 = prov2["d"]
    w2 = provodnici[prov2["tip"]][0]
    vetar2 = prov2["vetar"]
    param2 = prov2["p"]
    point1_2 = [prov2["x1"], 0, prov2["z1"]]
    point2_2 = [prov2["x2"], prov2["L"], prov2["z2"]]
    L_lanca_prov_2_tacka_1 = prov2["L_lanca_t1"]
    L_lanca_prov_2_tacka_2 = prov2["L_lanca_t2"]

    # PROVODNIK 1
    p = param1
    scale_x = 1
    scale_y = 1
    S = [point1[1], point2[1]]
    K = [point1[2], point2[2]]
    STACIONAZE = np.array(S)
    KOTE_VES = np.array(K)
    a1 = point2[1] - point1[1]
    ugao1 = ugao_otklona(diam1, vetar1, w1)
    raspon1 = point2[1] - point1[1]

    for jot in range(broj_stubova-1):
        s1 = STACIONAZE[jot]
        s2 = STACIONAZE[jot+1]
        k1 = KOTE_VES[jot]
        k2 = KOTE_VES[jot+1]
        a = s2 - s1
        h = k1 - k2
        if h >= 0:
            xoo_sis = s1 + a_tot(p, a, h)/2
            xx = np.linspace(-a_tot(p, a, h)/2, -a_tot(p, a, h)/2 + a, Rezolocija_tacaka)
            izlaz_y = f_x(xx, p) + k1 - f_max(p, a, h)
            x1 = -a_tot(p, a, h)/2
            x2 = x1 + a
            xx = xx + xoo_sis
            x2 = x2 + xoo_sis
            x1 = x1 + xoo_sis
            xx = xx * scale_x
            x2 = x2 * scale_x
            x1 = x1 * scale_x
            izlaz_y = izlaz_y * scale_y
        if h < 0:
            xoo_sis = s2 - a_tot(p, a, h)/2
            xx = np.linspace(a_tot(p, a, h)/2 - a, a_tot(p, a, h)/2, Rezolocija_tacaka)
            izlaz_y = f_x(xx, p) + k2 - f_max(p, a, h)
            x2 = a_tot(p, a, h)/2
            x1 = x2 - a
            xx = xx + xoo_sis
            x2 = x2 + xoo_sis
            x1 = x1 + xoo_sis
            xx = xx * scale_x
            x2 = x2 * scale_x
            x1 = x1 * scale_x
            izlaz_y = izlaz_y * scale_y
        xoo1 = xoo_sis
        y_low = np.min(izlaz_y)
        k1 = k1 * scale_y
        k2 = k2 * scale_y
        stubx1 = [x1, x1]
        stuby1 = [k1, k1-20]
        stubx2 = [x2, x2]
        stuby2 = [k2, k2-20]
        ugib_prov_1 = izlaz_y - param1
        plt.plot(xx, izlaz_y)
        plt.plot(stubx1, stuby1)
        plt.plot(stubx2, stuby2)
        plt.plot([xoo_sis*scale_x, xoo_sis*scale_x], [y_low+2, y_low-2])
        trans_koo_1 = (xx - point1[1]) / a1 * (point2[0] - point1[0]) + point1[0]
        xy_koordinate_lanc = np.vstack((xx, izlaz_y)).T
        xy_stubova = np.vstack((xx, izlaz_y)).T
        xyz_koordinate_prov1 = np.vstack((xx, trans_koo_1, izlaz_y,)).T

    # PROVODNIK 2
    p = param2
    scale_x = 1
    scale_y = 1
    S = [point1_2[1], point2_2[1]]
    K = [point1_2[2], point2_2[2]]
    STACIONAZE = np.array(S)
    KOTE_VES = np.array(K)
    a2 = point2_2[1] - point1_2[1]
    ugao2 = ugao_otklona(diam2, vetar2, w2)
    raspon2 = point2_2[1] - point1_2[1]

    for jot in range(broj_stubova-1):
        s1 = STACIONAZE[jot]
        s2 = STACIONAZE[jot+1]
        k1 = KOTE_VES[jot]
        k2 = KOTE_VES[jot+1]
        a = s2 - s1
        h = k1 - k2
        if h >= 0:
            xoo_sis = s1 + a_tot(p, a, h)/2
            xx = np.linspace(-a_tot(p, a, h)/2, -a_tot(p, a, h)/2 + a, Rezolocija_tacaka)
            izlaz_y_2 = f_x(xx, p) + k1 - f_max(p, a, h)
            x1 = -a_tot(p, a, h)/2
            x2 = x1 + a
            xx = xx + xoo_sis
            x2 = x2 + xoo_sis
            x1 = x1 + xoo_sis
            xx = xx * scale_x
            x2 = x2 * scale_x
            x1 = x1 * scale_x
            izlaz_y_2 = izlaz_y_2 * scale_y
        if h < 0:
            xoo_sis = s2 - a_tot(p, a, h)/2
            xx = np.linspace(a_tot(p, a, h)/2 - a, a_tot(p, a, h)/2, Rezolocija_tacaka)
            izlaz_y_2 = f_x(xx, p) + k2 - f_max(p, a, h)
            x2 = a_tot(p, a, h)/2
            x1 = x2 - a
            xx = xx + xoo_sis
            x2 = x2 + xoo_sis
            x1 = x1 + xoo_sis
            xx = xx * scale_x
            x2 = x2 * scale_x
            x1 = x1 * scale_x
            izlaz_y_2 = izlaz_y_2 * scale_y
        y_low = np.min(izlaz_y_2)
        k1 = k1 * scale_y
        k2 = k2 * scale_y
        stubx1 = [x1, x1]
        stuby1 = [k1, k1-20]
        stubx2 = [x2, x2]
        stuby2 = [k2, k2-20]
        xoo2 = xoo_sis
        ugib_prov_2 = izlaz_y_2 - param2
        plt.plot(xx, izlaz_y_2)
        plt.plot(stubx1, stuby1)
        plt.plot(stubx2, stuby2)
        plt.plot([xoo_sis*scale_x, xoo_sis*scale_x], [y_low+2, y_low-2])
        trans_koo_2 = (xx - point1_2[1]) / a2 * (point2_2[0] - point1_2[0]) + point1_2[0]
        xy_koordinate_lanc_2 = np.vstack((xx, izlaz_y_2)).T
        xy_stubova = np.vstack((xx, izlaz_y_2)).T
        xyz_koordinate_prov2 = np.vstack((xx, trans_koo_2, izlaz_y_2,)).T

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.step(xx, trans_koo_1, izlaz_y)
    ax.step(xx, trans_koo_2, izlaz_y_2)
    ax.set_xlabel('L pravac')
    ax.set_ylabel('T pravac')
    ax.set_zlabel('Z pravac')
    ax.set_title(naziv_provere)  # Dodaje se naziv_provere kao naslov 3D grafika

    LONGI_diff = xx - xx
    Z_diff = izlaz_y_2 - izlaz_y
    TRANS_diff = trans_koo_2 - trans_koo_1
    D_total = np.sqrt(TRANS_diff**2 + Z_diff**2)

    koef_provodnika_1, slovo_provodnika = koef_provodnika(TRANS_diff, Z_diff, ugao1, SR_za_sklopni_prenapon)
    koef_provodnika_2, slovo_provodnika = koef_provodnika(TRANS_diff, Z_diff, ugao2, SR_za_sklopni_prenapon)

    L_lanca1 = (xx - point1[1]) / a1 * (L_lanca_prov_1_tacka_2 - L_lanca_prov_1_tacka_1) + L_lanca_prov_1_tacka_1
    vert_koo_1 = (xx - point1[1]) / a1 * (point2[2] - point1[2]) + point1[2]
    prov1_y = param1 * np.cosh((xx - xoo1) / param1) - param1
    ugib1 = param1 * (np.cosh(raspon1/2/param1) - np.cosh((raspon1/2 - xx)/param1))
    D_potrebno_1 = D_potrebno(koef_provodnika_1, L_lanca1, ugib1*100, SR_za_sklopni_prenapon)

    L_lanca2 = (xx - point1_2[1]) / a2 * (L_lanca_prov_2_tacka_2 - L_lanca_prov_2_tacka_1) + L_lanca_prov_2_tacka_1
    vert_koo_2 = (xx - point1_2[1]) / a2 * (point2_2[2] - point1_2[2]) + point1_2[2]
    prov2_y = param2 * np.cosh((xx - xoo2) / param2) - param2
    ugib2 = param2 * (np.cosh(raspon2/2/param2) - np.cosh((raspon2/2 - xx)/param2))
    D_potrebno_2 = D_potrebno(koef_provodnika_2, L_lanca2, ugib2*100, SR_za_sklopni_prenapon)

    D_potrebno_max = np.maximum(D_potrebno_1, D_potrebno_2)
    Razlika1 = D_total*100 - D_potrebno_max

    Totalna_provera = " *** ok ***"
    Razlika2 = []
    for element in Razlika1:
        if element > 0:
            element = "OK"
        else:
            element = "NOT OK"
            Totalna_provera = "*** NOT OK ***"
        Razlika2 = np.append(Razlika2, element)

    df = DataFrame({"stac duž L ose": xx, "Rastojanje provodnika[cm]": (D_total)*100, "rastojanje Z[cm]": Z_diff*100, "rastojanje TRANS[cm]": TRANS_diff*100, "D potrebno[cm]": D_potrebno_max,
                   "L1 iz[cm]": L_lanca1, "L2 iz[cm]": L_lanca2, "ugib 1 [cm]": ugib1*100, "ugib 2 [cm]": ugib2*100, "Provera": Razlika2, "V/H/K ": slovo_provodnika,  "Rezerva Dtot-Dpot ": Razlika1})
    print(df)
    df.to_excel(Calc_name, sheet_name="rastojanja provodnika NOVO", float_format="%2.2f")
    print("Totalna_provera   ", Totalna_provera)

    fig2 = plt.figure()
    BX = fig2.add_subplot(211)
    BX.plot(xx, Razlika1)
    BX.set_xlabel('L pravac [m]')
    BX.set_ylabel('Rezerva D_stvarno - D_potrebno [cm]')

    RazlikaList = list(Razlika1)
    Index_min = RazlikaList.index(np.min(Razlika1))
    Index_min_MIN = Index_min
    Razlika1_min = np.round(np.min(Razlika1), 2)
    x_point = np.round(xx[Index_min], 2)
    y_point = np.round(np.min(Razlika1), 1)
    label_text = f"({x_point}, {y_point})"
    BX.scatter(x_point, y_point, color='red', zorder=5)
    BX.annotate(f"L={x_point} m, Δ={y_point} cm",
                (x_point, y_point),
                textcoords="offset points",
                xytext=(5, 5),
                ha='left',
                color='red')
    BX.grid()

    print("KRITICNA TACKA...........................Minimalno rastojanje [cm] ....", np.round(np.min(Razlika1_min),2), " na stacionazi [m]....", np.round(xx[Index_min],2), "@ index...",Index_min)
    print(f"KRITICNA TAČKA...........L={np.round(xx[Index_min],2)} [m], Dpotrebno={np.round(D_potrebno_max[Index_min],2)} [cm], Dstvarno={np.round(D_total[Index_min]*100,2)} [cm] , Raspored - {slovo_provodnika[Index_min]}"  )

    Rast_do_sred_raspona = np.abs(np.max(xx)/2 - xx)
    RazlikaList = list(Rast_do_sred_raspona)
    Index_min = RazlikaList.index(np.min(Rast_do_sred_raspona))
    print("Index sredine raspona (najbliza tacka)  ......", Index_min)
    print("SREDINA RASPONA", np.max(xx)/2, "\n Ukupno rastojanje ", D_total[Index_min]*100, "\n Potrebno rastojanje ", D_potrebno_max[Index_min], "\n Vrsta ", slovo_provodnika[Index_min], "izracunata sredina---", xx[Index_min])

    Ddiff = D_potrebno_max
    XX1 = raspon1/2
    idx_above = np.searchsorted(xx, XX1)
    idx_below = idx_above - 1
    x0, x1 = xx[idx_below], xx[idx_above]
    y0, y1 = D_total[idx_below], D_total[idx_above]
    DDIFF1 = y0 + (y1 - y0) * (XX1 - x0) / (x1 - x0)
    Inerpol_D_stvarno = np.round(DDIFF1*100, 1)
    x0, x1 = xx[idx_below], xx[idx_above]
    y0, y1 = D_potrebno_max[idx_below], D_potrebno_max[idx_above]
    DDIFF1 = y0 + (y1 - y0) * (XX1 - x0) / (x1 - x0)
    Inerpol_D_potrebno = np.round(DDIFF1, 1)
    print(f"INTERPOLACIJA SREDINA RASPONA........... L={XX1} [m], Dpotrebno={Inerpol_D_potrebno} [cm], Dstvarno={Inerpol_D_stvarno} [cm] , Raspored - {slovo_provodnika[idx_below]}"  )

    np.round(np.min(Razlika1),2)
    Inter_Raslika_1 = Inerpol_D_stvarno - Inerpol_D_potrebno
    x_point = np.round(XX1,2)
    y_point = np.round(np.min(Inter_Raslika_1),1)
    label_text = f"({x_point}, {y_point})"
    BX.scatter(x_point, y_point, color='blue', zorder=5)
    BX.annotate(f" ↓ L={x_point} m, Δ={y_point} cm",
                (x_point, y_point),
                textcoords="offset points",
                xytext=(-8,15),
                ha='left',
                color='blue')

    lines = []
    lines.append("\n" * 5)
    lines.append("-" * 100)
    lines.append("-" * 100)
    lines.append(f"{'OPIS TAČKE':<30} {'L [m]':<10} {'D potrebno [cm]':<18} {'D stvarno [cm]':<18} {'Raspored':<12} {'Ocena':<10}")
    lines.append("-" * 100)
    lines.append(f"{'KRITIČNA TAČKA':<30} "
                 f"{np.round(xx[Index_min_MIN],2):<10} "
                 f"{np.round(D_potrebno_max[Index_min_MIN],1):<18} "
                 f"{np.round(D_total[Index_min_MIN]*100,1):<18} "
                 f"{slovo_provodnika[Index_min_MIN]:<12}"
                 f"{Razlika2[Index_min_MIN]:<10}")
    lines.append("-" * 100)
    lines.append(f"{'I.SREDINA RASPONA':<30} "
                 f"{XX1:<10} "
                 f"{Inerpol_D_potrebno:<18} "
                 f"{Inerpol_D_stvarno:<18} "
                 f"{slovo_provodnika[idx_below]:<12}"
                 f"{Razlika2[Index_min]:<10}")
    lines.append("-" * 100)
    text_block = "\n".join(lines)
    CX = fig2.add_subplot(212)
    CX.axis("off")
    CX.text(0.0, 0.5, text_block, fontsize=7,
            ha="left", va="center", family="monospace")
    plt.savefig("BX.png", dpi=300, bbox_inches="tight")
    plt.show()







