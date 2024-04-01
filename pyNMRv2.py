import os
import tkinter as gui
from tkinter import ttk, font  # , Button, filedialog
import time
import webbrowser
from sys import platform as PLATFORM, executable as sys_executable
from math import ceil as Myceil
from pathlib import Path
from datetime import datetime
import threading
from PIL import Image as Image2, ImageTk
import numpy as np
from matplotlib.backend_bases import MouseButton
import keyboard
import json

# from queue import Queue

from BioSANS2020.test_data import test_data2 as test_data2
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.myglobal import proc_global
from BioSANS2020.gui_functs.prepare_canvas import prepare_frame_for_plot
from BioSANS2020.prepcodes.process import process
from BioSANS2020.analysis.numeric.transform_data import calc_cross_corr, \
    calc_covariance, prob_density_calc, \
    prob_density_calc2, prob_density_calc3, ave_traj_calc

from BioSANS2020.model import topology_view
from BioSANS2020.model.new_file import new_file
from BioSANS2020.gui_functs.draw_figure import draw_figure

import matplotlib.pylab as plt
import glob
import nmrglue as ng
import pandas as pd
from pandastable import Table, TableModel#, DataExplorer
from scipy.signal import find_peaks, detrend
from matplotlib.animation import FuncAnimation 
from matplotlib.widgets import Button as mpl_button

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
    NavigationToolbar2Tk as NavigationToolbar2TkAgg

if PLATFORM == "win32":
    from subprocess import Popen, CREATE_NEW_CONSOLE
elif PLATFORM == "darwin":
    try:
        from applescript import tell as my_tell_us
    except BaseException:
        pass
    from subprocess import check_output, call as my_call_us
elif PLATFORM == "linux":
    pass
else:
    from subprocess import Popen

try:
    import tempfile
    TEMPORARY_FOLDER = str(tempfile.gettempdir()).replace(
        "\\", "/") + "/BioSANS_temporary_folder"
except BaseException:
    TEMPORARY_FOLDER = "BioSANS_temporary_folder"

globals2.init(globals2)
if __name__ == '__main__':
    proc_global.init(proc_global)

TOP = gui.Tk()
TOP.title("pH-Intelligenceᴰᴱ")
TOP.geometry("1005x550")
#TOP.resizable(False, False)
#TOP.geometry(
#    "%dx%d+0+0" % (TOP.winfo_screenwidth(),  TOP.winfo_screenheight())
#)
TOP.state("zoomed")

HEADER = gui.Label(TOP, text="pH-Intelligenceᴰᴱ") #"ᴬᴮᶜᴰᴱᶠᴳ
HEADER.configure(
    bg="#f4bbff",
    fg="#a2006d",
    height=1,
    # width = 1005,
    font="Helvetica 18 bold italic"
)
HEADER.pack(fill="x")

FRAME = gui.Frame(TOP)
FRAME.configure(
    bg="light cyan",
    borderwidth=2,
    height=500,
    width=1005
)
FRAME.pack(fill="both", expand=True)

FOOTER = gui.Label(
    TOP, text="Powered by BioSANS2020")
FOOTER.configure(
    bg="#f4bbff",
    fg="#a2006d",
    # width = 1005,
    font="Helvetica 10 bold italic",
    anchor='w'
)
FOOTER.pack(fill="x")

FILE_NAME = {"current_folder": None}
CURRENT_DATA = None

from BioSANS2020.BioSANS import load_data 
from BioSANS2020.BioSANS import delete_this     
from BioSANS2020.BioSANS import canvas_update_widgets 
from BioSANS2020.BioSANS import TOP as topB
topB.withdraw()

def on_closing():
    topB.destroy()
    TOP.destroy()
        
def load_image(wdata=False):
    t_o = time.time()
    global CURRENT_DATA
    canvas, scroll_x, scroll_y = ITUPS
    file = gui.filedialog.askopenfilename(title="Select file")
    FILE_NAME["current_folder"] = file
    load = Image2.open(file)
    render = ImageTk.PhotoImage(load)
    img = gui.Label(canvas, image=render)
    img.image = render

    fframe = canvas.create_window(
        0, 426 * globals2.PLOT_I, anchor='nw', window=img)
    but_b = gui.Button(img, text=" X ", fg='red', highlightcolor='blue',
                       bg='white', height=1, relief='raised',
                       command=lambda: delete_this(fframe, canvas))
    but_b.place(rely=0.0, relx=1.0, x=-15, y=0, anchor="ne")

    canvas.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
    canvas.configure(scrollregion=canvas.bbox("all"))
    canvas.bind("<Configure>", lambda e: canvas_update_widgets(e, canvas))
    globals2.PLOT_I = globals2.PLOT_I + 1

    print(time.time() - t_o)
    canvas_update_widgets(None, canvas)    
    
def plot_pdata(data, slabels, items, plotted, mix_plot=True, logx=False,
              logy=False, normalize=False, si_ticked=None):

    miter = len(data)
    col = ['C' + str(i) for i in range(100)]
    if si_ticked is None:
        si_ticked = range(len(slabels))
    if mix_plot:
        plt.close()
        plt.figure(figsize=(9.68, 3.95))
        plt.xlabel("ppm")
        plt.ylabel("intensity")
        lines = []
        if logx:
            plt.xscale('log')
        if logy:
            plt.yscale('log')
        for j in range(miter):

            for i in si_ticked:
                if normalize:
                    line = plt.plot(
                        data[j][0], data[j][1]
                        / (npmax(data[j][1]) + 1.0e-30), col[j])
                else:
                    line = plt.plot(data[j][0], data[j][1], col[j])
        plt.legend([slabels[i] for i in si_ticked])
        plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        plt.tight_layout()
        plt.xlim(25, -5)
        fig = plt.gcf()
        canvas = items[0]
        #canvas_height = 0.75*TOP.winfo_screenmmheight()/10/2.54
        #canvas_width = 0.90*TOP.winfo_screenmmwidth()/10/2.54
        fig.set_size_inches(0.8*canvas.winfo_width()/80, 0.74*canvas.winfo_height()/80, forward=True)
        lines.append(line)
        plotted.append([plt.gca(), fig, lines])
        # fig_canvas_agg = draw_figure(items, fig)
        draw_figure(items, fig)
        #canvas = items[0]
        #canvas_update_widgets(None, canvas)
    else:
        lines = []
        figs = []
        for i in range(len(slabels)):
            plt.close()
            plt.figure(figsize=(9.68, 3.95))
            plt.gca().ticklabel_format(axis='y', style='sci')
            plt.xlabel("ppm")
            plt.ylabel("intensity")
            if logx:
                plt.xscale('log')
            if logy:
                plt.yscale('log')
            #for j in range(miter):
            j = i
            if normalize:
                line = plt.plot(data[j][0], data[j][1]/(npmax(data[j][1])
                                   + 1.0e-30), col[j])
            else:
                line = plt.plot(data[j][0], data[j][1], col[j])
            lines.append(line)

            plt.legend([slabels[i] for i in si_ticked])
            plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            plt.tight_layout()
            plt.xlim(25, -5)
            fig = plt.gcf()
            canvas = items[0]
            #canvas_height = 0.75*TOP.winfo_screenmmheight()/10/2.54
            #canvas_width = 0.90*TOP.winfo_screenmmwidth()/10/2.54
            fig.set_size_inches(0.8*canvas.winfo_width()/80, 0.74*canvas.winfo_height()/80, forward=True)
            figs.append(fig)
            plotted.append([plt.gca(), fig, lines])
            draw_figure(items, fig)
        canvas = items[0] 
        canvas_update_widgets(None, canvas)
    
def load_pdata(plot=False, def_data=None):
    global CURRENT_DATA
    CURRENT_DATA = None
    if not def_data:
        data = []
        slabels = []
        folder = gui.filedialog.askdirectory()
        nmr_files = list(Path(folder).rglob("pdata/1"))
        if not nmr_files:
            nmr_files = list(Path(folder).parent.rglob("pdata/1"))
        if not nmr_files:
            nmr_files = list(Path(folder).parent.parent.rglob("pdata/1"))
        c = 1
        for nmr_path in nmr_files:
            numr_file = str(nmr_path)
            try:
                key = "_".join(numr_file.split("\\")[-4:])
            except:
                try:
                    key = "_".join(numr_file.split("/")[-4:])
                except:
                    key = "None"
            label = key.split("-")[0]
            slabels.append(label)
            
            try:
                data_dir = numr_file
                dic, dta = ng.bruker.read_pdata(data_dir, scale_data=True)    
                udic = ng.bruker.guess_udic(dic, dta)
                uc = ng.fileiobase.uc_from_udic(udic)
                data.append([uc.ppm_scale(), dta])
                c = c + 1
            except:
                pass
    else:
        data, slabels = def_data
    if plot:
        plot_pdata(data, slabels, ITUPS, globals2.PLOTTED, mix_plot=True)
    CURRENT_DATA = (data, slabels)
    if not plot:
        gui.messagebox.showinfo("showinfo", "Trajectory loaded succesfully")
        #print(CURRENT_DATA)
        
def load_jcampdx(plot=False, def_data=None):
    global CURRENT_DATA
    data = []
    slabels = ["no-label"]
    
    file = gui.filedialog.askopenfilename(title="Select file")
    dic,dta = ng.jcampdx.read(file)
    udic = ng.jcampdx.guess_udic(dic, dta)
    uc = ng.fileiobase.uc_from_udic(udic)
    ref = int(dic['.SHIFTREFERENCE'][0].split(",")[2].split(".")[0])-1
    ppm_scale = np.array(uc.ppm_scale())
    ppm_scale = ppm_scale - ppm_scale[ref]
    data.append([ppm_scale, dta])
    if plot:
        plot_pdata(data, slabels, ITUPS, globals2.PLOTTED, mix_plot=True)
    CURRENT_DATA = (data, slabels)
    if not plot:
        gui.messagebox.showinfo("showinfo", "Trajectory loaded succesfully")
        #print(CURRENT_DATA)
        
def load_csv(ITUPS,show=True):
    file = gui.filedialog.askopenfilename(title="Select file")
    df = pd.read_csv(file,sep="\t|,")
    CURRENT_DATA = df
    if show:
        canvas, _, scroll_y = ITUPS
        frame = gui.Frame(canvas, height=455, width=940, borderwidth=10, bd=0)
        pt = Table(frame, dataframe=df)
        frame.pack(side="top", fill="both", expand=True)
        pt.show()
        #explorer = DataExplorer(frame)
        #explorer.show()
        #explorer.connect(pt)
        
        fframe = canvas.create_window(
            0, 450 * globals2.PLOT_I, anchor='nw', window=frame)
        bttn = gui.Button(frame, text=" X ", fg='red', highlightcolor='blue',
                      bg='white', height=1, relief='raised',
                      command=lambda: delete_this(fframe, canvas))
        bttn.place(rely=0.0, relx=1.0, x=-15, y=0, anchor="ne")
        #globals2.CONTAINER.append([canvas, frame])
        #canvas.configure(yscrollcommand=scroll_y.set,
                         #xscrollcommand=scroll_x.set)
        #canvas.configure(scrollregion=canvas.bbox("all"))
        #globals2.PLOT_I = globals2.PLOT_I + 1
        canvas.bind("<Configure>", lambda e: canvas_update_widgets(e, canvas))
        canvas_update_widgets(None, canvas)
    
def load_xls(ITUPS, show=True):
    file = gui.filedialog.askopenfilename(title="Select file")
    df = pd.read_excel(file,sheet_name=0,header=1)
    CURRENT_DATA = df
    if show:
        canvas, _, scroll_y = ITUPS
        frame = gui.Frame(canvas, height=455, width=940, borderwidth=10, bd=0)
        pt = Table(frame, dataframe=df,showtoolbar=True, showstatusbar=True)
        frame.pack(side="top", fill="both", expand=True)
        pt.show()
        #explorer = DataExplorer(frame)
        #explorer.show()
        #explorer.connect(pt)    
        
        fframe = canvas.create_window(
            0, 450 * globals2.PLOT_I, anchor='nw', window=frame)
        bttn = gui.Button(frame, text=" X ", fg='red', highlightcolor='blue',
                      bg='white', height=1, relief='raised',
                      command=lambda: delete_this(fframe, canvas))
        bttn.place(rely=0.0, relx=1.0, x=-15, y=0, anchor="ne")
        #globals2.CONTAINER.append([canvas, frame])
        #canvas.configure(yscrollcommand=scroll_y.set,
                         #xscrollcommand=scroll_x.set)
        #canvas.configure(scrollregion=canvas.bbox("all"))
        #globals2.PLOT_I = globals2.PLOT_I + 1
        canvas.bind("<Configure>", lambda e: canvas_update_widgets(e, canvas))
        canvas_update_widgets(None, canvas)
        
def is_float(x):
    try:
        float(x)
        return True
    except:
        return False

def read_pH_table_and_jdx_files():
    ph_tab = open("pH_and_conc_shift_Table.txt","r")
    metab = {}
    last = "None"
    for row in ph_tab:
        now = row.strip()
        if now:
            col = now.split("\t")
            if len(col) == 1:
                if last in metab:
                    metab[last] = np.array(metab[last]) 
                last = col[0].strip()
                metab[last] = []
            elif len(col) == 5:
                if is_float(col[0]):
                    metab[last].append([float(col[i]) for i in range(5)])
    metab[last] = np.array(metab[last])
    
    jdx_files = glob.glob("*.jdx")
    jdx = {}
    for row in jdx_files:
        col = row.split("_")
        jdx[col[0].strip()] = row

    return metab, jdx
    
metab_dict, jdx_dict = read_pH_table_and_jdx_files()
    
def process_metab_shifts_multiple(metabs, met):

    metab_table = metab_dict[met]
    xi = metab_table[:,2].mean()    
    xv = metab_table[:,2].var()**0.5


    if CURRENT_DATA:
        query = CURRENT_DATA[0][0]
        new_ppm, new_dta = query
        base_line = 0.85*np.mean(new_dta)
        new_peaks, _ = find_peaks(new_dta,  height=base_line*3)    
        idx = new_peaks[-1]
        q1 = 0
        while new_dta[idx-q1]>base_line:
            q1 = q1 + 1
        q2 = 0
        try:
            while new_dta[idx+q2]>base_line:
                q2 = q2 + 1    
        except:
            q2 = len(new_dta)-1
        
        new_dta = detrend(new_dta)
        area_std = np.trapz(new_dta[idx-q1:idx+q2])/12          
            
        scope = (new_ppm<xi+8*xv)*(new_ppm>xi-8*xv)
        new_ppm = new_ppm[scope]
        new_dta = new_dta[scope]
        #new_dta = detrend(new_dta)
        
        new_peaks, _ = find_peaks(new_dta,  height=0.85*max(new_dta))
    
    file = jdx_dict[metabs[met]]
    dic,dta = ng.jcampdx.read(file)
    udic = ng.jcampdx.guess_udic(dic, dta)
    uc = ng.fileiobase.uc_from_udic(udic)
    ref = int(dic['.SHIFTREFERENCE'][0].split(",")[2].split(".")[0])-1
    ppm_scale = np.array(uc.ppm_scale())
    ref_ppm = ppm_scale - ppm_scale[ref]
    ref_dta = dta
    ref_dta = detrend(ref_dta)
    
    scope = (ref_ppm<xi+8*xv)*(ref_ppm>xi-8*xv)
    ref_ppm = ref_ppm[scope]
    ref_dta = ref_dta[scope]
    ref_peaks, _ = find_peaks(ref_dta,  height=0.85*max(ref_dta))
    
    lines = []
    plt.close()
    plt.figure(figsize=(9.68, 3.95))
    plt.gca().ticklabel_format(axis='y', style='sci')
    plt.xlabel("ppm")
    plt.ylabel("Intensity(a.u.)")    
    
    if CURRENT_DATA:
        line = plt.plot(new_ppm, new_dta,label="new")
        lines.append(line)
        
    idx1 = (np.abs(ref_ppm[ref_peaks] - metab_table[3,2])).argmin()
    f = 1
    diffo = 1000
    diff = 0
    if CURRENT_DATA:
        for ih in range(6):
            ph_factor = metab_table[ih,3]
            idx2 = (np.abs(new_ppm[new_peaks] - metab_table[ih,2])).argmin()
            diff = new_ppm[new_peaks[idx2]] - metab_table[ih,2]
            if abs(diff)<abs(diffo):
                f = new_dta[new_peaks[idx2]] / ( ref_dta[ref_peaks[idx1]]*ph_factor )
                diffo = diff
        diff = diffo    

    for row in metab_table:
        corrected_ref = ref_dta*row[3]*f
        base_line = 0.99*solve_base_line(corrected_ref)
        idx = ref_peaks[idx1]
        q1 = 0
        while corrected_ref[idx-q1]>base_line:
            q1 = q1 + 1
        q2 = 0
        try:
            while corrected_ref[idx+q2]>base_line:
                q2 = q2 + 1
        except:
            q2 = len(corrected_ref)-1
        
        #corrected_ref = detrend(corrected_ref)
        area_ref = np.trapz(corrected_ref[idx-q1:idx+q2])/row[4]          
        corrected_ppm = ref_ppm-ref_ppm[ref_peaks[idx1]]+row[2]+diff
        
        if CURRENT_DATA:
            conc_ref = np.round(ref_standard_conc*area_ref/area_std,4)
            plt.text(corrected_ppm[ref_peaks[idx1]], corrected_ref[ref_peaks[idx1]]+base_line, str(conc_ref)+" mM")
        line = plt.plot(corrected_ppm, corrected_ref,label=str(round(row[0],2)),ls='--')
        lines.append(line)
    
    plt.xlim(xi+6*xv,xi-6*xv)
    plt.legend()
    plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.tight_layout()
    fig = plt.gcf()
    canvas = ITUPS[0]
    #canvas_height = 0.75*TOP.winfo_screenmmheight()/10/2.54
    #canvas_width = 0.90*TOP.winfo_screenmmwidth()/10/2.54
    fig.set_size_inches(0.8*canvas.winfo_width()/80, 0.74*canvas.winfo_height()/80, forward=True)
    lines.append(line)
    globals2.PLOTTED.append([plt.gca(), fig, lines])
    fig_canvas_agg = draw_figure(ITUPS, fig)
    #draw_figure(ITUPS, fig)
    
def process_metab_shifts_single(metabs, met):

    metab_table = metab_dict[met]
    xi = metab_table[:,2].mean()    
    xv = metab_table[:,2].var()**0.5


    if CURRENT_DATA:
        query = CURRENT_DATA[0][0]
        new_ppm, new_dta = query
        base_line = 0.85*np.mean(new_dta)
        #print(base_line)
        new_peaks, _ = find_peaks(new_dta,  height=base_line*3)    
        idx = new_peaks[-1]
        q1 = 0
        while new_dta[idx-q1]>base_line:
            q1 = q1 + 1
        q2 = 0
        try:
            while new_dta[idx+q2]>base_line:
                q2 = q2 + 1
        except:
            q2 = len(new_dta)-1
            
        new_dta = detrend(new_dta)
        area_std = np.trapz(new_dta[idx-q1:idx+q2])/12      
            
        scope = (new_ppm<xi+8*xv)*(new_ppm>xi-8*xv)
        new_ppm = new_ppm[scope]
        new_dta = new_dta[scope]
        #new_dta = detrend(new_dta)
        new_peaks, _ = find_peaks(new_dta,  height=0.85*max(new_dta))
    
    file = jdx_dict[metabs[met]]
    dic,dta = ng.jcampdx.read(file)
    udic = ng.jcampdx.guess_udic(dic, dta)
    uc = ng.fileiobase.uc_from_udic(udic)
    ref = int(dic['.SHIFTREFERENCE'][0].split(",")[2].split(".")[0])-1
    ppm_scale = np.array(uc.ppm_scale())
    ref_ppm = ppm_scale - ppm_scale[ref]
    ref_dta = dta
    ref_dta = detrend(ref_dta)
    
    scope = (ref_ppm<xi+8*xv)*(ref_ppm>xi-8*xv)
    ref_ppm = ref_ppm[scope]
    ref_dta = ref_dta[scope]
    ref_peaks, _ = find_peaks(ref_dta,  height=0.85*max(ref_dta))
    
    lines = []
    plt.close()
    plt.figure(figsize=(9.68, 3.95))
    plt.gca().ticklabel_format(axis='y', style='sci')
    plt.xlabel("ppm")
    plt.ylabel("Intensity(a.u.)")    
    
    if CURRENT_DATA:
        line = plt.plot(new_ppm, new_dta,label="new")
        lines.append(line)
        
    idx1 = (np.abs(ref_ppm[ref_peaks] - metab_table[3,2])).argmin()
    f = 1
    diffo = 1000
    diff = 0
    closest_ppm_idx = 0
    if CURRENT_DATA:
        for ih in range(6):
            ph_factor = metab_table[ih,3]
            idx2 = (np.abs(new_ppm[new_peaks] - metab_table[ih,2])).argmin()
            diff = new_ppm[new_peaks[idx2]] - metab_table[ih,2]
            if abs(diff)<abs(diffo):
                f = new_dta[new_peaks[idx2]] / ( ref_dta[ref_peaks[idx1]]*ph_factor )
                diffo = diff
                closest_ppm_idx = ih
        diff = diffo    

    row = metab_table[closest_ppm_idx]
    corrected_ref = ref_dta*row[3]*f
    base_line = 0.99*solve_base_line(corrected_ref)

    idx = ref_peaks[idx1]
    q1 = 0
    while corrected_ref[idx-q1]>base_line:
        q1 = q1 + 1
    q2 = 0
    try:
        while corrected_ref[idx+q2]>base_line:
            q2 = q2 + 1
    except:
        q2 = len(corrected_ref)-1
    
    #corrected_ref = detrend(corrected_ref)
    area_ref = np.trapz(corrected_ref[idx-q1:idx+q2])/row[4]       
    corrected_ppm = ref_ppm-ref_ppm[ref_peaks[idx1]]+row[2]+diff
    
    int_line1 = plt.axvline(corrected_ppm[idx-q1],ls='--',color='red',alpha=0.1)
    try:
        int_line2 = plt.axvline(corrected_ppm[idx+q2],ls='--',color='red',alpha=0.1)
    except:
        int_line2 = plt.axvline(corrected_ppm[-1],ls='--',color='red',alpha=0.1)
    
    if CURRENT_DATA:
        conc_ref = np.round(ref_standard_conc*area_ref/area_std,4)
        plt.text(corrected_ppm[ref_peaks[idx1]], corrected_ref[ref_peaks[idx1]]+base_line, str(conc_ref)+" mM")
    line = plt.plot(corrected_ppm, corrected_ref,label=str(round(row[0],2)),ls='--')
    lines.append(line)
    
    plt.xlim(xi+6*xv,xi-6*xv)
    plt.legend()
    plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.tight_layout()
    fig = plt.gcf()
    canvas = ITUPS[0]
    #canvas_height = 0.75*TOP.winfo_screenmmheight()/10/2.54
    #canvas_width = 0.90*TOP.winfo_screenmmwidth()/10/2.54
    fig.set_size_inches(0.8*canvas.winfo_width()/80, 0.74*canvas.winfo_height()/80, forward=True)
    lines.append(line)
    globals2.PLOTTED.append([plt.gca(), fig, lines])
    fig_canvas_agg = draw_figure(ITUPS, fig)
    #draw_figure(ITUPS, fig)
    
def solve_base_line(ref):
    ref_sort = np.sort(ref)
    pos = int(len(ref_sort)/2)
    k = 50
    while k>0:
        try:
            return ref_sort[pos-k:pos+k].mean()
        except:
            k = k - 1
    
def process_metab_shifts_manual(metabs, met):
    global idx_last, base_line, ref_peaks
    idx_last = None
    metab_table = metab_dict[met]
    xi = metab_table[:,2].mean()    
    xv = metab_table[:,2].var()**0.5


    if CURRENT_DATA:
        query = CURRENT_DATA[0][0]
        new_ppm, new_dta = query
        base_line = 0.85*np.mean(new_dta)
        new_peaks, _ = find_peaks(new_dta,  height=base_line*3)    
        idx = new_peaks[-1]
        q1 = 0
        while new_dta[idx-q1]>base_line:
            q1 = q1 + 1
        q2 = 0
        try:
            while new_dta[idx+q2]>base_line:
                q2 = q2 + 1
        except:
            q2 = len(new_dta)-1                
        new_dta = detrend(new_dta)
        area_std = np.trapz(new_dta[idx-q1:idx+q2])/12          
            
        scope = (new_ppm<xi+8*xv)*(new_ppm>xi-8*xv)
        new_ppm = new_ppm[scope]
        new_dta = new_dta[scope]
        #new_dta = detrend(new_dta)
        new_peaks, _ = find_peaks(new_dta,  height=0.85*max(new_dta))
    
    file = jdx_dict[metabs[met]]
    dic,dta = ng.jcampdx.read(file)
    udic = ng.jcampdx.guess_udic(dic, dta)
    uc = ng.fileiobase.uc_from_udic(udic)
    ref = int(dic['.SHIFTREFERENCE'][0].split(",")[2].split(".")[0])-1
    ppm_scale = np.array(uc.ppm_scale())
    ref_ppm = ppm_scale - ppm_scale[ref]
    ref_dta = dta
    ref_dta = detrend(ref_dta)
    
    scope = (ref_ppm<xi+8*xv)*(ref_ppm>xi-8*xv)
    ref_ppm = ref_ppm[scope]
    ref_dta = ref_dta[scope]
    ref_peaks, _ = find_peaks(ref_dta,  height=0.85*max(ref_dta))
    
    lines = []
    
    plt.close()
    plt.figure(figsize=(9.68, 3.95))
    plt.gca().ticklabel_format(axis='y', style='sci')
    plt.xlabel("ppm")
    plt.ylabel("Intensity(a.u.)")    
    
    if CURRENT_DATA:
        line_new, = plt.plot(new_ppm, new_dta,label="new")
        lines.append(line_new)
        
    idx1 = (np.abs(ref_ppm[ref_peaks] - metab_table[3,2])).argmin()
    f = 1
    diffo = 1000
    diff = 0
    closest_ppm_idx = 0
    if CURRENT_DATA:
        for ih in range(6):
            ph_factor = metab_table[ih,3]
            idx2 = (np.abs(new_ppm[new_peaks] - metab_table[ih,2])).argmin()
            diff = new_ppm[new_peaks[idx2]] - metab_table[ih,2]
            if abs(diff)<abs(diffo):
                f = new_dta[new_peaks[idx2]] / ( ref_dta[ref_peaks[idx1]]*ph_factor )
                diffo = diff
                closest_ppm_idx = ih
        diff = diffo    

    row = metab_table[closest_ppm_idx]
    corrected_ref = ref_dta*row[3]*f
    base_line = 0.99*solve_base_line(corrected_ref)
    idx = ref_peaks[idx1]
    q1 = 0
    while corrected_ref[idx-q1]>base_line:
        q1 = q1 + 1
    q2 = 0
    try:
        while corrected_ref[idx+q2]>base_line:
            q2 = q2 + 1
    except:
        q2 = len(corrected_ref)-1
    
    #corrected_ref = detrend(corrected_ref)
    area_ref = np.trapz(corrected_ref[idx-q1:idx+q2])/row[4]           
    corrected_ppm = ref_ppm-ref_ppm[ref_peaks[idx1]]+row[2]+diff
    base_line = 0.99*solve_base_line(corrected_ref)
    
    if CURRENT_DATA:
        conc_ref = np.round(ref_standard_conc*area_ref/area_std,4)
        textD = plt.text(corrected_ppm[ref_peaks[idx1]], corrected_ref[ref_peaks[idx1]]+base_line, "conc = "+str(conc_ref)+" mM")
        line, = plt.plot(corrected_ppm, corrected_ref,ls='--')
        lines.append(line)
        idx = ref_peaks[idx1]
        q1 = 0
        while corrected_ref[idx-q1]>base_line:
            q1 = q1 + 1
        q2 = 0
        try:
            while corrected_ref[idx+q2]>base_line:
                q2 = q2 + 1
        except:
            q2 = len(corrected_ref)-1

        int_line1 = plt.axvline(corrected_ppm[idx-q1],ls='--',color='red',alpha=0.1)
        try:
            int_line2 = plt.axvline(corrected_ppm[idx+q2],ls='--',color='red',alpha=0.1)
        except:
            int_line2 = plt.axvline(corrected_ppm[-1],ls='--',color='red',alpha=0.1)
        basis_line = plt.axhline(base_line,ls='--',color='red',alpha=0.1)
        line.set_picker(True)
    
    plt.xlim(xi+6*xv,xi-6*xv)
    plt.legend()
    plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.tight_layout()
    fig = plt.gcf()
    canvas = ITUPS[0]
    #canvas_height = 0.75*TOP.winfo_screenmmheight()/10/2.54
    #canvas_width = 0.90*TOP.winfo_screenmmwidth()/10/2.54
    fig.set_size_inches(0.8*canvas.winfo_width()/80, 0.74*canvas.winfo_height()/80, forward=True)
    
    def on_drag(event):
        global corrected_ppm, current_center, idx_last, base_line, ref_peaks
        #print(fig.canvas.manager.toolbar.mode)
        #print(fig.canvas.manager.toolbar.coordinates)
        #print(fig.canvas.manager.toolbar._zoom_info)
        if not fig.canvas.manager.toolbar.mode:
            corrected_ppm = line.get_xdata()
            corrected_ref = line.get_ydata()
            if event.button == 1:
                current_x = event.xdata            
                idx = 0
                while corrected_ppm[idx]>current_x:
                    idx += 1
                
                d = 1
                a = line.get_ydata()[idx+d:]
                b = line.get_ydata()[idx-d:idx+d]
                b = b + event.ydata - b
                c = line.get_ydata()[:idx-d]
                line.set_ydata([*c,*b,*a])
                
                corrected_ref = line.get_ydata()
                idx = ref_peaks[idx1]
                q1 = 0
                try:
                    while corrected_ref[idx-q1]>base_line:
                        q1 = q1 + 1
                except:
                    q1 = 0
                q2 = 0
                try:
                    while corrected_ref[idx+q2]>base_line:
                        q2 = q2 + 1
                except:
                    q2 = len(corrected_ref)-1
                try:
                    int_line1.set_xdata(corrected_ppm[idx-q1])
                except:
                    int_line1.set_xdata(corrected_ppm[0])
                try:
                    int_line2.set_xdata(corrected_ppm[idx+q2])
                except:
                    int_line2.set_xdata(corrected_ppm[-1])
                    
                area_ref = np.trapz(corrected_ref[idx-q1:idx+q2])/row[4]           

                #line.set_data(corrected_ppm, corrected_ref)
                conc_ref = np.round(ref_standard_conc*area_ref/area_std,4)
                textD.set_text("conc = "+str(conc_ref)+" mM")
                textD.set_x(corrected_ppm[ref_peaks[idx1]])
                textD.set_y(corrected_ref[ref_peaks[idx1]]+base_line)

                fig.canvas.draw()
            elif event.button == 2:
                current_yi = event.ydata    
                idx = 0
                try:
                    while corrected_ppm[idx]>current_yi:
                        idx += 1
                except:
                    idx = len(corrected_ppm)-1
                    
                if not idx_last:
                    idx_last = idx
                current_y = line.get_ydata()
                current_y = current_y + current_yi - current_y[idx_last]
                base_line = base_line + current_yi - current_y[idx_last]
                line.set_ydata(current_y)
                
                base_line = 0.99*solve_base_line(current_y)
                basis_line.set_ydata(base_line)
                fig.canvas.draw()                
                
            elif event.button == 3:
                current_xi = event.xdata    
                idx = 0
                while corrected_ppm[idx]>current_xi:
                    idx += 1
                if not idx_last:
                    idx_last = idx
                current_x = line.get_xdata()
                current_x = current_x + current_xi -  current_x[idx_last]
                corrected_ppm = current_x
                
                idx = ref_peaks[idx1]
                q1 = 0
                try:
                    while corrected_ref[idx-q1]>base_line:
                        q1 = q1 + 1
                except:
                    q1 = 0
                q2 = 0
                try:
                    while corrected_ref[idx+q2]>base_line:
                        q2 = q2 + 1
                except:
                    q2 = len(corrected_ref)-1

                try:
                    int_line1.set_xdata(corrected_ppm[idx-q1])
                except:
                    int_line1.set_xdata(corrected_ppm[0])
                try:
                    int_line2.set_xdata(corrected_ppm[idx+q2])
                except:
                    int_line2.set_xdata(corrected_ppm[i-1])

                line.set_xdata(current_x)
                fig.canvas.draw()
                
    def on_released(event):
        global idx_last
        idx_last = None
        
    def on_scrolled(event):
        global corrected_ppm, current_center, idx_last
        if not fig.canvas.manager.toolbar.mode:
            corrected_ppm = line.get_xdata()
            corrected_ref = line.get_ydata()
            if event.button == "up":
                corrected_ref = corrected_ref*1.05
            elif event.button == "down":
                try:
                    corrected_ref = corrected_ref*0.95
                except:
                    pass
            line.set_ydata(corrected_ref)

            idx = ref_peaks[idx1]
            q1 = 0
            try:
                while corrected_ref[idx-q1]>base_line:
                    q1 = q1 + 1
            except:
                q1 = 0
            q2 = 0
            try:
                while corrected_ref[idx+q2]>base_line:
                    q2 = q2 + 1
            except:
                q2 = len(corrected_ref)-1
            area_ref = np.trapz(corrected_ref[idx-q1:idx+q2])/row[4]           

            conc_ref = np.round(ref_standard_conc*area_ref/area_std,4)
            textD.set_text("conc = "+str(conc_ref)+" mM")
            textD.set_x(corrected_ppm[ref_peaks[idx1]])
            textD.set_y(corrected_ref[ref_peaks[idx1]]+base_line)
            fig.canvas.draw()            

    fig.canvas.mpl_connect('button_release_event', on_released)            
    fig.canvas.mpl_connect('motion_notify_event', on_drag)
    fig.canvas.mpl_connect('scroll_event', on_scrolled)
    
    def detrend_local(event):
        global idx1
        new_y = line_new.get_ydata()
        ref_y = line.get_ydata()
        new_y = new_y - solve_base_line(new_y)
        ref_y = ref_y - solve_base_line(ref_y)
        line_new.set_ydata(new_y)
        line.set_ydata(ref_y)
        idx1 = (np.abs(new_y[ref_peaks] - metab_table[3,2])).argmin()
        fig.canvas.draw()
    
    ax_det = fig.add_axes([0.1, 0.81, 0.05, 0.03])
    btn_detrend = mpl_button(ax_det, 'Detrend',color='#00bfff')
    btn_detrend.on_clicked(detrend_local)
    
    def overlay_local(event):
        global corrected_ppm, corrected_ref, idx1, ref_peaks
        new_y = line_new.get_ydata()
        new_y = new_y - solve_base_line(new_y)
        line_new.set_ydata(new_y)
        
        line.set_ydata(line_new.get_ydata())
        line.set_xdata(line_new.get_xdata())
        corrected_ppm = line.get_xdata()
        corrected_ref = line.get_ydata()
        
        ref_peaks, _ = find_peaks(corrected_ref,  height=0.85*max(corrected_ref))
        idx1 = (np.abs(corrected_ref[ref_peaks] - metab_table[3,2])).argmin()        
        fig.canvas.draw()
    
    ax_ove = fig.add_axes([0.1, 0.76, 0.05, 0.03])
    btn_overlay = mpl_button(ax_ove, 'Overlay',color='#00bfff')
    btn_overlay.on_clicked(overlay_local)
    
    #anim = FuncAnimation(fig, update_plot_step, frames = 6, interval = 300)
    lines.append(line)
    globals2.PLOTTED.append([plt.gca(), fig, lines])
    plt.show()
    #fig_canvas_agg = draw_figure(ITUPS, fig)
    #draw_figure(ITUPS, fig)    
    
def pH_from_Histidine_H17(ppm):
    return (7.8195 - ppm)/0.0922
    
def pH_from_Histidine_H18(ppm):
    return (9.5375 - ppm)/0.2047
    
def pH_from_Pseudouridine(ppm):
    return (ppm - 7.637)/0.0065
    
def predict_pH_using_Histidine_and_Pseudouridine(metab_indicator):
    if CURRENT_DATA:
        query = CURRENT_DATA[0][0]
        new_ppm, new_dta = query
        base_line = 0.85*np.mean(new_dta)
        new_peaks, _ = find_peaks(new_dta,  height=base_line*3)    
        idx = new_peaks[-1]
        q1 = 0
        while new_dta[idx-q1]>base_line:
            q1 = q1 + 1
        q2 = 0
        try:
            while new_dta[idx+q2]>base_line:
                q2 = q2 + 1
        except:
            q2 = len(new_dta)-1        
            
        new_dta = detrend(new_dta)    
        scope = (new_ppm<9)*(new_ppm>6)
        new_ppm = new_ppm[scope]
        new_dta = new_dta[scope]
        
        plt.close()
        plt.figure(figsize=(9.68, 3.95))
        plt.xlabel("ppm")
        plt.ylabel("Intensity(a.u.)")
        plt.plot(new_ppm, new_dta)
        plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        plt.tight_layout()
        fig = plt.gcf()
        ax = plt.gca()
        canvas = ITUPS[0]
        fig.set_size_inches(0.8*canvas.winfo_width()/80, 0.74*canvas.winfo_height()/80, forward=True)    
        pH_line = plt.axvline(7.2,lw=1,ls='--',color='red',alpha=1)
        
        ppm_val = 7.2
        if metab_indicator == "Histidine-H17":
            pH_val = pH_from_Histidine_H17(ppm_val)
        elif metab_indicator == "Histidine-H18":
            pH_val = pH_from_Histidine_H18(ppm_val)
        else:
            pH_val = pH_from_Pseudouridine(ppm_val)
        pH_val = " pH = "+str(round(pH_val,3))        
        pH_text = plt.text(ppm_val,0.95*ax.get_ylim()[1],pH_val)
        ppm_text = plt.text(ppm_val,0.90*ax.get_ylim()[1]," ppm = "+str(round(ppm_val,3)))
        
        def on_drag(event):
            if not fig.canvas.manager.toolbar.mode:
                ppm_val = event.xdata    
                pH_line.set_xdata(ppm_val)
                if metab_indicator == "Histidine-H17":
                    pH_val = pH_from_Histidine_H17(ppm_val)
                elif metab_indicator == "Histidine-H18":
                    pH_val = pH_from_Histidine_H18(ppm_val)
                else:
                    pH_val = pH_from_Pseudouridine(ppm_val)
                pH_val = " pH = "+str(round(pH_val,3))

                pH_text.set_text(pH_val)
                pH_text.set_x(ppm_val)
                pH_text.set_y(0.95*ax.get_ylim()[1])
                
                ppm_text.set_text(" ppm = "+str(round(ppm_val,3)))
                ppm_text.set_x(ppm_val)
                ppm_text.set_y(0.90*ax.get_ylim()[1])
                fig.canvas.draw()
                    
        fig.canvas.mpl_connect('motion_notify_event', on_drag)
        plt.gca().invert_xaxis()
        plt.show()
        
ref_standard_conc = 0.58 ################################################################################################  Change
if __name__ == "__main__":
###############################################################################
    VIEWMENU = gui.Menu(FRAME, tearoff=0)
    VIEWMENU.add_command(
        label="Text file", command=lambda: load_data(ITUPS),
        background="white", foreground="Blue")
    VIEWMENU.add_command(
        label="csv file", command=lambda: load_csv(ITUPS),
        background="white", foreground="Blue")        
    VIEWMENU.add_command(
        label="excel file", command=lambda: load_xls(ITUPS),
        background="white", foreground="Blue")        
    VIEWMENU.add_command(
        label="Image file", command=load_image,
        background="white", foreground="Blue")        
    VIEWMENU.add_command(
        label="NMR pdata", command=lambda: load_pdata(plot=True),
        background="white", foreground="Blue")
    VIEWMENU.add_command(
        label="jcampdx", command=lambda: load_jcampdx(plot=True),
        background="white", foreground="Blue")

    MENUBUT1 = gui.Menubutton(
        FRAME, text="View data", activebackground="#f2f20d",
        activeforeground="red", bg="#00bfff",
        fg="white" if PLATFORM.lower() != "darwin" else "green")
    MENUBUT1.menu = gui.Menu(MENUBUT1, tearoff=0)
    MENUBUT1["menu"] = MENUBUT1.menu    
    MENUBUT1.menu.add_cascade(label="Open", menu=VIEWMENU)
    MENUBUT1.place(x=2, y=5)
    
###############################################################################
    LOADMENU1 = gui.Menu(FRAME, tearoff=0)    
    LOADMENU1.add_command(
        label="csv file", command=lambda: load_csv(ITUPS,show=False),
        background="white", foreground="Blue")        
    LOADMENU1.add_command(
        label="excel file", command=lambda: load_xls(ITUPS,show=False),
        background="white", foreground="Blue")        
    LOADMENU1.add_command(
        label="NMR pdata", command=lambda: load_pdata(plot=False),
        background="white", foreground="Blue")
    LOADMENU1.add_command(
        label="jcampdx", command=lambda: load_jcampdx(plot=False),
        background="white", foreground="Blue")
    
    MENUBUT2 = gui.Menubutton(
        FRAME, text="Process data", activebackground="#f2f20d",
        activeforeground="red", bg="#00bfff",
        fg="white" if PLATFORM.lower() != "darwin" else "green")
    MENUBUT2.menu = gui.Menu(MENUBUT2, tearoff=0)
    MENUBUT2["menu"] = MENUBUT2.menu
    MENUBUT2.menu.add_cascade(label="load data in memory", menu=LOADMENU1)
    
    with open("metabs.json") as f:
        metabs = dict(json.load(f))
    
    LOADMENU2 = gui.Menu(FRAME, tearoff=0)
    for met in metabs:
        LOADMENU2.add_command(
            label=met, command=eval("lambda : process_metab_shifts_multiple(metabs, '"+met+"')"),
            background="white", foreground="Blue")    
            
    LOADMENU3 = gui.Menu(FRAME, tearoff=0)
    for met in metabs:
        LOADMENU3.add_command(
            label=met, command=eval("lambda : process_metab_shifts_single(metabs, '"+met+"')"),
            background="white", foreground="Blue")

    LOADMENU4 = gui.Menu(FRAME, tearoff=0)
    for met in metabs:
        LOADMENU4.add_command(
            label=met, command=eval("lambda : process_metab_shifts_manual(metabs, '"+met+"')"),
            background="white", foreground="Blue")    
            
    LOADMENU5 = gui.Menu(FRAME, tearoff=0)    
    LOADMENU5.add_command(
        label="Histidine-H17", command=lambda: predict_pH_using_Histidine_and_Pseudouridine("Histidine-H17"),
        background="white", foreground="Blue")        
    LOADMENU5.add_command(
        label="Histidine-H18", command=lambda: predict_pH_using_Histidine_and_Pseudouridine("Histidine-H18"),
        background="white", foreground="Blue")        
    LOADMENU5.add_command(
        label="Pseudouridine", command=lambda: predict_pH_using_Histidine_and_Pseudouridine("Pseudouridine"),
        background="white", foreground="Blue")
    
    LOADMENU6 = gui.Menu(FRAME, tearoff=0)
    LOADMENU6.add_command(
        label="Local", command=lambda: None,
        background="white", foreground="Blue")
    LOADMENU6.add_command(
        label="Global", command=lambda: None,
        background="white", foreground="Blue")    

    LOADMENU7 = gui.Menu(FRAME, tearoff=0)    
    LOADMENU7.add_command(
        label="excel file", command=lambda: load_xls(ITUPS,show=True),
        background="white", foreground="Blue")            

    MENUBUT2.menu.add_cascade(label="Urine pH detection", menu=LOADMENU5)
    MENUBUT2.menu.add_cascade(label="Single chemical shift profiling", menu=LOADMENU3)            
    MENUBUT2.menu.add_cascade(label="Multi-chemical shift profiling", menu=LOADMENU2)
    MENUBUT2.menu.add_cascade(label="Concentration determination", menu=LOADMENU4)
    MENUBUT2.menu.add_cascade(label="Statistical analysis", menu=LOADMENU7)

    #MENUBUT2.menu.add_cascade(label="Detrending", menu=LOADMENU6)
    MENUBUT2.place(x=79, y=5)
    
    #HMDB0000177
    
###############################################################################
    MENUBUT3 = gui.Menubutton(
        FRAME, text="Modelling", activebackground="#f2f20d",
        activeforeground="red", bg="#00bfff",
        fg="white" if PLATFORM.lower() != "darwin" else "green")
    MENUBUT3.menu = gui.Menu(MENUBUT3, tearoff=0)
    MENUBUT3["menu"] = MENUBUT3.menu
    MENUBUT3.place(x=171, y=5)
###############################################################################
    FRAME1 = gui.Frame(FRAME, height=435, width=972,
                       bg='#777696', borderwidth=2)
    FRAME1.place(x=0, y=35, relheight=0.93, relwidth=1.0)
    ITUPS = prepare_frame_for_plot(FRAME1, 972, 435)
    TOP.protocol("WM_DELETE_WINDOW", on_closing)
    #TOP.bind("<Map>", lambda e: canvas_update_widgets(e, ITUPS[0]))
    TOP.bind("<Configure>", lambda e: canvas_update_widgets(e, ITUPS[0]))
    TOP.mainloop()