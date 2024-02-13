#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os, sys
import seaborn as sns
import matplotlib.pyplot as plt 

def normalize_times_add_parameters(folder, OSym="even", run=2, protein="WT", IPTG="0.3 mM", replica =1, switching_time=0):
    """Reads the result from ImAnalysis, adds information on the experiment and time of the run.
    
        Parameters:        
            folder (str):the folder with the data.
            OSym (str): even or odd dependent on positions with the strain having even or odd numbers.
            run (int): number of the run to be analysed.
            protein (str): defined the protein mutant ("WT", "V52A" or "Q55N") contained in the cells.
            IPTG (str): defines the IPTG concentration used in the experiment.
            switching_time (int): time (s) of the media switch

        Returns:
            df (pandas DF): dataframe that contains the data from the image analysis with the time in s and specifics 
                about the experiment.

    """
    
    df = pd.read_excel(folder+'dots_run{}.xlsx'.format(run), engine='openpyxl')

    time1 = np.round(df.iloc[0]["time"]/1000)

    norm_time = np.round(df["time"]/1000)-time1-switching_time 
    df["Time (s)"] = norm_time
    
    max_pos = np.max(df["Position"].values)
    
    if OSym=="even":
        strains = np.where( df["Position"].isin(np.arange(1,max_pos+1,2)) , "term", "OSymL")
    if OSym=="odd":
        strains = np.where( df["Position"].isin(np.arange(1,max_pos+1,2)) , "OSymL", "term")

    df["Run"] = [run for x in range(len(df))]
    df["Strain"] = strains
    df["Run"] = [run for x in range(len(df))]
    df["LacI"] = [protein for x in range(len(df))]
    df["IPTG"] = [IPTG for x in range(len(df))]
    df["Replica"] = [replica for x in range(len(df))]
    
    return(df)

def average_dot_counts(df): 
    """Counts and averages the number of dots per cell in each position.
    
            Parameters:
                df (pandas DataFrame): from normalizeTimesAddStrain(OSym="even" or "odd", run = X)

            Returns:
                avgNumDotsCells (pandas DataFrame):contains the averaged number of dots, runs id, strain, time (s), 
                                        IPTG and LacI mutant.
    """
    run=df
    avgNumDotsCells = np.array(run.groupby("position").count()["numCells"])/ \
    np.array(run.groupby("position").mean(numeric_only=True)["numCells"])

    time = np.array(run.groupby("position").first()["Time (s)"])
    strain = np.array(run.groupby("position").first()["Strain"])

    runs = np.array(run.groupby("position").first()["Run"])
    IPTG = np.array(run.groupby("position").first()["IPTG"])
    LacI = np.array(run.groupby("position").first()["LacI"])
    replica = np.array(run.groupby("position").first()["Replica"])
    
    avgNumDotsCells = pd.DataFrame({"Run":runs, "Strain":strain, "avg. N(Dots)/Cell": \
             avgNumDotsCells, "Time (s)": time, "IPTG": IPTG, "LacI": LacI, "Replica":replica})
    return(avgNumDotsCells)

def all_runs(folders, ids_runs, locs_OSym, protein="V52A", switching_time=0):
    """Combines the average number of cells record for all runs into one DF.
    
            Parameters:
                ids_runs (list):the ids of the runs to sumarize.
                locs_OSym (list): "even" or "odd", is Osym imaged from the even or odd positions?
                protein (str): the protein (WT, V52A or Q55N)
                switching_time (int):time (s) of the media switch

            Returns:
                df_runs (pandas DataFrame): summarized data frame.

            Calls:
                normalize_times_add_parameters
                average_dot_counts
    """
    
    j=0
    df_list = []
    
    while j<len(folders):
        print(j)
        for i,m in zip(ids_runs[j], locs_OSym[j]):
                df = normalize_times_add_parameters(folders[j], OSym=m, run=i, replica=j+1,\
                                          protein=protein, switching_time=switching_time)
                runNum = average_dot_counts(df)
                df_list.append(runNum)
        j+=1
                
    df_runs = pd.concat(x for x in df_list)
    
    return(df_runs.reset_index())

def interpolate(df, j, i):
    """Get the linear interpolation between non-specific dots to substract from the mixed binding.
    
        Parameters:
            df (pandas DataFrame):
                contains data from microscopy
            j (int):
                replica (different days of the experiments, different chips, etc.)
            i (int):
                run (imaged on the same chip)
        
        Returns:
            lininp (numpy array):
                contains the interpolated values of non-specific bindings to match the time of the mixed one.
    """
    
    SymL = df[(df["Strain"] == "OSymL") & (df["Replica"] ==  j) ]
    term = df[(df["Strain"] == "term")& (df["Replica"] ==  j) ]
    
    lininp = np.interp(SymL[SymL["Run"]==i]["Time (s)"].values,\
                      term[term["Run"]==i]["Time (s)"].values,\
                      term[term["Run"]==i]["avg. N(Dots)/Cell"].values)
    return(lininp)

def calculate_delta_dots(df, LacI= "WT", IPTG="0.3 mM", ONPG=False ,runs = [1,2,3,4], replica = 1):
    """ Calculates the difference between number of dots in cells with and without operator --> specific binding.
    
    This works for a DataFrame that contains data from one or several experiments, where the relevant one is 
    extracted using the keywords. The time is adjusted by substracting the time of the media switch.

    The number of dots for the empty strain is substracted from the one with the operator. 
    As the data for the two strain cannot be recordet at the same time, the curves for the control strain
    are fitted using a cubic splines first. The curves serve to calculate the expected number of dots recorded 
    at the same time as the image for the operator containing strain was recorded.
    
    Finally the curves of the runs (Usually 4 runs per experiment) are averaged to obtain a mean curve and 
    standard deviation.
    
        Parameters:
            df (pandas DataFrame): contains the data from the microscopy experiment, where the strain and conditions are defined.
            LacI (str):defines the protein
            IPTG (str):defines the IPTG concentration in the media
            ONPG (bol, str):False when no ONPG has been used, otherwise the ONPG concentration
            runs (int):number of runs analysed for the df
            switching_time (int):time (s) of the media switch
    
        Return:
            dDots (pandas DataFrame): 
    """

    # get the experiments 
    if ONPG==True:
        rr_LacI = df[(df["IPTG"]==IPTG)&(df["LacI"]==LacI)&(df["ONPG"]==ONPG) & (df["Replica"]==replica) & (df["Run"].isin(runs))]
    else:
        rr_LacI = df[(df["IPTG"]==IPTG)&(df["LacI"]==LacI) & (df["Replica"]==replica) & (df["Run"].isin(runs))]
    
    # correct the time to 0 at the switch 
    rr_LacI_ = rr_LacI.copy()
    #rr_LacI_["Time (s)"] = rr_LacI["Time (s)"]-switching_time # Here you want a first point at 0
    
    #  exchange negative times with 0
    rr_LacI_.loc[rr_LacI_["Time (s)"] <0 , "Time (s)"] = 0
        
    # split the frame according to strain
    rr_LacI_OSymL = rr_LacI_[rr_LacI_["Strain"] == "OSymL"]
    rr_LacI_term = rr_LacI_[rr_LacI_["Strain"] == "term"]

    # calculate the difference between the average number of dots per cell in strains with and without operator
    relNdots = []

    for i in runs:
        # interpolate the number of non-specific dots for the specific binding time points    
        new_term_points = interpolate(rr_LacI_, replica, i)
        
        relNdots.append(np.array(rr_LacI_OSymL[rr_LacI_OSymL["Run"]==i]["avg. N(Dots)/Cell"])-new_term_points)
    rr_LacI_OSymL_rel = rr_LacI_OSymL.copy()
    rr_LacI_OSymL_rel["$\Delta$ avg. N(dots)/cell"] =  [x for i in relNdots for x in i]

    return(rr_LacI_OSymL_rel)
