#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import copy

def read_ms1(ms1_file1):
    ####----------Input: ms1 file
    ####----------Output: dictionary with scan number as key, 
    ####----------Output: return dictionary structure: dictionary[scan]{float(retention), list(tuple(m/z, intensity))}
    
    #-----open ms1 file
    with open(f'{ms1_file1}.ms1', 'r') as f:
        file1 = f.readlines()
    
    #-----declare variables
    temp_scan1 = 0
    temp_retention1 = 0
    temp_mz_itensity_tuple1 = ()
    temp_mz_intensity_list1 = []
    master_dict1 = {}
    
    #-----iterate through ms1 file
    for line1 in file1:
        if line1[0] == 'S':
            line2 = line1.split('\t')
            if int(line2[1]) > 1:
                master_dict1[temp_scan1] = [temp_retention1, temp_mz_intensity_list1]
                temp_scan1 = 0
                temp_retention1 = 0
                temp_intensity_tuple1 = ()
                temp_mz_intensity_list1 = []
            temp_scan1 = int(line2[1])
        elif line1[0] == 'I':
            line2 = line1.split('\t')
            if line2[1] == 'RetTime':
                temp_retention1 = round(float(line2[2]), 4)
        elif (line1[0] != 'I') and (line1[0] != 'S') and (line1[0] != 'H'):
            line2 = line1.split(' ')
            temp_mz_intensity_tuple1 = (float(line2[0]), float(line2[1]))
            temp_mz_intensity_list1.append(temp_mz_intensity_tuple1)
    
    master_dict1[temp_scan1] = [temp_retention1, temp_mz_intensity_list1]
    temp_scan1 = 0
    temp_retention1 = 0
    temp_intensity_tuple1 = ()
    temp_mz_intensity_list1 = []
    
    #-----return object_
    return master_dict1

def extract_ms1_tuples(input_ms1_dict1, input_mz1, input_ppm_tol1, baseline_interval1):
    ####----------Input: ms1 dictionary (input_ms1_dict1), m/z (input_mz1), ppm tolerance (input_ppm_tol1), time interval to add null values to baseline (baseline_interval1)
    ####----------Output: list of retention times, list of intensities
    
    #-----declare variables
    temp_tuple1 = ()
    temp_tuple_list1 = []
    temp_retention_list1 = []
    temp_intensity_list1 = []
    
    #-----add baseline values in periodic intervals according to baseline_interval1 variable
    for item1 in np.arange(0,10, baseline_interval1):
        temp_tuple1 = (float(item1), float(0))
        temp_tuple_list1.append(temp_tuple1)
    
    #-----iterate through dictionary and collect retention times and intensity
    for item1 in input_ms1_dict1.keys():
        for item2 in input_ms1_dict1[item1][1]:
            if (round((abs(float(item2[0]) - float(input_mz1))/input_mz1), 4) * 1000000) < float(input_ppm_tol1):
                temp_tuple1 = (input_ms1_dict1[item1][0], item2[1])
                temp_tuple_list1.append(temp_tuple1)
    
    #-----re-order tuple list
    temp_tuple_list1 = sorted(temp_tuple_list1, key=lambda x:x[0])
    
    #-----iterate through tuple list and remove baseline insertions that split peaks
    temp_removable_indices_list1 = []
    
    for item1 in range(1, len(temp_tuple_list1)-1):
        if temp_tuple_list1[item1][1] == 0:
            if (temp_tuple_list1[item1-1][1] > 0) and (temp_tuple_list1[item1+1][1] > 0):
                temp_removable_indices_list1.append(item1)
    
    temp_tuple_list2 = []
    for position1,item1 in enumerate(temp_tuple_list1):
        if position1 not in temp_removable_indices_list1:
            temp_tuple_list2.append(item1)
    
    #-----unzip tuple list
    temp_list1 = list(zip(*temp_tuple_list2))
    temp_retention_list1 = temp_list1[0]
    temp_intensity_list1 = temp_list1[1]
    
    #-----return
    return temp_retention_list1, temp_intensity_list1
    
def process_multiple_inputs(file_list1, input_mz, ppm_tol1, baseline1):
    ####----------Input: iteratively runs extract_ms1_tuples function; list of ms1 files(file_list1), mz (input_mz), ppm_tolerance (ppm_tol1), baseline intervals (baseline1)
    ####----------Output: dictionary with string keys (filenames) and list of lists as values (retention time and intensity)
    
    #-----declare variables
    locals_dict1 = locals()
    master_dict1 = {}
    
    #-----iterate through file list and fill variables
    for item1 in file_list1:
        item2 = item1.split('/')[-1].split('.')[0]
        temp_str1 = 'retention_list_' + str(item2)
        temp_str2 = 'intensity_list_' + str(item2)
        locals_dict1[temp_str1], locals_dict1[temp_str2] = extract_ms1_tuples(read_ms1(item1), input_mz, ppm_tol1, baseline1)
        master_dict1[item2] = [locals_dict1[temp_str1], locals_dict1[temp_str2]]
    
    #-----return
    return master_dict1

def plot_multiple_inputs(master_dict1, input_title1, save_fig1, fig_title1):
    ####----------Input: dictionary with string keys (filenames) and list of lists as values (retention time list and intensity list): master_dict1, plot title (input_title1)
    ####----------Output: matplotlib stacked plot with as many plots as input filenames
    
    #-----declare variables
    temp_max_intensity = 0
    max_intensity = 0
    
    #-----find max intensity from all file sets
    for position1, item1 in enumerate(master_dict1.keys()):
        temp_max_intensity = max(master_dict1[item1][1])
        if temp_max_intensity > max_intensity:
            max_intensity = copy.copy(temp_max_intensity)
    
    #-----plot
    fig, ax = plt.subplots(len(master_dict1), 1, sharex=True, sharey=True, figsize=(20,20), constrained_layout=True)
    fig.suptitle(input_title1, size=16)
    fig.text(-0.01, 0.5, 'Intensity', size=14, rotation=90, ha='center', va='center')

    for position1,item1 in enumerate(master_dict1.keys()):
        ax[position1].plot(master_dict1[item1][0], master_dict1[item1][1])
        temp_max_intensity = max(master_dict1[item1][1])
        
        #find retention time of max-intensity peak
        max_intensity_index1 = master_dict1[item1][1].index(temp_max_intensity)
        max_intensity_retention1 = master_dict1[item1][0][max_intensity_index1]
        
        ax[position1].set_ylim((max_intensity*-0.1), max_intensity*1.1)
        ax[position1].set_title(str(item1) + '\n' + '[max peak] intensity:' + str('{:.2E}'.format(temp_max_intensity)) + ', retention:' + str(round(max_intensity_retention1, 2)), loc='right', y=0.6)
    
    plt.xlim([0,10])
    plt.xticks(np.arange(0,10,0.5))
    plt.xlabel('Retention Time (min)', size=14)
    if save_fig1 == True:
        plt.savefig(f'{fig_title1}.jpg', bbox_inches='tight')
        plt.close()
    else:
        plt.show()

def find_Tallest_Peaks(retention_range1, input_ms1, mz_num1, input_ppmtol1):
    ####----------Input: retention time window (retention_range1), ms1 file (input_ms1), topN m/z to output (mz_num1), ppm tolerance (input_ppmtol1)
    ####----------Output: list of m/z with mz_num1 elements
    
    #-----declare variables
    ms1_dict1 = {}
    retention_min1 = 0
    retention_max1 = 0
    
    #-----read ms1 file
    ms1_dict1 = read_ms1(input_ms1)
    mz_intensity_list1 = []
    mz_intensity_list2 = []
    
    #-----define retention ranges
    retention_min1 = round(float(retention_range1.split('-')[0]), 4)
    retention_max1 = round(float(retention_range1.split('-')[1]), 4)
    
    #-----iterate through ms1 dictionary and select topN m/z by intensity (mz_num1) and put into list (mz_intensity_list1)
    for item1 in ms1_dict1.keys():
        if (ms1_dict1[item1][0] >= retention_min1) and (ms1_dict1[item1][0] <= retention_max1):
            item2 = sorted(ms1_dict1[item1][1], key=lambda x: x[1], reverse=True)[:mz_num1]
            mz_intensity_list1 += item2
    
    #-----sort m/z and intensity tuples list by m/z, combine entries based on m/z and ppm tolerance, keep highest intensity value
    mz_intensity_list1 = sorted(mz_intensity_list1, key=lambda x:x[0])
    temp_mz_intensity_tuple1 = (0,0)
    
    for item1 in range(1,len(mz_intensity_list1)):
        if (abs(mz_intensity_list1[item1][0] - mz_intensity_list1[item1-1][0])/(mz_intensity_list1[item1-1][0]))*1000000 < input_ppmtol1:
            if mz_intensity_list1[item1][1] > mz_intensity_list1[item1-1][1]:
                temp_mz_intensity_tuple1 = mz_intensity_list1[item1]
        else:
            mz_intensity_list2.append(temp_mz_intensity_tuple1)
            temp_mz_intensity_tuple1 = mz_intensity_list1[item1]
    mz_intensity_list2.append(temp_mz_intensity_tuple1)
    
    #-----select topN (based on mz_num1) highest intensity m/z
    mz_intensity_list2 = sorted(mz_intensity_list2, key=lambda x:x[1], reverse=True)[:mz_num1]
    
    #-----resort list by m/z
    mz_intensity_list2 = sorted(mz_intensity_list2, key=lambda x:x[0])
    
    #-----return
    return mz_intensity_list2

def create_Input_File_List(input_file_csv1):
    ####----------Input: csv file containing table of input files (input_file_csv1) where 'input_file_name' column contains names of input files 
    ####----------Output: list of strings containing all input file names (temp_input_file_list1)
    
    temp_df1 = pd.read_csv(input_file_csv1, sep=',')
    temp_input_file_list1 = list(temp_df1['input_file_name'])
    temp_input_label_tuple_list1 = zip(temp_input_file_list1, list(temp_df1['label']))
    return temp_input_file_list1, temp_input_label_tuple_list1




