#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt

def find_Nucleosides_From_MS2(high_res_ms2_file1, intens_cutoff1, ppm_tolerance1):
    ####----------Input: mgf file, ms2 intensity cutoff value (%), mass tolerance (ppm)
    ####----------Output: two dataframes, dataframe1 contains raw candidate nucleoside ms2 spectra, dataframe2 contains summarized nucleosides
    
    #-----Open mgf file
    with open(f'{high_res_ms2_file1}.mgf', 'r') as f:
        file1 = f.readlines()
    
    #-----ms2 relative intensity low % cut-off
    #--ruling out low intensity peaks that might be construed as coincidental noise
    #--%value
    intensity_percent_cutoff = intens_cutoff1

    #-----Neutral loss masses and tolerance settings
    #--stratified by charge state; since nucleosides are likely to be +1 and +2 charges, we won't consider ions of greater charge states
    #--DNA:
    neutral_loss_mass_list_charge1 = [116.0473, 132.04226]
    neutral_loss_mass_list_charge2 = [58.0237, 66.02113]
        #116.0473 = dRibose (+1; H+ or Na+; -H+ or -Na+)
        #58.0237 = dRibose(+2; 2H+ or 2Na+; -2H+ or -2Na+)
        #132.04226 = Ribose (+1; H+ or Na+)
        #66.02113 = Ribose (+2; 2H+ or 2Na+; -2H+ or -2Na+)

    #-----Mass tolerance setting
    #--this will be spectrometer and setting dependent; here we have lock-mass Orbitrap data collected at 100,000 resolution
    #--tolerance in ppm
    mass_tolerance1 = ppm_tolerance1

    #-----Declare variables
    #---Temporary variables
    temp_scan_num = 0
    temp_retention = 0
    temp_precursor = 0
    temp_precursor_intensity = 0
    temp_charge = 0
    temp_max_ms2_intensity = 0
    temp_ms2_intensity_list1 = []
    temp_ms2_mass = 0
    temp_ms2_intensity = 0
    temp_candidate_ms2_mass_list1 = []
    temp_candidate_ms2_intensity_list1 = []
    temp_tuple1 = ()
    #---Persisent variables
    candidate_list_of_tuples = []
    #---NEW
    thy_count1 = 0

    #-----Collect spectrum data iterating through mgf file one ms2 spectrum at a time
    #---ms2 masses that pass filter are stored as tuples in 'candidate_list_of_tuples list variable'

    for line1 in file1:
        if 'BEGIN IONS' in line1:
            temp_scan_num = 0
            temp_retention = 0
            temp_precursor = 0
            temp_precursor_intensity = 0
            temp_charge = 0
        elif 'TITLE=' in line1:
            line2 = line1.split(' ')
            for item1 in line2:
                if 'scan=' in item1:
                    temp_scan_num = int(item1.split('scan=')[1].split('"')[0])
        elif 'RTINSECONDS=' in line1:
            temp_retention = float(line1.split('=')[1].replace('\n', ''))
        elif 'PEPMASS=' in line1:
            temp_precursor = float(line1.split(' ')[0].split('=')[1])
            temp_precursor_intensity = float(line1.split(' ')[1].replace('\n', ''))
            #---NEW
            if ((round(abs(temp_precursor - 1.00783), 4) - 242.09027)/temp_precursor)*1000000 < mass_tolerance1 :
                thy_count1 += 1
        elif 'CHARGE=' in line1:
            temp_charge = int(line1.split('=')[1].split('+')[0])
        elif 'END ION' in line1:
            try:
                temp_max_ms2_intensity = max(temp_ms2_intensity_list1)
                temp_ms2_intensity_list1 = []
                for count, item1 in enumerate(temp_candidate_ms2_mass_list1):
                    temp_tuple1 = (temp_scan_num, temp_retention, temp_precursor, temp_precursor_intensity, temp_charge, item1, temp_candidate_ms2_intensity_list1[count], temp_max_ms2_intensity)
                    candidate_list_of_tuples.append(temp_tuple1)
                temp_tuple1 = ()
                temp_candidate_ms2_mass_list1 = []
                temp_candidate_ms2_intensity_list1 = []
            except:
                pass
        else:
            if temp_charge == 1:
                try:
                    temp_ms2_mass = float(line1.split(' ')[0])
                    temp_ms2_intensity = float(line1.split(' ')[1].replace('\n', ''))
                    temp_ms2_intensity_list1.append(temp_ms2_intensity)
                    for item1 in neutral_loss_mass_list_charge1:
                        if (abs(round(abs((temp_ms2_mass - temp_precursor)), 4) - item1)/temp_ms2_mass)*1000000 < mass_tolerance1 :
                            temp_candidate_ms2_mass_list1.append(temp_ms2_mass)
                            temp_candidate_ms2_intensity_list1.append(temp_ms2_intensity)
                except:
                    pass
            if temp_charge == 2:
                try:
                    temp_ms2_mass = float(line1.split(' ')[0])
                    temp_ms2_intensity = float(line1.split(' ')[1].replace('\n', ''))
                    temp_ms2_intensity_list1.append(temp_ms2_intensity)
                    for item1 in neutral_loss_mass_list_charge2:
                        if (abs(round(abs((temp_ms2_mass - temp_precursor)), 4) - item1)/temp_ms2_mass)*1000000 < mass_tolerance1 :
                            temp_candidate_ms2_mass_list1.append(temp_ms2_mass)
                            temp_candidate_ms2_intensity_list1.append(temp_ms2_intensity)
                except:
                    pass
    #---NEW
    print(thy_count1)
        
    #-----Create dataframe from tuples
    #---name columns accordingly
    df1 = pd.DataFrame(candidate_list_of_tuples, columns=['scan', 'retention_sec', 'precursor_mz', 'precursor_intensity', 'charge', 'ms2_mz', 'ms2_intensity', 'ms2_max_intensity'])

    #---Calculate ms2 peak relative intensities and insert in new column
    ms2_relative_intensity_list = []
    for row1 in df1.itertuples():
        temp_var1 = round(((row1.ms2_intensity/row1.ms2_max_intensity)*100), 1)
        ms2_relative_intensity_list.append(temp_var1)
    df1['ms2_relative_intensity'] = ms2_relative_intensity_list
        
    #---Calculate precursor and ms2_mz difference for each row to differentiate RNA and DNA nucleosides, create new column
    ms1_ms2_difference_list = []
    for row1 in df1.itertuples():
        temp_var1 = round((row1.precursor_mz - row1.ms2_mz), 1)
        ms1_ms2_difference_list.append(temp_var1)
    df1['ms1_ms2_diff'] = ms1_ms2_difference_list

    #---select only rows with relative intensities above user-defined cut-off value
    df2 = df1.loc[df1['ms2_relative_intensity'] > intensity_percent_cutoff].copy()

    #---create column with rounded precursor values for grouping
    rounded_precursor_list = []
    for row1 in df2.itertuples():
        temp_var1 = round(row1.precursor_mz, 3)
        rounded_precursor_list.append(temp_var1)
    df2['rounded_precursor'] = rounded_precursor_list

    #-----Create new dataframe containing mean aggregated values (simplified list of potential nucleosides)
    df3 = df2.groupby('rounded_precursor').mean()
    df3 = df3[['precursor_mz', 'precursor_intensity', 'charge', 'ms2_mz', 'ms2_intensity', 'ms2_max_intensity', 'ms2_relative_intensity', 'ms1_ms2_diff']]
    df3['charge'] = df3['charge'].astype(int)
    df3['ms2_relative_intensity'] = df3['ms2_relative_intensity'].round(1)
    
    #-----Add file column to denote file or origin
    df1['filename'] = high_res_ms2_file1
    df2['filename'] = high_res_ms2_file1
    df3['filename'] = high_res_ms2_file1

    #-----Dataframe key:
    #df1: raw dataframe
    #df2: dataframe with selected values (above relative intensity threshold)
    #df3: dataframe containing aggregate precursor values (simplified)
    
    return df1, df3

def assign_Nucleotides_to_Masses(query_mass1, query_tol1):
    ####----------Input: Nucleotide mass (query_mass1) and ppm tolerance (query_tol1)
    ####----------Output: list of strings describing outcome of mass query (nucleotide base or unknown)
    #-----Create reference lists
    base_tuple_list = [('adenine', 135.05450), ('guanine', 151.04941), ('thymine', 126.04293), 
                       ('uracil', 112.02728), ('hm_uracil', 142.03784), ('hmg_uracil', 304.09067), 
                       ('cytosine', 111.04326), ('hm_cytosine', 141.05383), ('hmg_cytosine', 303.10665)]
    charged_adduct_tuple_list = [('+H', 1.00783, 1), ('+Na', 22.98977, 1), ('+K', 38.96371, 1), ('+Ca', 39.96259, 2), ('+Mg,', 23.98504, 2), 
    #                             ('+2H', 2.01565, 2), ('+2Na', 45.97954, 2), ('+2K', 77.92741, 2), ('+H+Na', 23.99759, 2), ('+H+K', 39.97153, 2), ('+Na+K', 61.95348, 2)
                                ]
    neutral_adduct_tuple_list = [('', 0.0), ('+NH3', 17.02655), ('+CH3COOH', 58.00548), ('+H2O', 18.01056), ('+TFA', 113.99286), ('+MeCN', 41.02655), ('+MeOH', 32.02621), ('+PrOH', 60.05751)
                                ]

    #-----Declare variables
    list_of_hits1 = []
    #-----Iterate through mass combinations
    for item1 in charged_adduct_tuple_list:
        for item2 in neutral_adduct_tuple_list:
            for item3 in base_tuple_list:
                if (abs(query_mass1 - ((item3[1] + item2[1] + item1[1])/item1[2]))/query_mass1)*1000000 < query_tol1:
                    temp_identity1 = str('[' + item3[0] + item2[0] + item1[0] + ']')
                    list_of_hits1.append(temp_identity1)
    return list_of_hits1


    