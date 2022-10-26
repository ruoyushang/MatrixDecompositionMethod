
energy_bin_start_const=0
energy_bin_end_const=6
energy_bin_start=0
energy_bin_end=6
#analysis_type='off'
#analysis_type='on'
#analysis_type='fiction'
#analysis_type='imposter_on'
analysis_type='science_on'

source_name = []

#source_name += ['UrsaMajorII'] # 201.9 hrs # north LZA
#source_name += ['UrsaMinor'] # 129.2 hrs # north LZA
#source_name += ['RGB_J0710_p591'] # 118.4 hrs # north
#source_name += ['1ES0229'] # 191.7 hrs
#source_name += ['PKS1424'] # 187.5 hrs
#source_name += ['PG1553'] # 164.9 hrs
#source_name += ['3C273'] # 147.2 hrs
#source_name += ['Segue1'] # 111.7 hrs
#source_name += ['NGC1275'] # 97.6 hrs # north
#source_name += ['H1426'] # 91.1 hrs # north
#source_name += ['OJ287'] # 79.9 hrs
#source_name += ['Draco'] # 76.8 hrs # north LZA 
#source_name += ['BLLac'] # 76.7 hrs # north
#source_name += ['3C264'] # 66.3 hrs
#source_name += ['1ES0502'] # 65.8 hrs # north LZA
#source_name += ['M82'] # 59.6 hrs # north LZA
#source_name += ['1ES0414'] # 48 hrs 
#source_name += ['1ES1011'] # 41.0 hrs # north
#source_name += ['1ES0647'] # 29.3 hrs

source_name += ['MGRO_J1908'] # 127 hrs
#source_name += ['MGRO_J2019'] # north
#source_name += ['Boomerang']
#source_name += ['Geminga']
#source_name += ['CasA']
#source_name += ['HESS_J1825']
#source_name += ['IC443HotSpot']
#source_name += ['CrabRHV']
#source_name += ['Crab_Elev70']
#source_name += ['Crab_Elev60']
#source_name += ['Crab_Elev50']
#source_name += ['Crab_Elev30']
#source_name += ['Crab_Offset_1p0']
#source_name += ['Crab_Offset_1p5']
#source_name += ['Tycho']
#source_name += ['GammaCygni']
#source_name += ['WComae']
#source_name += ['LHAASO_J1843']
#source_name += ['LHAASO_J1929']
#source_name += ['LHAASO_J2108']
#source_name += ['LHAASO_J2032_Baseline']
#source_name += ['LHAASO_J2032_Fall2017']

#source_name += ['PSR_J1841_m0345']
#source_name += ['PSR_J1856_p0245']
#source_name += ['PSR_J1938_p2213']
#source_name += ['PSR_J2021_p3651']
#source_name += ['PSR_J2021_p4026']
#source_name += ['PSR_J2032_p4127_Baseline']
#source_name += ['PSR_J2032_p4127_Fall2017']

#source_name += ['V_V725_Tau']

#source_name += ['GalacticPlane_All_l30']
#source_name += ['GalacticPlane_All_l40']
#source_name += ['GalacticPlane_All_l50']
#source_name += ['GalacticPlane_All_l60']
#source_name += ['GalacticPlane_All_l70']
#source_name += ['GalacticPlane_All_l80']
#source_name += ['GalacticPlane_All_l90']
#source_name += ['GalacticPlane_All_l100']
#source_name += ['GalacticPlane_All_l110']
#source_name += ['GalacticPlane_All_l120']
#source_name += ['GalacticPlane_All_l130']

output_txt = open('run_plots.sh', 'w')
if analysis_type=='imposter_on':
    for src in range(0,len(source_name)):
        for imp in range(0,6):
            observation = source_name[src]
            observation += '_Imposter%s'%(imp+1)
            print ('python3 PlotAnalysisResults.py "%s" %s %s'%(observation,energy_bin_start_const,energy_bin_end_const))
            output_txt.write('python3 PlotAnalysisResults.py "%s" %s %s \n'%(observation,energy_bin_start_const,energy_bin_end_const))
        observation = source_name[src]
        observation += '_ON'
        print ('python3 PlotAnalysisResults.py "%s" %s %s'%(observation,energy_bin_start_const,energy_bin_end_const))
        output_txt.write('python3 PlotAnalysisResults.py "%s" %s %s \n'%(observation,energy_bin_start_const,energy_bin_end_const))
elif analysis_type=='on':
    for src in range(0,len(source_name)):
        observation = source_name[src]
        observation += '_ON'
        print ('python3 PlotAnalysisResults.py "%s" %s %s'%(observation,energy_bin_start_const,energy_bin_end_const))
        output_txt.write('python3 PlotAnalysisResults.py "%s" %s %s \n'%(observation,energy_bin_start_const,energy_bin_end_const))
elif analysis_type=='off':
    for src in range(0,len(source_name)):
        observation = source_name[src]
        observation += '_OFF'
        print ('python3 PlotAnalysisResults.py "%s" %s %s'%(observation,energy_bin_start_const,energy_bin_end_const))
        output_txt.write('python3 PlotAnalysisResults.py "%s" %s %s \n'%(observation,energy_bin_start_const,energy_bin_end_const))
elif analysis_type=='science_on':
    for src in range(0,len(source_name)):
        observation = source_name[src]
        observation += '_ON'
        print ('python3 ImposterAnalysis.py "%s" %s %s 1'%(observation,energy_bin_start,energy_bin_end))
        output_txt.write('python3 ImposterAnalysis.py "%s" %s %s 1 \n'%(observation,energy_bin_start,energy_bin_end))
elif analysis_type=='fiction':
    for src in range(0,len(source_name)):
        observation = source_name[src]
        observation += '_ON'
        print ('python3 ImposterAnalysis.py "%s" %s %s 0'%(observation,energy_bin_start,energy_bin_end))
        output_txt.write('python3 ImposterAnalysis.py "%s" %s %s 0 \n'%(observation,energy_bin_start,energy_bin_end))


