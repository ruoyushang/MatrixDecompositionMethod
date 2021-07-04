
import os
SMI_DIR = os.environ['SMI_DIR']
SMI_OUTPUT = os.environ['SMI_OUTPUT']
#CONDA_DIR = os.environ['CONDA_DIR']

source = []
mjd_cut = []
for_syst = []
for_raster = []
gamma_model = []

#epoch_bins = [53613,55074,56535,57996,59457]
epoch_bins = [0,0]

folder = '%s'%(SMI_OUTPUT)
rank = 3
MSCW_cut = 0.3
MSCL_cut = 0.3


elev = []
#elev_bins = [85,55,25]
elev_bins = [85,45]
#elev_bins = [85,75,65,55,45]
#elev_bins = [85,75]
theta2 = []

#source += ['SgrAV6_ON']
#elev += [[55,25]]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
source += ['3C273V6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]
source += ['3C273V5_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]

source += ['1ES0502V6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]
source += ['1ES0502V5_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]

source += ['DracoV6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]
source += ['DracoV5_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]

source += ['1ES0647V6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]

source += ['1ES1011V6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]

source += ['3C264V6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]

source += ['PKS1424V6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]

source += ['H1426V6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]
#
#source += ['NGC1275V6_OFF']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [0]
#
source += ['OJ287V6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]

source += ['BLLacV6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]
#source += ['BLLacV5_OFF']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['RGBJ0710V5_OFF']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [0]
#
source += ['1ES0229V6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]
#source += ['1ES0229V5_OFF']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['RBS0413V6_OFF']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [0]
#source += ['RBS0413V5_OFF']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [0]
#
source += ['PG1553V5_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]
#
source += ['Segue1V6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]
source += ['Segue1V5_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]
#
#source += ['Mrk421V5_OFF']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['Mrk421V5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
source += ['M82V6_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]
source += ['M82V5_OFF']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [True]
for_raster += [False]
gamma_model += [0]
#source += ['M82V4_OFF']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['CrabV6_OFF']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [0]
#source += ['CrabV5_OFF']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [0]
#source += ['CrabV4_OFF']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['CrabV6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['CrabV5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['CrabV4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#

source += ['3C273V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
source += ['3C273V5_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]

source += ['1ES0502V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
source += ['1ES0502V5_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]

source += ['DracoV6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
source += ['DracoV5_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]

source += ['1ES0647V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]

source += ['1ES1011V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]

source += ['3C264V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]

source += ['PKS1424V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]

source += ['H1426V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
#
#source += ['NGC1275V6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
source += ['OJ287V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]

source += ['CasAV6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
#
#source += ['2HWC_J1953V6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['2HWC_J1930V6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
source += ['SNR_G150p3Plus04p5_V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]

source += ['SNR_G150p3Plus04p5_Jamie_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]

source += ['BLLacV6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
#source += ['BLLacV5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['RGBJ0710V5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
source += ['1ES0229V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
#source += ['1ES0229V5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['RBS0413V6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['RBS0413V5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
source += ['PG1553V5_Raster']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [True]
gamma_model += [0]
source += ['PG1553V5_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
#
source += ['Segue1V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
source += ['Segue1V5_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
#
source += ['M82V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
source += ['M82V5_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
#source += ['M82V4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['CygnusV6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['CygnusV5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['GemingaV6_ON']
#elev += [elev_bins]
#mjd_cut += [epoch_bins]
#theta2 += [[0,4]]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['GemingaV5_ON']
#elev += [elev_bins]
#mjd_cut += [epoch_bins]
#theta2 += [[0,4]]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['M87V6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['M87V5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['BoomerangV6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['BoomerangV5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['BoomerangV4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
source += ['MGRO_J1908_V6_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
source += ['MGRO_J1908_V5_ON']
elev += [elev_bins]
theta2 += [[0,4]]
mjd_cut += [epoch_bins]
for_syst += [False]
for_raster += [False]
gamma_model += [0]
#source += ['MGRO_J1908_V4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['SS433_V6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['SS433_V5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['SS433_V4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['SS433Half1_V6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['SS433Half1_V5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['SS433Half1_V4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['SS433Half2_V6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['SS433Half2_V5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['SS433Half2_V4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
##
#source += ['MAGIC_J1857_V6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['MAGIC_J1857_V5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['MAGIC_J1857_V4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['IC443HotSpotV6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['IC443HotSpotV5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['IC443HotSpotV4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['WComaeV6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['WComaeV5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['WComaeV4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['MGRO_J2031_V6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['MGRO_J2031_V5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['MGRO_J2031_V4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['ComaV6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['ComaV4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['ComaV6_Model05']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [5]
#source += ['ComaV4_Model05']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [5]
#
#source += ['ComaV6_Model06']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [6]
#source += ['ComaV4_Model06']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [6]
#
#source += ['GammaCygniV6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['GammaCygniV5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['GammaCygniV4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['TychoV6_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['TychoV5_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#source += ['TychoV4_ON']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['HESS_J1825_V6_ON']
#elev += [[55,25]]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [False]
#for_raster += [False]
#gamma_model += [0]
#
#source += ['1ES1011V6_Model01']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [1]
#source += ['1ES1011V6_Model02']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [2]
#source += ['1ES1011V6_Model03']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [3]
#
#source += ['CrabV6_Model01']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [1]
#source += ['CrabV5_Model01']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [1]
#source += ['CrabV4_Model01']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [1]
#source += ['CrabV6_Model02']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [2]
#source += ['CrabV5_Model02']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [2]
#source += ['CrabV4_Model02']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [2]
#source += ['CrabV6_Model03']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [3]
#source += ['CrabV5_Model03']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [3]
#source += ['CrabV4_Model03']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [3]
#source += ['CrabV6_Model04']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [4]
#source += ['CrabV5_Model04']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [4]
#source += ['CrabV4_Model04']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [4]
#
#source += ['1ES0229V6_Model01']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [1]
#source += ['1ES0229V5_Model01']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [1]
#source += ['1ES0229V6_Model02']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [2]
#source += ['1ES0229V5_Model02']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [2]
#source += ['1ES0229V6_Model03']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [3]
#source += ['1ES0229V5_Model03']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [3]
#
#source += ['OJ287V6_Model01']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [1]
#source += ['OJ287V6_Model02']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [2]
#source += ['OJ287V6_Model03']
#elev += [elev_bins]
#theta2 += [[0,4]]
#mjd_cut += [epoch_bins]
#for_syst += [True]
#for_raster += [False]
#gamma_model += [3]


gamma = []
gamma += [0]
#gamma += [0,10,20,50,100]

job_counts = 0
for s in range(0,len(source)):
    job_counts += 1
    file = open("run/run_Netflix1_%s_MJD%sto%s.sh"%(source[s],mjd_cut[s][0],mjd_cut[s][1]),"w") 
    #file.write('#### submit_job.sh START ####\n')
    #file.write('#!/bin/bash\n')
    #file.write('#$ -cwd\n')
    #file.write('# error = Merged with joblog\n')
    #file.write('#$ -o joblog.$JOB_ID\n')
    #file.write('#$ -j y\n')
    #file.write('## Edit the line below as needed:\n')
    #file.write('#$ -l h_rt=2:00:00,h_data=1G\n')
    #file.write('## Modify the parallel environment\n')
    #file.write('## and the number of cores as needed:\n')
    #file.write('#$ -pe shared 1\n')
    #file.write('# Email address to notify\n')
    #file.write('#$ -M $USER@mail\n')
    #file.write('# Notify when\n')
    ##file.write('#$ -m bea\n')
    #file.write('#$ -m e\n')
    #file.write('. /u/local/Modules/default/init/modules.sh\n')
    #file.write('module load cern_root/6.12.06\n')
    ##file.write('module load anaconda3\n')
    ##file.write('source %s/etc/profile.d/conda.sh\n'%(CONDA_DIR))
    ##file.write('conda activate root_env\n')
    file.write('cd %s\n'%(SMI_DIR))
    file.write('source %s/setup_env.sh\n'%(SMI_DIR))
    file.write('rm -r %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('mkdir %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('cp GetRunList.h %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('cp NetflixParameters.h %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('cp PrepareDarkData.C %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('cp MakePrediction.C %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('cd %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('rm *_C*\n')
    elev_this = elev[s]
    theta2_this = theta2[s]
    mjd_this = mjd_cut[s]
    for e in range(0,len(elev_this)-1):
        for u in range(0,len(mjd_this)-1):
            if for_syst[s]:
                file.write("root -b -l -q 'PrepareDarkData.C+(\"%s\",%s,%s,%s,%s,%s,%s,false,false,%s)'\n"%(source[s],elev_this[e+1],elev_this[e],mjd_cut[s][u],mjd_cut[s][u+1],theta2_this[0],theta2_this[1],gamma_model[s])) 
            else:
                if for_raster[s]:
                    file.write("root -b -l -q 'PrepareDarkData.C+(\"%s\",%s,%s,%s,%s,%s,%s,true,true,%s)'\n"%(source[s],elev_this[e+1],elev_this[e],mjd_cut[s][u],mjd_cut[s][u+1],theta2_this[0],theta2_this[1],gamma_model[s])) 
                else:
                    file.write("root -b -l -q 'PrepareDarkData.C+(\"%s\",%s,%s,%s,%s,%s,%s,true,false,%s)'\n"%(source[s],elev_this[e+1],elev_this[e],mjd_cut[s][u],mjd_cut[s][u+1],theta2_this[0],theta2_this[1],gamma_model[s])) 
    #file.write('# echo job info on joblog:\n')
    #file.write('echo "Job $JOB_ID ended on:   " `hostname -s`\n')
    #file.write('echo "Job $JOB_ID ended on:   " `date `\n')
    #file.write('echo " "\n')
    #file.write('#### submit_job.sh STOP ####\n')
    file.close() 

job_counts = 0
qfile = open("run/qsub_Netflix1.sh","w") 
for s in range(0,len(source)):
    job_counts += 1
    qfile.write('qsub -V -N job_%s run_Netflix1_%s_MJD%sto%s.sh\n'%(source[s],source[s],mjd_cut[s][0],mjd_cut[s][1]))
    qfile.write('sleep 30s\n')
    #qfile.write('qsub -l nodes=1:gamma5:ppn=1 -V -N job_%s run_Netflix1_%s_MJD%sto%s_%s.sh\n'%(source[s],source[s],mjd_cut[s][0],mjd_cut[s][1],folder))
qfile.close() 


# 1ES 1215 flare
#source += ['WComaeV6']
#mjd_cut += [[56372,56464]]
#for_syst += [False]
#for_raster += [False]

# W Comae and 1ES 1218 flare
#source += ['WComaeV4']
#mjd_cut += [[54400,54700]]
#for_syst += [False]
#for_raster += [False]

# MGRO J2031
#source += ['MGRO_J2031_V5']
#mjd_cut += [[55000,57997]]
#for_syst += [False]
#for_raster += [False]
#source += ['MGRO_J2031_V6']
#mjd_cut += [[55000,57997]]
#for_syst += [False]
#for_raster += [False]
#source += ['MGRO_J2031_V6']
#mjd_cut += [[57997,59000]]
#for_syst += [False]
#for_raster += [False]

# 1ES 0229
#source += ['1ES0229V5']
#mjd_cut += [[55000,56500]]
#for_syst += [False]
#for_raster += [False]
#source += ['1ES0229V6']
#mjd_cut += [[55000,56500]]
#for_syst += [False]
#for_raster += [False]
#source += ['1ES0229V5']
#mjd_cut += [[56501,57500]]
#for_syst += [False]
#for_raster += [False]
#source += ['1ES0229V6']
#mjd_cut += [[56501,57500]]
#for_syst += [False]
#for_raster += [False]
#source += ['1ES0229V5']
#mjd_cut += [[57501,59500]]
#for_syst += [False]
#for_raster += [False]
#source += ['1ES0229V6']
#mjd_cut += [[57501,59500]]
#for_syst += [False]
#for_raster += [False]

job_counts = 0
for s in range(0,len(source)):
    job_counts += 1
    file = open("run/run_Netflix2_%s_MJD%sto%s.sh"%(source[s],mjd_cut[s][0],mjd_cut[s][1]),"w") 
    #file.write('#### submit_job.sh START ####\n')
    #file.write('#!/bin/bash\n')
    #file.write('#$ -cwd\n')
    #file.write('# error = Merged with joblog\n')
    #file.write('#$ -o joblog.$JOB_ID\n')
    #file.write('#$ -j y\n')
    #file.write('## Edit the line below as needed:\n')
    #file.write('#$ -l h_rt=0:30:00,h_data=1G\n')
    #file.write('## Modify the parallel environment\n')
    #file.write('## and the number of cores as needed:\n')
    #file.write('#$ -pe shared 1\n')
    #file.write('# Email address to notify\n')
    #file.write('#$ -M $USER@mail\n')
    #file.write('# Notify when\n')
    ##file.write('#$ -m bea\n')
    #file.write('#$ -m e\n')
    #file.write('. /u/local/Modules/default/init/modules.sh\n')
    #file.write('module load cern_root/6.12.06\n')
    ##file.write('module load anaconda3\n')
    ##file.write('source %s/etc/profile.d/conda.sh\n'%(CONDA_DIR))
    ##file.write('conda activate root_env\n')
    file.write('cd %s\n'%(SMI_DIR))
    file.write('source %s/setup_env.sh\n'%(SMI_DIR))
    file.write('rm -r %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('mkdir %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('cp GetRunList.h %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('cp NetflixParameters.h %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('cp PrepareDarkData.C %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('cp MakePrediction.C %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('cd %s/%s_MJD%sto%s\n'%(folder,source[s],mjd_cut[s][0],mjd_cut[s][1]))
    file.write('rm *_C*\n')
    elev_this = elev[s]
    theta2_this = theta2[s]
    mjd_this = mjd_cut[s]
    for e in range(0,len(elev_this)-1):
        for u in range(0,len(mjd_this)-1):
            if for_syst[s]:
                file.write("root -b -l -q 'MakePrediction.C+(\"%s\",%s,%s,%s,%s,%s,%s,false,%s)'\n"%(source[s],elev_this[e+1],elev_this[e],mjd_cut[s][u],mjd_cut[s][u+1],theta2_this[0],theta2_this[1],gamma_model[s]))
            else:
                file.write("root -b -l -q 'MakePrediction.C+(\"%s\",%s,%s,%s,%s,%s,%s,true,%s)'\n"%(source[s],elev_this[e+1],elev_this[e],mjd_cut[s][u],mjd_cut[s][u+1],theta2_this[0],theta2_this[1],gamma_model[s]))
    #file.write('# echo job info on joblog:')
    #file.write('echo "Job $JOB_ID ended on:   " `hostname -s`')
    #file.write('echo "Job $JOB_ID ended on:   " `date `')
    #file.write('echo " "')
    #file.write('#### submit_job.sh STOP ####')
    file.close() 

job_counts = 0
qfile = open("run/qsub_Netflix2.sh","w") 
for s in range(0,len(source)):
    job_counts += 1
    #if not for_syst[s]: continue
    qfile.write('qsub -V -N job_%s run_Netflix2_%s_MJD%sto%s.sh\n'%(source[s],source[s],mjd_cut[s][0],mjd_cut[s][1]))
    qfile.write('sleep 5s\n')
    #qfile.write('qsub -l nodes=1:gamma5:ppn=1 -V -N job_%s run_Netflix2_%s_MJD%sto%s_%s.sh\n'%(source[s],source[s],mjd_cut[s][0],mjd_cut[s][1],folder))
qfile.close() 

job_counts = 0
qfile = open("run/local_Netflix.sh","w") 
for s in range(0,len(source)):
    if job_counts==7:
        job_counts = 0
        qfile.write('wait\n')
    qfile.write('sh run_Netflix1_%s_MJD%sto%s.sh &\n'%(source[s],mjd_cut[s][0],mjd_cut[s][1]))
    job_counts += 1
qfile.close() 

file = open("run/run_ObservingEffect.sh","w") 
#file.write('#### submit_job.sh START ####\n')
#file.write('#!/bin/bash\n')
#file.write('#$ -cwd\n')
#file.write('# error = Merged with joblog\n')
#file.write('#$ -o joblog.$JOB_ID\n')
#file.write('#$ -j y\n')
#file.write('## Edit the line below as needed:\n')
#file.write('#$ -l h_rt=2:00:00,h_data=1G\n')
#file.write('## Modify the parallel environment\n')
#file.write('## and the number of cores as needed:\n')
#file.write('#$ -pe shared 1\n')
#file.write('# Email address to notify\n')
#file.write('#$ -M $USER@mail\n')
#file.write('# Notify when\n')
##file.write('#$ -m bea\n')
#file.write('#$ -m e\n')
#file.write('. /u/local/Modules/default/init/modules.sh\n')
#file.write('module load cern_root/6.12.06\n')
##file.write('module load anaconda3\n')
##file.write('source %s/etc/profile.d/conda.sh\n'%(CONDA_DIR))
##file.write('conda activate root_env\n')
file.write('cd %s\n'%(SMI_DIR))
file.write('source %s/setup_env.sh\n'%(SMI_DIR))
file.write('rm -r %s/ObservingEffect\n'%(folder))
file.write('mkdir %s/ObservingEffect\n'%(folder))
file.write('cp GetRunList.h %s/ObservingEffect\n'%(folder))
file.write('cp NetflixParameters.h %s/ObservingEffect\n'%(folder))
file.write('cp PrepareDarkData.C %s/ObservingEffect\n'%(folder))
file.write('cp MakePrediction.C %s/ObservingEffect\n'%(folder))
file.write('cp ObservingEffect.C %s/ObservingEffect\n'%(folder))
file.write('cd %s/ObservingEffect\n'%(folder))
file.write('rm *_C*\n')
file.write("root -b -l -q 'ObservingEffect.C+'\n")
#file.write('# echo job info on joblog:\n')
#file.write('echo "Job $JOB_ID ended on:   " `hostname -s`\n')
#file.write('echo "Job $JOB_ID ended on:   " `date `\n')
#file.write('echo " "\n')
#file.write('#### submit_job.sh STOP ####\n')
file.close() 

