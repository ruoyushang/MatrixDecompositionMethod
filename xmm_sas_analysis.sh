
stage=5

em_std_cut='(PATTERN <= 12)&&(PI in [200:12000])&&#XMMEA_EM'
ep_std_cut='(PATTERN <= 0)&&(PI in [200:12000])&&#XMMEA_EP'

# 1000 pix = 50 arcsec
test_region_cut='((X,Y) in CIRCLE(25112,22779,200))'
src_region_cut='((X,Y) in CIRCLE(24706.48, 24280.421, 340))'
bkg_region_cut='((X,Y) in ANNULUS(24706.48, 24280.421, 340,800))'
all_region_cut='((X,Y) in CIRCLE(24112,22779,800))'

em_filter=$src_region_cut' && (PI in [200:12000])&&#XMMEA_EM'
ep_filter=$src_region_cut' && (PI in [200:12000])&&#XMMEA_EP'

mytimebinsize=10

dir=$(pwd)

if test "$stage" -eq 0
then
    echo "Running stage 0"

    file=$(ls *.tar.gz)
    tar -xvf $file
    file=$(ls *.TAR)
    tar -xvf $file
    exit 0
fi

if test "$stage" -eq 1
then
    echo "Running stage 1"

    export SAS_ODF="$dir"
    cifbuild
    export SAS_CCF="$dir"/ccf.cif
    odfingest
    file=$(ls *.SAS)
    echo $file
    export SAS_ODF="$dir"/"$file"
    exit 0
fi

file=$(ls *.SAS)
export SAS_CCF="$dir"/ccf.cif
export SAS_ODF="$dir"/"$file"

#emchain and epchain creates the event files and are neccessary for timing analysis
if test "$stage" -eq 2
then
    echo "Running stage 2"
    emchain
    epchain
    exit 0
fi

mos1_fits=$(ls P*M1*EVL*)
mos2_fits=$(ls P*M2*EVL*)
pn_fits=$(ls P*PN*EVL*)

if test "$stage" -lt 4
then
    echo "Running stage 3"

    sh clean.sh
    #Filter Data By Energy. This depends on the particular viewing instrument. For MOS cameras the range is 0.2 - 12 keV while for PN it is 0.2 - 15 keV
    evselect table="$mos1_fits" withfilteredset=yes expression="$em_std_cut" filteredset=mos1_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes
    evselect table="$mos2_fits" withfilteredset=yes expression="$em_std_cut" filteredset=mos2_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes
    evselect table="$pn_fits" withfilteredset=yes expression="$ep_std_cut" filteredset=pn_filt.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

    #Create GTI files
    #mos-filter
    #pn-filter
    #exit 0

    cp mos1_filt.fits mos1_filt-clean.fits
    cp mos2_filt.fits mos2_filt-clean.fits
    cp pn_filt.fits pn_filt-clean.fits
    #em_rate_cut='RATE<=2.5e10'
    #ep_rate_cut='RATE<=2.5e10'
    em_rate_cut='RATE<=3.0'
    ep_rate_cut='RATE<=1.1'
    evselect table=mos1_filt.fits withrateset=Y rateset=mos1_ltcrv.fits maketimecolumn=Y timebinsize=$mytimebinsize makeratecolumn=Y
    tabgtigen table=mos1_ltcrv.fits expression="$em_rate_cut" gtiset=mos1_gti.fits
    evselect table=mos2_filt.fits withrateset=Y rateset=mos2_ltcrv.fits maketimecolumn=Y timebinsize=$mytimebinsize makeratecolumn=Y
    tabgtigen table=mos2_ltcrv.fits expression="$em_rate_cut" gtiset=mos2_gti.fits
    evselect table=pn_filt.fits withrateset=Y rateset=pn_ltcrv.fits maketimecolumn=Y timebinsize=$mytimebinsize makeratecolumn=Y
    tabgtigen table=pn_ltcrv.fits expression="$ep_rate_cut" gtiset=pn_gti.fits

    mos1_gti=$(ls mos1*gti.fits)
    mos2_gti=$(ls mos2*gti.fits)
    pn_gti=$(ls pn*gti.fits)
    mos1_time_cut='GTI('$mos1_gti', TIME)'
    mos2_time_cut='GTI('$mos2_gti', TIME)'
    pn_time_cut='GTI('$pn_gti', TIME)'
    mos1_clean=$(ls mos1*-clean.fits)
    mos2_clean=$(ls mos2*-clean.fits)
    pn_clean=$(ls pn*-clean.fits)

    #This filters the event files by GTI
    evselect table=$mos1_clean withfilteredset=yes expression="$mos1_time_cut" filteredset=mos1_filt_time.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes
    evselect table=$mos2_clean withfilteredset=yes expression="$mos2_time_cut" filteredset=mos2_filt_time.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes
    evselect table=$pn_clean withfilteredset=yes expression="$pn_time_cut" filteredset=pn_filt_time.fits filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes

    evselect table=mos1_filt_time.fits withrateset=Y rateset=mos1_filt_ltcrv.fits maketimecolumn=Y timebinsize=$mytimebinsize makeratecolumn=Y
    evselect table=mos2_filt_time.fits withrateset=Y rateset=mos2_filt_ltcrv.fits maketimecolumn=Y timebinsize=$mytimebinsize makeratecolumn=Y
    evselect table=pn_filt_time.fits withrateset=Y rateset=pn_filt_ltcrv.fits maketimecolumn=Y timebinsize=$mytimebinsize makeratecolumn=Y
    # dsplot table=mos1_filt_ltcrv.fits x=TIME y=RATE
    # dsplot table=mos2_filt_ltcrv.fits x=TIME y=RATE
    # dsplot table=pn_filt_ltcrv.fits x=TIME y=RATE

    #Apply barycenter correction to the GTI-filtered event Files
    barycen table=mos1_filt_time.fits:EVENTS
    barycen table=mos2_filt_time.fits:EVENTS
    barycen table=pn_filt_time.fits:EVENTS

    exit 0
fi

if test "$stage" -lt 5
then
    echo "Running stage 4"

    ## source dectection
    atthkgen atthkset=attitude.fits timestep=1
    evselect table=mos1_filt_time.fits withimageset=yes imageset=mos1image_filt.fits xcolumn=X ycolumn=Y imagebinning=imageSize ximagesize=600 yimagesize=600
    edetect_chain imagesets=mos1image_filt.fits eventsets=mos1_filt_time.fits attitudeset=attitude.fits pimin='300' pimax='10000' likemin=25 witheexpmap=yes ecf='0.878 0.220' eboxl_list=eboxlist_l_mos1.fits eboxm_list=eboxlist_m_mos1.fits eml_list=emllist_mos1.fits esp_withootset=no
    #
    ## display result
    #srcdisplay boxlistset=emllist_mos1.fits imageset=mos1image_filt.fits regionfile=regionfile_mos1.txt sourceradius=0.0166667 withregionfile=yes uselabel=yes
    # find region center and radius in physical coordinate
    # in ds9, click 'region -> list -> physical coordinate'

    #evselect table=mos1_filt_time.fits withimageset=yes imageset=mos1image_src.fits xcolumn=X ycolumn=Y imagebinning=imageSize ximagesize=600 yimagesize=600
    #echo "source location is at 245.4458333 -22.8862303"
    #echo "Src coord"
    #ecoordconv srcexp="$src_region_cut" imageset=mos1image_src.fits
    #echo "Test coord"
    #ecoordconv srcexp="$test_region_cut" imageset=mos1image_src.fits
    ## to view the image and make sure the region is defined correctly, do e.g.
    ## imgdisplay withimagefile=true imagefile=mos1image_src.fits

fi

if test "$stage" -lt 6
then
    echo "Running stage 5"

    infile_mos1=mos1_filt_time.fits
    infile_mos2=mos2_filt_time.fits
    infile_pn=pn_filt_time.fits

    rm mos1_src.fits
    rm mos2_src.fits
    rm pn_src.fits
    rm mos1_bkg.fits
    rm mos2_bkg.fits
    rm pn_bkg.fits
    rm mos1_src_pi.fits
    rm mos2_src_pi.fits
    rm pn_src_pi.fits
    rm mos1_bkg_pi.fits
    rm mos2_bkg_pi.fits
    rm pn_bkg_pi.fits

    #EXTRACT PROCESSED EVTS by region
    #with LC
    evselect table=$infile_mos1 energycolumn='PI' withfilteredset=yes filteredset='mos1_src.fits' keepfilteroutput=yes filtertype='expression' expression="$em_filter" withspectrumset=yes spectrumset='mos1_src_pi.fits' spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999
    evselect table=$infile_mos2 energycolumn='PI' withfilteredset=yes filteredset='mos2_src.fits' keepfilteroutput=yes filtertype='expression' expression="$em_filter" withspectrumset=yes spectrumset='mos2_src_pi.fits' spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999
    evselect table=$infile_pn energycolumn='PI' withfilteredset=yes filteredset='pn_src.fits' keepfilteroutput=yes filtertype='expression' expression="$ep_filter" withspectrumset=yes spectrumset='pn_src_pi.fits' spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479

    evselect table=$infile_mos1 energycolumn='PI' withfilteredset=yes filteredset='mos1_bkg.fits' keepfilteroutput=yes filtertype='expression' expression="$bkg_region_cut" withspectrumset=yes spectrumset='mos1_bkg_pi.fits' spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999
    evselect table=$infile_mos2 energycolumn='PI' withfilteredset=yes filteredset='mos2_bkg.fits' keepfilteroutput=yes filtertype='expression' expression="$bkg_region_cut" withspectrumset=yes spectrumset='mos2_bkg_pi.fits' spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999
    evselect table=$infile_pn energycolumn='PI' withfilteredset=yes filteredset='pn_bkg.fits' keepfilteroutput=yes filtertype='expression' expression="$bkg_region_cut" withspectrumset=yes spectrumset='pn_bkg_pi.fits' spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479

    evselect table=mos1_src.fits withimageset=yes imageset=mos1image_filt.fits xcolumn=X ycolumn=Y imagebinning=imageSize ximagesize=600 yimagesize=600
    # to view the image, do e.g.
    # imgdisplay withimagefile=true imagefile=mos1image_filt.fits

    exit 0
fi

if test "$stage" -eq 6
then
    echo "Running stage 6"
    #These are steps to create the rmf and arf files. These are correction matrices that calibrate data
    backscale spectrumset=mos1_src_pi.fits
    badpixlocation=mos1_filt_time.fits

    rmfgen rmfset=mos1_rmf.fits spectrumset=mos1_src_pi.fits
    arfgen arfset=mos1_arf.fits spectrumset=mos1_src_pi.fits withrmfset=yes rmfset=mos1_rmf.fits withbadpixcorr=yes badpixlocation=mos1_filt_time.fits

    backscale spectrumset=mos2_src_pi.fits
    badpixlocation=mos2_filt_time.fits

    rmfgen rmfset=mos2_rmf.fits spectrumset=mos2_src_pi.fits
    arfgen arfset=mos2_arf.fits spectrumset=mos2_src_pi.fits withrmfset=yes rmfset=mos2_rmf.fits withbadpixcorr=yes badpixlocation=mos2_filt_time.fits
fi

#if test "$stage" -eq 9
#then
#    echo "Running stage 9"
#    #HENreadevents reads the fits-formatted event files and converts them to .nc files to be processed by Hendrics
#    HENreadevents mos1_src.fits
#    HENreadevents mos2_src.fits
#    HENreadevents pn_src.fits
#
#    mos1_nc=$(ls mos1*nc)
#    mos2_nc=$(ls mos2*nc)
#    pn_nc=$(ls pn*nc)
#
#    #Finally, the Hendrics files are calibrated with the correctional matrices producing the final science products for timing analysis
#    #However this is not perfromed for XMM, you perform it for every other telescope I left this here commented it out to make you aware
#    #HENcalibrate mos1_nc -r mos1_rmf.fits
#    #HENcalibrate mos2_nc -r mos2_rmf.fits
#    #HENcalibrate pn_nc   -r pn_rmf.fits
#fi

echo "Done"
