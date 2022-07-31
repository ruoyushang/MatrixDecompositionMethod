# Base Image
FROM ubuntu:16.04

# Create software Directory
RUN mkdir software

# Update and Upgrade Apt Repo
RUN apt-get update && apt-get upgrade -y

# Install Some Tools
RUN apt-get install -y git htop vim
RUN apt-get install -y autoconf

# Install ROOT Dependencies
RUN apt-get install -y dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev \
libxft-dev libxext-dev python libssl-dev gfortran libpcre3-dev xlibmesa-glu-dev \
libglew1.5-dev libftgl-dev libmysqlclient-dev libfftw3-dev libcfitsio-dev \
graphviz-dev libavahi-compat-libdnssd-dev libldap2-dev python-dev libxml2-dev \
libkrb5-dev libgsl0-dev libqt4-dev

# Install MySQL Client for ROOT
RUN apt-get install -y mysql-client openssl

# Install ROOT Version 5.34/38
RUN mkdir /software/ROOT
WORKDIR /software/ROOT
RUN git clone https://github.com/root-project/root.git root-git && cd root-git && git checkout v5-34-38
RUN mkdir root-v5-34-38
WORKDIR /software/ROOT/root-v5-34-38
RUN cmake -D mysql=ON ../root-git
RUN make -j40

# Initialize ROOT Environment Variables
ENV ROOTSYS "/software/ROOT/root-v5-34-38"
ENV PATH "${PATH}:${ROOTSYS}/bin"
ENV LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${ROOTSYS}/lib"
ENV DYLD_LIBRARY_PATH "${DYLD_LIBRARY_PATH}:${ROOTSYS}/lib"

# # VBF Installation
WORKDIR /software
ADD software/VBF /software/VBF
WORKDIR /software/VBF
RUN ./autogen.sh
RUN ./configure --prefix=/software/VBF
RUN make clean && make && make install

# Required Environmental Variables for EventDisplay Installation
ENV CVSROOT ":pserver:cvsuser@romulus.ucsc.edu:/home/cvsuser/VERITAS"
ENV VBFSYS "/software/VBF"
ENV EDVERSION "v483"
ENV EVNDISPSYS "/software/EventDisplay_v4"
ENV VERITAS_EVNDISP_AUX_DIR "/software/Eventdisplay_AnalysisFiles_VTS"
ENV VERITAS_DATA_DIR "/scratch/rshang"
ENV VERITAS_USER_DATA_DIR "$VERITAS_DATA_DIR"
ENV VERITAS_USER_LOG_DIR "$VERITAS_DATA_DIR"
ENV PATH "$PATH:$VBFSYS/bin"
ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:$VBFSYS/lib"

# EventDisplay v483c Installation
WORKDIR /software
ADD software/EventDisplay_v4 /software/EventDisplay_v4
#ADD software/Eventdisplay_AnalysisFiles_VTS /software/Eventdisplay_AnalysisFiles_VTS
RUN cd EventDisplay_v4 && make clean && make -j40 VTS

# Add run scripts and data
ADD software/data /software/data
ADD software/script /software/script

