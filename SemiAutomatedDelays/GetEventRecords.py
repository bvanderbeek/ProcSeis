#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 09:53:29 2022

@author: bvanderbeek
"""

#%% Import Modules
from os import chdir as cd
from obspy.clients.fdsn import Client
from obspy import UTCDateTime

# Connect to event directory
cd('/Users/bvanderbeek/research/CASCADIA/Waveforms_S/EVT351')

#%% Steps for every trace:
# webService = Client("IRIS")
# LOC = '*'
# NW = 'CC'
# ST = 'PANH'
# AT = UTCDateTime('2012-02-03T04:10:11.583')
# CH = 'BHN'

# DT  = 300
# corners = (1/60/3,1/60,5,25)
# fq = 40.0

# TS = webService.get_waveforms(NW, ST, LOC, CH, AT - DT, AT + DT, attach_response=True)
# TS.remove_response(inventory=None, output='VEL',water_level=None, pre_filt=corners,zero_mean=True, taper=True, taper_fraction=0.1, plot=False)
# TS.resample(fq,window='hann')
# theFile = NW + '_' + ST + '_' + CH + '.mseed'
# TS.write(theFile, format='mseed')

# Filling data gaps via 
# Stream.merge(method=0, fill_value='interpolate')

#%% Input
DT  = 300 # Time before and after arrival of interest
LOC = '*' # Location code
# Specify channels to download in order of priority
channel_1 = ('BHE','HHE','BH1','HH1')
channel_2 = ('BHN','HHN','BH2','HH2')
channel_3 = ('BHZ','HHZ','BH3','HH3')
# Instrument response removal parameters
corners = (1/60/3,1/60,5,25)
# Desired sampling frequency (Hz)
fq = 40.0

#%% Main Routine

# Initialise web service client
webService = Client("IRIS") # FDSN web service client

# Open file to save missed traces
fid = open("missedTraces.txt",mode="a")

# Loop over arrivals and retrieve waveforms
NARR = sum(1 for line in open('ArrivalList.dat'))
k = 0
with open('ArrivalList.dat') as A:
    for line in A:
        k = k + 1
        print("Starting arrival " + str(k) + " of " + str(NARR))
        
        # Parse arrival file
        nline = line.split(' ')
        NW = nline[0]
        ST = nline[1]
        AT = UTCDateTime(nline[2] + 'T' + nline[3])
        
        # First Channel
        for CH1 in channel_1:
            # Request seismogram
            try:
                TS1 = webService.get_waveforms(NW, ST, LOC, CH1, AT - DT, AT + DT, attach_response=True)
                
                # Special treatment for COR. Use location 10.
                if (ST == 'COR') and (len(TS1) > 1):
                    TS1 = TS1.select(location='10')
                # Also NEW
                if (ST == 'NEW') and (len(TS1) > 1):
                    TS1 = TS1.select(location='00')
                
            except:
                TS1 = ()
                    
            
            # Data returned?
            N1 = len(TS1)
            if N1 == 1:
                # Remove instrument response
                TS1.remove_response(inventory=None, output='VEL',\
                                    water_level=None, pre_filt=corners,\
                                        zero_mean=True, taper=True, taper_fraction=0.1, plot=False)
                
                # Sometimes instrument response removal causes issues for write due to change in samples
                TS1.resample(fq,window='hann')
                
                # Write miniSEED file
                theFile = NW + '_' + ST + '_' + CH1 + '.mseed'
                TS1.write(theFile, format='mseed')
                
                # Exit loop once data is found
                break
            
            # Exit loop if multiple matches
            if N1 > 1:
                break
        
        # Display some results
        if N1 == 1:
            print("Returned " + CH1)
        elif N1 == 0:
            # print("No data returned for channel 1")
            print(str(k) + " No data returned for " + ST + " channel 1!",file=fid)
        else:
            # print("Multiple traces returned for channel 1")
            print(str(k) + " Multiple traces returned for " + ST + " channel 1!",file=fid)
        
        
        
        # Second Channel
        for CH2 in channel_2:
            # Request seismogram
            try:
                TS2 = webService.get_waveforms(NW, ST, LOC, CH2, AT - DT, AT + DT, attach_response=True)
                
                # Special treatment for COR. Use location 10.
                if (ST == 'COR') and (len(TS2) > 1):
                    TS2 = TS2.select(location='10')
                # Also NEW
                if (ST == 'NEW') and (len(TS2) > 1):
                    TS2 = TS2.select(location='00')
                    
            except:
                TS2 = ()
                    
            
            # Data returned?
            N2 = len(TS2)
            if N2 == 1:
                # Remove instrument response
                TS2.remove_response(inventory=None, output='VEL',\
                                    water_level=None, pre_filt=corners,\
                                        zero_mean=True, taper=True, taper_fraction=0.1, plot=False)
                # Resample
                TS2.resample(fq,window='hann')
                
                # Write miniSEED file
                theFile = NW + '_' + ST + '_' + CH2 + '.mseed'
                TS2.write(theFile, format='mseed')
                
                # Exit loop once data is found
                break
            
            # Exit loop if multiple matches
            if N2 > 1:
                break
        
        # Display some results
        if N2 == 1:
            print("Returned " + CH2)
        elif N2 == 0:
            # print("No data returned for channel 2")
            print(str(k) + " No data returned for " + ST + " channel 2!",file=fid)
        else:
            # print("Multiple traces returned for channel 2")
            print(str(k) + " Multiple traces returned for " + ST + " channel 2!",file=fid)
        
        
        
        # Third Channel
        for CH3 in channel_3:
            # Request seismogram
            try:
                TS3 = webService.get_waveforms(NW, ST, LOC, CH3, AT - DT, AT + DT, attach_response=True)
                
                # Special treatment for COR. Use location 10.
                if (ST == 'COR') and (len(TS3) > 1):
                    TS3 = TS3.select(location='10')
                # Also NEW
                if (ST == 'NEW') and (len(TS3) > 1):
                    TS3 = TS3.select(location='00')
                    
            except:
                TS3 = ()
                    
            
            # Data returned?
            N3 = len(TS3)
            if N3 == 1:
                # Remove instrument response
                TS3.remove_response(inventory=None, output='VEL',\
                                    water_level=None, pre_filt=corners,\
                                        zero_mean=True, taper=True, taper_fraction=0.1, plot=False)
                # Resample
                TS3.resample(fq,window='hann')
                
                # Write miniSEED file
                theFile = NW + '_' + ST + '_' + CH3 + '.mseed'
                TS3.write(theFile, format='mseed')
                
                # Exit loop once data is found
                break
            
            # Exit loop if multiple matches
            if N3 > 1:
                break
        
        # Display some results
        if N3 == 1:
            print("Returned " + CH3)
        elif N3 == 0:
            # print("No data returned for channel 3!")
            print(str(k) + " No data returned for " + ST + " channel 3!",file=fid)
        else:
            # print("Multiple traces returned for channel 3!")
            print(str(k) + " Multiple traces returned for " + ST + " channel 3!",file=fid)
            
# Close file
fid.close()