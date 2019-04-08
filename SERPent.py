#!/usr/bin/python

##############################################################################
# Scripted E-merlin Rfi-mitigation PypelinE for iNTerferometry - SERPent
#
# Original version 13/09/2011 - Luke Peck
#
# Updates by dmf, jmorford - numerous
#
# Update by dmf 02/2015 - modularised Lovell dropouts, zero-level dropouts and
#                         flagger functions. General de-bugging.
# Update by dmf 03/2015 - Re-parallelisation to use multiprocessing module, 
#                         parallelised over stokes as well.
# Update by dmf 03/2015 - Added ability to specify individual parameters per
#                         baseline
# Update by dmf 06/2015 - Updated the flagger function for more efficient
#                         processing. In-scan zero-level function added.
#                         Major changes to Lovell_dropout function.
# Update by dmf 09/2015 - Update the flag writing process to use the 
#                         ParselTongue table access
# Update by dmf 11/2015 - Written multi-source functionality - will now run on
#                         multiple sources using source dependent and baseline
#                         dependent parameters if requested
# Update by dmf 12/2015 - Will now load and apply previous flag tables if 
#                         requested. Can load multiple tables.
# Update by dmf 12/2016 - Will now write visibilities to file for processing
#                         to minimise memory usage. Also updated flag-writing
#                         code to group in time and frequency. This is quicker 
#                         and writes fewer flag entries.
# Update by dmf 04/2019 - Minor bugfixes. Temporary flags and data stored on
#                         disk during processing to limit memory usage.
# 
#
# Current version Aesculapian: 1.5 08/04/2019                
##############################################################################






##############################################################################
# AIPS and other Python Modules
##############################################################################


from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
## First three imports are starter modules but more might/ will be needed

from AIPSTV import AIPSTV
import Wizardry.AIPSData
from Wizardry.AIPSData import AIPSUVData, AIPSImage
## Wizardry is a useful module for clever manipulation of AIPS data

import collections
import LocalProxy, gc
from xmlrpclib import ServerProxy
import copy, optparse, os, sys
import re, string
import os.path
import multiprocessing as mp
#import numarray.ma
#import numarray.ieeespecial
#import numarray.nd_image
import inspect
import numpy as np
np.set_printoptions(threshold='nan')
import platform
import cPickle as pickle
import math, time, datetime
from time import gmtime, strftime, localtime
ti = time.time()    # To time the script

np.seterr(invalid='ignore')
import warnings
warnings.defaultaction = "always"
warnings.filterwarnings("ignore", 'Mean of empty slice.')
warnings.filterwarnings("ignore", 'invalid value encountered in double_scalar')
warnings.filterwarnings("ignore", 'Warning: converting a masked element to nan.')

#from recipe_sizeof import total_size as ts
'''
from guppy import hpy
h = hpy()
h.setref()

from pympler import muppy
all_objects = muppy.get_objects()
'''




###############################################################################
# Function definitions
################################################################################




 
#######################################################################
##### SERPent flagger multiprocess function definition ################
#######################################################################



def flagger_fits_process(send_q, rec_q, cpu):

    '''
    Parallelised flagging procedure. This reads inputs from the queue then loads selected data 
    from AIPS into python arrays then performs the Lovell dropout, zero_level dropout and main 
    flagger operations. Final flag tables and meta-data are written into pickle files in the 
    data area provided.
    '''

    # print 'Just checking'

    for value in iter(send_q.get, 'STOP'):


        ## retrieve inputs from queue ##
        #print value
        if value[0] == True:

            mem_low,flagsexist,jobcount,totaljobs,bline,nif,numif,source,nchan,ntimescans,times_array,polnames,path2folder,experiment,name,phasecal,lovell_bline,telescope,srcname,zero_level,do_lovell_cross,coadd_polarization_flags,pickle_list,flagging_options,parameters,uvname,uvklass,uvdisk,uvseq, singlesource, source_num, orderedpols, write_olds, timeflag_count, metadata, select_stokes,flag_all_stokes,coadd_zero_flags,coadd_lovell_flags=value[:]


        else:
            mem_low,flagsexist,jobcount,totaljobs,bline,nif,numif,source,nchan,ntimescans,times_array,polnames,path2folder,experiment,name,phasecal,lovell_bline,telescope,srcname,zero_level,do_lovell_cross,coadd_polarization_flags,pickle_list,flagging_options,parameters,uvname,uvklass,uvdisk,uvseq, singlesource, source_num, orderedpols, old_flags, write_olds, timeflag_count, metadata, select_stokes,flag_all_stokes,coadd_zero_flags,coadd_lovell_flags=value[:]

        
        print uvname,uvklass,uvdisk,uvseq
        uvdata = AIPSUVData(uvname,uvklass,uvdisk,uvseq)
        mycore = cpu
        JOBS = 'JOB %i of %i: ' % (jobcount, totaljobs)
        flag_count = 0
        ifnum = nif
        src = []
        subary = []
        frqid = [] 
        ifs = [] 
        chans = []
        pflagsarr = []
        antsarr = []
        time_range = []
        reasons = []

        ## Create empty arrays ready to read in the data ##

        if flagsexist and mem_low:
            fname = path2folder+str(name) + "__" + str(source) +"__"+ str(bline)+"--oldflags.npy"
            old_flags_load = np.load(fname)
            old_flags = copy.copy(old_flags_load[nif])
            del old_flags_load
        elif not flagsexist and mem_low:
            old_flags = np.zeros([npol,ntimescans,nchan],dtype=np.int8)
        else:
            old_flags = copy.copy(old_flags)
            del value[-2]

        # times_array = np.zeros([ntimescans], dtype='|S12')
        # print ntimescans[bline], bline
        # print len(times_array)
        times_array = np.array(times_array)
        appending_time = 0
        if 'RR' in polnames or 'XX' in polnames:
            if 'RR' in polnames:
                pol = 'RR'
            if 'XX' in polnames:
                pol = 'XX'
            # flag_RR = np.zeros([ntimescans, nchan])
            # print polnames[pol]
            stoke = polnames[pol]
            flag_RR = copy.deepcopy(old_flags[stoke])
            # print flag_RR.shape, amp_time_RR.shape
            appending_count_RR = 0
        if 'LL' in polnames or 'YY' in polnames:
            if 'LL' in polnames:
                pol = 'LL'
            if 'YY' in polnames:
                pol = 'YY'
            # print polnames[pol]
            stoke = polnames[pol]
            flag_LL = copy.deepcopy(old_flags[stoke])
            appending_count_LL = 0
            # print flag_LL.shape, amp_time_LL.shape
        if 'RL' in polnames or 'XY' in polnames:
            if 'RL' in polnames:
                pol = 'RL'
            if 'XY' in polnames:
                pol = 'XY'
            # flag_RR = np.zeros([ntimescans, nchan])
            # print polnames[pol]
            stoke = polnames[pol]
            flag_RL = copy.deepcopy(old_flags[stoke])
            appending_count_RL = 0
        if 'LR' in polnames or 'YX' in polnames:
            if 'LR' in polnames:
                pol = 'LR'
            if 'YX' in polnames:
                pol = 'YX'
            # flag_RR = np.zeros([ntimescans, nchan])
            # print polnames[pol]
            stoke = polnames[pol]
            flag_LR = copy.deepcopy(old_flags[stoke])
            appending_count_LR = 0
        del old_flags
        if flagsexist and mem_low:
            fname = path2folder+str(name) + "__" + str(source) +"__"+ str(bline)+"--oldflags.npy"
            os.remove(fname)



        # print stop
        ## Read in the data to the various arrays ##

        ap = time.time()
        memory_array = 0
        # print JOBS + " CPU", mycore, " Appending visibilities..."
        print JOBS + " Reading visibilities into memory..."
        # print singlesource

        #count_time = appending_time
        #print vis.time
        #times_array[count_time] = float(vis.time)
        #appending_time += 1
        visdata = []
        visinfile = open(path2folder + str(uvname)  +"__"+ str(srcname) + '__IF' + str(nif+1) + '__' + str(bline),"rb")
        checkcount = 0 
        while 1:
            try:
                visdata.append(np.load(visinfile))
                # if bline == '8-9':
                  #  checkcount += 1
                    #print 'loading', checkcount
            except:
                break
        visdata = np.array(visdata)

        # print 'loading', checkcount, JOBS, bline, visdata.shape, bline, nif
        '''
        if bline == '8-9':
            if nif == 1:
                print visdata
        '''
        #for k in xrange(len(polnames)):
        visdata = np.array(visdata)
        if 'RR' in polnames or 'XX' in polnames:
            if 'RR' in polnames:
                pol = 'RR'
            if 'XX' in polnames:
                pol = 'XX'
            amp_time_RR = copy.deepcopy(visdata[:,:,polnames[pol]])
            amp_freq_RR = copy.deepcopy(visdata[:,:,polnames[pol]])
            print JOBS, 'amp', amp_time_RR.shape
            print JOBS, 'flag', flag_RR.shape
        if 'LL' in polnames or 'YY' in polnames:
            if 'LL' in polnames:
                pol = 'LL'
            if 'YY' in polnames:
                pol = 'YY'
            amp_time_LL = copy.deepcopy(visdata[:,:,polnames[pol]])
            amp_freq_LL = copy.deepcopy(visdata[:,:,polnames[pol]])
            print JOBS, 'amp', amp_time_RR.shape
            print JOBS, 'flag', flag_RR.shape
        if 'RL' in polnames or 'XY' in polnames:
            if 'RL' in polnames:
                pol = 'RL'
            if 'XY' in polnames:
                pol = 'XY'
            amp_time_RL = copy.deepcopy(visdata[:,:,polnames[pol]])
            amp_freq_RL = copy.deepcopy(visdata[:,:,polnames[pol]])
            #print amp_time_RL
        if 'LR' in polnames or 'YX' in polnames:
            if 'LR' in polnames:
                pol = 'LR'
            if 'YX' in polnames:
                pol = 'YX'
            amp_time_LR = copy.deepcopy(visdata[:,:,polnames[pol]])
            amp_freq_LR = copy.deepcopy(visdata[:,:,polnames[pol]])
            '''
            testampinp = open('/local/python/RFI/SERPent/memory_tests_0316/testampinp_file.txt','ab')
            print >> testampinp, JOBS, 'amp', amp_time_LR
            testampinp.close()
            testflaginp = open('/local/python/RFI/SERPent/memory_tests_0316/testflaginp_file.txt','ab')
            print >> testflaginp, JOBS, 'flag', flag_LR
            testflaginp.close()
            '''
        del visdata
        visinfile.close()
        os.remove(path2folder + str(uvname)  +"__"+ str(srcname) + '__IF' + str(nif+1) + '__' + str(bline))


        ###  new code to apply previous flags from old_flags array

        if 'RR' in polnames or 'XX' in polnames:
            amp_time_RR[flag_RR!=0] = np.nan
            amp_freq_RR[flag_RR!=0] = np.nan
        if 'LL' in polnames or 'YY' in polnames:
            amp_time_LL[flag_LL!=0] = np.nan
            amp_freq_LL[flag_LL!=0] = np.nan
        if 'RL' in polnames or 'XY' in polnames:
            amp_time_RL[flag_RL!=0] = np.nan
            amp_freq_RL[flag_RL!=0] = np.nan
        if 'LR' in polnames or 'YX' in polnames:
            amp_time_LR[flag_LR!=0] = np.nan
            amp_freq_LR[flag_LR!=0] = np.nan





        # print JOBS+"CPU", mycore, "Append time (hh:mm:ss):", time2hms(time.time()-ap)
        print JOBS+ "Append time (hh:mm:ss):", time2hms(time.time()-ap)
        if 'RR' in polnames or 'XX' in polnames:
            memory_array += 3*amp_time_RR.nbytes
        if 'LL' in polnames or 'YY' in polnames:
            memory_array += 3*amp_time_LL.nbytes
        if 'RL' in polnames or 'XY' in polnames:
            memory_array += 3*amp_time_RL.nbytes
        if 'LR' in polnames or 'YX' in polnames:
            memory_array += 3*amp_time_LR.nbytes
        # print "Total memory size of arrays: %s" % array_size(memory_array)
        # print "AMPTIMELL",amp_time_LL
        # print "times", times_array

        

        # Test to see if entire IF is Nan's from the correlator for each stoke
        # print "\n"+JOBS+" CPU", mycore, 'Checking for NaN values in array...'
        print "\n"+JOBS+'Checking for NaN values in array...'
        percentage = 1.0    # 0 to 1.0
        if 'RR' in polnames or 'XX' in polnames:
            nan_count = 0
            nan_array_RR = 0
            nan_count = np.isnan(amp_time_RR).sum()
            if nan_count == (percentage*amp_time_RR.shape[0]*amp_time_RR.shape[1]):
                nan_array_RR = 1
        if 'LL' in polnames or 'YY' in polnames:
            nan_count = 0
            nan_array_LL = 0
            nan_count = np.isnan(amp_time_LL).sum()
            if nan_count == (percentage*amp_time_LL.shape[0]*amp_time_LL.shape[1]):
                nan_array_LL = 1
        if 'RL' in polnames or 'XY' in polnames:
            nan_count = 0
            nan_array_RL = 0
            nan_count = np.isnan(amp_time_RL).sum()
            if nan_count == (percentage*amp_time_RL.shape[0]*amp_time_RL.shape[1]):
                nan_array_RL = 1
        if 'LR' in polnames or 'YX' in polnames:
            nan_count = 0
            nan_array_LR = 0
            nan_count = np.isnan(amp_time_LR).sum()
            if nan_count == (percentage*amp_time_LR.shape[0]*amp_time_LR.shape[1]):
                nan_array_LR = 1


        
        keep_sep = 'no'
        keep_lov_sep = 'no'
        keep_zero_sep = 'no'
        if select_stokes == 'yes':
            if flag_all_stokes == 'yes':
                coadd_polarization_flags = 'yes'
                npflags = '1111'
                zeroflags = npflags
                lovellflags= npflags
            elif flag_all_stokes == 'no' and coadd_polarization_flags == 'yes':
                newstr=['0','0','0','0']
                for pol in polnames:
                    # print orderedpols[pol]
                    newstr[orderedpols[pol]]='1'
                    npflags = "".join(newstr)
                    zeroflags = npflags
                    lovellflags= npflags
            elif flag_all_stokes == 'no' and coadd_polarization_flags == 'no':
                keep_sep = 'yes'
                if coadd_zero_flags == 'yes':
                    newstr=['0','0','0','0']
                    for pol in polnames:
                        # print orderedpols[pol]
                        newstr[orderedpols[pol]]='1'
                        zeroflags = "".join(newstr)
                else:
                    keep_zero_sep = 'yes'
                if coadd_lovell_flags == 'yes':
                    newstr=['0','0','0','0']
                    for pol in polnames:
                        # print orderedpols[pol]
                        newstr[orderedpols[pol]]='1'
                        lovellflags = "".join(newstr)
                else:
                    keep_lov_sep = 'yes'
        elif select_stokes == 'no':
            if flag_all_stokes == 'yes':
                coadd_polarization_flags = 'yes'
                npflags = '1111'
                zeroflags = npflags
                lovellflags= npflags
            elif flag_all_stokes == 'no' and coadd_polarization_flags == 'yes':
                # coadd_polarization_flags = 'yes'
                npflags = '1111'
                zeroflags = npflags
                lovellflags= npflags




        ############################################################################
        # IN SCAN ZERO LEVEL DATA CODE. WHEN TELESCOPES DROPOUT FOR UNKNOWN REASONS#
        ############################################################################
        
        if zero_level == 'yes':
            #### print "\n CPU ", mycore, "Running Zero Level Dropout Passage on unprocessed baseline...", bline
            reason_def = 'ISZEROS ' + strftime("%d-%b-%y")
            con = np.zeros([len(times_array), 2], dtype='|S12')
            for t in xrange(len(times_array)):
                con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]

            if 'LL' in polnames or 'YY' in polnames:
                amp_time_LL, amp_freq_LL, flag_LL, zeros_LL = inscan_zero_level_func(amp_time_LL, amp_freq_LL, times_array, flag_LL,JOBS)
                # break
                if 'LL' in polnames:
                    pol = 'LL'
                elif 'YY' in polnames:
                    pol = 'YY'
                if len(zeros_LL) > 0:
                    if keep_zero_sep == 'yes':
                        newstr=['0','0','0','0']
                        # print orderedpols[pol]
                        newstr[orderedpols[pol]]='1'
                        zeroflags = "".join(newstr)
                    else:
                        zeroflags = zeroflags
                    if src:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_LL, bline, path2folder, timeflag_count, zero_flag_allif, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                    else:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_LL, bline, path2folder, timeflag_count, zero_flag_allif, flag_count)
                    '''
                    plname = str(name) + "__" + str(source) + "__"+str(bline) + "__IF" + str(nif+1)+ "__" + pol + "__inscan_zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_LL, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    '''
                    del zeros_LL
            if 'RR' in polnames or 'XX' in polnames:
                amp_time_RR, amp_freq_RR, flag_RR, zeros_RR = inscan_zero_level_func(amp_time_RR, amp_freq_RR, times_array, flag_RR,JOBS)
                # break
                if 'RR' in polnames:
                    pol = 'RR'
                elif 'XX' in polnames:
                    pol = 'XX'
                if len(zeros_RR) > 0:
                    if keep_zero_sep == 'yes':
                        newstr=['0','0','0','0']
                        # print orderedpols[pol]
                        newstr[orderedpols[pol]]='1'
                        zeroflags = "".join(newstr)
                    else:
                        zeroflags = zeroflags
                    if src:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_RR, bline, path2folder, timeflag_count, zero_flag_allif, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                    else:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_RR, bline, path2folder, timeflag_count, zero_flag_allif, flag_count)
                    '''
                    plname = str(name) + "__" + str(source) +"__"+ str(bline) + "__IF" + str(nif+1)+ "__" + pol + "__inscan_zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_RR, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    '''
                    del zeros_RR
            if 'RL' in polnames or 'XY' in polnames:
                amp_time_RL, amp_freq_RL, flag_RL, zeros_RL = inscan_zero_level_func(amp_time_RL, amp_freq_RL, times_array, flag_RL,JOBS)
                # break
                if 'RL' in polnames:
                    pol = 'RL'
                elif 'XY' in polnames:
                    pol = 'XY'
                if len(zeros_RL) > 0:
                    if keep_zero_sep == 'yes':
                        newstr=['0','0','0','0']
                        # print orderedpols[pol]
                        newstr[orderedpols[pol]]='1'
                        zeroflags = "".join(newstr)
                    else:
                        zeroflags = zeroflags
                    if src:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_RL, bline, path2folder, timeflag_count, zero_flag_allif, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                    else:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_RL, bline, path2folder, timeflag_count, zero_flag_allif, flag_count)
                    '''
                    plname = str(name) + "__" + str(source) + "__"+str(bline) + "__IF" + str(nif+1)+ "__" + pol + "__inscan_zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_RL, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    '''
                    del zeros_RL
            if 'LR' in polnames or 'YX' in polnames:
                amp_time_LR, amp_freq_LR, flag_LR, zeros_LR = inscan_zero_level_func(amp_time_LR, amp_freq_LR, times_array, flag_LR,JOBS)
                # break
                if 'LR' in polnames:
                    pol = 'LR'
                elif 'YX' in polnames:
                    pol = 'YX'
                if len(zeros_LR) > 0:
                    if keep_zero_sep == 'yes':
                        newstr=['0','0','0','0']
                        # print orderedpols[pol]
                        newstr[orderedpols[pol]]='1'
                        zeroflags = "".join(newstr)
                    else:
                        zeroflags = zeroflags
                    if src:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_LR, bline, path2folder, timeflag_count, zero_flag_allif, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                    else:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_LR, bline, path2folder, timeflag_count, zero_flag_allif, flag_count)
                    '''
                    plname = str(name) + "__"+ str(source) + "__"+ str(bline) + "__IF" + str(nif+1)+ "__" + pol + "__inscan_zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_LR, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    '''
                    del zeros_LR
            gc.collect()


        ##################################################################
        ########## End of in scan zero level dropouts run ################
        ##################################################################







        ########################################################
        #### LOVELL STATIONARY SCAN CODE. ONLY FOR E-MERLIN ####
        ########################################################
        
        
        ## Testing to see whether telescope is e-MERLIN and finding 
        ## Lovell antenna number:
        lovell_num = 0
        covbase = [] 
        reason_def = 'LOVELL ' + strftime("%d-%b-%y")
        if telescope in ('e-MERLIN', 'E-MERLIN', 'eMERLIN', 'EMERLIN'):
            #### print " Array is %s" % telescope
            for ant in xrange(len(uvdata.antennas)):
                if uvdata.antennas[ant].upper() in ('LO', 'LOV', 'LOVE', 'LOVEL', 'LOVELL'):
                    lovell_num = ant+1
                    #### print "Lovell number:", lovell_num
                    break
        # print lovell_num, bline, lovell_bline, phasecal, srcname


        ## Add pol check, only run on parallel hands unless specify do_lovell_cross

    
        if do_lovell_cross == 'no':
            if str(lovell_num) in bline and bline not in lovell_bline and srcname in phasecal:
                # print "\n"+JOBS+" CPU ", mycore, "Running Lovell Stationary Scan Passage on unprocessed Lovell baseline...", bline
                print "\n"+JOBS+"Running Lovell Stationary Scan Passage on unprocessed Lovell baseline...", bline
                if 'LL' in polnames or 'YY' in polnames:
                    amp_time_LL,amp_freq_LL,flag_LL,dropouts_LL,con = lovell_dropout_med(amp_time_LL,amp_freq_LL,times_array,flag_LL,nchan,JOBS)
                    print JOBS+'Running Lovell check on LL or YY' 
                    checkedpol = 1
                    if 'LL' in polnames:
                        pol = 'LL'
                    elif 'YY' in polnames:
                        pol = 'YY'
                    if len(dropouts_LL) > 0:
                        if keep_lov_sep == 'yes':
                            newstr=['0','0','0','0']
                            # print orderedpols[pol]
                            newstr[orderedpols[pol]]='1'
                            lovellflags = "".join(newstr)
                        else:
                            lovellflags = lovellflags
                        if src:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_LL, bline, path2folder, timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                        else:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_LL, bline, path2folder, timeflag_count, flag_count)
                        '''
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(nif+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_LL, open(path2folder + str(plname) + ".dropouts", "wb"))
                        '''
                        del dropouts_LL 
                        del con
                if 'RR' in polnames or 'XX' in polnames:
                    amp_time_RR,amp_freq_RR,flag_RR,dropouts_RR,con = lovell_dropout_med(amp_time_RR,amp_freq_RR,times_array,flag_RR,nchan,JOBS)
                    print JOBS+'Running Lovell check on RR or XX' 
                    checkedpol = 2
                    if 'RR' in polnames:
                        pol = 'RR'
                    elif 'YY' in polnames:
                        pol = 'XX'
                    if len(dropouts_RR) > 0:
                        if keep_lov_sep == 'yes':
                            newstr=['0','0','0','0']
                            # print orderedpols[pol]
                            newstr[orderedpols[pol]]='1'
                            lovellflags = "".join(newstr)
                        else:
                            lovellflags = lovellflags
                        if src:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_RR, bline, path2folder, timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                        else:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_RR, bline, path2folder, timeflag_count,flag_count)
                        '''
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(nif+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_RR, open(path2folder + str(plname) + ".dropouts", "wb"))
                        '''
                        del dropouts_RR
                        del con
                gc.collect()

                        ## below should run on all pols ##


        elif do_lovell_cross == 'yes':
            if str(lovell_num) in bline and bline not in lovell_bline and srcname in phasecal:
                print "\n"+JOBS+"Running Lovell Stationary Scan Passage on unprocessed Lovell baseline...", bline
                # print "\n"+JOBS+" CPU ", mycore, "Running Lovell Stationary Scan Passage on unprocessed Lovell baseline...", bline
                if 'LL' in polnames or 'YY' in polnames:
                    amp_time_LL,amp_freq_LL,flag_LL,dropouts_LL,con = lovell_dropout_med(amp_time_LL,amp_freq_LL,times_array,flag_LL,nchan,JOBS)
                    print JOBS+'Running Lovell check on LL or YY' 
                    checkedpol = 1
                    if 'LL' in polnames:
                        pol = 'LL'
                    elif 'YY' in polnames:
                        pol = 'YY'
                    if len(dropouts_LL) > 0:
                        if keep_lov_sep == 'yes':
                            newstr=['0','0','0','0']
                            # print orderedpols[pol]
                            newstr[orderedpols[pol]]='1'
                            lovellflags = "".join(newstr)
                        else:
                            lovellflags = lovellflags
                        if src:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_LL, bline, path2folder, timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                        else:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_LL, bline, path2folder, timeflag_count,flag_count)
                        '''
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(nif+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_LL, open(path2folder + str(plname) + ".dropouts", "wb"))
                        '''
                        del dropouts_LL
                        del con
                if 'RR' in polnames or 'XX' in polnames:
                    amp_time_RR,amp_freq_RR,flag_RR,dropouts_RR,con = lovell_dropout_med(amp_time_RR,amp_freq_RR,times_array,flag_RR,nchan,JOBS)
                    print JOBS+'Running Lovell check on RR or XX' 
                    checkedpol = 2
                    if 'RR' in polnames:
                        pol = 'RR'
                    elif 'YY' in polnames:
                        pol = 'XX'
                    if len(dropouts_RR) > 0:
                        if keep_lov_sep == 'yes':
                            newstr=['0','0','0','0']
                            # print orderedpols[pol]
                            newstr[orderedpols[pol]]='1'
                            lovellflags = "".join(newstr)
                        else:
                            lovellflags = lovellflags
                        if src:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_RR, bline, path2folder, timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                        else:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_RR, bline, path2folder, timeflag_count,flag_count )
                        '''
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(nif+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_RR, open(path2folder + str(plname) + ".dropouts", "wb"))
                        '''
                        del dropouts_RR
                        del con
                if 'RL' in polnames or 'XY' in polnames:
                    amp_time_RL,amp_freq_RL,flag_LL,dropouts_RL,con = lovell_dropout_med(amp_time_RL,amp_freq_RL,times_array,flag_RL,nchan,JOBS)
                    print JOBS+'Running Lovell check on LL or YY' 
                    checkedpol = 1
                    if 'RL' in polnames:
                        pol = 'RL'
                    elif 'XY' in polnames:
                        pol = 'XY'
                    if len(dropouts_RL) > 0:
                        if keep_lov_sep == 'yes':
                            newstr=['0','0','0','0']
                            # print orderedpols[pol]
                            newstr[orderedpols[pol]]='1'
                            lovellflags = "".join(newstr)
                        else:
                            lovellflags = lovellflags
                        if src:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_RL, bline, path2folder, timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                        else:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_RL, bline, path2folder, timeflag_count,flag_count)
                        '''
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(nif+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_RL, open(path2folder + str(plname) + ".dropouts", "wb"))
                        '''
                        del dropouts_RL
                        del con
                if 'LR' in polnames or 'YX' in polnames:
                    amp_time_LR,amp_freq_LR,flag_LR,dropouts_LR,con = lovell_dropout_med(amp_time_LR,amp_freq_LR,times_array,flag_LR,nchan,JOBS)
                    print JOBS+'Running Lovell check on RR or XX' 
                    checkedpol = 2
                    if 'RR' in polnames:
                        pol = 'RR'
                    elif 'YY' in polnames:
                        pol = 'XX'
                    if len(dropouts_LR) > 0:
                        if keep_lov_sep == 'yes':
                            newstr=['0','0','0','0']
                            # print orderedpols[pol]
                            newstr[orderedpols[pol]]='1'
                            lovellflags = "".join(newstr)
                        else:
                            lovellflags = lovellflags
                        if src:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_LR, bline, path2folder, timeflag_count,flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                        else:
                            timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_lovell_fgcols(source_num, lovellflags, ifnum, channum, reason_def, con, dropouts_LR, bline, path2folder, timeflag_count,flag_count)
                        '''
                        plname = str(name) + "__" + str(source) + "__" + str(bline) + "_IF" + str(nif+1) + '__' + pol + "__lovell_dummy"
                        # infolist = np.array([[source], [con], [dropouts]])
                        pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                        pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                        pickle.dump(dropouts_LR, open(path2folder + str(plname) + ".dropouts", "wb"))
                        '''
                        del dropouts_LR
                        del con
                gc.collect()



        #######################################################
        ################# End of Lovell dropout ###############
        #######################################################







        
        ####################################################################
        # ZERO LEVEL DATA CODE. WHEN TELESCOPES DROPOUT FOR UNKNOWN REASONS#
        ####################################################################
        
        if zero_level == 'yes':
            #### print "\n CPU ", mycore, "Running Zero Level Dropout Passage on unprocessed baseline...", bline
            con = np.zeros([len(times_array), 2], dtype='|S12')
            reason_def = 'ZEROS ' + strftime("%d-%b-%y")
            for t in xrange(len(times_array)):
                con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]
            

            if 'LL' in polnames or 'YY' in polnames:
                amp_time_LL, amp_freq_LL, flag_LL, zeros_LL = zero_level_func(amp_time_LL, amp_freq_LL, flag_LL,JOBS)
                # break
                if 'LL' in polnames:
                    pol = 'LL'
                elif 'YY' in polnames:
                    pol = 'YY'
                if len(zeros_LL) > 0:
                    if keep_zero_sep == 'yes':
                        newstr=['0','0','0','0']
                        # print orderedpols[pol]
                        newstr[orderedpols[pol]]='1'
                        zeroflags = "".join(newstr)
                    else:
                        zeroflags = zeroflags

                    if src:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_LL, bline, path2folder, timeflag_count, zero_flag_allif, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                    else:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_LL, bline, path2folder, timeflag_count, zero_flag_allif, flag_count)
                    '''
                    plname = str(name) + "__" + str(source) + "__" + str(bline) + "__IF" + str(nif+1)+ "__" + pol + "__zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_LL, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    '''
                    del zeros_LL
            if 'RR' in polnames or 'XX' in polnames:
                amp_time_RR, amp_freq_RR, flag_RR, zeros_RR = zero_level_func(amp_time_RR, amp_freq_RR, flag_RR,JOBS)
                # break
                if 'RR' in polnames:
                    pol = 'RR'
                elif 'XX' in polnames:
                    pol = 'XX'
                if len(zeros_RR) > 0:
                    if keep_zero_sep == 'yes':
                        newstr=['0','0','0','0']
                        # print orderedpols[pol]
                        newstr[orderedpols[pol]]='1'
                        zeroflags = "".join(newstr)
                    else:
                        zeroflags = zeroflags

                    if src:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_RR, bline, path2folder, timeflag_count, zero_flag_allif, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                    else:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_RR, bline, path2folder, timeflag_count, zero_flag_allif, flag_count)
                    '''
                    plname = str(name) + "__" + str(source) + "__" + str(bline) + "__IF" + str(nif+1)+ "__" + pol + "__zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_RR, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    '''
                    del zeros_RR
            if 'RL' in polnames or 'XY' in polnames:
                amp_time_RL, amp_freq_RL, flag_RL, zeros_RL = zero_level_func(amp_time_RL, amp_freq_RL, flag_RL,JOBS)
                # break
                if 'RL' in polnames:
                    pol = 'RL'
                elif 'XY' in polnames:
                    pol = 'XY'
                if len(zeros_RL) > 0:
                    if keep_zero_sep == 'yes':
                        newstr=['0','0','0','0']
                        # print orderedpols[pol]
                        newstr[orderedpols[pol]]='1'
                        zeroflags = "".join(newstr)
                    else:
                        zeroflags = zeroflags

                    if src:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_RL, bline, path2folder, timeflag_count, zero_flag_allif, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                    else:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_RL, bline, path2folder, timeflag_count, zero_flag_allif, flag_count)
                    '''
                    plname = str(name) + "__" + str(source) + "__" + str(bline) + "__IF" + str(nif+1)+ "__" + pol + "__zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_RL, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    '''
                    del zeros_RL
            if 'LR' in polnames or 'YX' in polnames:
                amp_time_LR, amp_freq_LR, flag_LR, zeros_LR = zero_level_func(amp_time_LR, amp_freq_LR, flag_LR,JOBS)
                # break
                if 'LR' in polnames:
                    pol = 'LR'
                elif 'YX' in polnames:
                    pol = 'YX'
                if len(zeros_LR) > 0:
                    if keep_zero_sep == 'yes':
                        newstr=['0','0','0','0']
                        # print orderedpols[pol]
                        newstr[orderedpols[pol]]='1'
                        zeroflags = "".join(newstr)
                    else:
                        zeroflags = zeroflags

                    if src:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_LR, bline, path2folder, timeflag_count, zero_flag_allif, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                    else:
                        timeflag_count, flag_count, src, subary, frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_zeros_fgcols(source_num, zeroflags, ifnum, numif, channum, reason_def, con, zeros_LR, bline, path2folder, timeflag_count, zero_flag_allif, flag_count)
                    '''
                    plname = str(name) + "__" + str(source) + "__" + str(bline) + "__IF" + str(nif+1)+ "__" + pol + "__zeros_dummy"
                    # infolist = np.array([[source], [con], [dropouts]])
                    pickle.dump(source, open(path2folder + str(plname) + ".info", "wb"))
                    pickle.dump(con, open(path2folder + str(plname) + ".con", "wb"))
                    pickle.dump(zeros_LR, open(path2folder + str(plname) + ".zeros", "wb"))
                    print JOBS+'checking pol:', pol
                    '''
                    del zeros_LR
            gc.collect()




        ##################################################################
        ############## End of zero level dropouts run ####################
        ##################################################################
        





        ##############################################################
        ############ Sum Threshold flagging sequence: ################
        ### Execute the Flagger for each polarisation and baseline ##


        reason_def = reason()
        if 'RR' in polnames or 'XX' in polnames:
            if 'RR' in polnames:
                pol = 'RR'
            if 'XX' in polnames:
                pol = 'XX'
            # print " \n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
            print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
            if amp_time_RR.shape[0] > 0 and amp_freq_RR.shape[0] > 0:
                # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and IF: %i" % (pol, bline, nif+1)
                flag_RR = flagger_new(amp_time_RR, amp_freq_RR, flag_RR,flagging_options,parameters,JOBS)
                if write_olds:
                    flag_RR[flag_RR==-9] = 9
                # only need one direction test as the other will also be equal to 0.
            del amp_time_RR
            del amp_freq_RR
            if nan_array_RR == 0:
                if keep_sep == 'yes':
                    newstr=['0','0','0','0']
                    # print orderedpols[pol]
                    newstr[orderedpols[pol]]='1'
                    pflags = "".join(newstr)
                else:
                    pflags = npflags
                if src:
                    print 'yes',src
                    timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_flag_columns(pflags, flag_RR, times_array, source_num, ifnum, reason_def, bline, timeflag_count, flag_count,src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                else:
                    print 'no'
                    timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_flag_columns(pflags, flag_RR, times_array, source_num, ifnum, reason_def, bline, timeflag_count, flag_count)
                print src
                '''
                #### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                pname = str(name) +"__"+ str(source) + "__IF" + str(nif+1) + "__" + str(bline) + "__" + pol
                pickle.dump(flag_RR, open(path2folder + str(pname)+".p", "wb"))
                pickle_list[pname] = [source, nif+1, bline, times_array]
                pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                #np.save(path2folder+str(pname)+".p",flag_RR)
                #pickle_list[pname] = [source, nif+1, bline, times_array]
                #np.save(path2folder+str(pname)+".info",times_array)
                # print len(flag)
                '''
                del flag_RR

        if 'LL' in polnames or 'YY' in polnames:
            if 'LL' in polnames:
                pol = 'LL'
            if 'YY' in polnames:
                pol = 'YY'
            print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
            # print "\n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
            if amp_time_LL.shape[0] > 0 and amp_freq_LL.shape[0] > 0:
                # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and IF: %i" % (pol, bline, nif+1)
                flag_LL = flagger_new(amp_time_LL, amp_freq_LL, flag_LL,flagging_options,parameters,JOBS)
                if write_olds:
                    flag_LL[flag_LL==-9] = 9
            del amp_time_LL
            del amp_freq_LL
            if nan_array_LL == 0:
                if keep_sep == 'yes':
                    newstr=['0','0','0','0']
                    # print orderedpols[pol]
                    newstr[orderedpols[pol]]='1'
                    pflags = "".join(newstr)
                else:
                    pflags = npflags
                if src:
                    timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_flag_columns(pflags, flag_LL, times_array, source_num, ifnum, reason_def, bline, timeflag_count, flag_count,src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                else:
                    timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_flag_columns(pflags, flag_LL, times_array, source_num, ifnum, reason_def, bline, timeflag_count, flag_count)
                '''
                #### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                pname = str(name)  +"__"+ str(source) + "__IF" + str(nif+1) + "__" + str(bline) + "__" + pol
                pickle.dump(flag_LL, open(path2folder + str(pname)+".p", "wb"))
                pickle_list[pname] = [source, nif+1, bline, times_array]
                pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                #np.save(path2folder+str(pname)+".p",flag_LL)
                #pickle_list[pname] = [source, nif+1, bline, times_array]
                #np.save(path2folder+str(pname)+".info",times_array)
                # print len(flag)
                # print len(flag)
                '''
                del flag_LL

        if 'RL' in polnames or 'XY' in polnames:
            if 'RL' in polnames:
                pol = 'RL'
            if 'XY' in polnames:
                pol = 'XY'
            print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
            # print "\n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
            if amp_time_RL.shape[0] > 0 and amp_freq_RL.shape[0] > 0:
                # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and IF: %i" % (pol, bline, nif+1)
                #flag_RL = flagger_new(amp_time_RL, amp_freq_RL, flag_RL,flagging_options,parameters,JOBS)
                if write_olds:
                    flag_RL[flag_RL==-9] = 9
            del amp_time_RL
            del amp_freq_RL
            if nan_array_RL == 0:
                if keep_sep == 'yes':
                    newstr=['0','0','0','0']
                    # print orderedpols[pol]
                    newstr[orderedpols[pol]]='1'
                    pflags = "".join(newstr)
                else:
                    pflags = npflags
                if src:
                    timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_flag_columns(pflags, flag_RL, times_array, source_num, ifnum, reason_def, bline, timeflag_count, flag_count,src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                else:
                    timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_flag_columns(pflags, flag_RL, times_array, source_num, ifnum, reason_def, bline, timeflag_count, flag_count)
                '''
                #### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                pname = str(name)  +"__"+ str(source) + "__IF" + str(nif+1) + "__" + str(bline) + "__" + pol
                pickle.dump(flag_RL, open(path2folder + str(pname)+".p", "wb"))
                pickle_list[pname] = [source, nif+1, bline, times_array]
                pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                # print len(flag)
                '''
                del flag_RL

        if 'LR' in polnames or 'YX' in polnames:
            if 'LR' in polnames:
                pol = 'LR'
            if 'YX' in polnames:
                pol = 'YX'
            # print "\n"+JOBS+" CPU %i: FLAGGING %s AMPLITUDES for source %s" % (mycore, pol, source)
            print " \n"+JOBS+"FLAGGING %s AMPLITUDES for source %s" % (pol, source)
            if amp_time_LR.shape[0] > 0 and amp_freq_LR.shape[0] > 0:
                # print "\n"+JOBS+" CPU %i: Flagging %s amplitudes for baseline: %s and IF: %i" % (mycore, pol, bline, nif+1)
                print "\n"+JOBS+"Flagging %s amplitudes for baseline: %s and IF: %i" % (pol, bline, nif+1)
                #flag_LR = flagger_new(amp_time_LR,amp_freq_LR,flag_LR,flagging_options,parameters,JOBS)
                if write_olds:
                    flag_LR[flag_LR==-9] = 9
            del amp_time_LR
            del amp_freq_LR
            if nan_array_LR == 0:
                if keep_sep == 'yes':
                    newstr=['0','0','0','0']
                    # print orderedpols[pol]
                    newstr[orderedpols[pol]]='1'
                    pflags = "".join(newstr)
                else:
                    pflags = npflags
                if src:
                    timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_flag_columns(pflags, flag_LR, times_array, source_num, ifnum, reason_def, bline, timeflag_count, flag_count,src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
                else:
                    timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons = write_flag_columns(pflags, flag_LR, times_array, source_num, ifnum, reason_def, bline, timeflag_count, flag_count)
                '''
                #### print "\n Source : %s, IF: %s, Baseline: %s, Stokes: %s" % (name, nif+1, bline, pol)
                pname = str(name) +"__"+ str(source)  + "__IF" + str(nif+1) + "__" + str(bline) + "__" + pol
                pickle.dump(flag_LR, open(path2folder + str(pname)+".p", "wb"))
                pickle_list[pname] = [source, nif+1, bline, times_array]
                pickle.dump(times_array, open(path2folder + str(pname)+".info", "wb"))
                # print len(flag)
                '''
                del flag_LR


        

        ## Now write all columns to fits file for this source,bline,if ##
        # need to pass metadata array #
        # write_fits_fg(experiment,metadata,path2folder,nif,bline,src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons)
        pname = str(name) +"__"+ str(source)  + "__IF" + str(nif+1) + "__" + str(bline)
        pickle.dump(timeflag_count, open(path2folder + str(pname)+".fgtimes", "wb"))
        pickle.dump(src, open(path2folder + str(pname)+".src", "wb"))
        pickle.dump(subary, open(path2folder + str(pname)+".subary", "wb"))
        pickle.dump(frqid, open(path2folder + str(pname)+".frqid", "wb"))
        pickle.dump(ifs, open(path2folder + str(pname)+".ifs", "wb"))
        pickle.dump(chans, open(path2folder + str(pname)+".chans", "wb"))
        pickle.dump(pflagsarr, open(path2folder + str(pname)+".pflag", "wb"))
        pickle.dump(antsarr, open(path2folder + str(pname)+".ants", "wb"))
        pickle.dump(time_range, open(path2folder + str(pname)+".trange", "wb"))
        pickle.dump(reasons, open(path2folder + str(pname)+".reasons", "wb"))
        
        del reasons


        # num_flagged_vis = (flag != 0).sum()
        #### print "Total number of flags: %i for IF %i and stokes %s" % (num_flagged_vis, nif+1,pol)
        # tsum = 0
        # tsum += ntimescans[bline]
        # total_tscans = 0
        # total_tscans += tsum
        # print "Total number of amplitudes: %i" % (tsum*nif*nchan*npol)
        # if num_flagged_vis > 0:
           #  print 'tsum is', tsum
            #print "Percentage of amplitudes flagged for IF %i, baseline %s: %f%%, stokes %s" % (nif+1, bline, (100.0*num_flagged_vis/(tsum*nchan)),pol)
                
        finish_str = "\n"+JOBS+" Finished flagging Source :" + source + ", IF:" +str(nif)+", Baseline:"+bline
        # print finish_str

        # return finish_str
        
        rec_q.put(finish_str, cpu)
        gc.collect()



#############################################################################
######### End of main parallel flagger_process definition ###################
#############################################################################



#############################################################################
################## Start of New SumThreshold flagger definition #################

def flagger_new(amp_time, amp_freq, flag, flagging_options, parameters,JOBS):

    '''
    SumThreshold flagging process. Takes in desired amplitude and flag data arrays 
    then uses the flagging_options and parameters to flag the data.
    Outputs updated amplitude and flag arrays.
    '''

    global_total = 0
    global_time = 0
    global_freq = 0



    if flagging_options == 'default':
        aggressiveness_first_run = 4
        max_subset_first_run = 2
        aggressiveness_second_run = 4
        max_subset_second_run = 2
        rho = 1.5
        kickout_sigma_level = 3.0
    else:
        aggressiveness_first_run = parameters[0]
        max_subset_first_run = parameters[1]
        aggressiveness_second_run = parameters[2]
        max_subset_second_run = parameters[3]
        rho = parameters[4]
        kickout_sigma_level = parameters[5]
        # print parameters[5]


    print 'Current parms are:'
    print aggressiveness_first_run, max_subset_first_run,aggressiveness_second_run,max_subset_second_run,rho,kickout_sigma_level
         
    # Initial run with subset size = 1 to remove v. strong RFI
    # Note only one direction needed window = 1

    # Need to consider any values the correlator has flagged. 
    # Written them in SERPent as 'nan'


    # Updated so that read_in_flags will turn old_flag data to nan's but maintain 
    # their positions in the flag array so they can still be written out if required
    
    flag[np.isnan(amp_time)] = np.where(flag[np.isnan(amp_time)]==-9, flag[np.isnan(amp_time)],-1)
    # flag[np.isnan(amp_time)] = -1.0

    # flag = np.where(amp_time!="nan",flag,-1.0)
    #print np.where(amp_time[0]!="nan")

    median = np.median(amp_time[~np.isnan(amp_time)])
    final_mad = mad(amp_time)
    chi_1 = median + (20*final_mad)

    i = 1
    if flagging_options == 'default':
        if freq_band(frequency) in ('L Band' or 'S Band'):
            kickout = median + 3.0*final_mad
            # print "Kickout aggressive for low frequency..."
        else:
            kickout = median + 3.0*final_mad
            # print "Kickout less aggressive for high frequency..."
            ### User defined kickout sigma level:
    else:
        # print kickout_sigma_level
        kickout = median + kickout_sigma_level*final_mad



    while i < amp_time.shape[0]+1 and i <= 1:
        timing = time.time()
        limit = chi_1/(rho**(math.log(i,2)))
        ## The below condition is to avoid flagging good data if the 
        ## thresholding gets too close to the mean
        if limit <= kickout:
            break
        # print amp_time, '\n \n \n \n'
        amp_time = np.where(flag!=1,amp_time,limit)
        # print amp_time
        amp_freq = np.where(flag!=1,amp_freq,limit)
        for row in xrange(1,amp_time.shape[0]-(i-2)):
            if i == 1:
                flag[row-1,:][np.where(amp_time[row-1,:]>limit)] = np.where(flag[row-1,:][np.where(amp_time[row-1,:]>limit)]<0,flag[row-1,:][np.where(amp_time[row-1,:]>limit)],1)
        i *= 2

    
    count_i = 0
    # count_i = np.count_nonzero(flag)
    count_i = (flag > 0).sum()
    old_count = (flag != 0).sum()
    # print "Flagged with initial run:", count_i, 'old_count', old_count

    ## Run the full Threshold algorithm for the first time.
    ## Flagged vis will be at lowest threshold for the following statistics:

    hm = time.time()
    median = np.median(amp_time[~np.isnan(amp_time)])
    final_mad = mad(amp_time)
    chi_1 = median + (aggressiveness_first_run*final_mad)

    ## Reset the amplitude arrays for the new statistics.
    ## Flagged visibilities get set to the first threshold level because 
    ## the one within the while loop is for setting limits within that 
    ## sumthreshold run
    
    amp_time = np.where(flag<=0,amp_time,chi_1)
    amp_freq = np.where(flag<=0,amp_freq,chi_1)

    ###amp_time = np.where(flag==0.0,amp_time,chi_1)
    ###amp_freq = np.where(flag==0.0,amp_freq,chi_1)

    '''
    testflag2out = open('/local/python/RFI/SERPent/memory_tests_0316/testflag2out_file.txt','ab')
    print >> testflag2out, JOBS , 'flag', flag
    testflag2out.close()

    testamp2out = open('/local/python/RFI/SERPent/memory_tests_0316/testamp2out_file.txt','ab')
    print >> testamp2out, JOBS , 'amp', amp_time
    testamp2out.close()
    
    print 'type', amp_time.dtype
    '''
    i = 1
    # print "\n FLAGGING IN TIME..."
    while i < amp_time.shape[0]+1 and i < max_subset_first_run+1:
        timing = time.time()
        limit = chi_1/(rho**(math.log(i,2)))
        ## The below condition is to avoid flagging good data if the 
        ## thresholding gets too close to the mean
        if limit <= kickout:
            break
        amp_time = np.where(flag!=2,amp_time,limit)
        amp_time = np.where(flag!=-2,amp_time,limit)
        for row in xrange(1, amp_time.shape[0]-(i-2)):
            if i == 1:
                flag[row-1,:][np.where(amp_time[row-1,:]>limit)] = np.where(flag[row-1,:][np.where(amp_time[row-1,:]>limit)]<0,flag[row-1,:][np.where(amp_time[row-1,:]>limit)],2)
                flag[row-1,:][np.where(amp_time[row-1,:]>limit)] = np.where(flag[row-1,:][np.where(amp_time[row-1,:]>limit)]>=0,flag[row-1,:][np.where(amp_time[row-1,:]>limit)],-2)
            else:
                nsum = amp_time[row-1:row-1+i,:].sum(0)
                flag[row-1:row-1+i,np.where(nsum/i>limit)] = np.where(flag[row-1:row-1+i,np.where(nsum/i>limit)]<0,flag[row-1:row-1+i,np.where(nsum/i>limit)],2)
                flag[row-1:row-1+i,np.where(nsum/i>limit)] = np.where(flag[row-1:row-1+i,np.where(nsum/i>limit)]>=0,flag[row-1:row-1+i,np.where(nsum/i>limit)],-2)
        i *= 2

    '''
    testflag2aftout = open('/local/python/RFI/SERPent/memory_tests_0316/testflag2aftout_file.txt','ab')
    print >> testflag2aftout, JOBS , 'flag', flag
    testflag2aftout.close()
    '''

    count_t = (flag == 2).sum()
    print "Flagged freq:", count_t

    i = 1
    # print "\n FLAGGING IN FREQUENCY..."
    while i < amp_freq.shape[1]+1 and i < max_subset_first_run+1:
        timing = time.time()
        limit = chi_1/(rho**(math.log(i,2)))
        ## The below condition is to avoid flagging good data if the 
        ## thresholding gets too close to the mean
        if limit <= kickout:
            break
        amp_freq = np.where(flag!=3,amp_freq,limit)
        amp_freq = np.where(flag!=-3,amp_freq,limit)
        for col in xrange(1, amp_freq.shape[1]-(i-2)):
            if i == 1:
                flag[:,col-1][np.where(amp_freq[:,col-1]>limit)] = np.where(flag[:,col-1][np.where(amp_freq[:,col-1]>limit)]<0,flag[:,col-1][np.where(amp_freq[:,col-1]>limit)],3)
                flag[:,col-1][np.where(amp_freq[:,col-1]>limit)] = np.where(flag[:,col-1][np.where(amp_freq[:,col-1]>limit)]>=0,flag[:,col-1][np.where(amp_freq[:,col-1]>limit)],-3)
            else:  
                nsum = amp_freq[:,col-1:col-1+i].sum(1)
                flag[np.where(nsum/i>limit),col-1:col-1+i] = np.where(flag[np.where(nsum/i>limit),col-1:col-1+i]<0,flag[np.where(nsum/i>limit),col-1:col-1+i],3)
                flag[np.where(nsum/i>limit),col-1:col-1+i] = np.where(flag[np.where(nsum/i>limit),col-1:col-1+i]>=0,flag[np.where(nsum/i>limit),col-1:col-1+i],-3)
        i *= 2


    count_f = (flag == 3).sum()
    # print "Flagged freq:", count_f

    
    ## second run of flagger removing already found strong RFI flags 
    ## and calculating thresholds with new stats
    ## if (really high max value a certain magnitude above average/median 
    ## then run flagger again removing the strong RFIs from the arrays for 
    ## the sumThreshold...)

    median1 = np.median(amp_time)
    median = np.median(amp_time[~np.isnan(amp_time)])
    final_mad = mad(amp_time)
    maxvalue = np.amax(amp_time)
    chi_1 = median + (aggressiveness_second_run*final_mad)


    #print median1, median
    
    amp_time = np.where(flag<=0,amp_time,chi_1)
    amp_freq = np.where(flag<=0,amp_freq,chi_1)


    ###amp_time = np.where(flag==0.0,amp_time,chi_1)
    ###amp_freq = np.where(flag==0.0,amp_freq,chi_1)

    i = 1
    # print "\n FLAGGING IN TIME..."
    if flagging_options == 'default':
        if freq_band(frequency) in ('L Band' or 'S Band'):
            kickout = median + 3.0*final_mad
            # print "Kickout aggressive for low frequency..."
        else:
            kickout = median + 3.0*final_mad
            # print "Kickout less aggressive for high frequency..."
    ## User defined kickout sigma level:
    else:
        kickout = median + kickout_sigma_level*final_mad
    while i < amp_time.shape[0]+1 and i < max_subset_second_run+1:
        if maxvalue > kickout:
            if count_t > 0 or count_f > 0:
                pass
            else:
                print JOBS+" Nothing flagged. Skipping run..."
                break
        else:
            print JOBS+" Max value < kickout. Skipping run..."
            break
        timing = time.time()
        limit = chi_1/(rho**(math.log(i,2)))
 
        ## The below condition is to avoid flagging good data if the 
        ## thresholding gets too close to the mean
        if limit <= kickout:
            #print 'breaking...'
            break
        amp_time = np.where(flag!=4,amp_time,limit)
        amp_time = np.where(flag!=-4,amp_time,limit)
        
        for row in xrange(1, amp_time.shape[0]-(i-2)):
            if i == 1:
                flag[row-1,:][np.where(amp_time[row-1,:]>limit)] = np.where(flag[row-1,:][np.where(amp_time[row-1,:]>limit)]<0,flag[row-1,:][np.where(amp_time[row-1,:]>limit)],4)
                flag[row-1,:][np.where(amp_time[row-1,:]>limit)] = np.where(flag[row-1,:][np.where(amp_time[row-1,:]>limit)]>=0,flag[row-1,:][np.where(amp_time[row-1,:]>limit)],-4)
            else:
                nsum = amp_time[row-1:row-1+i,:].sum(0)
                flag[row-1:row-1+i,np.where(nsum/i>limit)] = np.where(flag[row-1:row-1+i,np.where(nsum/i>limit)]<0,flag[row-1:row-1+i,np.where(nsum/i>limit)],4)
                flag[row-1:row-1+i,np.where(nsum/i>limit)] = np.where(flag[row-1:row-1+i,np.where(nsum/i>limit)]>=0,flag[row-1:row-1+i,np.where(nsum/i>limit)],-4)
        i *= 2


    count_t = (flag == 4).sum()
    #print "Flagged time:", count_t

    i = 1
    # print "\n FLAGGING IN FREQUENCY..."
    if flagging_options == 'default':
        if freq_band(frequency) in ('L Band' or 'S Band'):
            kickout = median + 3.0*final_mad
            # print "Kickout aggressive for low frequency..."
        else:
            kickout = median + 3.0*final_mad
            # print "Kickout less aggressive for high frequency..."
    ## User defined kickout sigma level:
    else:
        kickout = median + kickout_sigma_level*final_mad
    while i < amp_freq.shape[1]+1 and i < max_subset_second_run+1:
        if maxvalue > kickout:
            if count_t > 0 or count_f > 0:
                pass
            else:
                print JOBS+" Nothing flagged. Skipping run..."
                break
        else:
            print JOBS+" Max value < kickout. Skipping run..."
            break
        timing = time.time()
        limit = chi_1/(rho**(math.log(i,2)))
        ## The below condition is to avoid flagging good data if the 
        ## thresholding gets too close to the mean
        if limit <= kickout:
            break
        amp_freq = np.where(flag!=5,amp_freq,limit)
        amp_freq = np.where(flag!=-5,amp_freq,limit)
        for col in xrange(1, amp_freq.shape[1]-(i-2)):
            if i == 1:
                flag[:,col-1][np.where(amp_freq[:,col-1]>limit)] = np.where(flag[:,col-1][np.where(amp_freq[:,col-1]>limit)]<0,flag[:,col-1][np.where(amp_freq[:,col-1]>limit)],5)
                flag[:,col-1][np.where(amp_freq[:,col-1]>limit)] = np.where(flag[:,col-1][np.where(amp_freq[:,col-1]>limit)]>=0,flag[:,col-1][np.where(amp_freq[:,col-1]>limit)],-5)
            else:  
                nsum = amp_freq[:,col-1:col-1+i].sum(1)
                flag[np.where(nsum/i>limit),col-1:col-1+i] = np.where(flag[np.where(nsum/i>limit),col-1:col-1+i]<0,flag[np.where(nsum/i>limit),col-1:col-1+i],5)
                flag[np.where(nsum/i>limit),col-1:col-1+i] = np.where(flag[np.where(nsum/i>limit),col-1:col-1+i]>=0,flag[np.where(nsum/i>limit),col-1:col-1+i],-5)
        i *= 2


    count_f = (flag == 5).sum()
    # print "Flagged freq:", count_f
    

    # flagged_total = np.count_nonzero(flag)
    flagged_total = (flag > 0.0).sum()
    #### print "Total number of flagged amplitudes: %i" % flagged_total
    # global global_total
    global_total += flagged_total
    # global global_time
    global_time += count_t
    # global global_freq
    global_freq += count_f
    # print "CPU", mycore, "time (secs):", time.time() - hm
    #### print "time (secs):", time.time() - hm

    return (flag)


################ End of New SumThreshold flagger definition ##################
##########################################################################








##########################################################################
############ Start of in scan zero-level drop-out definition #############
##########################################################################


def inscan_zero_level_func(amp_time,amp_freq,times_array,flag_arr,JOBS):
    
    '''
    Inscan_zero_level dropout function. Takes as inputs the amplitude and 
    associated time and flag arrays. Performs a running mad stastical analysis 
    looking for data consistently less than 3*mad within a single scan.
    outputs updated amplitude and flag arrays.

    '''

    zeros = []

 
    con = np.zeros([len(times_array), 2], dtype='|S12')
    # print times_array
    test = len(times_array)
    # print test
        
    for t in xrange(len(times_array)):
        con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]
            
    int_time = np.zeros([len(times_array)-1])
    for t in xrange(1, len(times_array)):
        int_time[t-1] = float(times_array[t]) - float(times_array[t-1])
                        
    step = float(np.median(int_time))
    del int_time
    print JOBS+" Correlator/ Averaged integration time is:", time2dhms(step), 'seconds'
        
    dropouts = []

    scantimes = []
    start = 0
    for j in xrange(1,len(times_array)):
        if float(con[j][0]) - float(con[j-1][0]) > 5*step:
            end = j
            scantimes.append([start,end])
            start = end
        elif j == len(times_array)-1:
            end = j+1
            # print 'yes',end
            scantimes.append([start,end])
    
    scantimes = np.array(scantimes)
    # print scantimes

    for i in xrange(len(scantimes)):
        start = scantimes[i][0]
        end = scantimes[i][1]

    
        newamp_time = amp_time[start:end]
        mdat = np.ma.masked_array(newamp_time,np.isnan(newamp_time))
        # print 'MDAT is',len(mdat),mdat.shape, mdat
        single = np.ma.median(mdat, axis=1)*10**6

        #single = single

        zero_avg = np.ma.average(single[~np.isnan(single)])
        zero_std = np.ma.std(single[~np.isnan(single)])
        zero_median = np.ma.median(single[~np.isnan(single)])
        zero_mad = mad_1d(single[~np.isnan(single)])
        # print zero_avg,zero_std,zero_median,zero_mad

        # print single
        # print 'checking zeros',zeros
        for v in xrange(len(single)):
            # print v, v+start, single[v],zero_avg,zero_std,zero_median,zero_mad
            if single[v] < (0 + 3*zero_std) and single[v] < (zero_avg - 2*zero_std):
                # print 'time1:', times_array[v+start],v; alt_times(times_array[v+start])
                if v+start not in zeros:
                    # print v+start
                    zeros.append(v+start)
            if single[v] < (zero_median - 3*zero_mad):
                # print 'time2:', times_array[v+start],v; alt_times(times_array[v+start])
                if v+start not in zeros:
                    zeros.append(v+start)
            if single[v] < (zero_median - 3*zero_mad) and single[v] < (0 + 2*zero_mad):
                # print 'time3:', times_array[v+start],v; alt_times(times_array[v+start])
                if v not in zeros:
                    zeros.append(v+start)
            
        ## To get the Cambridge zero levels which affect 75% of observation?
        #if single[v] < (0 + 5*zero_std):
        # if single[v] < (zero_median - zero_mad):
        # if v not in zeros:
        # zeros.append(v)
        # print "Dropouts flagged:", pol, 'IF'+str(j+1), len(zeros)
        # print sorted(zeros), len(zeros)

    #for i in xrange(len(single)):
     #   if i in zeros:
      #      amp_time[i,:] = zero_median
       #     amp_freq[i,:] = zero_median
            # print  bline,str(j+1),pol+':',sorted(zeros), len(zeros)

    zeros = sorted(zeros)
    
    amp_time[zeros] = 'nan'
    amp_freq[zeros] = 'nan'
    flag_arr[zeros] = -6.
    #print UHOH
    # print len(amp_time)
    # print len(flag_arr)
    # return(amp_time, amp_freq, times, flag_arr, zeros)
    # print amp_time,amp_freq,flag_arr
    del con
    return(amp_time, amp_freq, flag_arr, zeros)



########## End of in-scan zero-level drop-out definition #############
######################################################################








##########################################################################
################ Start of zero-level drop-out definition #################
##########################################################################


def zero_level_func(amp_time,amp_freq,flag_arr,JOBS):

    '''
    Zero_level dropout function. Takes as inputs the amplitude and 
    associated time and flag arrays. Performs a running mad stastical analysis 
    looking for data consistently less than 3*mad over all scans.
    Outputs updated amplitude and flag arrays.

    '''

    
    zeros = []

    mdat = np.ma.masked_array(amp_time,np.isnan(amp_time))
    # print 'MDAT is',len(mdat),mdat.shape, mdat
    single = np.ma.median(mdat, axis=1)*10**6
    del mdat

    # con = np.zeros([len(times_array), 2], dtype='|S12')
    # for t in xrange(len(times_array)):
        # con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]

    #single = single

    zero_avg = np.ma.average(single[~np.isnan(single)])
    zero_std = np.ma.std(single[~np.isnan(single)])
    zero_median = np.ma.median(single[~np.isnan(single)])
    zero_mad = mad_1d(single[~np.isnan(single)])
    # print 'zero_avg,zero_std,zero_median,zero_mad', zero_avg,zero_std,zero_median,zero_mad,pol
    # print (0 + 3*zero_std), (zero_avg - 2*zero_std),pol
    # print (zero_median - 2*zero_mad),pol
    # print (zero_median - 2*zero_mad), (0 + 2*zero_mad),pol
    # print (zero_median - zero_mad)

    # print single
    # print 'checking zeros',zeros
    for v in xrange(len(single)):
        if single[v] < (0 + 3*zero_std) and single[v] < (zero_avg - 2*zero_std):
            #print 'time1:', times_array[v],v; alt_times(times_array[v])
            if v not in zeros:
                zeros.append(v)
        if single[v] < (zero_median - 3*zero_mad):
            #print 'time2:', times_array[v],v; alt_times(times_array[v])
            if v not in zeros:
                zeros.append(v)
        if single[v] < (zero_median - 3*zero_mad) and single[v] < (0 + 2*zero_mad):
            #print 'time3:', times_array[v],v; alt_times(times_array[v])
            if v not in zeros:
                zeros.append(v)
            
        ## To get the Cambridge zero levels which affect 75% of observation?
        #if single[v] < (0 + 5*zero_std):
        # if single[v] < (zero_median - zero_mad):
        # if v not in zeros:
        # zeros.append(v)

    # print "Dropouts flagged:", pol, 'IF'+str(j+1), len(zeros)
    # print sorted(zeros), len(zeros)

    for i in xrange(len(single)):
        if i in zeros:
            amp_time[i,:] = zero_median
            amp_freq[i,:] = zero_median
            # print  bline,str(j+1),pol+':',sorted(zeros), len(zeros)


    zeros = sorted(zeros)
    print zeros
    print len(zeros)
    print amp_time.shape
    print flag_arr.shape
    # print amp_time
    # print zeros
    # times_array = np.delete(times_array, zeros)
    # amp_time = np.delete(amp_time, zeros, axis=0)
    # amp_freq = np.delete(amp_freq, zeros, axis=0)
    # flag_arr = np.delete(flag_arr, zeros, axis=0)

    amp_time[zeros] = 'nan'
    amp_freq[zeros] = 'nan'
    flag_arr[zeros] = -6.


    #print len(amp_time)
    #print len(flag_arr)
    # return(amp_time, amp_freq, times, flag_arr, zeros)
    del single
    return(amp_time, amp_freq, flag_arr, zeros)



############## End of zero-level drop-out definition #################
######################################################################





#####################################################################
############# Start of mean Lovell drop-out definition ##############
#####################################################################


def lovell_dropout(amp_time,amp_freq,times_array,flag_arr,nchan,JOBS):

        '''
        Lovell dropout function. Takes as inputs the amplitude and 
        associated time and flag arrays. Performs a running mad stastical analysis 
        looking for data consistently low over all scans. Using the mean as reference.
        This is designed to pick up instances where the eMERLIN Lovell telescope skips 
        phase calibration observations becasue of slow-slewing times.
        Outputs updated amplitude and flag arrays.
        
        '''


        ## Select a single channel to perform test on. Choosing channel in the 
        ## middle of the IF for best quality amps...
        # test_chan = int(nchan/2)
        
        ## choose central 50% chans to average ##
        chan_start = int(nchan/4)
        chan_end = chan_start + int(nchan/2)
        print JOBS+'Checking channels:', chan_start, chan_end
        
        con = np.zeros([len(times_array), 2], dtype='|S12')
          
        for t in xrange(len(times_array)):
            # try:
                # print float(times_array[t])
            # except ValueError,e:
                # print "error",e,"on line",t,'value',times_array[t]
            con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]
            
        int_time = np.zeros([len(times_array)-1])
        for t in xrange(1, len(times_array)):
            int_time[t-1] = float(times_array[t]) - float(times_array[t-1])
                        
        # print 'Lovell polarisation is', pol, j+1 
        step = float(np.median(int_time))
        print JOBS+"Correlator/ Averaged integration time is:", time2dhms(step), 'seconds'
        del int_time

        checkedpol = 0
        dropouts = []
        
        onechans = np.ma.masked_array(amp_time[:,chan_start:chan_end],np.isnan(amp_time[:,chan_start:chan_end]))
        one = np.ma.median(onechans, axis=1)*10**6
        del onechans

        # print 'continue'
        # print 'checking pol', checkedpol
        
        levels = {'high': [0,0], 'low': [10**10,0]}
        
        ## Firstly cycle through baseline and find mean, std of scans and 
        ## set the lowest and highest levels.
        
        # print "ONE:",mycore, one
        count_drop = 0
        begin_step = 0
        end_step = 0
        for i in xrange(1,len(one)):
            # print i, one[i], con[i][0], con[i][1], con[i-1][0]
            # Some test level of 5 * integration time...
            if float(con[i][0]) - float(con[i-1][0]) > 5*step:
                # print "Scan has changed level"
                # print "loop1"
                end_step = i
                newo = one[begin_step:end_step]
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                if avg > levels['high'][0]:
                    levels['high'][:] = [avg, s]
                if avg < levels['low'][0]:
                    levels['low'][:] = [avg, s]
                begin_step = end_step
                count_drop += 1
                # print 'loop1 levels:', levels
            if i == len(one)-1:
                # print "Last scan"
                # print "loop2"
                end_step = i
                newo = one[begin_step:end_step+1]
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                if avg > levels['high'][0]:
                    levels['high'][:] = [avg, s]
                if avg < levels['low'][0]:
                    levels['low'][:] = [avg, s]
                begin_step = end_step
                # print 'loop2 levels:', levels
                # print "avg:", avg,s


        # now flag all scan levels within lowest avg + lowest s
        count_drop = 0
        begin_step = 0
        end_step = 0
        # print len(one)
        for i in xrange(1,len(one)):
            # print i, one[i], con[i][0], con[i][1], one[0],one[1]
            if float(con[i][0]) - float(con[i-1][0]) > 5*step:
                # print "loop3"
                end_step = i
                newo = one[begin_step:end_step]
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                newmed = np.ma.median(newo[~np.isnan(newo)])
                newmad = mad_1d(newo[~np.isnan(newo)])
                print JOBS+"Average: %f, std: %f" % (avg, s), newmed,newmad
                # print levels
                if levels['low'][0] + levels['low'][1] < levels['high'][0]-(2*levels['high'][1]):
                    if avg < levels['low'][0] + (6*levels['low'][1]) and avg > levels['low'][0] - (6*levels['low'][1]):
                        one[begin_step:end_step] = 0.0
                        # print 'flagging'
                        count_drop += 1
                begin_step = end_step
            if i == len(one)-1:
                # print "loop4"
                end_step = i
                newo = one[begin_step:end_step+1]
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                # print "Average: %f, std: %f" % (avg, s)
                # print levels
                #### print "Average: %f, std: %f" % (avg, s)
                if levels['low'][0] + levels['low'][1] < levels['high'][0]-(2*levels['high'][1]):
                    if avg < levels['low'][0] + (6*levels['low'][1]) and avg > levels['low'][0] - (6*levels['low'][1]):
                        one[begin_step:end_step+1] = 0.0
                        # print 'flagging'
                        count_drop += 1
                begin_step = end_step
                # print end_step, one
                
                # Now collect indexes for amp == 0.0
        
        # print count_drop
        if count_drop > 0:
            dropouts = []
            
            # print 'len one is', len(one)
            #dropoutsf = np.array(np.where(one==0.0))
            #dropouts = dropoutsf[0]
            #print dropoutsf

            for i in xrange(len(one)):
                #print i,one[i]
                # print i+1, one[i], con[i][0], con[i][1]
                if one[i] == 0.0:
                    dropouts.append(i)
                #if i == len(one)-1:
                 #   if one[i] == 0.0:
                  #      dropouts.append(i+1)
        # print dropouts

        del newo
        del one
        dropouts = sorted(dropouts)
        # print dropouts
        # print len(dropouts)
        # old_times = copy.copy(times_array)

        # times_array = np.delete(times_array, dropouts)
        # amp_time = np.delete(amp_time, dropouts, axis=0)
        # amp_freq = np.delete(amp_freq, dropouts, axis=0)
        # flag_arr = np.delete(flag_arr, dropouts, axis=0)
        # print len(amp_time), len(amp_freq), len(flag_arr),len(one)

        # old_times = np.delete(old_times, dropouts)
        amp_time[dropouts] = 'nan'
        amp_freq[dropouts] = 'nan'
        flag_arr[dropouts] = -7.
 
                

        #### print "Before Lovell:", test, "after Lovell:", len(old_times)
        return(amp_time,amp_freq,flag_arr,dropouts,con)
                                



############## End of mean Lovell droput definition ################
####################################################################








#####################################################################
########### Start of median Lovell drop-out definition ##############
#####################################################################


def lovell_dropout_med(amp_time,amp_freq,times_array,flag_arr,nchan,JOBS):


        '''
        Lovell dropout median function. Takes as inputs the amplitude and 
        associated time and flag arrays. Performs a running sad stastical analysis 
        looking for data consistently low over all scans. Using the median as reference.
        This is designed to pick up instances where the eMERLIN Lovell telescope skips 
        phase calibration observations becasue of slow-slewing times.
        Outputs updated amplitude and flag arrays.
        
        '''



        ## Select a single channel to perform test on. Choosing channel in the 
        ## middle of the IF for best quality amps...
        # test_chan = int(nchan/2)
        
        ## choose central 50% chans to average ##
        chan_start = int(nchan/4)
        chan_end = chan_start + int(nchan/2)
        print JOBS+'Checking channels:', chan_start, chan_end
        
        con = np.zeros([len(times_array), 2], dtype='|S12')
         
        for t in xrange(len(times_array)):
            # try:
                # print float(times_array[t])
            # except ValueError,e:
                # print "error",e,"on line",t,'value',times_array[t]
            con[t][:] = [float(times_array[t]), time2dhms(float(times_array[t]))]
            
        int_time = np.zeros([len(times_array)-1])
        for t in xrange(1, len(times_array)):
            int_time[t-1] = float(times_array[t]) - float(times_array[t-1])
                        
        # print 'Lovell polarisation is', pol, j+1 
        step = float(np.median(int_time))
        print JOBS+"Correlator/ Averaged integration time is:", time2dhms(step), 'seconds'
        del int_time

        checkedpol = 0
        dropouts = []
        
        onechans = np.ma.masked_array(amp_time[:,chan_start:chan_end],np.isnan(amp_time[:,chan_start:chan_end]))
        one = np.ma.median(onechans, axis=1)*10**6
        del onechans

        # print 'continue'
        # print 'checking pol', checkedpol
        
        levels = {'high': [0,0], 'low': [10**10,0]}
        
        ## Firstly cycle through baseline and find mean, std of scans and 
        ## set the lowest and highest levels.
        
        # print "ONE:",mycore, one
        count_drop = 0
        begin_step = 0
        end_step = 0
        for i in xrange(1,len(one)):
            # print i, one[i], con[i][0], con[i][1], con[i-1][0]
            # Some test level of 5 * integration time...
            if float(con[i][0]) - float(con[i-1][0]) > 5*step:
                # print "Scan has changed level"
                # print "loop1"
                end_step = i
                newo = one[begin_step:end_step]
                avg = np.ma.median(newo[~np.isnan(newo)])
                s = mad_1d(newo[~np.isnan(newo)])
                if avg > levels['high'][0]:
                    levels['high'][:] = [avg, s]
                if avg < levels['low'][0]:
                    levels['low'][:] = [avg, s]
                begin_step = end_step
                count_drop += 1
                # print 'loop1 levels:', levels
            if i == len(one)-1:
                # print "Last scan"
                # print "loop2"
                end_step = i
                newo = one[begin_step:end_step+1]
                avg = np.ma.median(newo[~np.isnan(newo)])
                s = mad_1d(newo[~np.isnan(newo)])
                if avg > levels['high'][0]:
                    levels['high'][:] = [avg, s]
                if avg < levels['low'][0]:
                    levels['low'][:] = [avg, s]
                begin_step = end_step
                # print 'loop2 levels:', levels
                # print "avg:", avg,s


        # now flag all scan levels within lowest avg + lowest s
        count_drop = 0
        begin_step = 0
        end_step = 0
        # print len(one)
        for i in xrange(1,len(one)):
            # print i, one[i], con[i][0], con[i][1], one[0],one[1]
            if float(con[i][0]) - float(con[i-1][0]) > 5*step:
                # print "loop3"
                end_step = i
                newo = one[begin_step:end_step]
                med = np.ma.median(newo[~np.isnan(newo)])
                mad = mad_1d(newo[~np.isnan(newo)])
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                # print "Average: %f, std: %f" % (avg, s)
                # print levels
                if levels['low'][0] + levels['low'][1] < levels['high'][0]-(2*levels['high'][1]):
                    # print levels['low'][0], levels['high'][0]-levels['high'][1]
                    # print levels['low'][0] + levels['low'][1], levels['low'][0] - levels['low'][1]
                    if avg < levels['low'][0] + (9*levels['low'][1]) and avg > levels['low'][0] - (9*levels['low'][1]):
                        # print avg, levels['low'][0] + levels['low'][1] , levels['low'][0] - levels['low'][1]
                        one[begin_step:end_step] = 0.0
                        count_drop += 1
                begin_step = end_step
            if i == len(one)-1:
                # print "loop4"
                end_step = i
                newo = one[begin_step:end_step+1]
                med = np.ma.median(newo[~np.isnan(newo)])
                mad = mad_1d(newo[~np.isnan(newo)])
                avg = np.ma.mean(newo[~np.isnan(newo)])
                s = np.ma.std(newo[~np.isnan(newo)])
                # print "Average: %f, std: %f" % (avg, s)
                # print levels
                #### print "Average: %f, std: %f" % (avg, s)
                if levels['low'][0] + levels['low'][1] < levels['high'][0]-(2*levels['high'][1]):
                    if avg < levels['low'][0] + (9*levels['low'][1]) and avg > levels['low'][0] - (9*levels['low'][1]):
                        one[begin_step:end_step+1] = 0.0
                        count_drop += 1
                begin_step = end_step
                # print end_step, one
                
                # Now collect indexes for amp == 0.0

        # print count_drop
        if count_drop > 0:
            dropouts = []
            
            # print 'len one is', len(one)
            for i in xrange(len(one)):
                # print i,one[i]
                # print i+1, one[i], con[i][0], con[i][1]
                if one[i] == 0.0:
                    dropouts.append(i)
                # if i == len(one)-1:
                #   if one[i] == 0.0:
                #    dropouts.append(i+1)
            # print dropouts

        del newo
        del one
        dropouts = sorted(dropouts)
        # old_times = copy.copy(times_array)

        # times_array = np.delete(times_array, dropouts)
        # amp_time = np.delete(amp_time, dropouts, axis=0)
        # amp_freq = np.delete(amp_freq, dropouts, axis=0)
        # flag_arr = np.delete(flag_arr, dropouts, axis=0)

        # old_times = np.delete(old_times, dropouts)
        amp_time[dropouts] = 'nan'
        amp_freq[dropouts] = 'nan'
        flag_arr[dropouts] = -7.
        # del onechans
                

        #### print "Before Lovell:", test, "after Lovell:", len(old_times)
        return(amp_time,amp_freq,flag_arr,dropouts,con)
                                



############ End of median Lovell droput definition ################
####################################################################












############################################################
### new code to write flags into fg table using wizardry ###
##### will replace each function call separately #####
############################################################

def create_blank_fg(data, path2folder, no_rows):
    '''
    Creates a blank AIPS compliant flag table
    '''
    # Need number of rows which we're going to enter.
    global fgfile
    # global fg_row
    fgfile = open(path2folder+'blank_table.fg', 'wr')#operation = 'wr'
    fgfile.flush()    # flush prints everything from the file onto screen
    # command below means print to this file
    print >> fgfile, "XTENSION= 'A3DTABLE'           / extension type\n",'BITPIX  =                    8 / printable ASCII codes\n','NAXIS   =                    2 / Table is a matrix\n','NAXIS1  =                  132 / Max. no. of characters/pass\n','NAXIS2  =             ','%7i' % no_rows,'/ Number of entries in table\n','PCOUNT  =                    0 / Random parameter count\n','GCOUNT  =                    1 / Group count\n','NOPASS  =                    1 / Number of passes thru table\n','TFIELDS =                    9 / Number of fields in each row\n',"EXTNAME = 'AIPS FG '           / AIPS table file\n",'EXTVER  =                    1 / Version Number of table\n','TBCOL1  =                    9 / Starting char. pos. of field\n',"TFORM1  = 'I11     '           / Fortran format of field  1\n",'TFDIM1  =                    1 / Dimension of field  1\n',"TTYPE1  = 'SOURCE  '           / type (heading) of field  1\n","TUNIT1  = '        '           / physical units of field  1\n",'TBCOL2  =                   20 / Starting char. pos. of field\n',"TFORM2  = 'I11     '           / Fortran format of field  2\n",'TFDIM2  =                    1 / Dimension of field  2\n',"TTYPE2  = 'SUBARRAY                '           / type (heading) of field  2\n","TUNIT2  = '        '           / physical units of field  2\n",'TBCOL3  =                   31 / Starting char. pos. of field\n',"TFORM3  = 'I11     '           / Fortran format of field  3\n",'TFDIM3  =                    1 / Dimension of field  3\n',"TTYPE3  = 'FREQ ID                 '           / type (heading) of field  3\n","TUNIT3  = '        '           / physical units of field  3\n",'TBCOL4  =                   42 / Starting char. pos. of field\n',"TFORM4  = 'I11     '           / Fortran format of field  4\n",'TFDIM4  =                    2 / Dimension of field  4\n',"TTYPE4  = 'ANTS                    '           / type (heading) of field  4\n","TUNIT4  = '        '           / physical units of field  4\n",'TBCOL5  =                   53 / Starting char. pos. of field\n',"TFORM5  = 'E18.8   '           / Fortran format of field  5\n",'TFDIM5  =                    2 / Dimension of field  5\n',"TTYPE5  = 'TIME RANGE              '           / type (heading) of field  5\n","TUNIT5  = 'DAYS    '           / physical units of field  5\n",'TBCOL6  =                   71 / Starting char. pos. of field\n',"TFORM6  = 'I9     '           / Fortran format of field  6\n",'TFDIM6  =                    2 / Dimension of field  6\n',"TTYPE6  = 'IFS                     '           / type (heading) of field  6\n","TUNIT6  = '        '           / physical units of field  6\n",'TBCOL7  =                   80 / Starting char. pos. of field\n',"TFORM7  = 'I9     '           / Fortran format of field  7\n",'TFDIM7  =                    2 / Dimension of field  7\n',"TTYPE7  = 'CHANS                   '           / type (heading) of field  7\n","TUNIT7  = '        '           / physical units of field  7\n",'TBCOL8  =                   89 / Starting char. pos. of field\n',"TFORM8  = 'X4      '           / Fortran format of field  8\n",'TFDIM8  =                    4 / Dimension of field  8\n',"TTYPE8  = 'PFLAGS                  '           / type (heading) of field  8\n","TUNIT8  = '        '           / physical units of field  8\n",'TBCOL9  =                  101 / Starting char. pos. of field\n',"TFORM9  = 'A24     '           / Fortran format of field  9\n",'TFDIM9  =                   26 / Dimension of field  9\n',"TTYPE9  = 'REASON                  '           / type (heading) of field  9\n","TUNIT9  = '        '           / physical units of field  9\n",'END\n','COL. NO.      1          2          3          4            5              6          7          8        9\n','     ROW   SOURCE     SUBARRAY   FREQ ID    ANTS         TIME RAN       IFS        CHANS      PFLAGS   REASON\n','  NUMBER                                                 DAYS\n','***BEGIN*PASS***'
    print >> fgfile, '***END*PASS***'
    fgfile.close()
    highfg = 0
    for tab in data.tables:
        if 'FG' in tab[1]:
            if tab[0] > highfg:
                highfg = tab[0]
    # print highfg


    tbin = AIPSTask('TBIN')
    tbin.outname = data.name
    tbin.outclass = data.klass
    tbin.outseq = data.seq
    tbin.outdisk = data.disk
    tbin.intext = path2folder + 'blank_table.fg'     # text file name
    tbin.bcount = 0
    tbin.ecount = 0
    #tbin.inp()
    tbin.go()

    print 'Adding empty flag table as FG',str(highfg+1)
    return (highfg+1)
    



def flag_lovellmk2(uvdata,fgtabnum,lovellmk2_flag, baseline_removal):

    '''
    Writes flag entry into selected AIPS table to flag Lovell-Mk2 baseline.
    Inputs are flag table number, data file, code switches and a list of
    baselines to write flags for.
    '''
    fg = uvdata.table('FG',fgtabnum)
    row = Wizardry.AIPSData.AIPSTableRow(fg)
    row['source'] = 0
    row['ifs'] = [1,0]
    row['pflags'] = [1,1,1,1]
    row['reason'] = [Lovell-Mk2]
    #row['ants'] = [int(bline[0]),int(bline[1])]
    row['chans'] = [1,0]
    row['time_range'] = [0.0,0.0]

    for bline in baseline_removal:
        print "Flagging Lovell-Mk2 baseline for all sources..."
        auto_bline = bline.rsplit('-')
        row['ants'] = [int(auto_bline[0]),int(auto_bline[1])]
        fg.append(row)

    lovellmk2_flag = 0
    fg.close()
    return lovellmk2_flag



def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)




def write_fg_table_times_alt(uvdata, fgtabnum, pflags, flag, times_array, source_num, ifnum, reason_def, bline, path2folder, timeflag_count, flag_count, no_rows):


    '''
    Writes flag entries into selected AIPS table using the provided amplitude and flag arrays
    Performs checks to ensure the number of entries per time and total does not exceed the limits
    set within AIPS. If so, a new flag table is appended and written.
    This writes consecutive flags in the array as a single entry and loops in frequency order to minimise 
    the number of entries written per time.
    '''

    #print flag
    ##flagarr = np.array(flag.reshape(flag.shape[0]*flag.shape[1]))
    #print flagarr
    numtimes = flag.shape[0]
    numchan = flag.shape[1]
    chancount = 0
    timecount = 0
    chanfirst = 0
    flag_count = 0
    tabledict = uvdata.tables
    fg = uvdata.table('FG',fgtabnum)
    row = Wizardry.AIPSData.AIPSTableRow(fg)
    # print source_num
    row['source'] = source_num
    row['ifs'] = [ifnum,ifnum]
    row['pflags'] = [int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])]
    row['reason'] = [reason_def]
    row['ants'] = [int(bline[0]),int(bline[1])]
    ti=time.time()
    #print len(flagarr)
    for i in xrange(numtimes):
        timecount = i
        chunks = consecutive(np.array(np.where(flag[i]>0)[0]))
        if np.array(chunks[0]).size > 0:
            for j in xrange(len(chunks)):
                #checkarr = flag[i:,chunks[j][0]:chunks[j][-1]]
                row['chans'] = [chunks[j][0]+1,chunks[j][-1]+1]
                t1 = float(times_array[timecount])-0.00000116
                t2 = float(times_array[timecount])+0.00000116
                #t1 = float(timecount)-0.00000116
                #t2 = float(timecount)+0.00000116
                row['time_range'] = [t1,t2]
                fg.append(row)
                flag[timecount, chunks[j][0]:chunks[j][-1]+1] = 0
                flag_count +=1
                if times_array[timecount] in timeflag_count:
                    timeflag_count[times_array[timecount]] += 1
                if timeflag_count[times_array[timecount]] >= 600000:
                    print 'Maximum flags per time reached, proceeding with new flag table'
                    fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                    fg.close() #close old table
                    fg = uvdata.table('FG',fgtabnum) #access new table
                    timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                    flag_count = 0
                ### Think following only applies to TBIN ###
                if flag_count >= 100000000:
                    print 'Maximum flag entries per table reached, proceeding with new flag table'
                    fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                    fg.close() #close old table
                    fg = uvdata.table('FG',fgtabnum)
                    timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                    flag_count = 0
    print flag_count
    del flag
    print 'finished alt writing flag table in '+ str(time.time()-ti)
    #del flag
    fg.close() #close table, finished with for now
    return fgtabnum, timeflag_count, flag_count






def write_fg_table_times_alt4(uvdata, fgtabnum, pflags, flag, times_array, source_num, ifnum, reason_def, bline, path2folder, timeflag_count, flag_count, no_rows):

    '''
    Writes flag entries into selected AIPS table using the provided amplitude and flag arrays
    Performs checks to ensure the number of entries per time and total does not exceed the limits
    set within AIPS. If so, a new flag table is appended and written.
    This writes consecutive flags in the array as a single entry and loops in frequency order to minimise 
    the number of entries written per time.
    '''

    #print flag
    ##flagarr = np.array(flag.reshape(flag.shape[0]*flag.shape[1]))
    #print flagarr
    numtimes = flag.shape[0]
    numchan = flag.shape[1]
    chancount = 0
    timecount = 0
    chanfirst = 0
    flag_count = 0
    tabledict = uvdata.tables
    fg = uvdata.table('FG',fgtabnum)
    row = Wizardry.AIPSData.AIPSTableRow(fg)
    # print source_num
    row['source'] = source_num
    row['ifs'] = [ifnum,ifnum]
    row['pflags'] = [int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])]
    row['reason'] = [reason_def]
    row['ants'] = [int(bline[0]),int(bline[1])]
    ti=time.time()
    totaltinew = 0.
    totaltiold = 0.
    #print len(flagarr)
    for i in xrange(numtimes):
        timecount = i
        #print np.where(flag[i]>0)[0]
        chunks = consecutive(np.array(np.where(flag[i]>0)[0]))
        if np.array(chunks[0]).size > 0:
            for j in xrange(len(chunks)):
                checkarr = np.array(flag[i:,chunks[j][0]:chunks[j][-1]+1])
                #tin3 = time.time()
                k=0
                while (checkarr[k]>0).all(axis=0):
                    #print k,'k'
                    k += 1
                    if k == len(checkarr):
                        #print k,'k2'
                        break
                row['chans'] = [chunks[j][0]+1,chunks[j][-1]+1]
                t1 = float(times_array[timecount])-0.00000116
                t2 = float(times_array[timecount+k-1])+0.00000116
                #t1 = float(timecount)-0.00000116
                #t2 = float(timecount)+0.00000116
                row['time_range'] = [t1,t2]
                fg.append(row)
                flag[timecount:timecount+k, chunks[j][0]:chunks[j][-1]+1] = 0
                flag_count +=1
                #print flag
                '''
                timeflag_count[timecount:timecount+k,1] += 1
                if np.where(timeflag_count[timecount:timecount+k,1] >= 600000)[0].size > 0:
                    #if timeflag_count[timecount:timecount+k,1].any() >= 600000:
                    print 'Maximum flags per time reached, proceeding with new flag table'
                    ## fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                    ## fg.close() #close old table
                    ## fg = uvdata.table('FG',fgtabnum) #access new table
                    timeflag_count[timecount:timecount+k,1] = 0
                    #timeflag_count = dict.fromkeys(timeflag_count+times, 0) # set all timecounts back to 0
                    flag_count = 0
                
                for times in xrange(k):
                    timeflag_count[timecount+times,1] += 1
                    if timeflag_count[timecount+times,1] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        ## fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                        ## fg.close() #close old table
                        ## fg = uvdata.table('FG',fgtabnum) #access new table
                        timeflag_count[timecount+times,1] = 0
                        #timeflag_count = dict.fromkeys(timeflag_count+times, 0) # set all timecounts back to 0
                        flag_count = 0
                '''
                for times in xrange(k):
                    if times_array[timecount+times] in timeflag_count:
                        timeflag_count[times_array[timecount+times]] += 1
                    if timeflag_count[times_array[timecount+times]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                        fg.close() #close old table
                        fg = uvdata.table('FG',fgtabnum) #access new table
                        timeflag_count = dict.fromkeys(timeflag_count+times, 0) # set all timecounts back to 0
                        flag_count = 0
                
                ### Think following only applies to TBIN ###
                if flag_count >= 100000000:
                    print 'Maximum flag entries per table reached, proceeding with new flag table'
                    fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                    fg.close() #close old table
                    fg = uvdata.table('FG',fgtabnum)
                    timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                    flag_count = 0
    # print flag_count
    del flag
    print 'finished alt writing flag table in '+ str(time.time()-ti)
    #del flag
    fg.close() #close table, finished with for now
    return fgtabnum, timeflag_count, flag_count
    










def write_fg_table_times(uvdata, fgtabnum, pflags, flag, times_array, source_num, ifnum, reason_def, bline, path2folder, timeflag_count, flag_count, no_rows):

    '''
    Writes flag entries into selected AIPS table using the provided amplitude and flag arrays
    Performs checks to ensure the number of entries per time and total does not exceed the limits
    set within AIPS. If so, a new flag table is appended and written.
    This writes consecutive flags in the array as a single entry and loops in frequency order to minimise 
    the number of entries written per time.
    '''

    #print flag
    flagarr = np.array(flag.reshape(flag.shape[0]*flag.shape[1]))
    #print flagarr
    numtimes = flag.shape[0]
    numchan = flag.shape[1]
    chancount = 0
    timecount = 0
    chanfirst = 0
    tabledict = uvdata.tables
    fg = uvdata.table('FG',fgtabnum)
    row = Wizardry.AIPSData.AIPSTableRow(fg)
    # print source_num
    row['source'] = source_num
    row['ifs'] = [ifnum,ifnum]
    row['pflags'] = [int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])]
    row['reason'] = [reason_def]
    row['ants'] = [int(bline[0]),int(bline[1])]
    #print len(flagarr)
    for i in xrange(len(flagarr)):
        #print chancount, numchan-1
        if chancount < numchan-1: 
            #print 'yes', chanfirst, chancount
            #print flagarr[i],flagarr[i+1], chancount
            if flagarr[i] <= 0.0:
                chanfirst = chancount+1
                #print flagarr[i],flagarr[i+1], chancount
                #chancount +=1
            if flagarr[i] > 0.0 and flagarr[i+1] <= 0.0:
                #print flagarr[i],flagarr[i+1], chancount
                #print 'yes2', chanfirst, chancount
                row['chans'] = [chanfirst+1,chancount+1]
                #print row['chans']
                t1 = float(times_array[timecount])-0.00000116
                t2 = float(times_array[timecount])+0.00000116
                row['time_range'] = [t1,t2]
                fg.append(row)
                flag[timecount, chanfirst:chancount+1] = 0.0
                flag_count +=1
                #chancount +=1 
                chanfirst = chancount+1
                if times_array[timecount] in timeflag_count:
                    timeflag_count[times_array[timecount]] += 1
                    if timeflag_count[times_array[timecount]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                        fg.close() #close old table
                        fg = uvdata.table('FG',fgtabnum) #access new table
                        timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                        flag_count = 0
                ### Think following only applies to TBIN ###
                if flag_count >= 100000000:
                    print 'Maximum flag entries per table reached, proceeding with new flag table'
                    fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                    fg.close() #close old table
                    fg = uvdata.table('FG',fgtabnum)
                    timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                    flag_count = 0
            #chanfirst = chancount+1
            chancount += 1
        elif chancount == numchan-1:
            #print 'no'
            if flagarr[i] > 0.0:
                #print 'no2'
                #print chanfirst, chancount+1
                row['chans'] = [chanfirst+1,chancount+1]
                #print row['chans']
                t1 = float(times_array[timecount])-0.00000116
                t2 = float(times_array[timecount])+0.00000116
                row['time_range'] = [t1,t2]
                fg.append(row)
                flag[timecount, chanfirst:chancount+1] = 0.0
                flag_count +=1
                #chanfirst = 0
                if times_array[timecount] in timeflag_count:
                    timeflag_count[times_array[timecount]] += 1
                    if timeflag_count[times_array[timecount]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        fgtabnum = create_blank_fg(uvdata, path2folder, no_rows) #load new dummy table
                        fg.close() #close old table
                        fg = uvdata.table('FG',fgtabnum)
                        timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                        flag_count = 0
                ### Think following only applies to TBIN ###
                if flag_count >= 100000000:
                    print 'Maximum flag entries per table reached, proceeding with new flag table'
                    fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                    fg.close() #close old table
                    fg = uvdata.table('FG',fgtabnum)
                    timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                    flag_count = 0
            chanfirst = 0
            chancount = 0
            timecount += 1

    del flagarr
    #del flag
    fg.close() #close table, finished with for now
    return fgtabnum, timeflag_count, flag_count



def write_lovell_times(uvdata, fgtabnum, source_num, pflags, ifnum, channum, reason_def, times_array, dropouts, bline, path2folder, timeflag_count, flag_count, no_rows):

    '''
    Writes flag entries from Lovell dropout code into selected AIPS table using the provided times
    and index (dropouts) lists.
    Performs checks to ensure the number of entries per time and total does not exceed the limits
    set within AIPS. If so, a new flag table is appended and written.
    '''


    uvdata.tables
    fg = uvdata.table('FG',fgtabnum)
    row = Wizardry.AIPSData.AIPSTableRow(fg)
    row['source'] = source_num
    row['ifs'] = [ifnum,ifnum]
    row['pflags'] = [int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])]
    row['reason'] = [reason_def]
    row['ants'] = [int(bline[0]),int(bline[1])]
    row['chans'] = [1,channum]

    start = dropouts[0]
    tag = 0
    for i in xrange(1,len(dropouts)):
        if dropouts[i] - dropouts[i-1] > 1:
            t1 = float(times_array[start][0])-0.00000116
            t2 = float(times_array[dropouts[i-1]][0])+0.00000116
            row['time_range'] = [t1,t2]
            fg.append(row)
            flag_count +=1
            ### Think following only applies to TBIN ###
            if flag_count >= 100000000:
                print 'Maximum flag entries per table reached, proceeding with new flag table'
                fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                fg.close() #close old table
                fg = uvdata.table('FG',fgtabnum)
                timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                flag_count = 0
            for j in xrange(tag, i):
                timecount = dropouts[j]
                if times_array[timecount][0] in timeflag_count:
                    timeflag_count[times_array[timecount][0]] += 1
                    if timeflag_count[times_array[timecount][0]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        fgtabnum = create_blank_fg(uvdata, path2folder, no_rows) #load new dummy table
                        fg.close() #close old table
                        fg = uvdata.table('FG',fgtabnum)
                        timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                        flag_count = 0
            start = dropouts[i]
            tag = i
            if i == len(dropouts)-1:
                t1 = float(times_array[start][0])-0.00000116
                t2 = float(times_array[dropouts[i]][0])+0.00000116
                row['time_range'] = [t1,t2]
                fg.append(row)
                flag_count +=1
                if flag_count >= 100000000:
                    print 'Maximum flag entries per table reached, proceeding with new flag table'
                    fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                    fg.close() #close old table
                    fg = uvdata.table('FG',fgtabnum)
                    timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                    flag_count = 0
                for j in xrange(tag, i+1):
                    timecount = dropouts[j]
                    if times_array[timecount][0] in timeflag_count:
                        timeflag_count[times_array[timecount][0]] += 1
                        if timeflag_count[times_array[timecount][0]] >= 600000:
                            print 'Maximum flags per time reached, proceeding with new flag table'
                            fgtabnum = create_blank_fg(uvdata, path2folder, no_rows) #load new dummy table
                            fg.close() #close old table
                            fg = uvdata.table('FG',fgtabnum)
                            timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                            flag_count = 0
        elif i == len(dropouts)-1:
            t1 = float(times_array[start][0])-0.00000116
            t2 = float(times_array[dropouts[i]][0])+0.00000116
            row['time_range'] = [t1,t2]
            fg.append(row)
            flag_count +=1
            if flag_count >= 100000000:
                print 'Maximum flag entries per table reached, proceeding with new flag table'
                fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                fg.close() #close old table
                fg = uvdata.table('FG',fgtabnum)
                timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                flag_count = 0
            for j in xrange(tag, i+1):
                timecount = dropouts[j]
                if times_array[timecount][0] in timeflag_count:
                    timeflag_count[times_array[timecount][0]] += 1
                    if timeflag_count[times_array[timecount][0]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        fgtabnum = create_blank_fg(uvdata, path2folder, no_rows) #load new dummy table
                        fg.close() #close old table
                        fg = uvdata.table('FG',fgtabnum)
                        timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                        flag_count = 0
    #catch single entries
    if len(dropouts) == 1:
        t1 = float(times_array[start][0])-0.00000116
        t2 = float(times_array[start][0])+0.00000116
        row['time_range'] = [t1,t2]
        fg.append(row)
        flag_count +=1
        if flag_count >= 100000000:
            print 'Maximum flag entries per table reached, proceeding with new flag table'
            fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
            fg.close() #close old table
            fg = uvdata.table('FG',fgtabnum)
            timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
            flag_count = 0
        timecount = dropouts[0]
        if times_array[timecount][0] in timeflag_count:
            timeflag_count[times_array[timecount][0]] += 1
            if timeflag_count[times_array[timecount][0]] >= 600000:
                print 'Maximum flags per time reached, proceeding with new flag table'
                fgtabnum = create_blank_fg(uvdata, path2folder, no_rows) #load new dummy table
                fg.close() #close old table
                fg = uvdata.table('FG',fgtabnum)
                timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                flag_count = 0
    
    
    fg.close() #done with for now
    return fgtabnum, timeflag_count, flag_count





def write_zeros_times(uvdata, fgtabnum, source_num, pflags, ifnum, channum, reason_def, times_array, zeros, bline, path2folder, timeflag_count, flag_count, zero_flag_allif, no_rows):

    '''
    Writes flag entries from Inscan_zero and Zero dropout code into selected AIPS table using the provided times
    and array index (zeros) lists.
    Performs checks to ensure the number of entries per time and total does not exceed the limits
    set within AIPS. If so, a new flag table is appended and written.
    '''


    if zero_flag_allif == 'yes':
        startif = 1
    elif zero_flag_allif == 'no':
        startif = ifnum
    
    uvdata.tables
    fg = uvdata.table('FG',fgtabnum)
    row = Wizardry.AIPSData.AIPSTableRow(fg)
    row['source'] = source_num
    row['ifs'] = [startif,ifnum]
    row['pflags'] = [int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])]
    row['reason'] = [reason_def]
    row['ants'] = [int(bline[0]),int(bline[1])]
    row['chans'] = [1,channum]

    start = zeros[0]
    tag = 0

    for i in xrange(1,len(zeros)):
        if zeros[i] - zeros[i-1] > 1:
            #print '1'
            t1 = float(times_array[start][0])-0.00000116
            t2 = float(times_array[zeros[i-1]][0])+0.00000116
            row['time_range'] = [t1,t2]
            fg.append(row)
            flag_count +=1
            if flag_count >= 100000000:
                print 'Maximum flag entries per table reached, proceeding with new flag table'
                fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                fg.close() #close old table
                fg = uvdata.table('FG',fgtabnum)
                timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                flag_count = 0
            for j in xrange(tag, i):
                timecount = zeros[j]
                if times_array[timecount][0] in timeflag_count:
                    timeflag_count[times_array[timecount][0]] += 1
                    if timeflag_count[times_array[timecount][0]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        fgtabnum = create_blank_fg(uvdata, path2folder, no_rows) #load new dummy table
                        fg.close() #close old table
                        fg = uvdata.table('FG',fgtabnum)
                        timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                        flag_count = 0
            start = zeros[i]
            tag = i
            if i == len(zeros)-1:
                #print '2'
                t1 = float(times_array[start][0])-0.00000116
                t2 = float(times_array[zeros[i]][0])+0.00000116
                row['time_range'] = [t1,t2]
                fg.append(row)
                flag_count +=1
                if flag_count >= 100000000:
                    print 'Maximum flag entries per table reached, proceeding with new flag table'
                    fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                    fg.close() #close old table
                    fg = uvdata.table('FG',fgtabnum)
                    timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                    flag_count = 0
                for j in xrange(tag, i+1):
                    timecount = zeros[j]
                    if times_array[timecount][0] in timeflag_count:
                        timeflag_count[times_array[timecount][0]] += 1
                        if timeflag_count[times_array[timecount][0]] >= 600000:
                            print 'Maximum flags per time reached, proceeding with new flag table'
                            fgtabnum = create_blank_fg(uvdata, path2folder, no_rows) #load new dummy table
                            fg.close() #close old table
                            fg = uvdata.table('FG',fgtabnum)
                            timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                            flag_count = 0
        elif i == len(zeros)-1:
            #print '3'
            t1 = float(times_array[start][0])-0.00000116
            t2 = float(times_array[zeros[i]][0])+0.00000116
            row['time_range'] = [t1,t2]
            fg.append(row)
            flag_count +=1
            if flag_count >= 100000000:
                print 'Maximum flag entries per table reached, proceeding with new flag table'
                fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                fg.close() #close old table
                fg = uvdata.table('FG',fgtabnum)
                timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                flag_count = 0
            for j in xrange(tag, i+1):
                timecount = zeros[j]
                if times_array[timecount][0] in timeflag_count:
                    timeflag_count[times_array[timecount][0]] += 1
                    if timeflag_count[times_array[timecount][0]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        fgtabnum = create_blank_fg(uvdata, path2folder, no_rows) #load new dummy table
                        fg.close() #close old table
                        fg = uvdata.table('FG',fgtabnum)
                        timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                        flag_count = 0
    #catch single entries
    if len(zeros) == 1:
        #print source_num, 'source'
        #print '4'
        #print zeros, times_array
        #print times_array[start]
        #row['reason'] = ['ISZEROS TEST']
        t1 = float(times_array[start][0])-0.00000116
        t2 = float(times_array[start][0])+0.00000116
        row['time_range'] = [t1,t2]
        #print row['time_range']
        fg.append(row)
        flag_count +=1
        #for nrow in fg:
        #    print nrow
        if flag_count >= 100000000:
            print 'Maximum flag entries per table reached, proceeding with new flag table'
            fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
            fg.close() #close old table
            fg = uvdata.table('FG',fgtabnum)
            timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
            flag_count = 0
        timecount = zeros[0]
        if times_array[timecount][0] in timeflag_count:
            timeflag_count[times_array[timecount][0]] += 1
            if timeflag_count[times_array[timecount][0]] >= 600000:
                print 'Maximum flags per time reached, proceeding with new flag table'
                fgtabnum = create_blank_fg(uvdata, path2folder, no_rows) #load new dummy table
                fg.close() #close old table
                fg = uvdata.table('FG',fgtabnum)
                timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                flag_count = 0
        
    
    fg.close() #done with for now
    return fgtabnum, timeflag_count, flag_count






#### need to update the following chunk of code ###

def write_fg_table_chans(uvdata, fgtabnum, no_rows, tab_row, pflags, flag, times_array, source_num, ifnum, bline, timeflag_count):
    '''Not being used'''

    numtimes = flag.shape[0]
    numchan = flag.shape[1]
    flagarr = np.array(flag.transpose().reshape(numchan*numtimes))            
    chancount = 0
    timecount = 0
    first = 0
    fg = uvdata.table('FG',fgtabnum)
    row = Wizardry.AIPSData.AIPSTableRow(fg)
    row['source'] = source_num
    row['ifs'] = [ifnum,ifnum]
    row['pflags'] = [int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])]
    row['reason'] = [reason_def]
    row['ants'] = [bline[0],bline[1]]
    for i in xrange(len(flagarr)):
        if timecount < numtimes-1: 
            if flagarr[i] > 0.0 and flagarr[i+1] <= 0.0:
                row['chans'] = [chancount+1,chancount+1]
                t1 = float(times_array[first])-0.00000116
                t2 = float(times_array[timecount])+0.00000116
                row['time_range'] = [t1,t2]
                first = timecount+1
                fg.append(row)
                flag[first:timecount+1, chancount+1] = 0.0
            timecount += 1
        elif timecount == numtimes-1:
            if flagarr[i] > 0.0:
                row['chans'] = [chancount+1,chancount+1]
                t1 = float(times_array[first])-0.00000116
                t2 = float(times_array[timecount])+0.00000116
                row['time_range'] = [t1,t2]
                first = timecount+1
                fg.append(row)
                flag[first:timecount+1, chancount+1] = 0.0
            timecount = 0
            chancount += 1
            
 
def write_flag_columns(pflags, flag, times_array, source_num, ifnum, reason_def, bline, timeflag_count, flag_count=0,src=[],subary=[],frqid=[],ifs=[],chans=[],pflagsarr=[],antsarr=[],time_range=[],reasons=[]):

    numtimes = flag.shape[0]
    numchan = flag.shape[1]
    chancount = 0
    timecount = 0
    chanfirst = 0
    flag_count = 0
    print 'flag cols'
    nbline = bline.split('-')

    ti=time.time()
    totaltinew = 0.
    totaltiold = 0.

    for i in xrange(numtimes):
        timecount = i
        chunks = consecutive(np.array(np.where(flag[i]>0)[0]))
        if np.array(chunks[0]).size > 0:
            for j in xrange(len(chunks)):
                checkarr = np.array(flag[i:,chunks[j][0]:chunks[j][-1]+1])
                k=0
                while (checkarr[k]>0).all(axis=0):
                    k += 1
                    if k == len(checkarr):
                        break
                chans.append([chunks[j][0]+1,chunks[j][-1]+1])
                t1 = float(times_array[timecount])-0.00000116
                t2 = float(times_array[timecount+k-1])+0.00000116
                time_range.append([t1,t2])
                flag[timecount:timecount+k, chunks[j][0]:chunks[j][-1]+1] = 0
                flag_count +=1
                #print flag
                src.append(source_num)
                subary.append(0)
                frqid.append(-1)
                ifs.append([ifnum,ifnum])
                pflagsarr.append([int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])])
                reasons.append(reason_def)
                antsarr.append([int(nbline[0]),int(nbline[1])])

                for times in xrange(k):
                    if times_array[timecount+times] in timeflag_count:
                        timeflag_count[times_array[timecount+times]] += 1
                    if timeflag_count[times_array[timecount+times]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        '''
                        fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                        fg.close() #close old table
                        fg = uvdata.table('FG',fgtabnum) #access new table
                        timeflag_count = dict.fromkeys(timeflag_count+times, 0) # set all timecounts back to 0
                        flag_count = 0
                        '''
                ### Think following only applies to TBIN ###
                if flag_count >= 100000000:
                    print 'Maximum flag entries per table reached, proceeding with new flag table'
                    '''
                    fgtabnum = create_blank_fg(uvdata, path2folder, no_rows=0) #load new dummy table
                    fg.close() #close old table
                    fg = uvdata.table('FG',fgtabnum)
                    timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                    flag_count = 0
                    '''
    # print flag_count
    # del flag
    print 'finished alt writing flag table in '+ str(time.time()-ti)
    #del flag
    # fg.close() #close table, finished with for now
    print 'src is', src
    return timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons
    


def aipsfitld(path2folder,experiment,num,userno,disk):

    infile = path2folder+experiment+'_'+str(num)+'.fits'
    fitld = AIPSTask('FITLD')
    fitld.datain = infile
    fitld.outname = experiment
    fitld.outclass = 'SERP'
    fitld.outdisk = disk
    fitld.outseq = num
    fitld.dotable = 1
    fitld.douvcomp = -1
    fitld.doconcat = -1
    fitld.inp()
    fitld.go()

    try:
        tablefits = AIPSImage(experiment,'SERP',disk,num)
    except:
        print "Can't find file"

    return experiment,'SERP',disk,num


def aipsfgcop(data1name,data1klass,data1disk,data1seq,uvdata,innum,outnum=None):

    '''
    Uses AIPs task TACOP to copy a flag table (version innum) from data1 to uvdata. An input
    for the table outversion can be used or aipsfgcop will work out the next highest number 
    available, ignoring any gaps in numbers for the attached file. E.g. if uvdata has 3 FG tables
    FG 1,2 and 4, unless told otherwise aipsfgcop will copy to FG 5 and not 3.
    '''

    tacop = AIPSTask('TACOP')
    tacop.inname = data1name
    tacop.inclass = data1klass
    tacop.inseq = data1seq
    tacop.indisk = data1disk
    tacop.inext = 'FG'
    tacop.invers = num
    tacop.ncount = 1
    tacop.outname = uvdata.name
    tacop.outclass = uvdata.klass
    tacop.outseq = uvdata.seq
    tacop.outdisk = uvdata.disk
    highfg = 0
    for tab in uvdata.tables:
        if 'FG' in tab[1]:
            if tab[0] > highfg:
                highfg = tab[0]
    # print highfg
    if outnum == None:
        outnum = highfg+1
    tacop.outvers = outnum
    tacop.go()

    return outnum


def write_fits_fg(sourcename, metadata,path2folder,nif,bline,src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons):


    '''
    Writes a dummy image fits file with an attached AIPS style flag (FG) table made from the arrays src,subary,frqid,ifs,chans,
    pflagsarr,antsarr,time_range,reasons
    '''

    obsJD,obs_date_str,obs_RA,obs_DEC,freq,dfrq,pfrq = metadata

    imdata = np.zeros((1,1,1,8,512,4,3))
    pnames = ['UU---SIN','VV---SIN','WW---SIN',
              'DATE','DATE','BASELINE','SOURCE',
              'FREQSEL','INTTIM','CORR-ID']
    pdata = [0.,0.,0.,obsJD,0.,0.,0.,0.,0.,0.]
    
    gdata = fits.GroupData(imdata,parnames=pnames,
                           pardata=pdata,bitpix=-32)
    
    prihdu = fits.GroupsHDU(gdata)
    
    prihdu.header.set('BLOCKED',value=True,comment='Tape may be blocked')
    prihdu.header.set('OBJECT',value=sourcename,comment='Source name')
    prihdu.header.set('TELESCOP',value='e-MERLIN',comment=' ')
    prihdu.header.set('INSTRUME',value='VLBA',comment=' ')
    prihdu.header.set('OBSERVER',value='Calibrat',comment=' ')
    prihdu.header.set('DATE-OBS',value=obs_date_str,
                      comment='Obs start date YYYY-MM-DD')
    
    prihdu.header.set('BSCALE',value=1.0E0,
                      comment='REAL = TAPE * BSCALE + BZERO')
    prihdu.header.set('BZERO',value=0.0E0,comment=' ')
    prihdu.header.set('BUNIT',value='UNCALIB',comment='Units of flux')
    
    prihdu.header.set('EQUINOX',value=2.0E3,comment='Epoch of RA DEC')
    prihdu.header.set('ALTRPIX',value=1.0E+0,
                      comment='Altenate FREQ/VEL ref pixel')
    
    prihdu.header.set('OBSRA',value=obs_RA,
                      comment='Antenna pointing RA')
    prihdu.header.set('OBSDEC',value=obs_DEC,
                      comment='Antenna pointing DEC')
    
    prihdu.header.set('CRVAL2',value=1.0E+0,comment=' ') 
    prihdu.header.set('CDELT2',value=1.0E+0,comment=' ') 
    prihdu.header.set('CRPIX2',value=1.0E+0,comment=' ') 
    prihdu.header.set('CROTA2',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('CRVAL3',value=-1.0E+0,comment=' ') 
    prihdu.header.set('CDELT3',value=-1.0E+0,comment=' ') 
    prihdu.header.set('CRPIX3',value=1.0E+0,comment=' ') 
    prihdu.header.set('CROTA3',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('CRVAL4',value=freq,comment=' ') 
    prihdu.header.set('CDELT4',value=dfrq,comment=' ') 
    prihdu.header.set('CRPIX4',value=pfrq,comment=' ') 
    prihdu.header.set('CROTA4',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('CRVAL5',value=1.0E+0,comment=' ') 
    prihdu.header.set('CDELT5',value=1.0E+0,comment=' ') 
    prihdu.header.set('CRPIX5',value=1.0E+0,comment=' ') 
    prihdu.header.set('CROTA5',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('CRVAL6',value=obs_RA,comment=' ') 
    prihdu.header.set('CDELT6',value=1.0E+0,comment=' ') 
    prihdu.header.set('CRPIX6',value=1.0E+0,comment=' ') 
    prihdu.header.set('CROTA6',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('CRVAL7',value=obs_DEC,comment=' ') 
    prihdu.header.set('CDELT7',value=1.0E+0,comment=' ') 
    prihdu.header.set('CRPIX7',value=1.0E+0,comment=' ') 
    prihdu.header.set('CROTA7',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL1',value=1./freq,comment=' ') 
    prihdu.header.set('PZERO1',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL2',value=1./freq,comment=' ') 
    prihdu.header.set('PZERO2',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL3',value=1./freq,comment=' ') 
    prihdu.header.set('PZERO2',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL4',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO3',value=obsJD,comment=' ') 
    
    prihdu.header.set('PSCAL5',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO4',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL6',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO6',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL7',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO7',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL8',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO8',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL9',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO9',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL10',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO10',value=0.0E+0,comment=' ') 
    
    
    col1 = fits.Column(name='SOURCE',format='1J',unit=' ',array=src)
    col2 = fits.Column(name='SUBARRAY',format='1J',unit=' ',array=subary)
    col3 = fits.Column(name='FREQ ID',format='1J',unit=' ',array=frqid)
    col4 = fits.Column(name='ANTS',format='2J',unit=' ',array=antsarr)
    col5 = fits.Column(name='TIME RANGE',format='2E',unit='DAYS',array=time_range)
    col6 = fits.Column(name='IFS',format='2J',unit=' ',array=ifs)
    col7 = fits.Column(name='CHANS',format='2J',unit=' ',array=chans)
    col8 = fits.Column(name='PFLAGS',format='4X',unit=' ',array=pflagsarr)
    col9 = fits.Column(name='REASON',format='24A',unit=' ',array=reasons)
    
    cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9])
    
    try:
        fg_hdu = fits.BinTableHDU.from_columns(cols)
    except AttributeError:
        fg_hdu = fits.new_table(cols)
    fg_hdu.header.set('EXTNAME',value='AIPS FG',
                      comment='AIPS table file')
    fg_hdu.header.set('EXTVER',value=1,
                      comment='Version number of table')
    
    hdulist = fits.HDUList([prihdu,fg_hdu])
    outfile = path2folder+sourcename+'__'+bline+'__IF'+str(nif+1)+'.fits'
    hdulist.writeto(outfile)


def write_full_fits_fg(sourcename, metadata,path2folder,nif,bline,src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons,num):


    obsJD,obs_date_str,obs_RA,obs_DEC,freq,dfrq,pfrq = metadata

    imdata = np.zeros((1,1,1,8,512,4,3))
    pnames = ['UU---SIN','VV---SIN','WW---SIN',
              'DATE','DATE','BASELINE','SOURCE',
              'FREQSEL','INTTIM','CORR-ID']
    pdata = [0.,0.,0.,obsJD,0.,0.,0.,0.,0.,0.]
    
    gdata = fits.GroupData(imdata,parnames=pnames,
                           pardata=pdata,bitpix=-32)
    
    prihdu = fits.GroupsHDU(gdata)
    
    prihdu.header.set('BLOCKED',value=True,comment='Tape may be blocked')
    prihdu.header.set('OBJECT',value=sourcename,comment='Source name')
    prihdu.header.set('TELESCOP',value='e-MERLIN',comment=' ')
    prihdu.header.set('INSTRUME',value='VLBA',comment=' ')
    prihdu.header.set('OBSERVER',value='Calibrat',comment=' ')
    prihdu.header.set('DATE-OBS',value=obs_date_str,
                      comment='Obs start date YYYY-MM-DD')
    
    prihdu.header.set('BSCALE',value=1.0E0,
                      comment='REAL = TAPE * BSCALE + BZERO')
    prihdu.header.set('BZERO',value=0.0E0,comment=' ')
    prihdu.header.set('BUNIT',value='UNCALIB',comment='Units of flux')
    
    prihdu.header.set('EQUINOX',value=2.0E3,comment='Epoch of RA DEC')
    prihdu.header.set('ALTRPIX',value=1.0E+0,
                      comment='Altenate FREQ/VEL ref pixel')
    
    prihdu.header.set('OBSRA',value=obs_RA,
                      comment='Antenna pointing RA')
    prihdu.header.set('OBSDEC',value=obs_DEC,
                      comment='Antenna pointing DEC')
    
    prihdu.header.set('CRVAL2',value=1.0E+0,comment=' ') 
    prihdu.header.set('CDELT2',value=1.0E+0,comment=' ') 
    prihdu.header.set('CRPIX2',value=1.0E+0,comment=' ') 
    prihdu.header.set('CROTA2',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('CRVAL3',value=-1.0E+0,comment=' ') 
    prihdu.header.set('CDELT3',value=-1.0E+0,comment=' ') 
    prihdu.header.set('CRPIX3',value=1.0E+0,comment=' ') 
    prihdu.header.set('CROTA3',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('CRVAL4',value=freq,comment=' ') 
    prihdu.header.set('CDELT4',value=dfrq,comment=' ') 
    prihdu.header.set('CRPIX4',value=pfrq,comment=' ') 
    prihdu.header.set('CROTA4',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('CRVAL5',value=1.0E+0,comment=' ') 
    prihdu.header.set('CDELT5',value=1.0E+0,comment=' ') 
    prihdu.header.set('CRPIX5',value=1.0E+0,comment=' ') 
    prihdu.header.set('CROTA5',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('CRVAL6',value=obs_RA,comment=' ') 
    prihdu.header.set('CDELT6',value=1.0E+0,comment=' ') 
    prihdu.header.set('CRPIX6',value=1.0E+0,comment=' ') 
    prihdu.header.set('CROTA6',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('CRVAL7',value=obs_DEC,comment=' ') 
    prihdu.header.set('CDELT7',value=1.0E+0,comment=' ') 
    prihdu.header.set('CRPIX7',value=1.0E+0,comment=' ') 
    prihdu.header.set('CROTA7',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL1',value=1./freq,comment=' ') 
    prihdu.header.set('PZERO1',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL2',value=1./freq,comment=' ') 
    prihdu.header.set('PZERO2',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL3',value=1./freq,comment=' ') 
    prihdu.header.set('PZERO2',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL4',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO3',value=obsJD,comment=' ') 
    
    prihdu.header.set('PSCAL5',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO4',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL6',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO6',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL7',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO7',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL8',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO8',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL9',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO9',value=0.0E+0,comment=' ') 
    
    prihdu.header.set('PSCAL10',value=1.0E+0,comment=' ') 
    prihdu.header.set('PZERO10',value=0.0E+0,comment=' ') 
    
    
    col1 = fits.Column(name='SOURCE',format='1J',unit=' ',array=src)
    col2 = fits.Column(name='SUBARRAY',format='1J',unit=' ',array=subary)
    col3 = fits.Column(name='FREQ ID',format='1J',unit=' ',array=frqid)
    col4 = fits.Column(name='ANTS',format='2J',unit=' ',array=antsarr)
    col5 = fits.Column(name='TIME RANGE',format='2E',unit='DAYS',array=time_range)
    col6 = fits.Column(name='IFS',format='2J',unit=' ',array=ifs)
    col7 = fits.Column(name='CHANS',format='2J',unit=' ',array=chans)
    col8 = fits.Column(name='PFLAGS',format='4X',unit=' ',array=pflagsarr)
    col9 = fits.Column(name='REASON',format='24A',unit=' ',array=reasons)
    
    cols = fits.ColDefs([col1,col2,col3,col4,col5,col6,col7,col8,col9])
    
    try:
        fg_hdu = fits.BinTableHDU.from_columns(cols)
    except AttributeError:
        fg_hdu = fits.new_table(cols)
    fg_hdu.header.set('EXTNAME',value='AIPS FG',
                      comment='AIPS table file')
    fg_hdu.header.set('EXTVER',value=1,
                      comment='Version number of table')
    
    hdulist = fits.HDUList([prihdu,fg_hdu])
    outfile = path2folder+experiment+'_'+str(num)+'.fits'
    hdulist.writeto(outfile)



def write_lovell_fgcols(uvdata, fgtabnum, source_num, pflags, ifnum, channum, reason_def, times_array, dropouts, bline, path2folder, timeflag_count, flag_count=0,src=[],subary=[],frqid=[],ifs=[],chans=[],pflagsarr=[],antsarr=[],time_range=[],reasons=[]):

    '''
    Writes flag entries from Lovell dropout code into fits file using the provided times and index (dropouts) lists.
    Performs checks to ensure the number of entries per time and total does not exceed the limits
    set within AIPS. 
    '''
    nbline = bline.split('-')
    start = dropouts[0]
    tag = 0
    for i in xrange(1,len(dropouts)):
        if dropouts[i] - dropouts[i-1] > 1:
            t1 = float(times_array[start][0])-0.00000116
            t2 = float(times_array[dropouts[i-1]][0])+0.00000116
            src.append(source_num)
            subary.append(0)
            frqid.append(-1)
            ifs.append([ifnum,ifnum])
            pflagsarr.append([int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])])
            reasons.append(reason_def)
            antsarr.append([int(nbline[0]),int(nbline[1])])
            chans.append([1,channum])
            flag_count +=1
            ### Think following only applies to TBIN ###
            if flag_count >= 100000000:
                print 'Maximum flag entries per table reached, proceeding with new flag table'
                #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                #flag_count = 0
            for j in xrange(tag, i):
                timecount = dropouts[j]
                if times_array[timecount][0] in timeflag_count:
                    timeflag_count[times_array[timecount][0]] += 1
                    if timeflag_count[times_array[timecount][0]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                        #flag_count = 0
            start = dropouts[i]
            tag = i
            if i == len(dropouts)-1:
                t1 = float(times_array[start][0])-0.00000116
                t2 = float(times_array[dropouts[i]][0])+0.00000116
                src.append(source_num)
                subary.append(0)
                frqid.append(-1)
                ifs.append([ifnum,ifnum])
                pflagsarr.append([int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])])
                reasons.append(reason_def)
                antsarr.append([int(nbline[0]),int(nbline[1])])
                chans.append([1,channum])
                flag_count +=1
                if flag_count >= 100000000:
                    print 'Maximum flag entries per table reached, proceeding with new flag table'
                    #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                    #flag_count = 0
                for j in xrange(tag, i+1):
                    timecount = dropouts[j]
                    if times_array[timecount][0] in timeflag_count:
                        timeflag_count[times_array[timecount][0]] += 1
                        if timeflag_count[times_array[timecount][0]] >= 600000:
                            print 'Maximum flags per time reached, proceeding with new flag table'
                            #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                            #flag_count = 0
        elif i == len(dropouts)-1:
            t1 = float(times_array[start][0])-0.00000116
            t2 = float(times_array[dropouts[i]][0])+0.00000116
            src.append(source_num)
            subary.append(0)
            frqid.append(-1)
            ifs.append([ifnum,ifnum])
            pflagsarr.append([int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])])
            reasons.append(reason_def)
            antsarr.append([int(nbline[0]),int(nbline[1])])
            chans.append([1,channum])
            flag_count +=1
            if flag_count >= 100000000:
                print 'Maximum flag entries per table reached, proceeding with new flag table'
                #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                #flag_count = 0
            for j in xrange(tag, i+1):
                timecount = dropouts[j]
                if times_array[timecount][0] in timeflag_count:
                    timeflag_count[times_array[timecount][0]] += 1
                    if timeflag_count[times_array[timecount][0]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                        #flag_count = 0
    #catch single entries
    if len(dropouts) == 1:
        t1 = float(times_array[start][0])-0.00000116
        t2 = float(times_array[start][0])+0.00000116
        src.append(source_num)
        subary.append(0)
        frqid.append(-1)
        ifs.append([ifnum,ifnum])
        pflagsarr.append([int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])])
        reasons.append(reason_def)
        antsarr.append([int(nbline[0]),int(nbline[1])])
        chans.append([1,channum])
        flag_count +=1
        if flag_count >= 100000000:
            print 'Maximum flag entries per table reached, proceeding with new flag table'
            #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
            #flag_count = 0
        timecount = dropouts[0]
        if times_array[timecount][0] in timeflag_count:
            timeflag_count[times_array[timecount][0]] += 1
            if timeflag_count[times_array[timecount][0]] >= 600000:
                print 'Maximum flags per time reached, proceeding with new flag table'
                #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                #flag_count = 0
    
    return timeflag_count, flag_count, src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons







def write_zeros_fgcols(source_num, pflags, ifnum, nif, channum, reason_def, times_array, zeros, bline, path2folder, timeflag_count, zero_flag_allif, flag_count=0,  src=[],subary=[],frqid=[],ifs=[],chans=[],pflagsarr=[],antsarr=[],time_range=[],reasons=[]):

    '''
    Writes flag entries from Inscan_zero and Zero dropout code into selected AIPS table using the provided times
    and array index (zeros) lists.
    Performs checks to ensure the number of entries per time and total does not exceed the limits
    set within AIPS. If so, a new flag table is appended and written.
    '''

    if zero_flag_allif == 'yes':
        startif = 1
        ifnum = nif
    elif zero_flag_allif == 'no':
        startif = ifnum

    nbline = bline.split('-')
    start = zeros[0]
    tag = 0

    for i in xrange(1,len(zeros)):
        if zeros[i] - zeros[i-1] > 1:
            #print '1'
            t1 = float(times_array[start][0])-0.00000116
            t2 = float(times_array[zeros[i-1]][0])+0.00000116
            src.append(source_num)
            subary.append(0)
            frqid.append(-1)
            ifs.append([startif,ifnum])
            pflagsarr.append([int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])])
            reasons.append(reason_def)
            antsarr.append([int(nbline[0]),int(nbline[1])])
            chans.append([1,channum])
            flag_count +=1
            if flag_count >= 100000000:
                print 'Maximum flag entries per table reached, proceeding with new flag table'
                #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                #flag_count = 0
            for j in xrange(tag, i):
                timecount = zeros[j]
                if times_array[timecount][0] in timeflag_count:
                    timeflag_count[times_array[timecount][0]] += 1
                    if timeflag_count[times_array[timecount][0]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                        #flag_count = 0
            start = zeros[i]
            tag = i
            if i == len(zeros)-1:
                #print '2'
                t1 = float(times_array[start][0])-0.00000116
                t2 = float(times_array[zeros[i]][0])+0.00000116
                src.append(source_num)
                subary.append(0)
                frqid.append(-1)
                ifs.append([startif,ifnum])
                pflagsarr.append([int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])])
                reasons.append(reason_def)
                antsarr.append([int(nbline[0]),int(nbline[1])])
                chans.append([1,channum])
                flag_count +=1
                if flag_count >= 100000000:
                    print 'Maximum flag entries per table reached, proceeding with new flag table'
                    #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                    #flag_count = 0
                for j in xrange(tag, i+1):
                    timecount = zeros[j]
                    if times_array[timecount][0] in timeflag_count:
                        timeflag_count[times_array[timecount][0]] += 1
                        if timeflag_count[times_array[timecount][0]] >= 600000:
                            print 'Maximum flags per time reached, proceeding with new flag table'
                            #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                            #flag_count = 0
        elif i == len(zeros)-1:
            #print '3'
            t1 = float(times_array[start][0])-0.00000116
            t2 = float(times_array[zeros[i]][0])+0.00000116
            src.append(source_num)
            subary.append(0)
            frqid.append(-1)
            ifs.append([startif,ifnum])
            pflagsarr.append([int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])])
            reasons.append(reason_def)
            antsarr.append([int(nbline[0]),int(nbline[1])])
            chans.append([1,channum])
            flag_count +=1
            if flag_count >= 100000000:
                print 'Maximum flag entries per table reached, proceeding with new flag table'
                #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                #flag_count = 0
            for j in xrange(tag, i+1):
                timecount = zeros[j]
                if times_array[timecount][0] in timeflag_count:
                    timeflag_count[times_array[timecount][0]] += 1
                    if timeflag_count[times_array[timecount][0]] >= 600000:
                        print 'Maximum flags per time reached, proceeding with new flag table'
                        #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                        #flag_count = 0
    #catch single entries
    if len(zeros) == 1:
        t1 = float(times_array[start][0])-0.00000116
        t2 = float(times_array[start][0])+0.00000116
        src.append(source_num)
        subary.append(0)
        frqid.append(-1)
        ifs.append([startif,ifnum])
        pflagsarr.append([int(pflags[0]),int(pflags[1]),int(pflags[2]),int(pflags[3])])
        reasons.append(reason_def)
        antsarr.append([int(nbline[0]),int(nbline[1])])
        chans.append([1,channum])
        flag_count +=1
        if flag_count >= 100000000:
            print 'Maximum flag entries per table reached, proceeding with new flag table'
            #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
            #flag_count = 0
        timecount = zeros[0]
        if times_array[timecount][0] in timeflag_count:
            timeflag_count[times_array[timecount][0]] += 1
            if timeflag_count[times_array[timecount][0]] >= 600000:
                print 'Maximum flags per time reached, proceeding with new flag table'
                #timeflag_count = dict.fromkeys(timeflag_count, 0) # set all timecounts back to 0
                #flag_count = 0
        
    
    return timeflag_count, flag_count,src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons


       

     
################ End of new direct FG table writing functions ##############
############################################################################





#######################################################################
# Read in previous flags multiprocess function
#######################################################################


def read_in_flags(send_q, rec_q, cpu):

    '''
    Parallel process to read in from an existing AIPS FG table and output appropriate flags 
    in a numpy flag array.
    Information about data and flag tables to select are pulled from the queue. The output
    flag array along with source and baseline information are sent to the output queue.

    '''

    for value in iter(send_q.get, 'STOP'):

        time0 = time.time()
        bline,nif,source,nchan,npol,ntimescans,times_array,orderedpols,uvname,uvklass,uvdisk,uvseq,old_flags,write_old_flags,flag_tables,mem_low,path2folder,srcname,job,totaljobs = value
        # print 'jobnum is',job
        # print fg,bline,nif,source,nchan,npol,ntimescans,orderedpols,uvname,uvklass,uvdisk,uvseq
        uvdata = AIPSUVData(uvname,uvklass,uvdisk,uvseq)
        #fgtab = uvdata.table('FG',fg)
        times_array=np.array(times_array)
        # print times_array.nbytes
        old_flags = copy.copy(old_flags)
        JOB = 'JOB %i of %i: ' % (job, totaljobs)

        if write_old_flags == 'yes':
            out_val = -9
            print JOB+'Including old flag information in new tables'
        elif write_old_flags == 'no':
            out_val = -9
            print JOB+'Not writing old flag information in new tables'
            
        
        if not type(old_flags).__module__ == np.__name__:
            # if not old_flags:
            old_flags = np.zeros([nif,npol,ntimescans,nchan],dtype=np.int8)
            #old_flags = np.zeros([nif,npol,ntimescans,nchan],dtype=np.float)
            # print old_flags.nbytes, nif, npol, ntimescans, nchan
            print JOB+'Creating empty flag array'
        #else:
            #print 'Using existing flag array'
        
        #old_flags2 = copy.copy(old_flags)
        
        for fg in flag_tables:
            #print fg
            # fgtab = uvdata.table('FG',fg)
            #check if flag table exists#
            #try:
            fgtab = uvdata.table('FG',fg)
            #except IOError:
             #   print '\nERROR: Flag table %i selected for loading, but it does not seem to exist' % fg
              #  print "ERROR: Either remove this from the flag_tables list or set do_loadflag = 'no'\n"
               # sys.exit(0)

            print '\n'+JOB+'Loading in previous flag table %i \n' % fg
            
            for row in fgtab:
                #print sys.getsizeof(row)
                if row.source == source or row.source == 0:
                    if (row.ants[0] == int(bline[0]) or row.ants[0] == 0) and (row.ants[1] == int(bline[-1]) or row.ants[1] == 0):
                        # print 'row.pflags is', row.pflags
                        tempol = np.array(row.pflags)
                        # print 'these are sorted', np.array(np.sort(orderedpols.values()))
                        tempol = tempol[np.array(np.sort(orderedpols.values()))]
                        # print tempol
                        polarr = np.where(tempol>0)
                        # print polarr
                        schan = row.chans[0]-1
                        if row.chans[1] == 0:
                            echan = nchan
                        else:
                            echan = row.chans[1]
                        sif = row.ifs[0]-1
                        if row.ifs[1] == 0:
                            eif = nif
                        else:
                            eif = row.ifs[1]
                        stime = row.time_range[0]
                        etime = row.time_range[1]
                        # print times_array
                        #print len(times_array), stime, etime
                        temptime = np.where((times_array>=stime) & (times_array<=etime))
                        #print polarr, temptime
                        #print npol, old_flags.shape, len(times_array)
                        #old_flags[sif:eif,polarr[0],temptime[0],schan:echan] = 10
                        for pol in polarr[0]:
                            # print 'pol is', pol
                            old_flags[sif:eif,pol,temptime,schan:echan] = out_val
                    

        time2 = time.time()
        print JOB+'Creating flags on CPU', cpu, 'took:', time2hms(time2-time0)
        # print mem_low
        if mem_low == True:
            t1=time.time()
            fname = path2folder+str(name) + "__" + str(srcname) +"__"+ str(bline)+"--oldflags"
            np.save(fname,old_flags)
            del old_flags
            rec_q.put((source, bline, nif, npol, cpu))
            # print time2hms(time.time()-t1), source, bline, nif
        # else:
          #  rec_q.put((source, bline, old_flags, nif, npol, cpu))

        del source 
        del bline
        #del old_flags
        del nif
        del orderedpols
        del times_array
        #print 'deleted old_flags'

############ End of read_in_flags multiprocess function #############
#######################################################################








###############################################################################
# Global Functions
################################################################################



## Task to find the source names and numbers if an SU table doesn't exist.
def indxr(uvdata):
    '''
    Performs AIPS task INDXR on the selected data
    '''
    indxr = AIPSTask('INDXR')
    indxr.inname = uvdata.name
    indxr.inclass = uvdata.klass
    indxr.inseq = uvdata.seq
    indxr.indisk = uvdata.disk
    #indxr.cparm = # sets length of scan, begin and end times of scan
    #indxr.bparm = # opacity and gain curve control
    indxr.inp()
    indxr.go()



## Definition to find computer memory stats (works only for csh and bash 
## at the moment)...
def computer_memory(path2folder):
    '''
    Will find the available computer memory information where possible
    '''
    import os
    memfile = open(path2folder + 'mem_stats.txt', 'wr')
    memfile.flush()
    file = path2folder + 'mem_stats.txt'
    os.system('free -k >' + file)    # need to test this!
    #os.system('free -k > mem_stats.txt')    # need to put in file path
    memfile = open(path2folder + 'mem_stats.txt', 'r')
    stats = []
    for line in memfile:
        string = ''
        for char in xrange(len(line)):
            if line[char].isdigit():
                string += line[char]
                if line[char+1] == ' ' or line[char+1] == '\n':
                    stats.append(int(string))
                    string = ''
    global mem_stats
    mem_stats = {}
    mem_stats['mem_total'] = stats[0]
    mem_stats['mem_used'] = stats[1]
    mem_stats['mem_free'] = stats[2]
    mem_stats['mem_shared'] = stats[3]
    mem_stats['mem_buffers'] = stats[4]
    mem_stats['mem_cached'] = stats[5]
    mem_stats['buf/cach_used'] = stats[6]
    mem_stats['buf/cach_free'] = stats[7]
    mem_stats['swap_total'] = stats[8]
    mem_stats['swap_used'] = stats[9]
    mem_stats['swap_free'] = stats[10]
    memfile.close()
    return mem_stats




## Frequency band calculator. Only works for freq between 1.2GHz - 26.5GHz 
## but easy to add on more...
def freq_band(freq):
    '''
    Will determine eMERLIN frequency band for information.
    '''
    if freq > 1.2E+9 and freq < 2.0E+9:
        return "L Band"
    elif freq > 2.0E+9 and freq < 4.0E+9:
        return "S Band"
    elif freq > 4.0E+9 and freq < 8.0E+9:
        return "C Band"
    elif freq > 8.0E+9 and freq < 12.0E+9:
        return "X Band"
    elif freq > 12.0E+9 and freq < 18.0E+9:
        return "Ku Band"
    elif freq > 18.0E+9 and freq < 26.5E+9:
        return "K Band"
    else:
        return "Out of range for calculator"


## Frequency prefix calculator
def freq_prefix(num):
    '''
    Will determine an appropriate formatting 
    for the data frequency.
    '''
    if num/(10**12) > 1:
        return "%.5f THz" % (num / 10.0**12)
    elif num/(10**9) > 1:
        return "%.5f GHz" % (num / 10.0**9)
    elif num/(10**6) > 1:
        return "%.5f MHz" % (num / 10.0**6)
    elif num/(10**3) > 1:
        return "%.5f KHz" % (num / 10.0**3)





## Determines the array size
def array_size(data_mem):
    '''
    Determines the memory required to store the specified
    numpy array
    '''
    ## Returns the size of the array with the appropriate suffix
    if len(str(data_mem)) <= 6:
        return "%.3f KB" % (data_mem / 10.0**3)
    elif len(str(data_mem)) <= 9:
        return "%.3f MB" % (data_mem / 10.0**6)
    elif len(str(data_mem)) <= 12:
        return "%.3f GB" % (data_mem / 10.0**9)
    elif len(str(data_mem)) <= 15:
        return "%.3f TB" % (data_mem / 10.0**12)


## Definition to give REASON string to find the date and give REASON in the FG table in 
## the format like IBLED does i.e. name of program + date + time...
def reason():
    '''
    Will give REASON string for use when writing entried in an AIPS flag table
    '''
    reason_str = 'SERPent ' + strftime("%d-%b-%y") + ' ' + strftime("%H:%M", localtime())
    return reason_str.upper()





## Function to convert time to d/hh:mm:ss.s format which is what is required for FG tables
def time2dhms(seconds):
    '''
    Function to convert time to d/hh:mm:ss.s format
    '''
    conv = seconds * 86400
    d = int(conv/ 86400)
    h=int(conv % 86400)/ 3600
    m=int(conv % 3600)/60
    s=conv-(d*86400)-(h*3600)-(m*60)
    d=`d`
    h=`h`
    m=`m`
    s="%3.1f" % s
    dhms= d.zfill(1) + "/" + h.zfill(2) + ":" + m.zfill(2) + ":" + s.zfill(4)
    return dhms



## Function to convert seconds to hh:mm:ss.ss format, returns a string
def time2hms(seconds):
    '''
    Function to convert seconds to hh:mm:ss.ss format, returns a string
    '''
    h=int(seconds/3600)
    m=int(seconds % 3600)/60
    s=seconds-(h*3600)-(m*60)
    h=`h`
    m=`m`
    s="%4.2f" % s
    hms=h.zfill(2)+":"+m.zfill(2)+":"+s.zfill(4)
    return hms







def mad(array):
    '''
    Calculates the median absolute deviation of the input array
    '''
    ab_dev = np.zeros_like(array)
    coef = 1.4286
    # median_i = float(np.median(array, axis=None))
    median_i = float(np.median(array[~np.isnan(array)], axis=None))
    ab_dev[:,:] = abs(array[:,:] - median_i)
    # median_j = float(np.median(ab_dev, axis=None))
    median_j = float(np.median(ab_dev[~np.isnan(array)], axis=None))
    final_mad = coef*median_j
    # print "The Median Absolute Deviation for this baseline is %f" % final_mad
    return final_mad






def mad_1d(array):
    '''
    Calculates the median absolute deviation of the input array in one dimension
    '''
    ab_dev = np.zeros_like(array)
    coef = 1.4286
    median_i = float(np.ma.median(array[~np.isnan(array)]))
    ab_dev[:] = abs(array[:] - median_i)
    median_j = float(np.ma.median(ab_dev[~np.isnan(ab_dev)]))
    final_mad = coef*median_j
    return final_mad




def alt_times(time):
    '''
    Function to print time information in human-readable format.
    '''

    time=float(time)
    rhrs = (time*24.0)
    hrs = math.floor(rhrs)
    rmins = (rhrs-hrs)*60.0
    mins = math.floor(rmins)
    rsecs = (rmins-mins)*60.0
    secs = math.floor(rsecs)
    # print hrs,mins,rsecs
    # print '%.2f:%.2f:%.2f' %(hrs,mins,rsecs)





## On the new version of AIPS a verb called REFLG is available which concatenates 
## rows in the FG table more efficiently. This will be very useful when flagging 
## L band data as the RFI is no longer '1 dimensional' but 2D, but it also has a 
## function where you can flag the entire time scan or channel if x% of the time 
## scan or channel is flagged... Very useful as some RFI in a channel (e.g.) might 
## be too weak to be picked up, but most of the rest of the channel is flagged... 
## Makes the flagger more robust and efficient with its rows in the FG table...

def reflg(cparm_1, cparm_2, cparm_3, cparm_4, cparm_5, cparm_6, cparm_7,path2folder,experiment,experiment_klass,experiment_seq,experiment_disk,flagging_options):

    '''
    Function to perform AIPS task REFLG on input data and flag table.
    '''
    
    reflg = AIPSTask('REFLG')
    ## changed from uvdata.name etc...
    reflg.inname = experiment
    reflg.inclass = experiment_klass
    reflg.inseq = experiment_seq
    reflg.indisk = experiment_disk
    reflg.flagver = 0
    #reflg.cparm[1:] = [0, 2, 0.70, 0.99, 0.70, 0.70, 2]
    if flagging_options == 'default':
        reflg.cparm[1:] = [0, 2, 0.70, 0.99, 0.70, 0.70, 2]
    else:
        reflg.cparm[1:] = [cparm_1, cparm_2, cparm_3, cparm_4, cparm_5, cparm_6, cparm_7]
    ## These control the flagging conditions:
    ## 1: interval
    ## 2: flags channels if between cparm[2] flagged channels
    ## 3: fraction of channels flagged at specific time - flag all channels
    ## 4: fraction of times flagged at specific channel - flag all times
    ## 5: fraction of baselines flagged at specific channel and time - flag all bline
    ## 6: fraction of baselines flagged at specific antenna - flag all bline with ant
    ## 7: cparm[7] = 1, flag cross hand pols, cparm[7] = 2, flag all pols
    reflg.inp()
    reflg.go()




## Definition to write text file to FG table.
def tbin(path2folder,experiment,experiment_klass,experiment_seq,experiment_disk):
    
    '''
    Performs AIPS task TBIN to append a new table to the selected data.
    '''
    tbin = AIPSTask('TBIN')
    tbin.outname = experiment
    tbin.outclass = experiment_klass
    tbin.outseq = experiment_seq
    tbin.outdisk = experiment_disk
    tbin.intext = path2folder + experiment + '_combined.fg'     # text file name
    tbin.bcount = 0
    tbin.ecount = 0
    tbin.inp()
    tbin.go()







## Write out the REFLG'd FG table to a text file...
def tbout(path2folder,experiment,experiment_klass,experiment_seq,experiment_disk):
    '''
    Performs AIPS task TBOUT to write a copy of the selected table to harddisk.
    '''
    tbout = AIPSTask('TBOUT')
    tbout.inname = experiment
    tbout.inclass = experiment_klass
    tbout.inseq = experiment_seq
    tbout.indisk = experiment_disk
    tbout.inext = 'FG'
    tbout.invers = 0
    tbout.bcount = 0
    tbout.ecount = 0
    tbout.outtext = path2folder + experiment + "_r.fg"    # reflg'd file name
    tbout.inp()
    tbout.go()








def coinc_chan_flags_alt(flagnew,coincs):

    '''
    Function to add flags in input flag array for instances where 'coincs' indices are encompassed
    by True flag entries.
    Outputs an updated flag array.
    '''
    import copy
    ncoincs = coincs+1
    # print 'Flagging %i coincident channels' % coincs
    # new_flag = copy.copy(flagnew)
    
    for col in xrange(flagnew.shape[1]-1):
        if flagnew[np.where((flagnew[:,col] > 0) & (flagnew[:,col+1] <= 0))].size > 0:
            short_arr = flagnew[np.where((flagnew[:,col] > 0) & (flagnew[:,col+1] <= 0)),col+1:col+ncoincs+1]
            if np.any(np.where(short_arr>0)):
                test_arr = np.where(short_arr>0,flagnew[np.where((flagnew[:,col] > 0) & (flagnew[:,col+1] <= 0)),col+1:col+ncoincs+1],20)
                farr = np.argmax(np.any(short_arr>0,axis=0),axis=1)
                for k in xrange(len(farr)):
                    short_arr[0,k,0:farr[k]] = 50
                flagnew[np.where((flagnew[:,col] > 0) & (flagnew[:,col+1] <= 0)),col+1:col+ncoincs+1]=short_arr
    return(flagnew)



def coinc_flags(flagnew,time_coincs,chan_coincs):
    
    '''
    Function to add flags in input flag array for instances where 'coincs' indices are encompassed
    by True flag entries. Performs the process iteratively over each dimension, updating across channels
    first.
    Outputs an updated flag array.
    '''

    import copy
    # print 'Flagging %i coincident channels and %i coincident times' % (chan_coincs,  time_coincs)
    new_flag = copy.copy(flagnew)
    if time_coincs > 0:
        new_flag = coinc_chan_flags_alt(new_flag,chan_coincs)
    if chan_coincs > 0:
        new_flag = new_flag.transpose()
        new_flag = coinc_chan_flags_alt(new_flag,time_coincs)
        new_flag = new_flag.transpose()
    return(new_flag)

def coinc_flags_noncumu(flagnew,time_coincs,chan_coincs):
    
    '''
    Function to add flags in input flag array for instances where 'coincs' indices are encompassed
    by True flag entries. Performs the process iteratively over each dimension but in a non-cumulative
    fashion, updating across channels first.
    Outputs an updated flag array.
    '''


    import copy
    # print 'Flagging %i coincident channels and %i coincident times' % (chan_coincs,  time_coincs)
    new_flag = copy.copy(flagnew)
    if time_coincs > 0:
        new_flag = coinc_chan_flags_alt(new_flag,chan_coincs)
    if chan_coincs > 0:
        new_flag2 = copy.copy(flagnew.transpose())
        new_flag2 = coinc_chan_flags_alt(new_flag2,time_coincs)
        new_flag2 = new_flag2.transpose()
        new_flag[np.where(new_flag2==50)]=50
    return(new_flag)





# Code to flatten lists into single dimension. Code from #
# http://rightfootin.blogspot.co.uk/2006/09/more-on-python-flatten.html #
def flatten(l, ltypes=(list, tuple)):
    '''
    Function to flatten lists into single dimension. Code from 
    http://rightfootin.blogspot.co.uk/2006/09/more-on-python-flatten.html
    '''

    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)





def srce_num_lookup(uvdata,tnum):

    '''
    Function to get the source information for the selected AIPS data file
    using the selected AIPS SU table.
    '''
    source_list = {}
    source_lookup = []
    sutab = uvdata.table('SU',tnum)
    nsource = len(uvdata.sources)
    print " Total number of sources (nsource): %i" % nsource
    #print uvdata.sources
    for s in uvdata.sources:
        for row in sutab:
            #print s,row.source.strip(' ')
            if s == row.source.strip(' '):
                #print 'yes',s
                #print row.id__no
                source_list[row.id__no] = s
                source_lookup.append(s)
    # print source_list, source_lookup
    return source_list, source_lookup



def get_ordered_pols(uvdata):
    '''
    Function to determine the AIPS flag entry appropriate index for 
    the present stokes parameters.
    '''
    ## Figure out number of Stokes and polarizations
    npol = len(uvdata.stokes)
    # print " Number of Polarizations (npol): %i" % npol
    orderedpols = {}
    for x in xrange(npol):
        if uvdata.stokes[x] == 'RR' or uvdata.stokes[x] == 'XX':
            orderedpols[str(uvdata.stokes[x])] = 0
        elif uvdata.stokes[x] == 'LL' or uvdata.stokes[x] == 'YY':
            orderedpols[str(uvdata.stokes[x])] = 1
        elif uvdata.stokes[x] == 'RL' or uvdata.stokes[x] == 'XY':
            orderedpols[str(uvdata.stokes[x])] = 2
        elif uvdata.stokes[x] == 'LR' or uvdata.stokes[x] == 'YX':
            orderedpols[str(uvdata.stokes[x])] = 3
    return orderedpols




###############################################################################
# End of Global Functions
################################################################################









###############################################################################
###############################################################################
## Main SERPent Script 
###############################################################################
###############################################################################

version_name = "Aesculapian"
version_number = '1.2'
version_date = "12/12/16"

#print sys.argv
if len(sys.argv) > 1:
    if sys.argv[1] == '--v' or sys.argv[1] == '--version':
        print 'SERPent version is %s: %s dated %s' % (version_name, version_number, version_date)
        sys.exit(0)
    else:
        print 'ERROR: unrecognised option'
        sys.exit(0)


print '\n Started running SERPent version %s: %s' % (version_name, version_number), 'on %s' % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()),'\n'

#print '\n Started running SERPent version: %s' % version, 'on %s' % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()),'\n'


## Execute the SERPent_input.py file with all the observation information 
## and flagging parameters...
try:
    execfile("SERPent_input.py")
    
except:
    print "\n Could not find SERPent_input.py!"
    print " Make sure input file is in the same folder as this script (SERPent.py)"
    print " Aborting!"
    sys.exit(0)



## Test to see whether all variables have been set:
try:
    AIPS_user_number
    Name
    Klass
    Disk
    Seq
    NCPU
    path2folder
    path2folder = os.path.join(path2folder, '')
    if os.path.exists(path2folder) == False:
        print " Folder for outputs does not exist! Please check inputs and try again."
        print " Aborting SERPent!\n"
        sys.exit(0)
    do_Lovell_scan
    phasecal
    if do_Lovell_scan == 'no':
        phasecal = []
    do_lovell_cross
    zero_level
    which_baselines
    flagging_options
    flag_sources
    if flag_sources == 'choose':
        flag_list
    if flagging_options == 'choose' or flagging_options == 'source':
        aggressiveness_first_run
        max_subset_first_run
        aggressiveness_second_run
        max_subset_second_run
        rho
        kickout_sigma_level
    if which_baselines == 'choose':
        baselines
    dobline_specific
    dosource_specific
    select_stokes
    if select_stokes == 'yes':
        flag_all_stokes 
        stokes
    coadd_polarization_flags
    if coadd_polarization_flags == 'no':
        coadd_zero_flags
        coadd_lovell_flags
    flag_coinc_chans
    flag_coinc_times
    zero_flag_allif
    use_fits
    if do_loadflag:
        flag_tables
        write_old_flags
except NameError:
    print " Please check you\'ve set all the variables in the input file."
    print " ABORTING SERPent!\n"
    sys.exit(0)





### Check to see if pyfits is available for use when writing flag entries
### If not will default to parseltongue table entry in AIPS
### If it is, flags will be written out to dummy fits files, concatenated 
### and loaded into AIPS using either pyfits or astropy

if use_fits == 'yes':
    try:
        import pyfits as fits
        use_fits = 'yes'
        print 'Using pyfits to write flag table fits files'
    except ImportError:
        try:
            import astropy.io.fits as fits
            use_fits = 'yes'
            print 'Using astropy to write flag table fits files'
        except ImportError:
            use_fits = 'no'
            print 'Using parseltongue'

else:
    print 'Using parseltongue to write flag tables in AIPS'





## Read parameters from input file to local variable
# parameters = [aggressiveness_first_run, max_subset_first_run, aggressiveness_second_run, max_subset_second_run, rho, kickout_sigma_level]
# print parameters


AIPS.userno = AIPS_user_number
experiment = Name
experiment_klass = Klass
experiment_disk = Disk
experiment_seq = Seq
# print AIPS.userno
# print experiment, experiment_klass, experiment_disk, experiment_seq
print Disk
cat = AIPSCat(Disk)
uvdata = AIPSUVData(Name, Klass, Disk, Seq)
no_rows = 0
pickle_list = {}




## Destroy the fg table for trial runs (and mem_stats file)...
if os.path.exists(path2folder+experiment+'.fg') == True:
    os.remove(path2folder+experiment+'.fg')
if os.path.exists(path2folder+experiment+'_r.fg') == True:
    os.remove(path2folder+experiment+'_r.fg')
if os.path.exists(path2folder+'dummy.fg') == True:
    os.remove(path2folder+'dummy.fg')
if os.path.exists(path2folder+'mem_stats.txt') == True:
    os.remove(path2folder+'mem_stats.txt')
if os.path.exists(path2folder+experiment+'_lovell_dropouts.fg') == True:
    os.remove(path2folder+experiment+'_lovell_dropouts.fg')
if os.path.exists(path2folder+experiment+'_lovell_dummy.fg') == True:
    os.remove(path2folder+experiment+'_lovell_dummy.fg')
if os.path.exists(path2folder+experiment+'_combined.fg') == True:
    os.remove(path2folder+experiment+'_combined.fg')
if os.path.exists(path2folder+experiment+'_zeros_dropouts.fg') == True:
    os.remove(path2folder+experiment+'_zeros_dropouts.fg')
if os.path.exists(path2folder+experiment+'_zeros_dummy.fg') == True:
    os.remove(path2folder+experiment+'_zeros_dummy.fg')
if os.path.exists(path2folder+'blank_table.fg') == True:
    os.remove(path2folder+'blank_table.fg')



''' Don't need this any longer, now counts tables'''
## Destroy AIPS FG extensions for trial runs...
#uvdata.zap_table('FG', -1)

print flag_list

## new code to find proper source numbers ##
## First try and see whether there are any SU tables, and get source numbers and 
## names from that.
try:
    source_list, source_lookup = srce_num_lookup(uvdata,1)
    print source_list,'1'
    singlesource = False
    srces2remove=[]
    if flag_sources == 'choose':
        # print 'flag_list1:', flag_list
        for srce in source_list:
            # print srce
            sdi = source_lookup.index(source_list[srce])
            print sdi, source_list[srce]
            if source_list[srce] not in flag_list:
                print 'nope'
                srces2remove.append(srce)
                #source_lookup.remove(source_list[srce])
                #del source_list[srce]
            else:
                continue
                print 'yep'
            # srces2remove=[]
        for srce in srces2remove:
            # print srce
            source_lookup.remove(source_list[srce])
            del source_list[srce]
        source_lookup = flag_list
    # print 'sources are', source_lookup
    # print source_list, source_lookup
    # check against source list from inputs #
    print " Names of sources selected for flagging (with dictionary index):"
    nsource = len(source_lookup)
    for s in source_list:
        print " %i: %s" % (s, source_list[s])
except:
    print "\n No SU table found... "
    print " Assuming single source..."
    singlesource = True
    source_lookup = []
    source_list = {}
    # nsource = 1
    if flag_sources == 'choose':
        source_list[1] = flag_list[0]
        srce = flag_list[0]
        nsource = len(source_lookup)
    else:
        print "\n ERROR: Single-source file but no sources specified in flag_list\n Please check flag_sources='choose' and the correct source is entered.\n"
        sys.exit(0)
    # print 'source_list, source_lookup:', source_list, source_lookup
    # source_list[1] = Name.upper()
    source_lookup.append(flag_list[0])
    print " %i: %s" % (1, source_list[1])



## If basic flagging options used, load in parameters ##

source_base_parms = []
if flagging_options == 'choose':
    if len(aggressiveness_first_run) > 1 or len(max_subset_first_run) > 1 or len(aggressiveness_second_run) > 1 or len(max_subset_second_run) > 1 or len(rho) > 1 or len(kickout_sigma_level) > 1:
        print "\n ERROR: flagging_options = 'choose' set but multiple instead of single flagging \n ERROR: parameters given, please check inputs. \n ERROR: If separate options required for multiple sources, \n ERROR: set flagging_options = 'source'\n"
        sys.exit(0)
    else:
        source_parms = [aggressiveness_first_run[0], max_subset_first_run[0], aggressiveness_second_run[0], max_subset_second_run[0], rho[0], kickout_sigma_level[0]]
        for i in xrange(len(source_lookup)):
            source_base_parms.append(source_parms)
        
elif flagging_options == 'source':
    noptions = len(source_lookup)
    if len(aggressiveness_first_run) != noptions or len(max_subset_first_run) != noptions or len(aggressiveness_second_run) != noptions or len(max_subset_second_run) != noptions or len(rho) != noptions or len(kickout_sigma_level) != noptions:
        print "\n ERROR: flagging_options = 'source' set but number of parameters given does not \n ERROR: match the number of sources specified for flagging. \n ERROR: If flag_sources = 'all', the number of parameters should match \n ERROR: the number of sources in the data. \n ERROR: If flag_sources = 'choose', the number of parameters should match \n ERROR: the number of sources given in flag_list. \n"
        sys.exit(0)
    else:
        source_base_parms = []
        for i in xrange(noptions):
            newparm = [aggressiveness_first_run[i],max_subset_first_run[i],aggressiveness_second_run[i],max_subset_second_run[i],rho[i],kickout_sigma_level[i]]
            source_base_parms.append(newparm)


## If using source and/or baseline specific options, load in the parameters ###

if dobline_specific == 'no':
    src_baselines = [[] for i in range(len(source_list))]
    for srcething in source_list:
        sdi = source_lookup.index(source_list[srcething])
        source_num = srcething
        src_baselines[sdi] = baselines
        #    for i in xrange(len(source_lookup)):
        #       src_baselines.append(baselines)
    # print src_baselines


elif dobline_specific == 'yes' and dosource_specific == 'no':
    print source_list
    source_base_parms = [[] for i in range(len(source_list))]
    src_baselines = [[] for i in range(len(source_list))]
    for srcething in source_list:
        sdi = source_lookup.index(source_list[srcething])
        lbaselines = []
        full_parms_dict = {}
        for i in xrange(1000):
            try:
                lbaselines.extend(globals()['baselines_'+str(i)]) 
                for base in globals()['baselines_'+str(i)]:
                    full_parms_dict[base] = globals()['parameters_'+str(i)]
            except: continue
        if full_parms_dict:
            source_base_parms[sdi]=full_parms_dict
            src_baselines[sdi]=lbaselines
    # print src_baselines, source_base_parms



elif dobline_specific == 'yes' and dosource_specific == 'yes':
    source_base_parms = [[] for i in range(len(source_list))]
    src_baselines = [[] for i in range(len(source_list))]
    for srcething in source_list:
        sdi = source_lookup.index(source_list[srcething])
        base_parms_dict = {}
        lbaselines=[]
        # print sdi
        for i in xrange(1000):
            try:
                lbaselines.extend(globals()['source_'+str(sdi)+'_baselines_'+str(i)])
                for base in globals()['source_'+str(sdi)+'_baselines_'+str(i)]:
                    base_parms_dict[base] = globals()['source_'+str(sdi)+'_parameters_'+str(i)] 
                # print base_parms_dict
            except:
                continue
        if base_parms_dict:
            source_base_parms[sdi]=base_parms_dict
            src_baselines[sdi]=lbaselines




## Information about the observations:
name = uvdata.name
print "\n RUNNING SCRIPT ON DATA: %s" % name

## Name of Telescope:
telescope = uvdata.header['telescop']
print " Name of Telescope: %s" % telescope

# sets up local variable
nvis = len(uvdata)    ## prints the number of visibilities
print " Number of visibilities (nvis): %i" % nvis

## Find the Frequency of the observation
frequency = uvdata.header['crval'][uvdata.header['ctype'].index('FREQ')]
print " Frequency of the Observation: %s" % freq_prefix(frequency)
print " Frequency band of the Observation: %s" % freq_band(frequency)

## Figure out number of Stokes and polarizations
npol = len(uvdata.stokes)
print " Number of Polarizations (npol): %i" % npol
polnames = {}
for x in xrange(npol):
    polnames[str(uvdata.stokes[x])] = x
## Below is just to print out the stokes in a nice way :)
string = ""
for x in xrange(npol):
    if x != npol-1:
        string += (str(uvdata.stokes[x])+ ", ")
    else:
        string += str(uvdata.stokes[x])
print " Polarizations: %s" % string


# Get the stokes in AIPS array order #
orderedpols = get_ordered_pols(uvdata)

## Figure out number of IFs
if 'IF' in uvdata.header['ctype']:
    nif = uvdata.header['naxis'][uvdata.header['ctype'].index('IF')]
    print " Number of IFs (nif): %i" % nif
else:
    print " Keyword IF not found in header, assuming number of IFs is 1"
    nif = 1


## Figure out number of channels
nchan = uvdata.header['naxis'][uvdata.header['ctype'].index('FREQ')]
print " Number of Channels (nchan): %i" % nchan


## Figure out number of baselines
nant = len(uvdata.antennas)
nbase = ((nant*(nant-1))/2)
nacbase = nbase + nant
print " Number of Possible Baselines (nbase): %i" % nacbase



obs_date_str = uvdata.header.date_obs
obs_RA = uvdata.header['crval'][4]
obs_DEC = uvdata.header['crval'][5]
freq = frequency
dfrq = uvdata.header['cdelt'][uvdata.header['ctype'].index('FREQ')]
pfrq = uvdata.header['crpix'][uvdata.header['ctype'].index('FREQ')]
yr,mn,dy = obs_date_str.split('-')
obsJD = datetime.date(int(yr),int(mn),int(dy)).toordinal() + 1721424.5  # time = 00:00:00

metadata=[obsJD,obs_date_str,obs_RA,obs_DEC,freq,dfrq,pfrq]

## If selected stokes requested, remove extra from polnames
if select_stokes == 'yes':
    if not stokes:
        print ' No stokes specified, processing all'
    else:
        for pol in polnames.keys():
            if pol not in stokes:
                # print pol
                # print 'no'
                del polnames[pol]

# print 'polnames', polnames



## List for baselines to remove from SERPent job list and flag at end of script
baseline_removal = {}

## Removes Lovell-Mk2 baselines as this baseline is pretty much useless.
antable = uvdata.table('AN', 0)
lovellmk2 = 0
lovellmk2flag = 0
# print baselines
if telescope in ('e-MERLIN', 'E-MERLIN', 'eMERLIN', 'EMERLIN'):
    for ant in xrange(len(uvdata.antennas)):
        if uvdata.antennas[ant].upper() in ('LO', 'LOV', 'LOVE', 'LOVEL', 'LOVELL'):
            lo = antable[ant]['nosta']
        if uvdata.antennas[ant].upper() in ('MK2', 'M2', 'MARK2', 'MK'):
            mk = antable[ant]['nosta']
            #print 'MK is', mk
    try:
        lo
        mk
    except:
        lo = 0
        mk = 0
    # src_baselines = [[] for i in range(len(source_list))]
    for srcething in source_list:
        sdi = source_lookup.index(source_list[srcething])
        for bline in src_baselines[sdi]:
            if str(lo) in bline and str(mk) in bline:
                print "\n Removing Lovell-Mk2 baseline from baseline list."
                lovellmk2 = bline
                baseline_removal[bline] = 1
                src_baselines[sdi].remove(bline)
                lovellmk2flag = 1
            break


## Remove auto-correlations...
antable = uvdata.table('AN', 0)
for ant in xrange(len(uvdata.antennas)):
    ant_num = antable[ant]['nosta']
    auto_bline = str(ant_num) + '-' + str(ant_num)
    #print ant_num, auto_bline
    for srcething in source_list:
        sdi = source_lookup.index(source_list[srcething])
        for bline in src_baselines[sdi]:
            if auto_bline == bline:
                print " Removing Auto-correlation baseline %s from baseline list." % (auto_bline)
                baseline_removal[bline] = 1
                baselines.remove(bline)
            break

# print 'baseline_removal is:', baseline_removal


## Dictionary containing number of time scans for each baseline.
if dobline_specific == 'no' and which_baselines == 'all':
    
    print "\n ALL baselines have been selected for flagging"
    print "\n Writing visibilities to file..."
    
    ## List of baselines. Note that this uses the actual visibility data to find the 
    ## antenna numbers and doesn't assume they start at 1 and increase linearly by 
    ## the same amount. This has been written because the VLA and early commissioning 
    ## e-MERLIN data have a unique (i.e. stupid) numbering system.
    ## Note that this is only in the 'all' baseline choice as it then needs to find 
    ## all the baselines but doesn't if they've been stated before...

    # print 'flag_list is:', flag_list


    ## New 11/2015 code, assume have removed any sources from source_list and source_lookup
    ## that are not in the requested source list from the inputs file.

    sbline = [[] for i in range(len(source_list))]
    nbases = [0 for i in range(len(source_list))]
    baseline_dict = [{} for i in range(len(source_list))]
    dtimescans = [[] for i in range(len(source_list))]
    dvisnums = [[] for i in range(len(source_list))]
    viscount = 0
    checkcount=0
    n=0

    vistime1 = time.time()

    for vis in uvdata:
        #print vis.visibility.shape
        if singlesource:
            srce = 1
        else: 
            srce = vis.source
        if srce in source_list:
            sdi = source_lookup.index(source_list[srce])
            bline = "%i-%i" % (vis.baseline[0], vis.baseline[1])
            if bline not in sbline[sdi]:
                sbline[sdi].append(bline)
                nbases[sdi] += 1
                baseline_dict[sdi][bline] = 1
                dtimescans[sdi].append([])
                dvisnums[sdi].append([])
                dtimescans[sdi][sbline[sdi].index(bline)].append(float(vis.time))
                dvisnums[sdi][sbline[sdi].index(bline)].append(int(viscount))
            else:
                baseline_dict[sdi][bline] += 1
                dtimescans[sdi][sbline[sdi].index(bline)].append(float(vis.time))
                dvisnums[sdi][sbline[sdi].index(bline)].append(int(viscount))

            for j in xrange(nif):

                # if os.path.exists(path2folder+ str(name)  +"__"+ str(source_list[srce]) + '__IF' + str(j+1) + '__' + str(bline)) == True:
                    # os.remove(path2folder+ str(name)  +"__"+ str(source_list[srce]) + '__IF' + str(j+1) + '__' + str(bline))
                    
                visoutfile = open(path2folder + str(name)  +"__"+ str(source_list[srce]) + '__IF' + str(j+1) + '__' + str(bline),"ab")
                #visarray = np.array(np.where(np.array(vis.visibility)[j,:,:,2] > 0.0, (np.sqrt(np.array(vis.visibility)[j,:,:,0]**2 + np.array(vis.visibility)[j,:,:,1]**2)), float('NaN')))
                visarray = np.array(np.where(np.asfarray(vis.visibility,dtype=np.float64)[j,:,:,2] > 0.0, (np.sqrt(np.asfarray(vis.visibility,dtype=np.float64)[j,:,:,0]**2 + np.asfarray(vis.visibility,dtype=np.float64)[j,:,:,1]**2)), float('NaN')))
                np.save(visoutfile,visarray)
                #print visarray
                visoutfile.close()
                '''
                if bline == '8-9' and j == 1:
                    checkcount += 1
                    # print 'appending amps', checkcount
                '''
        viscount += 1
	#if viscount/nvis > n*0.10*nvis:
	#	print n*10, '%'
	#	n +=1
    #print "List of ALL baselines with corresponding number of time scans:", ntimescans
    
    newvisti=time.time()
    print 'Finished loading visies in:', time2hms(newvisti-vistime1)

    


## Testing for specified baselines.
elif which_baselines == 'choose' or dobline_specific == 'yes':
    print "\n SPECIFIC baselines selected for flagging..."
    print "\n Writing visibilities to file..."

    sbline = [[] for i in range(len(source_list))]
    nbases = [0 for i in range(len(source_list))]
    baseline_dict = [{} for i in range(len(source_list))]
    dtimescans = [[] for i in range(len(source_list))]
    dvisnums = [[] for i in range(len(source_list))]
    viscount = 0
    n=0

    vistime1 = time.time()
    print src_baselines, 'src'

    for vis in uvdata:
        if singlesource:
            srce = 1
        else: 
            srce = vis.source
        if srce in source_list:
            sdi = source_lookup.index(source_list[srce])
            bline = "%i-%i" % (vis.baseline[0], vis.baseline[1])
            # if (bline not in sbline[sdi]) and (bline in src_baselines[sdi]):
            if (bline in src_baselines[sdi]):
                if (bline not in sbline[sdi]):
                    sbline[sdi].append(bline)
                    nbases[sdi] += 1
                    baseline_dict[sdi][bline] = 1
                    dtimescans[sdi].append([])
                    dvisnums[sdi].append([])
                    dtimescans[sdi][sbline[sdi].index(bline)].append(float(vis.time))
                    dvisnums[sdi][sbline[sdi].index(bline)].append(int(viscount))
                elif bline in src_baselines[sdi]:
                    baseline_dict[sdi][bline] += 1
                    dtimescans[sdi][sbline[sdi].index(bline)].append(float(vis.time))
                    dvisnums[sdi][sbline[sdi].index(bline)].append(int(viscount))
                for j in xrange(nif):
                    # if os.path.exists(path2folder+ str(name)  +"__"+ str(source_list[srce]) + '__IF' + str(j+1) + '__' + str(bline)) == True:
                       #  os.remove(path2folder+ str(name)  +"__"+ str(source_list[srce]) + '__IF' + str(j+1) + '__' + str(bline))
                    visoutfile = open(path2folder + str(name)  +"__"+ str(source_list[srce]) + '__IF' + str(j+1) + '__' + str(bline),"ab")
                    # visarray = np.array(np.where(np.array(vis.visibility)[j,:,:,2] > 0.0, (np.sqrt(np.array(vis.visibility)[j,:,:,0]**2 + np.array(vis.visibility)[j,:,:,1]**2)), float('NaN')))
                    visarray = np.array(np.where(np.asfarray(vis.visibility,dtype=np.float64)[j,:,:,2] > 0.0, (np.sqrt(np.asfarray(vis.visibility,dtype=np.float64)[j,:,:,0]**2 + np.asfarray(vis.visibility,dtype=np.float64)[j,:,:,1]**2)), float('NaN')))
                    #print visarray
                    np.save(visoutfile,visarray)  
                    visoutfile.close()
        viscount += 1
	#if viscount/nvis > n*0.10*nvis:
#		print n*10, '%'
#		n +=1
    newvisti=time.time()
    print 'Finished loading visies in:', time2hms(newvisti-vistime1)


print "List of ALL baselines with corresponding number of time scans:", baseline_dict #ntimescans

# print '\n order is sbline,baseline_dict,dtimescans,dvisnums,nbases'
# print sys.getsizeof(sbline),sys.getsizeof(baseline_dict), sys.getsizeof(dtimescans),sys.getsizeof(dvisnums),sys.getsizeof(nbases)
#print ts(sbline,verbose=True) 
# print ts(sbline),ts(baseline_dict),ts(dtimescans),ts(dvisnums),ts(nbases)
# print '\n'
# print sbline

# Check data exists for all sources #
nodata_list = []
for srce in source_list:
    sdi = source_lookup.index(source_list[srce])
    if not baseline_dict[sdi]:
        print "\n WARNING: No data found for source %s, it will not be flagged!\n" % source_list[srce]
        nodata_list.append(source_list[srce])
# print baseline_dict




'''

# We already have:
# 1. a dictionary that contains the baseline and number of timescans associated with it...(ntimescans)
# 2. a dictionary that contains the baseline and each timescan associated with it that we want to plot (visnumdict)
# 1. vis_out[label] = file pointer to pytables visibility file
# 2. vis_array[label] = pointer to numpy array to be saved within pytables file
# 3. times_out=[label] = file pointer to pytables time file
# 4. times_array[label] = pointer to numpy times array to be saved within pytables file
# 5. countdic[label] = to count each visibility of each page (i.e. baseline/time)
# 6. label_list = list to keep track of saved file names to return at end of visappend function
vis_array = collections.OrderedDict()
vis_out = collections.OrderedDict()
times_array_new = collections.OrderedDict()
times_out = collections.OrderedDict()
countdict = collections.OrderedDict()
label_list = []

if os.path.exists(path2folder+'/timesfile')
if os.path.exists(path2folder+'/visfile')

timesfile = tables.openFile(path2folder+'/timesfile', mode='w')
visfile = tables.openFile(path2folder+'/visfile', mode='w')


for srce in source_list:
    sdi = source_lookup.index(source_list[srce])
    for bline in src_baselines[sdi]:
        labeltime = '{}_{}'.format(sdi, bline)
        atom = tables.Float64Atom()
        timesarrayidentifier = 'timesarray'+str(timecount)
        timeidentitydic[labeltime] = 'times_file.root.'+str(timesarrayidentifier)
        timecount += 1
        times_array[labeltime] = timesfile.createEArray(timesfile.root, timesarrayidentifier, atom, (0,), labeltime)
        for st in stokes:
            for IF in xrange(IF_start-1,IF_end):
                label = '{}_{}_{}-{}_{}'.format(bline, st, tblock[0], tblock[1], IF+1)
                visarrayname = label+'.hdf5'
                label_list.append(visarrayname)
                visarrayname = label+'.hdf5'
                label_list.append(visarrayname)
							
                visarrayidentifier = 'visarray'+str(viscount)
                visidentitydic[label] = 'visfile.root.'+str(visarrayidentifier)
                viscount += 1
                vis_array[label] = visfile.createEArray(visfile.root, visarrayidentifier, atom, (0,nchan), label)

		countdict[bline] = countlist
	

'''











'''
## Checks for Linux-based OS then runs memory definition
try:
    mem_stats = computer_memory(path2folder)
except:
    print "\n Sorry, computer_memory() definition does not work on your system!"
    
    
if 'mem_total' in mem_stats:
    print "System Memory Information:"
    print "Total Memory  :    %s" % array_size(mem_stats['mem_total']*1000)
    print "Used Memory   :    %s" % array_size(mem_stats['mem_used']*1000)
    print "Free Memory   :    %s" % array_size(mem_stats['mem_free']*1000)
    print "Total Swap    :    %s" % array_size(mem_stats['swap_total']*1000)
    print "Used Swap     :    %s" % array_size(mem_stats['swap_used']*1000)
    print "Free Swap     :    %s" % array_size(mem_stats['swap_free']*1000)
'''

## The predicted array sizes, these predict the arrays when flagging all baselines 
## simultaneously
pred_array = 0
srce_times = [[] for i in range(len(source_list))]
for srce in source_list:
    timeflag_count = {}
    sdi = source_lookup.index(source_list[srce])
    pred_array += sum(baseline_dict[sdi].itervalues())*nif*nchan*npol*3*8

    ## use dtimescans to create dictionary of total timescans 
    ## per source to keep count when writing flags 
    print 'Creating timecount dictionary'
    # sdi = source_lookup.index(source_list[srce])
    temparr = flatten(dtimescans[sdi])
    for i in xrange(len(temparr)):
        if temparr[i] not in timeflag_count:
            timeflag_count[temparr[i]] = 0
        else: 
            timeflag_count[temparr[i]] = 0
            #timeflag_count[temparr[i]] += 1
    del temparr
    srce_times[sdi] = timeflag_count
print "Total predicted memory usage: %s" % array_size(int(pred_array))



# print stop

#########################################################################
# 
# Reading in previous flag tables. Will read into arrays for each source
# nif,npol,ntimes,nchan
# Is parallelised for each source
#
#########################################################################


if do_loadflag == 'yes':

    # print '\nReading in previous flag tables %s\n' % str(', '.join(flag_tables)) 
    # mem_low = True 

    #### check flag tables exist ####
    notabs = []
    for fg in flag_tables:
        # print fg
        # fgtab = uvdata.table('FG',fg)
            #check if flag table exists#
        try:
            fgtab = uvdata.table('FG',fg)
        except IOError:
            print '\nERROR: Flag table %i selected for loading, but it does not seem to exist' % fg
            print "ERROR: You should either remove this from the flag_tables list or set do_loadflag = 'no'"
            print "ERROR: Removing table %i from flag table list" % fg
            # sys.exit(0)
            notabs.append(fg)
    
    for tab in notabs:
        flag_tables.remove(tab)

    if not flag_tables:
        print '\nERROR: No flag tables selected or no flag tables found.'
        print 'ERROR: Aborting loading previous flag tables and continuing with flagging. \n'
        flagsexist = False

    else:  
        flagsexist = True
        load_results = [{} for i in range(len(source_list))]
        old_flags = []
        for srce in source_list:
            sdi = source_lookup.index(source_list[srce])
            for i in xrange(len(sbline[sdi])):
                base = sbline[sdi][i]
                load_results[sdi][base] = old_flags

        joblist3 = []
        total_jobs = 0
        for s in sbline:
            total_jobs += len(s)
        jobcount = 1

        for source in source_list:
            sdin = source_lookup.index(source_list[source])
            srcname = source_list[source]
            for base in sbline[sdin]:
                ntimescans = baseline_dict[sdin][base]
                times_array = dtimescans[sdin][sbline[sdin].index(base)]
                old_flags = load_results[sdin][base]
                joblist3.append([base,nif,source,nchan,npol,ntimescans,times_array,orderedpols,uvdata.name,uvdata.klass,uvdata.disk,uvdata.seq,old_flags, write_old_flags,flag_tables,mem_low,path2folder,srcname,jobcount,total_jobs])
                jobcount +=1


        sendfg_q = mp.Queue()
        recfg_q = mp.Queue()
        ncpus = mp.cpu_count()
        if NCPU == 0:
            NCPU = ncpus
        elif ncpus > NCPU and not NCPU == 0:
            ncpus = NCPU
        # print ncpus
    

        for i in xrange(len(joblist3)):
            print 'Sending to queue job', i+1, 'of', len(joblist3)
            # print "CPU", c, "will flag baselines and IFs: \n"
            sendfg_q.put(joblist3[i])

        if len(joblist3)<ncpus:
            jobnum = len(joblist3)
        else:
            jobnum = ncpus
    

        for i in xrange(jobnum):
            proc = mp.Process(target=read_in_flags, args=(sendfg_q,recfg_q,i))
            proc.start()
            print 'Starting process on CPU:', i



        for i in xrange(len(joblist3)):
            results = recfg_q.get()
            srce = results[0]
            base = results[1]
            sdin = source_lookup.index(source_list[srce])
            if not mem_low == True:
                load_results[sdin][base] = results[2]
                # print srce
                # print 'doing this'

        for i in xrange(ncpus):
            sendfg_q.put('STOP')
        

        newti=time.time()
        print 'Finished in:', time2hms(newti-ti)

        del joblist3


else:
    flagsexist = False
#print stop
#del sendfg_q
#del recfg_q
#thefile = open('testflags_m82only.txt','wb')



###############################################################################
#
#   Flagging algorithm split into baselines and IFs (for each baseline work 
#   through the IFs individually then move on to the next baseline)
#   This now uses the process/queue approach to multiprocessing from the 
#   multiprocessing python module.
#
###############################################################################

# Begin the parallelization block:

total_tscans = 0
#global pickle_list
pickle_list = {}
lovell_bline = []
global lovell_row
lovell_row = 1
zeros_bline = []
global zeros_row
zeros_row = 1

nbline = len(baselines)    ## need to change
total_jobs = 0
for s in sbline:
    total_jobs += len(s)
total_jobs = total_jobs * nif
#print total_jobs,'yes'
# total_jobs = nbline*nif*nsource
size_total = 0
# mem_low = True
jobcount = 1

### Create array of jobs to do ####
joblist = []
for srcething in source_list:
    sdi = source_lookup.index(source_list[srcething])
    source_num = srcething
    source = source_list[srcething]
    srcname = source_list[srcething]
    #print source_list, source_lookup
    print sorted(srce_times)
    #print dtimescans[sdi]
    #timeflag_count = srce_times[source_list.index(srcething)]
    for i in xrange(len(sbline[sdi])):
        # print sdi,bline
        bline = sbline[sdi][i]
        ntimescans = baseline_dict[sdi][bline]
        print ntimescans
        times_array = dtimescans[sdi][sbline[sdi].index(bline)]
        # print ntimescans
        # print baseline_dict[sdi]
        if dobline_specific == 'yes':
            parameters = source_base_parms[sdi][bline]
            # print 'using parms', srcname, parameters, bline
        else:
            # parameters = parameters
            parameters = source_base_parms[sdi]
            # print 'using parms', srcname, parameters, bline
        if do_loadflag == 'yes':
            if write_old_flags == 'yes':
                write_olds = True
            else:
                write_olds = False
            if mem_low:
                print 'Low memory: loading old flags from file'
                # continue
            else:
                if not type(load_results[sdi][bline]).__module__ == np.__name__:
                    # if not old_flags:
                    old_flags = np.zeros([nif,npol,ntimescans,nchan],dtype=np.int8)
                    print 'Creating empty flag array'
                    # print srcname, bline
                else:
                    print 'Using existing flag array'
                    # print srcname, bline
                    old_flags = load_results[sdi][bline]
                    # print old_flags.shape,nif,npol,ntimescans,nchan
             
        else:
            write_olds = False
            # old_flags = np.zeros([nif,npol,ntimescans,nchan],dtype=int)
            # print 'Creating empty flag array'
            # print srcname, bline 
        # print 'write_olds is:', write_olds 
        for j in xrange(nif):
            if mem_low:
                # print jobcount, total_jobs
                joblist.append([mem_low,flagsexist,jobcount,total_jobs,bline,j,nif,source,nchan,ntimescans,times_array,polnames,path2folder,experiment,name,phasecal,lovell_bline,telescope,srcname,zero_level,do_lovell_cross,coadd_polarization_flags,pickle_list,flagging_options,parameters,uvdata.name,uvdata.klass,uvdata.disk,uvdata.seq,singlesource,source_num,orderedpols, write_olds, timeflag_count, metadata, select_stokes,flag_all_stokes,coadd_zero_flags,coadd_lovell_flags])
            else:
                joblist.append([mem_low,flagsexist,jobcount,total_jobs,bline,j,nif,source,nchan,ntimescans,times_array,polnames,path2folder,experiment,name,phasecal,lovell_bline,telescope,srcname,zero_level,do_lovell_cross,coadd_polarization_flags,pickle_list,flagging_options,parameters,uvdata.name,uvdata.klass,uvdata.disk,uvdata.seq,singlesource,source_num,orderedpols,old_flags[j], write_olds, timeflag_count, metadata, select_stokes,flag_all_stokes,coadd_zero_flags,coadd_lovell_flags])
                # print old_flags[j].nbytes
                # size_total += old_flags[j].nbytes
            jobcount +=1    
    
#print ts(sbline,verbose=True) 
#print ts(joblist[0],verbose=True)
#print size_total





########################################################################
### Parallelisation block for full flagger process, calls function 
### which runs lovell dropouts, zero-level dropouts and flagger on 
### each source, baseline and IF separately
########################################################################


send_q = mp.Queue()
rec_q = mp.Queue()
ncpus = mp.cpu_count()

if NCPU == 0:
    NCPU = ncpus
elif ncpus > NCPU and not NCPU == 0:
    ncpus = NCPU

for i in xrange(len(joblist)):
    print 'Sending to queue job', i+1, 'of', len(joblist)
    # print "CPU", c, "will flag baselines and IFs: \n"
    send_q.put(joblist[i])


if use_fits == 'yes':
    for i in xrange(ncpus):
        proc = mp.Process(target=flagger_fits_process, args=(send_q,rec_q,i))
        proc.start()
        print 'Starting process on CPU:', i
else:
    for i in xrange(ncpus):
        proc = mp.Process(target=flagger_process, args=(send_q,rec_q,i))
        proc.start()
        print 'Starting process on CPU:', i




## where does this list need to be? before or after the for bline in blist[mycore]?
lovell_bline = []

# print "Source:", name, ", Baseline:", bline, ", IF:", j+1, " on CPU", mycore, ", Job #:", relsc+1, " of", mynbline
# print "\n CPU", mycore, "has finished all jobs...\n"


#### End Parallelization and wait for all processes to finish before 
#### reading everything into fg file.


results = []
for i in xrange(len(joblist)):
    results.append(rec_q.get())

for i in xrange(ncpus):
    send_q.put('STOP')

# print(results)

print "\nTime taken to append visibilities (hh:mm:ss):", time2hms(time.time()-ti)





##################################################################
#### Now write all of the final flag table and add to the data ###
##################################################################



try:
    time_limit = flag_coinc_times
except NameError:
    time_limit = 0

try:
    count_limit = flag_coinc_chans
except NameError:
    count_limit = 0
    

try:
    zero_flag_allif
except NameError:
    zero_flag_allif = 'yes'

if select_stokes == 'yes':
    if flag_all_stokes == 'yes':
        coadd_polarization_flags = 'yes'
        npflags = '1111'
    elif flag_all_stokes == 'no' and coadd_polarization_flags == 'yes':
        newstr=['0','0','0','0']
        for pol in polnames:
            # print orderedpols[pol]
            newstr[orderedpols[pol]]='1'
            npflags = "".join(newstr)



print " Creating first FG text file from pickle files..."

startwrite = time.time()
tcount = 0
new_fg_tables = {}
#timeflag_count = 0
#flag_count = 0
####

#### Need to create timeflag_count dictionary, should be dictionary containing full list of times for each source###
#### Will need to combine list of time_arrays for all blines

if use_fits == 'no':


    if coadd_polarization_flags == 'yes':
        print " Coadding polarizations for FG table..."
        # tabnum = create_blank_fg(data, path2folder, no_rows=0)

        lovell_row = 1
        name = uvdata.name
        for srce in source_list:
            timeflag_count = {}
            source = source_list[srce]
            if source in nodata_list:
                continue
            source_num = srce

            ## use dtimescans to create dictionary of total timescans 
            ## per source to keep count when writing flags 
            print 'Creating timecount dictionary'
            sdi = source_lookup.index(source_list[srce])
            temparr = flatten(dtimescans[sdi])
            for i in xrange(len(temparr)):
                if temparr[i] not in timeflag_count:
                    timeflag_count[temparr[i]] = 0
                else: 
                    timeflag_count[temparr[i]] += 1
            del temparr
            flag_count = 0
            tabnum = create_blank_fg(uvdata, path2folder, no_rows=0)
            new_fg_tables[source] = tabnum 
            orig_tabnum = tabnum
            if lovellmk2flag == 1:
                flag_lovellmk2(uvdata, tabnum,lovellmk2flag,baseline_removal)
            # if nsource > 1:
            #   name = source_list[source]
            for bline in sbline[sdi]:
                bline_name = bline
                bline = bline.rsplit('-')
                for j in xrange(nif):
                    atleastone = False
                    flag = []
                    if not flag:
                        del flag
                    for pol in polnames:
                        pname = str(name)  +"__"+ str(source) + '__IF' + str(j+1) + '__' + str(bline_name) +  '__' + pol
                        # print "Reading file:", pname
                        if os.path.exists(path2folder+str(pname)+'.p') == True:
                            times_array = pickle.load(open(path2folder + str(pname) +".info", "rb"))
                            atleastone = True
                            print "Reading file:", pname
                            flag_one = pickle.load(open(path2folder + str(pname) +".p", "rb"))
                            try:
                                flag = np.where(flag>0,flag,flag_one)
                            except NameError:
                                flag = flag_one
                    if atleastone:
                        if count_limit > 0 or time_limit > 0:
                            print '\nExtra flagging requested'
                            print 'Flagging every %i coincident channels and every %i coincident times\n' % (count_limit, time_limit)
                            flag = coinc_flags_noncumu(flag,time_limit, count_limit)
                        times_array = pickle.load(open(path2folder + str(pname) +".info", "rb"))
                        ifnum = j+1
                        pflags = '1111'
                        totalchan = nchan
                        totaltime = len(times_array)
                        reason_def = reason()
                        # if totalchan > totaltime:
                        # only writing in time order currently #
                        tabnum, timeflag_count, flag_count =  write_fg_table_times_alt4(uvdata, tabnum, pflags, flag, times_array, source_num, ifnum, reason_def, bline, path2folder, timeflag_count, flag_count, no_rows=0)
                        del times_array
                        del flag
                        del flag_one
                        for pol in polnames:
                            pname = str(name)  +"__"+ str(source) + '__IF' + str(j+1) + '__' + str(bline_name) +  '__' + pol
                            os.remove(path2folder+str(pname) +".p")
                            os.remove(path2folder+str(pname) +".info")
                    for pol in polnames:
                        plname = str(name) + "__" +str(source) + "__" + str(bline_name) + "_IF" + str(j+1) + '__' + pol +"__lovell_dummy"
                        # print "Reading Lovell file:", pname
                        if os.path.exists(path2folder + str(plname) + ".con") == True:
                            print "Reading Lovell file:", plname
                            con = pickle.load(open(path2folder + str(plname) + ".con", "rb"))
                            dropouts = pickle.load(open(path2folder + str(plname) + ".dropouts", "rb"))
                            pflags = '1111'
                            ifnum = nif
                            channum = nchan
                            reason_def = 'LOVELL ' + strftime("%d-%b-%y")
                            tabnum, timeflag_count, flag_count = write_lovell_times(uvdata, tabnum, source_num, pflags, ifnum, channum, reason_def, con, dropouts, bline, path2folder, timeflag_count, flag_count, no_rows=0)
                            os.remove(path2folder + str(plname) + ".con")
                            os.remove(path2folder + str(plname) + ".info")
                            os.remove(path2folder + str(plname) + ".dropouts")
                    for pol in polnames:
                        zname = str(name) + "__" +str(source)+"__" + str(bline_name) + "__IF" + str(j+1) + "__" + pol + "__zeros_dummy"
                        if os.path.exists(path2folder + str(zname) + ".con") == True:
                            con = pickle.load(open(path2folder + str(zname) + ".con", "rb"))
                            zeros = pickle.load(open(path2folder + str(zname) + ".zeros", "rb"))
                            pflags = '1111'
                            if zero_flag_allif == 'yes':
                                ifnum = nif
                            elif zero_flag_allif == 'no':
                                ifnum = j+1
                            channum = nchan
                            reason_def = 'ZEROS ' + strftime("%d-%b-%y")
                            tabnum, timeflag_count, flag_count = write_zeros_times(uvdata, tabnum, source_num, pflags, ifnum, channum, reason_def, con, zeros, bline, path2folder, timeflag_count, flag_count, zero_flag_allif, no_rows=0)
                            os.remove(path2folder + str(zname) + ".con")
                            os.remove(path2folder + str(zname) + ".info")
                            os.remove(path2folder + str(zname) + ".zeros")
                    for pol in polnames:
                        zisname = str(name)  +"__"+ str(source) + "__" + str(bline_name) + "__IF" + str(j+1) + "__" + pol + "__inscan_zeros_dummy"
                        if os.path.exists(path2folder + str(zisname) + ".con") == True:
                            con = pickle.load(open(path2folder + str(zisname) + ".con", "rb"))
                            zeros = pickle.load(open(path2folder + str(zisname) + ".zeros", "rb"))
                            pflags = '1111'
                            if zero_flag_allif == 'yes':
                                ifnum = nif
                            elif zero_flag_allif == 'no':
                                ifnum = j+1
                            channum = nchan
                            reason_def = 'ISZEROS ' + strftime("%d-%b-%y")
                            tabnum, timeflag_count, flag_count = write_zeros_times(uvdata, tabnum, source_num, pflags, ifnum, channum, reason_def, con, zeros, bline, path2folder, timeflag_count, flag_count, zero_flag_allif, no_rows=0)
                            os.remove(path2folder + str(zisname) + ".con")
                            os.remove(path2folder + str(zisname) + ".info")
                            os.remove(path2folder + str(zisname) + ".zeros")
                if tabnum > orig_tabnum:
                    change = orig_tabnum-tabnum
                    for i in xrange(1,change+1):
                        new_fg_tables[source].append(tabnum+1)


    # if coadd_polarization_flags == 'no':

    else:
        print " Creating separate polarizations flags for FG table..."

        lovell_row = 1
        name=uvdata.name
        for srce in source_list:
            timeflag_count = {}
            source = source_list[srce]
            if source in nodata_list:
                continue
            source_num = srce
            # print srce, source, source_num

            ## use dtimescans to create dictionary of total timescans 
            ## per source to keep count when writing flags 
            print 'Creating timecount dictionary'
            sdi = source_lookup.index(source_list[srce])
            temparr = flatten(dtimescans[sdi])
            for i in xrange(len(temparr)):
                if temparr[i] not in timeflag_count:
                    timeflag_count[temparr[i]] = 0
                else: 
                    timeflag_count[temparr[i]] += 1
            del temparr

            flag_count = 0
            tabnum = create_blank_fg(uvdata, path2folder, no_rows=0)
            new_fg_tables[source] = tabnum 
            orig_tabnum = tabnum
            # print 'lovell is', lovellmk2flag
            if lovellmk2flag == 1:
                flag_lovellmk2(uvdata,tabnum,lovellmk2flag,baseline_removal)
            # if nsource > 1:
            # name = source_list[srce]
            for bline in sbline[sdi]:
                bline_name = bline
                bline = bline.rsplit('-')
                for j in xrange(nif):
                    for pol in polnames:
                        pname = str(name)  +"__"+ str(source) + '__IF' + str(j+1) + '__' + str(bline_name) + '__' + pol
                        print "Reading file:", pname
                        if os.path.exists(path2folder+str(pname)+'.p') == True:
                            flag = pickle.load(open(path2folder + str(pname) +".p", "rb"))
                            times_array = pickle.load(open(path2folder + str(pname) +".info", "rb"))
                            # flag = np.load(path2folder + str(pname) +".p.npy")
                            # times_array = np.load(path2folder + str(pname) +".info.npy")
                            # print times_array
                            ifnum = j+1
                            reason_def = reason()
                            if count_limit > 0 or time_limit > 0:
                                print '\nExtra flagging requested'
                                print 'Flagging every %i coincident channels and every %i coincident times\n' % (count_limit, time_limit)
                                flag = coinc_flags_noncumu(flag,time_limit, count_limit)
                            if pol == 'RR' or pol == 'XX':
                                pflags = '1000'
                            if pol == 'LL' or pol == 'YY':
                                pflags = '0100'
                            if pol == 'RL' or pol == 'XY':
                                pflags = '0010'
                            if pol == 'LR' or pol == 'YX':
                                pflags = '0001'
                            # pflags = '1111'
                            # print time_limit,count_limit
                            totalchan = nchan
                            totaltime = len(times_array)
                            # if totalchan > totaltime:
                            testflagout = open('/local/python/RFI/SERPent/memory_tests_0316/testflagout_file.txt','ab')
                            print >> testflagout, pol, ifnum, 'flag', flag
                            testflagout.close()
                            tabnum, timeflag_count, flag_count =  write_fg_table_times_alt4(uvdata, tabnum, pflags, flag, times_array, source_num, ifnum, reason_def, bline, path2folder, timeflag_count, flag_count, no_rows=0)
                            os.remove(path2folder+str(pname) +".p")
                            os.remove(path2folder+str(pname) +".info")
                            del flag
                            del times_array
                            zname = str(name) + "__"+ str(source) +"__" + str(bline_name) + "__IF" + str(j+1) + "__" + pol + "__zeros_dummy"
                            if os.path.exists(path2folder + str(zname) + ".con") == True:
                                con = pickle.load(open(path2folder + str(zname) + ".con", "rb"))
                                zeros = pickle.load(open(path2folder + str(zname) + ".zeros", "rb"))
                                if coadd_zero_flags == 'no':
                                    if pol == 'RR' or pol == 'XX':
                                        pflags = '1000'
                                    if pol == 'LL' or pol == 'YY':
                                        pflags = '0100'
                                    if pol == 'RL' or pol == 'XY':
                                        pflags = '0010'
                                    if pol == 'LR' or pol == 'YX':
                                        pflags = '0001'
                                else:
                                    pflags = '1111'
                                if zero_flag_allif == 'yes':
                                    ifnum = nif
                                elif zero_flag_allif == 'no':
                                    ifnum = j+1
                                channum = nchan
                                reason_def = 'ZEROS ' + strftime("%d-%b-%y")
                                tabnum, timeflag_count, flag_count = write_zeros_times(uvdata, tabnum, source_num, pflags, ifnum, channum, reason_def, con, zeros, bline, path2folder, timeflag_count, flag_count, zero_flag_allif, no_rows)
                                os.remove(path2folder + str(zname) + ".con")
                                os.remove(path2folder + str(zname) + ".info")
                                os.remove(path2folder + str(zname) + ".zeros")
                            zisname = str(name) +"__"+ str(source) + "__" + str(bline_name) + "__IF" + str(j+1) + "__" + pol + "__inscan_zeros_dummy"
                            if os.path.exists(path2folder + str(zisname) + ".con") == True:
                                con = pickle.load(open(path2folder + str(zisname) + ".con", "rb"))
                                zeros = pickle.load(open(path2folder + str(zisname) + ".zeros", "rb"))
                                if coadd_zero_flags == 'no':
                                    if pol == 'RR' or pol == 'XX':
                                        pflags = '1000'
                                    if pol == 'LL' or pol == 'YY':
                                        pflags = '0100'
                                    if pol == 'RL' or pol == 'XY':
                                        pflags = '0010'
                                    if pol == 'LR' or pol == 'YX':
                                        pflags = '0001'
                                else:
                                    pflags = '1111'
                                if zero_flag_allif == 'yes':
                                    ifnum = nif
                                elif zero_flag_allif == 'no':
                                    ifnum = j+1
                                channum = nchan
                                reason_def = 'ISZEROS ' + strftime("%d-%b-%y")
                                tabnum, timeflag_count, flag_count = write_zeros_times(uvdata, tabnum, source_num, pflags, ifnum, channum, reason_def, con, zeros, bline, path2folder, timeflag_count, flag_count, zero_flag_allif, no_rows)
                                os.remove(path2folder + str(zisname) + ".con")
                                os.remove(path2folder + str(zisname) + ".info")
                                os.remove(path2folder + str(zisname) + ".zeros")
                            plname = str(name)+"__"+ str(source) + "__" + str(bline_name) + "_IF" + str(j+1) + "__" + pol + "__lovell_dummy"
                            if os.path.exists(path2folder + str(plname) + ".con") == True:
                                con = pickle.load(open(path2folder + str(plname) + ".con", "rb"))
                                dropouts = pickle.load(open(path2folder + str(plname) + ".dropouts", "rb"))
                                if coadd_lovell_flags == 'no':
                                    if pol == 'RR' or pol == 'XX':
                                        pflags = '1000'
                                    if pol == 'LL' or pol == 'YY':
                                        pflags = '0100'
                                    if pol == 'RL' or pol == 'XY':
                                        pflags = '0010'
                                    if pol == 'LR' or pol == 'YX':
                                        pflags = '0001'
                                else:
                                    pflags = '1111'
                                ifnum = nif
                                channum = nchan
                                reason_def = 'LOVELL ' + strftime("%d-%b-%y")
                                tabnum, timeflag_count, flag_count = write_lovell_times(uvdata, tabnum, source_num, pflags, ifnum, channum, reason_def, con, dropouts, bline, path2folder, timeflag_count, flag_count, no_rows=0)

                            readtime=time.time()
                            if tcount == 0:
                                print 'readtime is',time2hms(readtime - startwrite)
                            else:
                                print 'readtime is',time2hms(readtime - oldreadtime)
                            oldreadtime = readtime
                            tcount +=1
                if tabnum > orig_tabnum:
                    change = orig_tabnum-tabnum
                    for i in xrange(1,change+1):
                        new_fg_tables[source].append(tabnum+1)



    # Calculate time taken and write log file
    print '\n---------------------------------------------\n'
    print 'Finished on %s' % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()),''
    print "Time taken (hh:mm:ss):", time2hms(time.time()-ti),"\n"
    print 'Final flag tables are:'#, new_fg_tables
    print "Source: Flag tables"
    for srce in source_list:
        sdi = source_lookup.index(source_list[srce])
        source = source_list[srce]
        if source in nodata_list:
            # print 'Source %s not flagged, no data found' % source
            continue
        # print srce,sdi,source
        print "%s: %s" % (source, new_fg_tables[source])
    if nodata_list:
        print 'Source(s) %s not flagged, no data found' % ', '.join(nodata_list)
    
    print '\n---------------------------------------------\n'




else:


    # Read in files to write all into single fits file ###
    src = []
    subary = []
    frqid = []
    ifs = []
    chans = []
    pflagsarr = []
    antsarr = []
    time_range = []
    reasons = []
    numtables = 1
    num = 1
    writetime1 = time.time()


    while numtables > 0:
        for srce in source_list:
            sdi = source_lookup.index(source_list[srce])
            source = source_list[srce]
            if source in nodata_list:
                continue
            for bline in sbline[sdi]:
                for j in xrange(nif):
                    pname = str(name) +"__"+ str(source)  + "__IF" + str(j+1) + "__" + str(bline)
                    print pname
                    numtables = 0
                    startload = time.time()
                    if os.path.exists(path2folder + str(pname) + ".src") == True:
                        numtables = 0
                        timeflags = pickle.load(open(path2folder + str(pname)+".fgtimes", "rb"))
                        for key in timeflags.keys():
                            try:
                                if key in timeflag_count:
                                    if timeflag_count[key] >= 600000:
                                        print "We're going to need more than one FG file"
                                        numtables += 1
                                        write_full_fits_fg(experiment,metadata,path2folder,nif,bline,src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons,num)
                                        timeflag_count = {}
                                        num += 1
                                    elif key not in timeflag_count:
                                        timeflag_count[key] = timeflags[key]
                                    else:
                                        timeflag_count[key] = timeflag_count[key]+timeflags[key]
                            except NameError:
                                timeflag_count = timeflags
                                # extra check that it's not all in one file
                                if timeflag_count[key] >= 600000:
                                    print "Single Source/baseline/IF requires multiple flag tables"
                                    print "Seek help" 
                                    break
                        src += pickle.load(open(path2folder + str(pname)+".src", "rb"))
                        # print src
                        subary += pickle.load(open(path2folder + str(pname)+".subary", "rb"))
                        frqid += pickle.load(open(path2folder + str(pname)+".frqid", "rb"))
                        ifs += pickle.load(open(path2folder + str(pname)+".ifs", "rb"))
                        chans += pickle.load(open(path2folder + str(pname)+".chans", "rb"))
                        pflagsarr += pickle.load(open(path2folder + str(pname)+".pflag", "rb"))
                        antsarr += pickle.load(open(path2folder + str(pname)+".ants", "rb"))
                        time_range += pickle.load(open(path2folder + str(pname)+".trange", "rb"))
                        reasons += pickle.load(open(path2folder + str(pname)+".reasons", "rb"))


                        os.remove(path2folder + str(pname) +".subary")
                        os.remove(path2folder + str(pname) + ".src")
                        os.remove(path2folder + str(pname) + ".frqid")
                        os.remove(path2folder + str(pname) + ".ifs")
                        os.remove(path2folder + str(pname) + ".chans")
                        os.remove(path2folder + str(pname) + ".pflag")
                        os.remove(path2folder + str(pname) + ".ants")
                        os.remove(path2folder + str(pname) + ".trange")
                        os.remove(path2folder + str(pname) + ".reasons")
                        os.remove(path2folder + str(pname) + ".fgtimes")
                        time2load = time.time()
                        print "File loading (hh:mm:ss):", time2hms(time2load-startload),"\n"
        timeflag_count = {} #Assume no two sources can be observed at same time so reset time flag count for each source
         

    write_full_fits_fg(experiment,metadata,path2folder,nif,bline,src,subary,frqid,ifs,chans,pflagsarr,antsarr,time_range,reasons,num)      
    print "Total fits file creation (hh:mm:ss):", time2hms(time.time()-writetime1),"\n"

    new_fg_tables = []
    # Now load fitsfiles and copy to original datafile
    for i in xrange(1,num+1):
        if os.path.exists(path2folder+experiment+'_'+str(num)+'.fits'):
            print 'file exists'
        tablefilename,tablefileklass,tablefiledisk,tablefileseq = aipsfitld(path2folder,experiment,i,AIPS.userno,uvdata.disk)
        print tablefilename
        outnum = aipsfgcop(tablefilename,tablefileklass,tablefiledisk,tablefileseq,uvdata,innum=i)
        new_fg_tables.append(outnum)



    # Calculate time taken and write log file
    print '\n---------------------------------------------\n'
    print 'Finished on %s' % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()),''
    print "Time taken (hh:mm:ss):", time2hms(time.time()-ti),"\n"
    if len(new_fg_tables) > 1:
        print 'Final flag tables are: %s' % ', '.join(new_fg_tables)
    else:
        print 'Final flag table is: %s' % new_fg_tables[0]
    if nodata_list:
        print 'Source(s) %s not flagged, no data found' % ', '.join(nodata_list)
    
    print '\n---------------------------------------------\n'








## On the new version of AIPS a verb called REFLG is available which concatenates 
## rows in the FG table more efficiently. This will be very useful when flagging 
## L band data as the RFI is no longer '1 dimensional' but 2D, but it also has a 
## function where you can flag the entire time scan or channel if x% of the time 
## scan or channel is flagged... Very useful as some RFI in a channel (e.g.) might 
## be too weak to be picked up, but most of the rest of the channel is flagged... 
## Makes the flagger more robust and efficient with its rows in the FG table...


## Write out the REFLG'd FG table to a text file...
## Try and use REFLG, if AIPS version does not contain verb then skip over it.

#print 'flaggin options are', flagging_options
'''
try:
    reflg(cparm_1, cparm_2, cparm_3, cparm_4, cparm_5, cparm_6, cparm_7,path2folder,experiment,experiment_klass,experiment_seq,experiment_disk,flagging_options)
    tbout(path2folder,experiment,experiment_klass,experiment_seq,experiment_disk)
except:
    print "\n Current AIPS Version does not have REFLG verb. Need a newer version."
    print " Passing over REFLG function..."
'''



#for src in source_lookup:
 #   print " %i: %s" % (s, source_list[s])


if write2log == 'yes':
    try:
        log = open(path2folder + experiment + "_run_log", 'a')
        log.flush()
        t = time2hms(time.time()-ti)
        d = strftime("%d-%b-%y")
        dt = strftime("%H:%M:%S", localtime())
        print >> log, "NCPUs:", NCPU, ", Experiment Name:", experiment, ", Run #", runs, ", Time taken (hh:mm:ss):", time2hms(time.time()-ti), ", Date and Time: %s" % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()), "SERPent Version:", version, '\n'
        log.close()
    except:
        pass



#############################################################################
#############################################################################
########### END OF MAIN SERPENT CODE
#############################################################################
#############################################################################


'''
print h.heap()


from pympler import summary
sum1 = summary.summarize(all_objects)
summary.print_(sum1)             
'''
