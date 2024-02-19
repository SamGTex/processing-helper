#
# Set of common definitions used by all filters.  Name definitions mostly
#
from icecube import icetray, dataclasses

class Domain(object):

    SPTS = 1
    SPS = 2
    WISC = 3
    DESY = 4
    UNKOWN = 5
    ACCESS = 6

def getDomain():
    import socket
    host = socket.gethostname()
    domain = host[host.find(".")+1:]

    # ease filterscript testing at access.spts where tables 
    #   are in different locations
    if host[0:6] == 'access':
        return Domain.ACCESS
    elif 'sps' in domain:
        return Domain.SPS
    elif 'spts' in domain:
        return Domain.SPTS
    else:
        return None


#
# functions to conditionally execute based on split
#
InIceSplitter = 'InIceSplit'
IceTopSplitter = 'IceTopSplit'
NullSplitter = 'NullSplit'
SLOPSplitter = 'SLOPSplit'

def print_split(frame):
    if frame.Stop == icetray.I3Frame.Physics:
        print((frame['I3EventHeader'].sub_event_stream))


#
#DOMLaunches
#

#InIce Bad dom cleaned list:
CleanInIceRawData = 'CleanInIceRawData'

#IceTop Bad dom cleaned list:
CleanIceTopRawData = 'CleanIceTopRawData'

#
# RecoPulseSeries/Maps
#
UncleanedInIcePulses       = 'UncleanedInIcePulses'
SplitUncleanedInIcePulses  = 'SplitUncleanedInIcePulses'
SplitRTCleanedInIcePulses  = 'SplitRTCleanedInIcePulses'
SplitUncleanedITPulses  = 'SplitUncleanedITPulses'
SLOPPulses = 'SLOPPulseMask'

IceTopPulses_HLC = 'IceTopPulses_HLC'
IceTopPulses_SLC = 'IceTopPulses_SLC'
HLCTankPulses = 'HLCTankPulses'
CleanedHLCTankPulses = 'CleanedHLCTankPulses'
IceTopBadDoms = 'IceTopBadDOMs'
IcetopBadTanks = 'IceTopBadTanks'
TankPulseMergerExcludedTanks = 'TankPulseMergerExcludedTanks'
ClusterCleaningExcludedTanks = 'ClusterCleaningExcludedTanks'
IsSmallShower = 'IsSmallShower'

DSTPulses = 'DSTPulses'
InIceDSTPulses = 'InIce' + DSTPulses
IceTopDSTPulses = 'IceTop' + DSTPulses

CleanedMuonPulses = 'CleanedMuonPulses'


# trigger I3Bools from TriggerChecker
inicesmttriggered = 'InIceSMTTriggered'
icetopsmttriggered = 'IceTopSMTTriggered'
inicestringtriggered = 'InIceStringTriggered'
deepcoresmttriggered = 'DeepCoreSMTTriggered'
volumetrigtriggered = 'VolumeTriggerFlag'
slowparticletriggered = 'SlowParticleTriggered'
fixedratetriggered = 'FixedRateTriggered'
icetopvolumetriggered = 'IceTopVolumeTriggered'
scintminbiastriggered = 'ScintMinBiasTriggered'
iceactsmttriggered = 'IceActSMTTriggered'

#TODO: make sure configIDs have not changed in 20XX, 
# commit any changes at last minute for test run
# The triggerConfigIDs that are specified.
deepcoreconfigid = 1011
inicesmtconfigid = 1006
inicestringconfigid = 1007
slowparticleconfigid = 24002
volumetriggerconfigid = 21001
icetopsmtconfigid = 102
fixedrateconfigid = 23050
scintminbiasconfigid = 210
iceactsmtconfigid = 30002

#
# Shared reconstructions:
#
muon_linefit = 'PoleMuonLinefit'
muon_llhfit = 'PoleMuonLlhFit'

## I3Doubles for the cascaders -- does this belong here?
cascade_toi = 'PoleCascadeFilter_ToiVal'
cascade_lfvel = 'PoleCascadeFilter_LFVel'
cascade_cscdllh = 'PoleCascadeFilter_CscdLlh'


## standard items
rawdaqdata = 'I3DAQData'
seatbeltdata = 'I3DAQDataTrimmed'
icetoprawdata = 'I3DAQDataIceTop'
superdst = 'I3SuperDST'
eventheader = 'I3EventHeader'
triggerhierarchy = 'I3TriggerHierarchy'
qtriggerhierarchy = 'QTriggerHierarchy'
smalltriggerhierarchy = 'DSTTriggers'
dst  = 'I3DST'
dstheader  = 'I3DSTHeader'

filter_mask = 'FilterMask'
qfilter_mask = 'QFilterMask'
icetop_vemcal = 'I3VEMCalData'

client_start_time = 'ClientStartTime'
homogenized_qtot = 'Homogenized_QTot'

estres_fit_name = 'SplineMPE_Estres'
estres_homogenized_qtot = 'Estres_Homogenized_QTot'

## follow-up frame items
alert_candidate_list = 'AlertNamesPassed'
alert_candidate_short_message = 'AlertShortFollowupMsg'
alert_candidate_full_message = 'AlertFullFollowupMsg'

greco_candidate_list = 'GRECONames'
greco_short_followup_message = 'GRECO_valuesdict'
greco_full_followup_message = 'GRECOFullFollowupMsg'

## Moni Frame items from other modules
flaring_dom_message = 'FlaringDOMs'


#
# FilterNames
#
FilterMinBias = 'FilterMinBias_13'
MuonFilter = 'MuonFilter_13'
CosmicRayFilter = 'CosmicRayFilter_13'
SunFilter = 'SunFilter_13'
MoonFilter = 'MoonFilter_13'
CascadeFilter = 'CascadeFilter_13'
EHEFilter = 'EHEFilter_13'
OnlineL2Filter = 'OnlineL2Filter_17'
GFUFilter = 'GFUFilter_17'
DeepCoreFilter = 'DeepCoreFilter_13'
FSSFilter = 'FSSFilter_13'
FSSCandidate = 'FSSCandidate_13'
VEFFilter = 'VEF_13'
LowUpFilter = 'LowUp_13'
SlopFilter = 'SlopFilter_13'
FixedRateFilterName = 'FixedRateFilter_13'
HighQFilter = 'HighQFilter_17'
HESEFilter = 'HESEFilter_15'
MESEFilter = 'MESEFilter_15'
EHEAlertFilter = 'EHEAlertFilter_15'
EHEAlertFilterHB = 'EHEAlertFilterHB_15'
MonopoleFilter = 'MonopoleFilter_16'
ScintMinBiasFilter = 'ScintMinBias_16'
IceTopTwoStationFilter = 'IceTop_InFill_STA2_17'
EstresAlertFilter = 'EstresAlertFilter_18'
IceActFilter = 'IceActTrigFilter_18'
GRECO_OnlineFilter = 'GRECOOnlineFilter_19'

filter_pairs = [(FilterMinBias,1000),
                (MuonFilter,1),
                (MoonFilter,1),
                (SunFilter,1),
                (OnlineL2Filter,1),
                (GFUFilter,1),
                (CascadeFilter,1),
                ('IceTopSTA3_13',10),
                ('IceTopSTA5_13',1),
                ('IceTop_InFill_STA3_13',10),	                
                ('InIceSMT_IceTopCoincidence_13',100),
                (DeepCoreFilter,1),
                (EstresAlertFilter,1),
                (GRECO_OnlineFilter,1),
                (FSSFilter,1),
                (FSSCandidate,0),
                (VEFFilter,1),
                (LowUpFilter,1),
                (SlopFilter, 1),
                (FixedRateFilterName,1),
                (HighQFilter, 1),
                (HESEFilter, 1),
                (MESEFilter, 1),
                ('I3DAQDecodeException', 1),
                (EHEAlertFilter,1),
                (EHEAlertFilterHB,1),
                (MonopoleFilter,1),
                (ScintMinBiasFilter,1),
                (IceTopTwoStationFilter,1),
                ('SDST_InIceSMT_IceTopCoincidence_13',1),
                ('SDST_IceTopSTA3_13',1),
                ('SDST_IceTop_InFill_STA3_13',1),
                (IceActFilter,1)
                ]

filter_streams,filter_prescales = zip(*filter_pairs)

## Filters that select events for "SuperDST only" transmission
sdst_pairs = []

sdst_streams = ()
sdst_prescales = ()

## List of filters that want *ALL* waveforms as well as SuperDST
filters_keeping_allraw = [FilterMinBias,
                          'IceTopSTA3_13',
                          'IceTopSTA5_13',
                          'IceTop_InFill_STA3_13',
                          'InIceSMT_IceTopCoincidence_13',
                          HighQFilter,
                          'I3DAQDecodeException',
                          HESEFilter,
                          FixedRateFilterName,
                          EHEAlertFilter,
                          EHEAlertFilterHB]

q_frame_keeps = [smalltriggerhierarchy,
                 seatbeltdata,
                 icetoprawdata,
                 superdst,
                 eventheader,
                 'JEBClientInfo',
                 'PFClientInfo',
                 qfilter_mask,
                 dst,
                 dstheader,
                 icetop_vemcal]

sdstall_keeps = [superdst,
                 seatbeltdata,
                 smalltriggerhierarchy,
                 icetoprawdata]

null_split_keeps = [filter_mask]

# objects the HESE filter wants to keep in case of a positive decision
hese_reco_keeps = [
    "HESE_HomogenizedQTot",
    "HESE_CausalQTot",
    "HESE_VHESelfVeto",
    "HESE_VHESelfVetoVertexPos",
    "HESE_VHESelfVetoVertexTime",
    "HESE_CascadeLlhVertexFit",
    "HESE_CascadeLlhVertexFitParams",
    "HESE_MuonImprovedLineFit",
    "HESE_MuonImprovedLineFitParams",
    "HESE_SPEFitSingle",
    "HESE_SPEFitSingleFitParams",
    "HESE_SPEFit2",
    "HESE_SPEFit2FitParams",
    "HESE_llhratio",
]

estres_reco_keeps = [
   	estres_homogenized_qtot,
    "Estres_CausalQTot",
    "Estres_p_miss",
    "Estres_direct_length",
    estres_fit_name,
    estres_fit_name+"FitParams",
]

ehe_reco_keeps = ["PoleEHEOpheliaParticle_ImpLF",
                  "PoleEHEOphelia_ImpLF"]

greco_reco_keeps = ['GRECO_SPEFit11',
                    'GRECO_Variables',
                    'GRECO_bdt_2019',
                    'IC2018_LE_L3_bools']

inice_split_keeps = [filter_mask,
        muon_linefit,
        muon_llhfit,
        muon_llhfit+'FitParams',
        cascade_toi,
        cascade_lfvel,
        cascade_cscdllh,
        'CorsikaMoonMJD', # only affects simulation: saves an I3Double with the generated fake MJD
        'CorsikaSunMJD', # only affects simulation: saves an I3Double with the generated fake MJD
        dst,
        dstheader,
        flaring_dom_message,
        'PoleEHESummaryPulseInfo', # EHE info
        homogenized_qtot,  # high Qtot filter info
        'MPInfoDict', #Monopole keeps             
    ] + hese_reco_keeps + ehe_reco_keeps + estres_reco_keeps + greco_reco_keeps

icetop_split_keeps = [filter_mask,
                      #IceTopRecoStandard,
                      #IceTopRecoSmallShowers,
                      #IceTopRecoStandard + 'Params',
                      #IceTopRecoSmallShowers + 'Params',
                      IsSmallShower,]

slop_split_keeps = [filter_mask,
                    'SLOPTuples_CosAlpha_Launches',
                    'SLOPTuples_RelV_Launches',
                    'SLOPLaunchMapTuples']

keep_nofilterpass = [smalltriggerhierarchy,
                     superdst,
                     eventheader,
                     'JEBClientInfo',
                     'PFClientInfo',
                     qfilter_mask,
                     dst,
                     dstheader,
                     icetop_vemcal]

keep_dst_only = ['JEBClientInfo',
                 'PFClientInfo',
                 dst,
                 dstheader,
                 icetop_vemcal]

alert_followup_keeps = [alert_candidate_list,
                        alert_candidate_short_message,
                        alert_candidate_full_message,
                        greco_candidate_list,
                        greco_short_followup_message,
                        greco_full_followup_message]

onlinel2filter_keeps = []  ## added to in l2
ofufilter_keeps = []  ## not used, but kept for compat..
gfufilter_keeps = []  ## added to in gfu

online_bad_doms = [icetray.OMKey(1,46), icetray.OMKey(5,21), icetray.OMKey(5,22),
                   icetray.OMKey(7,35), icetray.OMKey(7,36), icetray.OMKey(7,37),
                   icetray.OMKey(7,38), icetray.OMKey(7,39), icetray.OMKey(7,40),
                   icetray.OMKey(7,41), icetray.OMKey(7,42), icetray.OMKey(7,44),
                   icetray.OMKey(7,45), icetray.OMKey(7,46), icetray.OMKey(8,59),
                   icetray.OMKey(8,60), icetray.OMKey(9,15), icetray.OMKey(9,16),
                   icetray.OMKey(11,1), icetray.OMKey(11,2), icetray.OMKey(18,45),
                   icetray.OMKey(18,46), icetray.OMKey(22,49), icetray.OMKey(24,56),
                   icetray.OMKey(26,46), icetray.OMKey(27,7), icetray.OMKey(27,8),
                   icetray.OMKey(28,41), icetray.OMKey(28,42), icetray.OMKey(29,59),
                   icetray.OMKey(29,60), icetray.OMKey(30,23), icetray.OMKey(30,60),
                   icetray.OMKey(32,57), icetray.OMKey(32,58), icetray.OMKey(33,6),
                   icetray.OMKey(34,11), icetray.OMKey(34,12), icetray.OMKey(34,13),
                   icetray.OMKey(34,14), icetray.OMKey(34,15), icetray.OMKey(34,17),
                   icetray.OMKey(34,22), icetray.OMKey(38,59), icetray.OMKey(39,22),
                   icetray.OMKey(39,61), icetray.OMKey(40,33), icetray.OMKey(40,34),
                   icetray.OMKey(40,51), icetray.OMKey(40,52), icetray.OMKey(42,47),
                   icetray.OMKey(42,48), icetray.OMKey(44,25), icetray.OMKey(44,26),
                   icetray.OMKey(44,47), icetray.OMKey(44,48), icetray.OMKey(45,19),
                   icetray.OMKey(47,55), icetray.OMKey(47,56), icetray.OMKey(49,55),
                   icetray.OMKey(49,56), icetray.OMKey(50,36), icetray.OMKey(53,17),
                   icetray.OMKey(53,18), icetray.OMKey(54,47), icetray.OMKey(59,45),
                   icetray.OMKey(59,46), icetray.OMKey(60,55), icetray.OMKey(66,33),
                   icetray.OMKey(66,34), icetray.OMKey(66,45), icetray.OMKey(66,46),
                   icetray.OMKey(68,42), icetray.OMKey(69,23), icetray.OMKey(69,24),
                   icetray.OMKey(69,44), icetray.OMKey(69,47), icetray.OMKey(69,48),
                   icetray.OMKey(71,39), icetray.OMKey(74,9), icetray.OMKey(74,17),
                   icetray.OMKey(74,18), icetray.OMKey(79,55), icetray.OMKey(79,56),
                   icetray.OMKey(86,27), icetray.OMKey(86,28) 
                   ]

# enable/disable calibration harvesting moni items
do_harvesting = False
