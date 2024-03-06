from . import filter_globals
import copy
from icecube import dataclasses
from icecube import icetray

@icetray.traysegment
def OnlineL2Filter(tray, name, pulses=filter_globals.CleanedMuonPulses,
                   linefit_name = filter_globals.muon_linefit,
                   llhfit_name  = filter_globals.muon_llhfit,
                   SplineRecoAmplitudeTable=None, SplineRecoTimingTable=None,
                   PathToCramerRaoTable=None, forceOnlineL2BadDOMList=None,
                   If = lambda f: True):

    from math import cos, radians
    import os.path

    from icecube.icetray import I3Units
    from icecube import filterscripts
    from icecube import linefit, gulliver, gulliver_modules, lilliput, spline_reco
    from icecube import bayesian_priors
    from icecube import cramer_rao
    #from icecube import double_muon
    from icecube import mue, truncated_energy
    from icecube.common_variables import direct_hits, hit_multiplicity, hit_statistics, track_characteristics
    import icecube.lilliput.segments

    # set path to Spline tables
    if (SplineRecoAmplitudeTable is None) or (SplineRecoTimingTable is None):
        domain = filter_globals.getDomain()
        if (domain == filter_globals.Domain.SPS) or (domain == filter_globals.Domain.SPTS):
            SplineDir                = '/opt/cvmfs/current/data/photon-tables/splines/'
            SplineRecoAmplitudeTable = os.path.join(SplineDir, 'InfBareMu_mie_abs_z20a10.fits')
            SplineRecoTimingTable    = os.path.join(SplineDir, 'InfBareMu_mie_prob_z20a10.fits')
        elif domain == filter_globals.Domain.ACCESS:
            SplineDir                = '/opt/cvmfs/current/data/photon-tables/splines/'
            SplineRecoAmplitudeTable = os.path.join(SplineDir, 'InfBareMu_mie_abs_z20a10.fits')
            SplineRecoTimingTable    = os.path.join(SplineDir, 'InfBareMu_mie_prob_z20a10.fits')
        else:
            SplineDir                = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/"
            SplineRecoAmplitudeTable = os.path.join(SplineDir, 'InfBareMu_mie_abs_z20a10_V2.fits')
            SplineRecoTimingTable    = os.path.join(SplineDir, 'InfBareMu_mie_prob_z20a10_V2.fits')

    # set path to Cramer-Rao tables
    if (PathToCramerRaoTable is None):
        domain = filter_globals.getDomain()
        if (domain == filter_globals.Domain.SPS) or (domain == filter_globals.Domain.SPTS):
            PathToCramerRaoTable = '/usr/local/pnf/jeb/cramer-rao/resources/data/'
        else:
            PathToCramerRaoTable = os.path.expandvars("$I3_BUILD/cramer-rao/resources/data/")


    def DoBasicReco(frame):
        return True
        '''
        Which events should get the basic OnlineL2 recos? (e.g. SPE2it, MPE)
        '''
        if frame.Has('FilterMask'):
            _filter = frame['FilterMask'].get(filter_globals.HighQFilter)
            if _filter.condition_passed:
                return True
            else:
                return False
        
        # filter not need because qfilter is already applied
    
        #passed_Muon = False
        #_filter = frame['FilterMask'].get(filter_globals.MuonFilter)
        #if _filter.condition_passed:
        #    passed_Muon = True
        #    print('MuonFilter passed')

        #return passed_Muon


    def DoAdvancedReco(frame):
        return True
        # same as basic reco
        if frame.Has('FilterMask'):
            _filter = frame['FilterMask'].get(filter_globals.HighQFilter)
            if _filter.condition_passed:
                return True
            else:
                return False


    def FindBestTrack(frame, tracks, output):
        '''
        Store first succesful fit from a given list of tracks under a new name.
        '''
        for track in tracks:
            if frame[track].fit_status == dataclasses.I3Particle.OK:
                frame[output] = frame[track]
                frame[output+'_Name'] = dataclasses.I3String(track)
                if frame.Has(track+'FitParams'):
                    frame[output+'FitParams'] = frame[track+'FitParams']
                return True
        # when we end up here, all fits failed
        return True


    def EventHasGoodTrack(frame, fitname=name+'_BestFit'):
        '''
        Test whether a track fit has has converged.
        '''
        if frame.Has(fitname):
            return (frame[fitname].fit_status == 0)

        return False

    def SelectOnlyIceCubePulses(frame, pulses):
        '''
        Create a masked pulsemap which contains only IceCube DOMs.
        '''
        mask = dataclasses.I3RecoPulseSeriesMapMask(frame, pulses,
               lambda om, idx, pulse: om.string <= 78)
        frame[pulses+'IC'] = mask
        return True


    # direct hits definitions (default IceCube direct hits definitions: A, B, C, D)
    dh_defs = direct_hits.get_default_definitions()
    # append 'non-standard' time window E
    if 'E' not in [ dh_def.name for dh_def in dh_defs ]:
        dh_defs.append(direct_hits.I3DirectHitsDefinition("E", -15., 250.))



    #########################################################
    # Begin calculating input variables for OnlineL2 filter #
    #########################################################


    ###############################
    # Perform SPE2it and MPE fits #
    ###############################

    # iterative SPE reco; seeded with single SPE LLH fit
    tray.AddSegment(lilliput.segments.I3IterativePandelFitter, name+'_SPE2itFit',
            pulses       = pulses,
            n_iterations = 2,
            seeds        = [ llhfit_name ],
            domllh       = 'SPE1st',
            tstype       = 'TFirst',
            fitname      = name+'_SPE2itFit',
            If = lambda f: If(f) and DoBasicReco(f) )   

    # seed MPE with SPE2it
    tray.AddSegment(lilliput.segments.I3SinglePandelFitter, name+'_MPEFit',
            pulses       = pulses,
            domllh       = 'MPE',
            seeds        = [ name+'_SPE2itFit' ],
            tstype       = 'TNone', # OnlineL2_16: TNone, OfflineL2: unset (=TFirst)
            fitname      = name+'_MPEFit',
            If = lambda f: If(f) and DoBasicReco(f) )

    # Find best track so far
    tray.Add(FindBestTrack,
             tracks = [ name+'_MPEFit', name+'_SPE2itFit', llhfit_name, linefit_name ],
             output = name+'_BestFit',
             If = lambda f: If(f) and DoBasicReco(f) )

    # Calculate track parameters for best track
    tray.AddSegment(direct_hits.I3DirectHitsCalculatorSegment, name+'_BestFit_DirectHits',
        DirectHitsDefinitionSeries       = dh_defs,
        PulseSeriesMapName               = pulses,
        ParticleName                     = name+'_BestFit',
        OutputI3DirectHitsValuesBaseName = name+'_BestFit_DirectHits',
        BookIt                           = True,
        If = lambda f: If(f) and DoBasicReco(f) )


    # Calculate basic event parameters with all pulses
    tray.AddSegment(hit_multiplicity.I3HitMultiplicityCalculatorSegment, name+'_HitMultiplicityValues',
            PulseSeriesMapName                = pulses,
            OutputI3HitMultiplicityValuesName = name+'_HitMultiplicityValues',
            BookIt                            = True,
            If = lambda f: If(f) and DoBasicReco(f) )

    tray.AddSegment(hit_statistics.I3HitStatisticsCalculatorSegment, name+'_HitStatisticsValues',
            PulseSeriesMapName              = pulses,
            OutputI3HitStatisticsValuesName = name+'_HitStatisticsValues',
            BookIt                          = True,
            If = lambda f: If(f) and DoBasicReco(f) )


    #######################
    # Apply OnlineL2 cuts #
    #######################

    #tray.AddModule("I3FilterModule<I3OnlineL2Filter_13>", name+"_OnlineL2Filter2013",
    #    # input variable names:
    #    PriParticleKey        = name+'_BestFit',
    #    DirectHitValues       = name+'_BestFit_DirectHitsC',
    #    HitMultiplicityValues = name+'_HitMultiplicityValues',
    #    HitStatisticsValues   = name+'_HitStatisticsValues',
    #    LLHFitParamsKey       = llhfit_name+'FitParams',
    #    # cut zones:
    #    CosZenithZone1 = [-1., cos(radians(82.))],
    #    CosZenithZone2 = [cos(radians(82.)), cos(radians(66.))],
    #    CosZenithZone3 = [cos(radians(66.)), 1.],
    #    # cut values:
    #    LDirCZone1 = 160.,
    #    NDirCZone1 = 9,
    #    PLogLParamZone1 = 4.5,
    #    PLogLZone1 = 8.3,
    #    QTotZone1 = 2.7,
    #    QTotSlopeZone2 = 3.3,
    #    QTotInterceptZone2 = 1.05 + 0.08,
    #    PLogLParamZone2 = 4.5,
    #    PLogLZone2 = 8.3,
    #    QTotOrCutZone2 = 2.5,
    #    QTotKinkZone3 = 0.5,
    #    QTotSlopeZone3 = 0.6,
    #    QTotInterceptZone3 = 2.4 + 0.13,
    #    PLogLParamZone3 = 4.5,
    #    PLogLZone3 = 10.5,
    #    QTotOrCutZone3 = 3.,
    #    DiscardEvents = False,
    #    DecisionName = filter_globals.OnlineL2Filter,
    #    If = lambda f: If(f) and DoCut(f)
    #    )


    ##########################
    # Advanced reco's follow #
    ##########################

    #################################
    # Cramer-Rao for BestFit so far #
    #################################

    tray.AddModule("CramerRao", name+'_BestFit_CramerRao',
            InputResponse = pulses,
            InputTrack    = name+'_BestFit',
            OutputResult  = name+'_BestFit_CramerRao',
            AllHits       = True,
            DoubleOutput  = True,
            PathToTable   = PathToCramerRaoTable,
            z_dependent_scatter = True,
            If = lambda f: If(f) and DoAdvancedReco(f) )

    
    ############################################
    # Calculate some variables only on IC DOMs #
    ############################################

    tray.Add(SelectOnlyIceCubePulses, name+'_SelectICPulses', pulses=pulses,
             If = lambda f: If(f) and DoAdvancedReco(f))

    tray.AddSegment(hit_multiplicity.I3HitMultiplicityCalculatorSegment, name+'_HitMultiplicityValuesIC',
            PulseSeriesMapName                = pulses+'IC',
            OutputI3HitMultiplicityValuesName = name+'_HitMultiplicityValuesIC',
            BookIt                            = True,
            If = lambda f: If(f) and DoAdvancedReco(f)
            )

    tray.AddSegment(hit_statistics.I3HitStatisticsCalculatorSegment, name+'_HitStatisticsValuesIC',
            PulseSeriesMapName              = pulses+'IC',
            OutputI3HitStatisticsValuesName = name+'_HitStatisticsValuesIC',
            BookIt                          = True,
            If = lambda f: If(f) and DoAdvancedReco(f)
            )

    ##############################################################
    # To improve SplineMPE, first run a MuEx reco on the MPE fit #
    ##############################################################

    tray.AddModule("muex", name+'_BestFit_MuEx',
            pulses = pulses,
            rectrk = name+'_BestFit',
            result = name+'_BestFit_MuEx',
            energy = True,
            detail = True,
            compat = False,
            lcspan = 0,
            If = lambda f: If(f) and DoAdvancedReco(f) )

    #################
    # Run SplineMPE #
    #################

    tray.AddSegment(spline_reco.SplineMPE, name+'_SplineMPE',
                    fitname               = name+'_SplineMPE',
                    configuration         = 'fast',
                    PulsesName            = pulses,
                    TrackSeedList         = [ name+'_BestFit' ],
                    BareMuTimingSpline    = SplineRecoTimingTable,
                    BareMuAmplitudeSpline = SplineRecoAmplitudeTable,
                    EnergyEstimators      = [ name+'_BestFit_MuEx' ],
                    If = lambda f: If(f) and DoAdvancedReco(f) )

    tray.AddModule("CramerRao", name+'_SplineMPE_CramerRao',
            InputResponse = pulses,
            InputTrack    = name+'_SplineMPE',
            OutputResult  = name+'_SplineMPE_CramerRao',
            AllHits       = True,
            DoubleOutput  = True,
            PathToTable   = PathToCramerRaoTable,
            z_dependent_scatter = True,
            If = lambda f: If(f) and DoAdvancedReco(f) )

    # calculate track properties with all DOMs, and with only IceCube DOMs
    for suffix in [ '', 'IC' ]:
        tray.AddSegment(direct_hits.I3DirectHitsCalculatorSegment, name+'_SplineMPE_DirectHits'+suffix,
                DirectHitsDefinitionSeries       = dh_defs,
                PulseSeriesMapName               = pulses+suffix,
                ParticleName                     = name+'_SplineMPE',
                OutputI3DirectHitsValuesBaseName = name+'_SplineMPE_DirectHits'+suffix,
                BookIt                           = True,
                If = lambda f: If(f) and DoAdvancedReco(f) )

        tray.AddSegment(track_characteristics.I3TrackCharacteristicsCalculatorSegment, name+'_SplineMPE_Characteristics'+suffix,
                PulseSeriesMapName                     = pulses+suffix,
                ParticleName                           = name+'_SplineMPE',
                OutputI3TrackCharacteristicsValuesName = name+'_SplineMPE_Characteristics'+suffix,
                TrackCylinderRadius                    = 150.*I3Units.m,
                BookIt                                 = True,
                If = lambda f: If(f) and DoAdvancedReco(f) )

        tray.AddSegment(track_characteristics.I3TrackCharacteristicsCalculatorSegment, name+'_SplineMPE_CharacteristicsNoRCut'+suffix,
                PulseSeriesMapName                     = pulses+suffix,
                ParticleName                           = name+'_SplineMPE',
                OutputI3TrackCharacteristicsValuesName = name+'_SplineMPE_CharacteristicsNoRCut'+suffix,
                TrackCylinderRadius                    = 4000.*I3Units.m,
                BookIt                                 = True,
                If = lambda f: If(f) and DoAdvancedReco(f) )


    ############################################
    # Run energy estimators on SplineMPE track #
    ############################################
    def PutBadDomListInFrame(frame):
        if forceOnlineL2BadDOMList is None:
            frame[name+'_BadDomList'] = dataclasses.I3VectorOMKey(filter_globals.online_bad_doms)
        else:
            # If requested, store a copy of the forced onlineL2 bad DOM list. This is used for pass2 
            # where the hard-coded bad DOM list is not applicable to old data. We made the choice 
            # to re-use the L2 bad DOM list we have access to in pass2 for this instead. 
            frame[name+'_BadDomList'] = copy.copy(frame[forceOnlineL2BadDOMList])
    tray.AddModule(PutBadDomListInFrame, name+"_BadDOMsToFrame",)

    # find spline table service name
    splineServiceName = None
    for key in tray.context.keys():
        if key.startswith('BareMuSpline'):
            assert splineServiceName is None
            splineServiceName = key

    tray.AddModule("I3TruncatedEnergy", name+'_SplineMPE_TruncatedEnergy',
            RecoPulsesName =          pulses,
            RecoParticleName =        name+'_SplineMPE',
            ResultParticleName =      name+'_SplineMPE_TruncatedEnergy',
            I3PhotonicsServiceName =  splineServiceName,
            UseRDE =                  True, # default - correct for HQE DOMs
            BadDomListName =          name+'_BadDomList',
            If = lambda f: If(f) and DoAdvancedReco(f) )

    tray.AddModule("I3mue", name+'_SplineMPE_MuE',
            RecoPulseSeriesNames = [ pulses ],
            RecoResult           = name+'_SplineMPE',
            RecoIntr             = 1,
            OutputParticle       = name+'_SplineMPE_MuE',
            If = lambda f: If(f) and DoAdvancedReco(f) )

    tray.AddModule("muex", name+'_SplineMPE_MuEx',
            pulses = pulses,
            rectrk = name+'_SplineMPE',
            result = name+'_SplineMPE_MuEx',
            energy = True,
            detail = True,
            compat = False,
            lcspan = 0,
            If = lambda f: If(f) and DoAdvancedReco(f) )


    ##################################
    # Variables to keep in the frame #
    ##################################

    filter_globals.onlinel2filter_keeps += [
        name+'_SPE2itFit',
        name+'_SPE2itFitFitParams',
        name+'_MPEFit',
        name+'_MPEFitFitParams',
        name+'_BestFit',
        name+'_BestFitFitParams',
        name+'_BestFit_Name',
        name+'_BestFit_DirectHitsA',
        name+'_BestFit_DirectHitsB',
        name+'_BestFit_DirectHitsC',
        name+'_BestFit_DirectHitsD',
        name+'_BestFit_DirectHitsE',
        name+'_HitMultiplicityValues',
        name+'_HitStatisticsValues',

        name+'_'+pulses,

        name+'_BestFit_MuEx',
        name+'_BestFit_MuEx_r',
        name+'_BestFit_CramerRao',
        name+'_BestFit_CramerRaoParams',
        name+'_BestFit_CramerRao_cr_azimuth',
        name+'_BestFit_CramerRao_cr_zenith',

        name+'_HitMultiplicityValuesIC',
        name+'_HitStatisticsValuesIC',

        name+'_SplineMPE',
        name+'_SplineMPEFitParams',
        name+'_SplineMPE_CramerRao',
        name+'_SplineMPE_CramerRaoParams',
        name+'_SplineMPE_CramerRao_cr_azimuth',
        name+'_SplineMPE_CramerRao_cr_zenith',
        name+'_SplineMPE_DirectHitsA',
        name+'_SplineMPE_DirectHitsB',
        name+'_SplineMPE_DirectHitsC',
        name+'_SplineMPE_DirectHitsD',
        name+'_SplineMPE_DirectHitsE',
        name+'_SplineMPE_Characteristics',
        name+'_SplineMPE_CharacteristicsNoRCut',
        name+'_SplineMPE_DirectHitsICA',
        name+'_SplineMPE_DirectHitsICB',
        name+'_SplineMPE_DirectHitsICC',
        name+'_SplineMPE_DirectHitsICD',
        name+'_SplineMPE_DirectHitsICE',
        name+'_SplineMPE_CharacteristicsIC',
        name+'_SplineMPE_CharacteristicsNoRCutIC',
        name+'_SplineMPE_TruncatedEnergy_AllDOMS_Muon',
        name+'_SplineMPE_TruncatedEnergy_AllDOMS_MuEres',
        name+'_SplineMPE_TruncatedEnergy_AllDOMS_Neutrino',
        name+'_SplineMPE_TruncatedEnergy_DOMS_Muon',
        name+'_SplineMPE_TruncatedEnergy_DOMS_MuEres',
        name+'_SplineMPE_TruncatedEnergy_DOMS_Neutrino',
        name+'_SplineMPE_TruncatedEnergy_AllBINS_Muon',
        name+'_SplineMPE_TruncatedEnergy_AllBINS_MuEres',
        name+'_SplineMPE_TruncatedEnergy_AllBINS_Neutrino',
        name+'_SplineMPE_TruncatedEnergy_BINS_Muon',
        name+'_SplineMPE_TruncatedEnergy_BINS_MuEres',
        name+'_SplineMPE_TruncatedEnergy_BINS_Neutrino',
        name+'_SplineMPE_TruncatedEnergy_ORIG_Muon',
        name+'_SplineMPE_TruncatedEnergy_ORIG_Neutrino',
        name+'_SplineMPE_TruncatedEnergy_ORIG_dEdX',
        name+'_SplineMPE_MuE',
        name+'_SplineMPE_MuEx',
        name+'_SplineMPE_MuEx_r',
    ]

    # remove leftover reconstructions for events that did not pass OnlineL2
    #tray.Add('Delete', name+'_Delete_Leftovers',
    #        Keys = filter_globals.onlinel2filter_keeps,
    #        If = lambda f: If(f) and (not DoAdvancedReco(f)) )