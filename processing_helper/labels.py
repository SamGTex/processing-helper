from icecube import icetray
from icecube import dataclasses

import numpy as np

class labels_on_cleanedpulses(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        #super(labels_on_cleanedpulses, self).__init__(context)

    #def Configure(self):
    #    self.qcut = self.GetParameter("qcut")

    def Physics(self, frame):
        if frame['I3EventHeader'].sub_event_stream == "InIceSplit":
            #print(frame['InIcePulses'])
            DomPulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'CleanedMuonPulses')
            charges = []
            ndoms = 0
            for omkey, pulses in DomPulses:
                charge_per_dom  = 0.
                ndoms += 1
                for j in pulses:
                    charge_per_dom += j.charge
                charges.append(charge_per_dom)       
            #find max charge at a DOM
            Qmax = max(charges)
            Qtot = sum(charges)

            frame.Put('Qmax_CleanedPulses', dataclasses.I3Double(Qmax))
            frame.Put('Qtot_CleanedPulses', dataclasses.I3Double(Qtot))
            frame.Put('NDoms_CleanedPulses', icetray.I3Int(ndoms))

            self.PushFrame(frame)


# -----------------------------------------------

#def findHEMuon(frame):
#     if frame['I3EventHeader'].sub_event_stream == "InIceSplit":
#         multi = len(frame['MMCTrackList']) 
#         Emu = []
#         Emu_in = []
#         for i in range(len(frame['MMCTrackList'])):
#             Emu.append(frame['MMCTrackList'][i].Ec)#.particle
#             Emu_in.append(frame['MMCTrackList'][i].Ei)    
#         #search for max E muon
#         Emax = max(Emu)
#         Emax_in = max(Emu_in)
#         index = Emu.index(Emax)
#         muonlist = np.array(Emu)
#         muons_Ec = muonlist[muonlist>0.]
#         multi_Ec = len(muons_Ec)
# 
#         # multiplicty at entry in detector
#         index_in = Emu_in.index(Emax_in)
#         muonlist_in = np.array(Emu_in)
#         muons_Ei = muonlist_in[muonlist_in>0.]
#         multi_Ei = len(muons_Ei)
#
#         # calculate total bundle energy
#         Ebundle = sum(Emu)
#         Ebundle_in = sum(Emu_in)
#                  
#         # create new frames
#         # muon at detector center
#         HEmuon = frame["MMCTrackList"][index].particle
#         newParticle = dataclasses.I3Particle(HEmuon)
#         frame['HEMuon'] = newParticle
#
#         # muon at entry to simulation volume
#         HEmuon_in = frame["MMCTrackList"][index_in].particle
#         newParticle_in = dataclasses.I3Particle(HEmuon_in)
#         frame['HEMuon_in'] = newParticle_in
#
#        
#
#         m = icetray.I3Int(multi)
#         frame['Multiplicity'] = m
#         frame['Multiplicity_1950'] = icetray.I3Int(multi_Ec) # was called _1500 before
#         energy = dataclasses.I3Double(Ebundle)
#         frame['Etot_bundle'] = energy
#         frame['HEMuon_Ec'] = dataclasses.I3Double(Emax)
#         
#         # muons at entry to sim volume
#         frame['Etot_bundle_in'] = dataclasses.I3Double(Ebundle_in)
#         frame['HEMuon_Ei'] = dataclasses.I3Double(Emax_in)
#         #print(frame['I3MCTree_preMuonProp'][0])
#
#        
#
#def CountHits(frame, shieldhits, tmin, tmax, OutputName):
#    CountHits = 0
#    if not frame.Has(shieldhits):
#        frame.Put(OutputName, icetray.I3Int(-1))
#    else:
#        Hits = frame[shieldhits]
#        for item in Hits:
#             if tmin < item.time_residual < tmax:
#                 CountHits+=1
#        frame.Put(OutputName, icetray.I3Int(CountHits))
#
#def MakeSurface(gcdName):
#    file = dataio.I3File(gcdName, "r")
#    geometry = file.pop_frame()["I3Geometry"]
#
#    xyList = []
#    zmax = -1e100
#    zmin = 1e100
#    step = int(len(geometry.omgeo.keys())/10)
#    for i, key in enumerate(geometry.omgeo.keys()):
#        #if i % step == 0:
#	#    print( "{0}/{1} = {2}%".format(i,len(geometry.omgeo.keys()), int(round(i/len(geometry.omgeo.keys())*100))))
#
#        if key.om in [61, 62, 63, 64] and key.string <= 81: #Remove IT...
#            continue
#
#        pos = geometry.omgeo[key].position
#        xyList.append(pos)
#        i+=1
#
#    return MuonGun.ExtrudedPolygon(xyList, 60.*I3Units.meter)
#
#
#def EnergyLoss(frame):
#    if frame['I3EventHeader'].sub_event_stream == "InIceSplit":
#        eLossInIce = 0.
#        totalEmuInIce = 0.
#        surface = MakeSurface('/data/ana/SterileNeutrino/IC86/HighEnergy/MC/scripts/jobs/files/GCD/GeoCalibDetectorStatus_2014.56784_V0.i3.gz')
#        
#        # make sure that corresponding MCTree is opnend
#        for track in MuonGun.Track.harvest(frame['I3MCTree'], frame['MMCTrackList']):
#            intersections = surface.intersection(track.pos, track.dir)
#            e0, e1 = track.get_energy(intersections.first), track.get_energy(intersections.second)
#            totalEmuInIce += e0
#            #if e0 > 1.0 * I3Units.TeV:
#            #  eMu1TeVInIce += e0
#            eLossInIce += (e0-e1)
#    #print(eLossInIce, totalEmuInIce)
#    frame['InIceEMu'] = dataclasses.I3Double(totalEmuInIce)
#    frame['InIceELoss'] = dataclasses.I3Double(eLossInIce)
#
#def getMultiplicity(frame):
#    if frame['I3EventHeader'].sub_event_stream == "InIceSplit":
#        surface = MakeSurface('/data/ana/SterileNeutrino/IC86/HighEnergy/MC/scripts/jobs/files/GCD/GeoCalibDetectorStatus_2014.56784_V0.i3.gz')
#        muonList = MuonGun.muons_at_surface(frame, surface)
#        Emu_entry = []
#        for i in range(len(muonList)):
#            Emu_entry.append(muonList[i].energy)
#        
#        # multiplicity in detector at entry
#        nMu = len(muonList)
#        frame['MultiplicityAtDetectorEntry'] =  icetray.I3Int(nMu)
#
#        # find highest energy muon within detector
#        Emax_entry = max(Emu_entry)
#        index_entry = Emu_entry.index(Emax_entry)
#        Ebundle_entry = sum(Emu_entry)
#
#        HEmuon_entry = muonList[index_entry]#.particle
#        frame['HEMuon_entry'] = dataclasses.I3Particle(HEmuon_entry)
#        frame['EMumax_entry'] = dataclasses.I3Double(Emax_entry)
#        frame['EMuBundle_entry'] = dataclasses.I3Double(Ebundle_entry)
#
