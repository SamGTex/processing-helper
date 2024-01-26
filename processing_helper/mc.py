from icecube import icetray
from icecube import dataclasses
from icecube import MuonGun
from icecube import dataio
from I3Tray import I3Units

from ic3_labels.labels.utils.high_level import get_muon_bundle_information
import numpy as np

class EnergyLoss(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
    
    def Physics(self, frame):
        if frame['I3EventHeader'].sub_event_stream == "InIceSplit":
            eLossInIce = 0.
            totalEmuInIce = 0.
            surface = self.MakeSurface('/data/ana/SterileNeutrino/IC86/HighEnergy/MC/scripts/jobs/files/GCD/GeoCalibDetectorStatus_2014.56784_V0.i3.gz')
            
            # make sure that corresponding MCTree is opnend
            for track in MuonGun.Track.harvest(frame['I3MCTree'], frame['MMCTrackList']):
                intersections = surface.intersection(track.pos, track.dir)
                e0, e1 = track.get_energy(intersections.first), track.get_energy(intersections.second)
                totalEmuInIce += e0
                #if e0 > 1.0 * I3Units.TeV:
                #  eMu1TeVInIce += e0
                eLossInIce += (e0-e1)
            
            frame.Put('InIceEMu', dataclasses.I3Double(totalEmuInIce))
            frame.Put('InIceELoss', dataclasses.I3Double(eLossInIce))
            
            self.PushFrame(frame)
            

    def MakeSurface(self, gcdName):
        file = dataio.I3File(gcdName, "r")
        geometry = file.pop_frame()["I3Geometry"]

        xyList = []
        zmax = -1e100
        zmin = 1e100
        step = int(len(geometry.omgeo.keys())/10)
        for i, key in enumerate(geometry.omgeo.keys()):
            #if i % step == 0:
        #    print( "{0}/{1} = {2}%".format(i,len(geometry.omgeo.keys()), int(round(i/len(geometry.omgeo.keys())*100))))

            if key.om in [61, 62, 63, 64] and key.string <= 81: #Remove IT...
                continue

            pos = geometry.omgeo[key].position
            xyList.append(pos)
            i+=1

        return MuonGun.ExtrudedPolygon(xyList, 60.*I3Units.meter)




class muon_corsika(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("energy_threshold", "energy threshold for muon bundle", 20)

    def Configure(self):
        self.energy_threshold = self.GetParameter("energy_threshold")
    
    def Physics(self, frame):
        if frame['I3EventHeader'].sub_event_stream == "InIceSplit":

            # get bundle information: dict
            labels = get_muon_bundle_information(frame, energy_threshold=self.energy_threshold)

            # check if coinc event
            status_coinc = self._check_coinc(frame)
            
            
            frame.Put('muon_labels', dataclasses.I3MapStringDouble(labels))


    # check if coincident event: False if mcBackgroundTree is empty
    def _check_coinc(self, frame):
        for parent in frame('BackgroundI3MCTree'):
            print(parent)
        