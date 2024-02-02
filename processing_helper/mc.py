from icecube import icetray
from icecube import dataclasses
from icecube import MuonGun
from icecube import dataio
from icecube import phys_services
from I3Tray import I3Units
from ic3_labels.labels.base_module import MCLabelsBase

#from ic3_labels.labels.utils.high_level import get_muon_bundle_information
from ic3_labels.labels.utils import muon as mu_utils
from ic3_labels.labels.utils import detector

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

class MergeTrees(icetray.I3ConditionalModule):#icetray.I3ConditionalModule
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("OutputKey", "Save labels to this frame key.", 'CoincLabel')
        self.AddParameter("Tree1", "Name of the I3MCTree.", 'I3MCTree_copy')
        self.AddParameter("Tree2", "Name of the BackgroundI3MCTree.", 'BackgroundI3MCTree')
        self.AddParameter("OutputTree", "Name of the output tree.", 'I3MCTree')

    def Configure(self):
        self._output_key = self.GetParameter("OutputKey")
        self._tree1 = self.GetParameter("Tree1")
        self._tree2 = self.GetParameter("Tree2")
        self._output_tree_name = self.GetParameter("OutputTree")

    def Physics(self, frame):
        # create dictonary
        labels = {}

        # check if coinc event
        labels['coinc'] = self._check_coinc(frame)

        # merge trees if coinc
        if labels['coinc']:
            frame.Put(self._output_tree_name, self._merge_trees(frame))
        else:
            frame.Put(self._output_tree_name, frame[self._tree1])
        
        frame.Put(self._output_key, dataclasses.I3MapStringDouble(labels))
        self.PushFrame(frame)

    # check if coincident event: False if mcBackgroundTree is empty
    def _check_coinc(self, frame):
        # if particles in background tree exist -> coincident event
        if frame[self._tree2].size() > 0:
            print('Coincident event found.')
            return 1
        else:
            return 0

    def _merge_trees(self, frame):
        # merge trees
        print('Merging trees.')
        _tree = frame[self._tree1]
        _tree.merge(frame[self._tree2])
        #MuonGun.MergeMCTree(frame['I3MCTree'], frame[self._background_tree_name])
        return _tree