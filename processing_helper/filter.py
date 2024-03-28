from icecube import icetray
from icecube import dataclasses

from ic3_labels.labels.utils import muon as mu_utils
from ic3_labels.labels.utils.high_level import get_muon_bundle_information
from dnn_selections.selections.atmospheric_muon_stopping.ic3.stopping_muon_label import inf_muon_track_hit_detector

class drop_p(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        #super(labels_on_cleanedpulses, self).__init__(context)
        self.AddParameter("column_name", "Drop p frame if column is not in frame.", 'CleanedMuonPulses')
    def Configure(self):
        self.column_name = self.GetParameter("column_name")

    def Physics(self, frame):
        if frame.Has(self.column_name):
            self.PushFrame(frame)
        else:
            return False        

# level 0: add bundle energy at surf + inice, same for leading energy, zenith
class filter_detectorhits(icetray.I3ConditionalModule):#icetray.I3ConditionalModule
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("MMCTrackList_name", "Name of MMC TrackList.", 'MMCTrackList')

    def Configure(self):
        self.MMC_name = self.GetParameter("MMCTrackList_name")

    def Physics(self, frame):
        hit_detector = False
        
        # loop over muons in MMC tracklist and check if they hit the detector
        for particle in frame[self.MMC_name]:
            muon = particle.particle
            if mu_utils.is_muon(muon) and inf_muon_track_hit_detector(muon):
                hit_detector = True
                break
        
        # hold frame if any muon hits the detector
        if hit_detector:
            self.PushFrame(frame)
        else:
            return False