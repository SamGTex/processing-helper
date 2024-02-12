from icecube import icetray
from icecube import dataclasses

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
