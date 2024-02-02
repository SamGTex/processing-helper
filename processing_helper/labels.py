from icecube import icetray
from icecube import dataclasses

# General imports for icecube modules: shieldhits
from icecube import shield

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


class ShieldhitsLabels(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("tmin", "minimum time", -1000)
        self.AddParameter("tmax", "maximum time", 1000)
        self.AddParameter("InputShieldHits", "Name of ShieldHits variable", 'shieldhits')
        self.AddParameter("OutputName", "Name of the output", 'ShieldHits')

    def Configure(self):
        self.tmin = self.GetParameter("tmin")
        self.tmax = self.GetParameter("tmax")
        self.shieldhits = self.GetParameter("InputShieldHits")
        self.outputname = self.GetParameter("OutputName")

    def Physics(self, frame):
        if frame['I3EventHeader'].sub_event_stream == "InIceSplit":
            CountHits = 0
            if not frame.Has(self.shieldhits):
                frame.Put(self.outputname, icetray.I3Int(-1))
            else:
                Hits = frame[self.shieldhits]
                for item in Hits:
                    if self.tmin < item.time_residual < self.tmax:
                        CountHits+=1

            frame.Put(self.outputname, icetray.I3Int(CountHits))
            self.PushFrame(frame)