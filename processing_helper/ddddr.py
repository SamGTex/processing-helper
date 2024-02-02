from icecube import icetray
from icecube import dataclasses
from icecube import ddddr
from icecube import MuonGun

import numpy as np


class I3Vectorize(icetray.I3ConditionalModule):
    '''
    Extend a given frame object to match the length of another vector-like object.
    Written by Stefan Coenders for the 7-year PS sample.
    '''
    def __init__(self, context):
        super(I3Vectorize, self).__init__(context)
        self.AddParameter("Input",  "Name of the input object","")
        self.AddParameter("Func",   "Optional function to apply on input values",  lambda x: x)
        self.AddParameter("Vector", "Corresponding vector object to match length", "")
        self.AddParameter("Output", "Name of the output object","")
        self.AddOutBox("OutBox")

    def Configure(self):
        self.input  = self.GetParameter("Input")
        self.func   = self.GetParameter("Func")
        self.vector = self.GetParameter("Vector")
        self.output = self.GetParameter("Output")

    def Physics(self, frame):
        if (not frame.Has(self.input)) or (not frame.Has(self.vector)):
            self.PushFrame(frame)
            return
        # create output vector with the same
        vec = np.repeat(self.func(frame[self.input]), len(frame[self.vector]))
        frame[self.output] = dataclasses.I3VectorDouble(vec)
        self.PushFrame(frame)

class ddddr_labels(icetray.I3ConditionalModule):
    '''
    Calculate following one-dimensional variables for the DDDDR variables:
    - hits
    - max
    - mean
    - first hit
    - last hit
    - sum of first x% hits
    - sum of last x% hits

    Parameters
    ----------
    frame : I3Frame
        frame object
    percentile : float
        percentile to calculate the first/last x% of the track
    '''


    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("percentile", "percentile to calculate the first/last x% of the track", 0.1)

    def Configure(self):
        self.percentile = self.GetParameter("percentile")

    def Physics(self, frame):
        if frame['I3EventHeader'].sub_event_stream == "InIceSplit":
            
            dE_dX = 'OnlineL2_SplineMPE_DDDDR_150_dEdX'
            eloss = 'OnlineL2_SplineMPE_DDDDR_150_CascadeEnergyLosses'
            domQ = 'OnlineL2_SplineMPE_DDDDR_150_CascadeDomCharges'

            # drop frame
            if not frame.Has(dE_dX) and not frame.Has(eloss) and not frame.Has(domQ):
                return False
        
            # dE/dX variables
            if frame.Has(dE_dX):
                labels_dEdX = self._create_dict(frame[dE_dX], percentile=self.percentile)
                frame.Put('labels_DDDDR_150_dEdX', dataclasses.I3MapStringDouble(labels_dEdX))

            # cascade energy loss variables
            if frame.Has(eloss):
                labels_eloss = self._create_dict(frame[eloss], percentile=self.percentile)
                frame.Put('labels_DDDDR_150_CascadeEloss', dataclasses.I3MapStringDouble(labels_eloss))

            # cascade dom charge variables
            if frame.Has(domQ):
                labels_domQ = self._create_dict(frame[domQ], percentile=self.percentile)
                frame.Put('labels_DDDDR_150_CascadeDomQ', dataclasses.I3MapStringDouble(labels_domQ))
            
            self.PushFrame(frame)


    def _create_dict(self, values, percentile=0.1):
        '''
        Create a dictionary with the following keys:
        - hits
        - max
        - mean
        - first hit
        - last hit
        - first 10%
        - last 10%

        Parameters
        ----------
        values : list
            list of values
        percentile : float
            percentile to calculate the first/last 10% of the track

        Returns
        -------
        labels : dict
        '''
        
        labels = {}
        
        # number of doms hit
        labels['hits'] = len(values)

        # max energy loss
        labels['max'] = max(values)

        # mean energy loss
        labels['mean'] = sum(values)/len(values)

        # first(early)/last(late) energy lost on DOM
        labels['first'] = values[0]
        labels['last'] = values[-1]

        # energy loss on first/last 10% of track
        n_010 = int(len(values)*percentile)
        labels[f'first_{percentile}'] = sum(values[:n_010])
        labels[f'last_{percentile}'] = sum(values[-n_010:])

        return labels


@icetray.traysegment
def MuonReco(tray, name,
        pulses = 'CleanedMuonPulses',
        If = lambda f: True):

    
    splineMPE_name = 'OnlineL2_SplineMPE'

    def ModifyDDDDROutput(frame, prefix):
        # adjust the slant such that the first bin is always at zero
        slant = prefix+'Slantbinned'
        if slant in frame:
            frame[prefix+'dXbinned'] = dataclasses.I3VectorDouble(frame[slant] - np.amin(frame[slant]))
        return True


    tray.Add('I3MuonEnergy', 'GFU_SplineMPE_DDDDR_150',
            BinWidth       = 50.,
            InputPulses    = pulses,
            MaxImpact      = 150.,
            Method         = -1, # fit not required
            Seed           = splineMPE_name,
            Prefix         = splineMPE_name+'_DDDDR_150_',
            SaveDomResults = True)#,
            #If = lambda f: If(f) and EventHasBasicReco(f) and EventIsDowngoing(f) )
	
    tray.Add(ModifyDDDDROutput, 'GFU_FixDDDDR_DDDDR_150',
            prefix = splineMPE_name+'_DDDDR_150_')#,
             #If = lambda f: If(f) and EventHasBasicReco(f) and EventIsDowngoing(f) )
	
    tray.Add(I3Vectorize, 'GFU_SplineMPE_DDDDR_150_vecE',
            Input  = splineMPE_name+'_MuEx',
            Func   = lambda x: np.log10(x.energy),
            Vector = splineMPE_name+'_DDDDR_150_dEdXbinned',
            Output = splineMPE_name+'_DDDDR_150_logE')