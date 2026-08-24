import numpy as np
import awkward as ak
from analysis.filesets.utils import get_nano_version
from analysis.working_points.utils import get_ctag_mask

    
class WorkingPoints:

    # -----------------------------------------------------------------
    # Electrons
    # -----------------------------------------------------------------
    def electron_id(self, events, wp):
        wps = {
            "wp80iso": (
                events.Electron.mvaFall17V2Iso_WP80
                if hasattr(events.Electron, "mvaFall17V2Iso_WP80")
                else events.Electron.mvaIso_WP80
            ),
            "wp90iso": (
                events.Electron.mvaFall17V2Iso_WP90
                if hasattr(events.Electron, "mvaFall17V2Iso_WP90")
                else events.Electron.mvaIso_WP90
            ),
            "Fail": events.Electron.cutBased >= 0,
            "Veto": events.Electron.cutBased >= 1,
            "Loose": events.Electron.cutBased >= 2,
            "Medium": events.Electron.cutBased >= 3,
            "Tight": events.Electron.cutBased >= 4,
        }

        # add electron BDT ID wp
        if hasattr(events.Electron, "mvaHZZIso"):
            bdt_id = (
                (
                    (np.abs(events.Electron.eta) < 0.8)
                    & (
                        (
                            (events.Electron.pt > 5)
                            & (events.Electron.pt < 10)
                            & (events.Electron.mvaHZZIso > 0.9044)
                        )
                        | (
                            (events.Electron.pt > 10)
                            & (events.Electron.mvaHZZIso > 0.1969)
                        )
                    )
                )
                | (
                    (np.abs(events.Electron.eta) > 0.8)
                    & (np.abs(events.Electron.eta) < 1.479)
                    & (
                        (
                            (events.Electron.pt > 5)
                            & (events.Electron.pt < 10)
                            & (events.Electron.mvaHZZIso > 0.9094)
                        )
                        | (
                            (events.Electron.pt > 10)
                            & (events.Electron.mvaHZZIso > 0.0759)
                        )
                    )
                )
                | (
                    (np.abs(events.Electron.eta) > 1.479)
                    & (
                        (
                            (events.Electron.pt > 5)
                            & (events.Electron.pt < 10)
                            & (events.Electron.mvaHZZIso > 0.9444)
                        )
                        | (
                            (events.Electron.pt > 10)
                            & (events.Electron.mvaHZZIso > -0.5169)
                        )
                    )
                )
            )
            wps["bdt"] = bdt_id

        if wp not in wps:
            raise ValueError(
                f"Invalid value {wp} for electron ID working point. Please specify {list(wps.keys())}"
            )

        return wps[wp]

    def electron_promptMVA(self, events, wp):
        wps = {
            "ttHMVA_Run3": (
                (events.Electron.mvaTTH > 0.5)
                & (
                    ((np.abs(events.Electron.eta) <= 1.479) & (np.abs(events.Electron.dxy) < 0.05))
                    |
                    ((np.abs(events.Electron.eta) > 1.479) & (np.abs(events.Electron.dxy) < 0.1))
                )
                & (
                    ((np.abs(events.Electron.eta) <= 1.479) & (np.abs(events.Electron.dz) < 0.1))
                    |
                    ((np.abs(events.Electron.eta) > 1.479) & (np.abs(events.Electron.dz) < 0.2))
                )
                & (events.Electron.convVeto == True)
            ),
            "ttHMVA_HWW": (
                (
                    ((events.Electron.pt <= 20) & (events.Electron.mvaTTH > 0.35))
                    |
                    ((events.Electron.pt > 20) & (events.Electron.mvaTTH > 0.90))
                )
                & (
                    ((np.abs(events.Electron.eta) <= 1.479) & (np.abs(events.Electron.dxy) < 0.05))
                    |
                    ((np.abs(events.Electron.eta) > 1.479) & (np.abs(events.Electron.dxy) < 0.1))
                )
                & (
                    ((np.abs(events.Electron.eta) <= 1.479) & (np.abs(events.Electron.dz) < 0.1))
                    |
                    ((np.abs(events.Electron.eta) > 1.479) & (np.abs(events.Electron.dz) < 0.2))
                )
                & (events.Electron.convVeto == True)
            ),
        }
        if wp not in wps:
            raise ValueError(
                f"Invalid value {wp} for electron ID working point. Please specify {list(wps.keys())}"
            )

        return wps[wp]

    # -----------------------------------------------------------------
    # Muons
    # -----------------------------------------------------------------
    def muon_id(self, events, wp):
        wps = {
            "loose": events.Muon.looseId,
            "medium": events.Muon.mediumId,
            "tight": events.Muon.tightId,
            "tight_HWW": ( 
                (events.Muon.pt > 10) 
                & (np.abs(events.Muon.eta) < 2.4)
                & (events.Muon.tightId)
                & (np.abs(events.Muon.dz) < 0.1)
                & (
                    ((events.Muon.pt < 20) & (np.abs(events.Muon.dxy) < 0.01))
                    |
                    ((events.Muon.pt >= 20) & (np.abs(events.Muon.dxy) < 0.02))
                )
            ),
        }
        if wp not in wps:
            raise ValueError(
                f"Invalid value for muon ID working point. Please specify {list(wps.keys())}"
            )
        return wps[wp]

    def muon_iso(self, events, wp):
        wps = {
            "loose": (
                events.Muon.pfRelIso04_all < 0.25
                if hasattr(events.Muon, "pfRelIso04_all")
                else events.Muon.pfRelIso03_all < 0.25
            ),
            "medium": (
                events.Muon.pfRelIso04_all < 0.20
                if hasattr(events.Muon, "pfRelIso04_all")
                else events.Muon.pfRelIso03_all < 0.20
            ),
            "tight": (
                events.Muon.pfRelIso04_all < 0.15
                if hasattr(events.Muon, "pfRelIso04_all")
                else events.Muon.pfRelIso03_all < 0.15
            ),
            "tight_HWW": (
                events.Muon.pfRelIso04_all < 0.15
                if hasattr(events.Muon, "pfRelIso04_all")
                else events.Muon.pfRelIso03_all < 0.15
            ),
        }
        if wp not in wps:
            raise ValueError(
                f"Invalid value {wp} for muon ISO working point. Please specify {list(wps.keys())}"
            )
        return wps[wp]

    def muon_promptMVA(self, events, year, wp):
        if year in ['2016preVFP', '2016postVFP', '2017', '2018']:
            wps = {
                "tight_HWW":  events.Muon.mvaTTH > 0.67, # this value to use as a cut is not verified
            }
        elif year in ['2022preEE', '2022postEE', '2023preBPix', '2023postBPix']:
            wps = {
                "tight_HWW":  events.Muon.tthMVA > 0.67, 
            }
        else:
            # for other years in run 3 starting from 2024, using nanoAODv15 -> for other versions, check whether promptMVA is available (replaces tthMVA from Run 2)
            wps = {
                "tight_HWW":  events.Muon.promptMVA > 0.67, # this value to use as a cut is not verified
            }            
        return wps[wp]

    # -----------------------------------------------------------------
    # Jets
    # -----------------------------------------------------------------
    def jet_id(self, events, year, wp):
        # Run 3 NanoAODs have a bug in jetId
        # Implement fix from:
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13p6TeV#nanoAOD_Flags
        nano_version = get_nano_version(year)
        if nano_version == "15":
            barrel = (
                (events.Jet.neHEF < 0.99)
                & (events.Jet.neEmEF < 0.9)
                & (events.Jet.chMultiplicity + events.Jet.neMultiplicity > 1)
                & (events.Jet.chHEF > 0.01)
                & (events.Jet.chMultiplicity > 0)
            )
            t1 = (events.Jet.neHEF < 0.9) & (events.Jet.neEmEF < 0.99)
            t2 = events.Jet.neHEF < 0.99
            endcap = (events.Jet.neMultiplicity >= 2) & (events.Jet.neEmEF < 0.4)

            jetid_tight = ak.where(
                abs(events.Jet.eta) <= 2.6,
                barrel,
                ak.where(
                    (abs(events.Jet.eta) > 2.6) & (abs(events.Jet.eta) <= 2.7),
                    t1,
                    ak.where(
                        (abs(events.Jet.eta) > 2.7) & (abs(events.Jet.eta) <= 3.0),
                        t2,
                        ak.where(
                            (abs(events.Jet.eta) > 3.0),
                            endcap,
                            ak.zeros_like(events.Jet.pt, dtype=bool),
                        ),
                    ),
                ),
            )
            jetid_tightlepveto = ak.where(
                np.abs(events.Jet.eta) <= 2.7,
                jetid_tight & (events.Jet.muEF < 0.8) & (events.Jet.chEmEF < 0.8),
                jetid_tight,
            )
        elif nano_version == "12":
            jetid_tight = ak.where(
                np.abs(events.Jet.eta) <= 2.7,
                (events.Jet.jetId >= 2)
                & (events.Jet.muEF < 0.8)
                & (events.Jet.chEmEF < 0.8),
                ak.where(
                    (abs(events.Jet.eta) > 2.7) & (abs(events.Jet.eta) <= 3.0),
                    (events.Jet.jetId >= 2) & (events.Jet.neHEF < 0.99),
                    ak.where(
                        (abs(events.Jet.eta) > 3.0),
                        (events.Jet.jetId & (1 << 1)) & (events.Jet.neEmEF < 0.4),
                        ak.zeros_like(events.Jet.pt, dtype=bool),
                    ),
                ),
            )
            jetid_tightlepveto = ak.where(
                np.abs(events.Jet.eta) <= 2.7,
                jetid_tight & (events.Jet.muEF < 0.8) & (events.Jet.chEmEF < 0.8),
                jetid_tight,
            )
        elif nano_version == "9":
            jetid_tight = events.Jet.jetId == 2
            jetid_tightlepveto = events.Jet.jetId == 6

        wps = {
            "tight": ak.values_astype(jetid_tight, bool),
            "tightlepveto": ak.values_astype(jetid_tightlepveto, bool),
        }
        if wp not in wps:
            raise ValueError(
                f"Invalid value {wp} for jet ID working point. Please specify {list(wps.keys())}"
            )
        return wps[wp]

    def jet_pileup_id(self, events, year, wp):
        nano_version = get_nano_version(year)
        if nano_version == "9":
            wps = {
                "2016preVFP": {
                    "loose": events.Jet.puId == 1,
                    "medium": events.Jet.puId == 3,
                    "tight": events.Jet.puId == 7,
                },
                "2016postVFP": {
                    "loose": events.Jet.puId == 1,
                    "medium": events.Jet.puId == 3,
                    "tight": events.Jet.puId == 7,
                },
                "2017": {
                    "loose": events.Jet.puId == 4,
                    "medium": events.Jet.puId == 6,
                    "tight": events.Jet.puId == 7,
                },
                "2018": {
                    "loose": events.Jet.puId == 4,
                    "medium": events.Jet.puId == 6,
                    "tight": events.Jet.puId == 7,
                },
            }
            if wp not in wps[year]:
                raise ValueError(
                    f"Invalid value {wp} for jet pileup ID working point. Please specify {list(wps[year].keys())}"
                )
            # break up selection for low and high pT jets
            # to apply jets_pileup only to jets with pT < 50 GeV
            return ak.where(
                events.Jet.pt < 50,
                wps[year][wp],
                events.Jet.pt >= 50,
            )
        else:
            return np.ones_like(events.Jet.pt, dtype=bool)

    def jet_ctagging(self, events, wp, year):
        return get_ctag_mask(events.Jet, year, wp)
