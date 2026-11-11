import json
from pathlib import Path
import correctionlib
import numpy as np
import awkward as ak
from typing import Type
from coffea.analysis_tools import Weights
from analysis.filesets.utils import get_nano_version
from analysis.corrections.utils import correction_files, unflat_sf
from analysis.corrections.trigger import *
from IPython.display import display, JSON

class MuonWeights:
    """
    Muon weights class

    Parameters:
    -----------
        events:
            pruned events
        weights:
            Weights container
        year:
            Year of the dataset {2016preVFP, 2016postVFP, 2017, 2018, 2022preEE, 2022postEE, 2023preBPix, 2023postBPix, 2024}
        variation:
            syst variation
        id_wp:
            ID working point {loose, medium, tight}
        iso_wp:
            Iso working point {loose, medium, tight}

    more info at: https://cms-analysis-corrections.docs.cern.ch/corrections/MUO/
    """

    def __init__(
        self,
        events: ak.Array,
        weights: Type[Weights],
        year: str,
        variation: str = "nominal",
    ) -> None:
        self.events = events
        self.muons = events.selected_muons
        self.weights = weights
        self.year = year
        self.variation = variation
        self.nano_version = get_nano_version(year)

        self.flat_muons = ak.flatten(self.muons)
        self.muons_counts = ak.num(self.muons)
        # get muon id/iso/HLT correction set
        # get muon id/iso/HLT correction set from muon POG
        self.cset = correctionlib.CorrectionSet.from_file(
            correction_files["muon"][year]
        )
        # get muon HLT eff correction set
        self.cset_eff = {}
        self.cset_eff["mc"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon_hlt"][year]
        )
        self.cset_eff["data"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon_hlt"][year]
        )
        self.cset_eff["sf"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon"][year]
        )
        self.cset_eff_hlt = correctionlib.CorrectionSet.from_file(
            correction_files["muon_hlt"][year]
        )
        self.cset_eff_data = self.cset_eff_hlt["Muon-HLT-DataEff"]
        self.cset_eff_mc = self.cset_eff_hlt["Muon-HLT-McEff"]

    def add_id_weights(self, id_wp):
        """
        add muon ID weights to weights container
        """
        nominal_weights = self.get_id_weights(id_wp, variation="nominal")
        if self.variation == "nominal":
            # get 'up' and 'down' weights
            up_weights = self.get_id_weights(id_wp, variation="systup")
            down_weights = self.get_id_weights(id_wp, variation="systdown")
            # add scale factors to weights container
            self.weights.add(
                name=f"CMS_eff_m_id_{self.year[:4]}",
                weight=nominal_weights,
                weightUp=up_weights,
                weightDown=down_weights,
            )
        else:
            self.weights.add(
                name=f"CMS_eff_m_id_{self.year[:4]}",
                weight=nominal_weights,
            )

    def add_iso_weights(self, id_wp, iso_wp):
        """
        add muon iso weights to weights container
        """
        # get nominal scale factors
        nominal_weights = self.get_iso_weights(id_wp, iso_wp, variation="nominal")
        if self.variation == "nominal":
            # get 'up' and 'down' weights
            up_weights = self.get_iso_weights(id_wp, iso_wp, variation="systup")
            down_weights = self.get_iso_weights(id_wp, iso_wp, variation="systdown")
            # add nominal, up and down weights to weights container
            self.weights.add(
                name=f"CMS_eff_m_iso_{self.year[:4]}",
                weight=nominal_weights,
                weightUp=up_weights,
                weightDown=down_weights,
            )
        else:
            # add nominal weights to weights container
            self.weights.add(
                name=f"CMS_eff_m_iso_{self.year[:4]}",
                weight=nominal_weights,
            )
            
    def add_promptMVA_weights(self, id_wp, iso_wp, promptMVA_wp):
        """
        add muon promptMVA weights to weights container
        """
        # get nominal scale factors
        nominal_weights = self.get_promptMVA_weights(id_wp, iso_wp, promptMVA_wp, variation="nominal")
        if self.variation == "nominal":
            # get 'up' and 'down' weights
            up_weights = self.get_promptMVA_weights(id_wp, iso_wp, promptMVA_wp, variation="systup")
            down_weights = self.get_promptMVA_weights(id_wp, iso_wp, promptMVA_wp, variation="systdown")
            # add nominal, up and down weights to weights container
            self.weights.add(
                name=f"CMS_eff_m_promptMVA_{self.year[:4]}",
                weight=nominal_weights,
                weightUp=up_weights,
                weightDown=down_weights,
            )
        else:
            # add nominal weights to weights container
            self.weights.add(
                name=f"CMS_eff_m_promptMVA_{self.year[:4]}",
                weight=nominal_weights,
            )

    def get_id_weights(self, id_wp, variation):
        """Compute muon ID weights"""
        useHWW_sf = False
        if self.year in ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"] and id_wp == "tight_HWW":
            # get muon id/iso/HLT correction set from HWW main analysis
            self.cset = correctionlib.CorrectionSet.from_file(
                correction_files["muon_HWW"][self.year]
            )
            useHWW_sf = True
        id_corrections = {
            "loose": "NUM_LooseID_DEN_TrackerMuons",
            "medium": "NUM_MediumID_DEN_TrackerMuons",
            "tight": "NUM_TightID_DEN_TrackerMuons",
            "tight_HWW": "NUM_TightID_HWW_DEN_TrackerMuons",
        }
        # get muons that pass the id wp, and within SF binning
        if useHWW_sf:
            pt_upper_limit = 1000.0 
            pt_lower_limit = 10.0
        else:
            pt_upper_limit = 500.0 
            pt_lower_limit = 10.0 if self.nano_version == "15" else 15.0
            
        muon_pt_mask = (self.flat_muons.pt > pt_lower_limit) & (
            self.flat_muons.pt < pt_upper_limit
        )
        muon_eta_mask = np.abs(self.flat_muons.eta) < 2.399
        in_muon_mask = muon_pt_mask & muon_eta_mask
        in_muons = self.flat_muons.mask[in_muon_mask]

        # get muons pT and abseta (replace None values with some 'in-limit' value)
        muon_pt = ak.fill_none(in_muons.pt, pt_lower_limit)
        muon_eta = np.abs(ak.fill_none(in_muons.eta, 0.0))

        sfs = self.cset[id_corrections[id_wp]].evaluate(muon_eta, muon_pt, variation)
        weights = unflat_sf(
            sfs,
            in_muon_mask,
            self.muons_counts,
        )
        # set back to default: get muon id/iso/HLT correction set from muon POG
        self.cset = correctionlib.CorrectionSet.from_file(
            correction_files["muon"][self.year]
        )
        return weights

    def get_iso_weights(self, id_wp, iso_wp, variation):
        """Compute muon iso weights"""
        useHWW_sf = False
        if self.year in ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"] and (iso_wp == "tight_HWW" or id_wp == "tight_HWW"):
            # get muon id/iso/HLT correction set from HWW main analysis
            self.cset = correctionlib.CorrectionSet.from_file(
                correction_files["muon_HWW"][self.year]
            )
            useHWW_sf = True
        iso_corrections = {
            "9": {
                "loose": {
                    "loose": "NUM_LooseRelIso_DEN_LooseID",
                    "medium": None,
                    "tight": None,
                },
                "medium": {
                    "loose": "NUM_LooseRelIso_DEN_MediumID",
                    "medium": None,
                    "tight": "NUM_TightRelIso_DEN_MediumID",
                },
                "tight": {
                    "loose": "NUM_LooseRelIso_DEN_TightIDandIPCut",
                    "medium": None,
                    "tight": "NUM_TightRelIso_DEN_TightIDandIPCut",
                },
            },
            "12": {
                "loose": {
                    "loose": "NUM_LoosePFIso_DEN_LooseID",
                    "medium": "NUM_LoosePFIso_DEN_MediumID",
                    "tight": "NUM_LoosePFIso_DEN_TightID",
                },
                "medium": {
                    "loose": None,
                    "medium": None,
                    "tight": None,
                },
                "tight": {
                    "loose": None,
                    "medium": "NUM_TightPFIso_DEN_MediumID",
                    "tight": "NUM_TightPFIso_DEN_TightID",
                },
                "tight_HWW": {
                    "loose": None,
                    "medium": "NUM_TightPFIso_DEN_MediumID",
                    "tight": "NUM_TightPFIso_DEN_TightID",
                    "tight_HWW": "NUM_TightPFIso_DEN_TightID_HWW",
                },
            },
            "15": {
                "loose": {
                    "loose": "NUM_LoosePFIso_DEN_LooseID",
                    "medium": "NUM_LoosePFIso_DEN_MediumID",
                    "tight": "NUM_LoosePFIso_DEN_TightID",
                },
                "medium": {
                    "loose": None,
                    "medium": None,
                    "tight": None,
                },
                "tight": {
                    "loose": None,
                    "medium": "NUM_TightPFIso_DEN_MediumID",
                    "tight": "NUM_TightPFIso_DEN_TightID",
                },
                "tight_HWW": {
                    "loose": None,
                    "medium": "NUM_TightPFIso_DEN_MediumID",
                    "tight": "NUM_TightPFIso_DEN_TightID",
                    "tight_HWW": "NUM_TightPFIso_DEN_TightID_HWW",
                },
            },
        }
        correction_name = iso_corrections[self.nano_version][iso_wp][id_wp]
        if correction_name is None:
            raise ValueError(
                f"There are no muon ISO weights for id wp '{id_wp}' and iso wp '{iso_wp}' combination"
            )
        # get 'in-limits' muons
        # get muons that pass the id wp, and within SF binning
        if useHWW_sf:
            pt_upper_limit = 1000.0 
            pt_lower_limit = 10.0
        else:
            pt_upper_limit = 500.0 
            pt_lower_limit = 10.0 if self.nano_version == "15" else 15.0
        muon_pt_mask = (self.flat_muons.pt > pt_lower_limit) & (
            self.flat_muons.pt < pt_upper_limit
        )
        muon_eta_mask = np.abs(self.flat_muons.eta) < 2.399
        in_muon_mask = muon_pt_mask & muon_eta_mask
        in_muons = self.flat_muons.mask[in_muon_mask]

        # get muons pT and abseta (replace None values with some 'in-limit' value)
        muon_pt = ak.fill_none(in_muons.pt, pt_lower_limit)
        muon_eta = ak.fill_none(in_muons.eta, 0.0)
        if self.nano_version == "9":
            muon_eta = np.abs(muon_eta)

        weights = unflat_sf(
            self.cset[correction_name].evaluate(
                muon_eta,
                muon_pt,
                variation,
            ),
            in_muon_mask,
            self.muons_counts,
        )
        # set back to default: get muon id/iso/HLT correction set from muon POG
        self.cset = correctionlib.CorrectionSet.from_file(
            correction_files["muon"][self.year]
        )
        return weights

    def get_promptMVA_weights(self, id_wp, iso_wp, promptMVA_wp, variation):
        """Compute muon promptMVA weights"""
        useHWW_sf = False
        if self.year in ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"] and (iso_wp == "tight_HWW" or id_wp == "tight_HWW" or promptMVA_wp == "tight_HWW"):
            # get muon id/iso/HLT correction set from HWW main analysis
            self.cset = correctionlib.CorrectionSet.from_file(
                correction_files["muon_HWW"][self.year]
            )
            useHWW_sf = True
        # zeroth level in the dictionary is the nanoAOD version
        # first level in the dictionary is the promptMVA wp
        # second level is iso wp
        # third level is ID wp
        promptMVA_corrections = {
            "9": {
                "tight_HWW": {
                    "loose": {
                        "loose": None,
                        "medium": None,
                        "tight": None,
                    },
                    "medium": {
                        "loose": None,
                        "medium": None,
                        "tight": None,
                    },
                    "tight": {
                        "loose": None,
                        "medium": None,
                        "tight": None,
                    },
                },
            },
            "12": {
                "tight_HWW": {
                    "loose": {
                        "loose": None,
                        "medium": None,
                        "tight": None,
                    },
                    "medium": {
                        "loose": None,
                        "medium": None,
                        "tight": None,
                    },
                    "tight": {
                        "loose": None,
                        "medium": None,
                        "tight": None,
                    },
                    "tight_HWW": {
                        "loose": None,
                        "medium": None,
                        "tight": None,
                        "tight_HWW": "NUM_TightID_HWW_TightIso_tthMVA_DEN_TightPFIso",
                    },
                },
            },
            "15": {
                "tight_HWW": {
                    "loose": {
                        "loose": None,
                        "medium": None,
                        "tight": None,
                    },
                    "medium": {
                        "loose": None,
                        "medium": None,
                        "tight": None,
                    },
                    "tight": {
                        "loose": None,
                        "medium": None,
                        "tight": None,
                    },
                    "tight_HWW": {
                        "loose": None,
                        "medium": None,
                        "tight": None,
                        "tight_HWW": "NUM_TightID_HWW_TightIso_tthMVA_DEN_TightPFIso",
                    },
                },
            },
        }
        correction_name = promptMVA_corrections[self.nano_version][promptMVA_wp][iso_wp][id_wp]
        if correction_name is None:
            raise ValueError(
                f"There are no muon promptMVA weights for id wp '{id_wp}', iso wp '{iso_wp}' and promptMVA '{promptMVA_wp}' combination for this nano version: {nano_version}."
            )
        # get 'in-limits' muons
                # get muons that pass the id wp, and within SF binning
        if useHWW_sf:
            pt_upper_limit = 1000.0 
            pt_lower_limit = 10.0
        else:
            pt_upper_limit = 500.0 
            pt_lower_limit = 10.0 if self.nano_version == "15" else 15.0
        muon_pt_mask = (self.flat_muons.pt > pt_lower_limit) & (
            self.flat_muons.pt < pt_upper_limit
        )
        muon_eta_mask = np.abs(self.flat_muons.eta) < 2.399
        in_muon_mask = muon_pt_mask & muon_eta_mask
        in_muons = self.flat_muons.mask[in_muon_mask]

        # get muons pT and abseta (replace None values with some 'in-limit' value)
        muon_pt = ak.fill_none(in_muons.pt, pt_lower_limit)
        muon_eta = ak.fill_none(in_muons.eta, 0.0)
        if self.nano_version == "9":
            muon_eta = np.abs(muon_eta)

        weights = unflat_sf(
            self.cset[correction_name].evaluate(
                muon_eta,
                muon_pt,
                variation,
            ),
            in_muon_mask,
            self.muons_counts,
        )
        # set back to default: get muon id/iso/HLT correction set from muon POG
        self.cset = correctionlib.CorrectionSet.from_file(
            correction_files["muon"][self.year]
        )
        return weights

