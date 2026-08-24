import json
import correctionlib
import numpy as np
import awkward as ak
from typing import Type
from coffea.analysis_tools import Weights
from analysis.filesets.utils import get_nano_version
from analysis.corrections.utils import correction_files, unflat_sf


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
        self.cset_eff = correctionlib.CorrectionSet.from_file(
            correction_files["muon_hlt"][year]
        )

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
            
    def add_trigger_weights(self, id_wp, iso_wp, hlt_paths, dataset):
        """
        add muon iso weights to weights container
        """
        # get nominal scale factors
        nominal_weights = self.get_hlt_weights(
            id_wp=id_wp,
            iso_wp=iso_wp,
            variation="nominal",
        )
        if self.variation == "nominal":
            """
            # get 'up' and 'down' weights
            up_weights = self.get_hlt_weights(
                id_wp=id_wp,
                iso_wp=iso_wp,
                variation="systup",
            )
            down_weights = self.get_hlt_weights(
                id_wp=id_wp,
                iso_wp=iso_wp,
                variation="systdown",
            )
            """
            # add nominal, up and down weights to weights container
            self.weights.add(
                name=f"CMS_eff_m_trigger_{self.year[:4]}",
                weight=nominal_weights,
                # weightUp=up_weights,
                # weightDown=down_weights,
            )
        else:
            # add nominal weights to weights container
            self.weights.add(
                name=f"CMS_eff_m_trigger_{self.year[:4]}",
                weight=nominal_weights,
            )

    def get_id_weights(self, id_wp, variation):
        """Compute muon ID weights"""
        useHWW_sf = False
        if self.year in ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix"] and id_wp == "tight_HWW":
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
        if self.year in ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix"] and (iso_wp == "tight_HWW" or id_wp == "tight_HWW"):
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
        if self.year in ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix"] and (iso_wp == "tight_HWW" or id_wp == "tight_HWW" or promptMVA_wp == "tight_HWW"):
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

    def get_hlt_weights(self, id_wp, iso_wp, variation):
        """Compute muon HLT weights"""

        hlt_path_id_map = {
            "2016preVFP": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2016postVFP": "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2017": "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
            "2018": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2022preEE": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2022postEE": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2023preBPix": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2023postBPix": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
            "2024": "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
        }

        if not ((id_wp == "tight") & (iso_wp == "tight")):
            raise ValueError(
                f"There are no muon HLT weights for id wp '{self.id_wp}' and iso wp '{self.iso_wp}' combination"
            )

        # get muons within SF binning
        pt_lower_limit = 29.0 if self.year == "2017" else 26.0
        muon_pt_mask = (self.flat_muons.pt > pt_lower_limit) & (
            self.flat_muons.pt < 199.99
        )
        muon_eta_mask = np.abs(self.flat_muons.eta) < 2.399

        in_muons_mask = muon_pt_mask & muon_eta_mask
        in_muons = self.flat_muons.mask[in_muons_mask]

        # get muons pT and abseta (replace None values with some 'in-limit' value)
        muon_pt = ak.fill_none(in_muons.pt, 30.0)
        muon_eta = ak.fill_none(in_muons.eta, 0)
        muon_abseta = ak.fill_none(np.abs(in_muons.eta), 0)

        # for single muon events, compute SF from POG SF
        # for double muon events, compute SF from data/mc muon hlt efficiencies
        kind = "single" if ak.all(ak.num(self.muons) == 1) else "double"

        if kind == "single":
            sf = self.cset[hlt_path_id_map[self.year]].evaluate(
                muon_abseta, muon_pt, variation
            )
            nominal_sf = unflat_sf(sf, in_muons_mask, self.muons_counts)

        elif kind == "double":
            data_eff = self.cset_eff["Muon-HLT-DataEff"].evaluate(
                variation,
                hlt_path_id_map[self.year],
                muon_abseta if self.nano_version == "9" else muon_eta,
                muon_pt,
            )
            data_eff = ak.where(in_muons_mask, data_eff, ak.ones_like(data_eff))
            data_eff = ak.unflatten(data_eff, self.muons_counts)
            data_eff_leading = ak.firsts(data_eff)
            data_eff_subleading = ak.pad_none(data_eff, target=2)[:, 1]
            full_data_eff = (
                data_eff_leading
                + data_eff_subleading
                - data_eff_leading * data_eff_subleading
            )
            full_data_eff = ak.fill_none(full_data_eff, 1)

            mc_eff = self.cset_eff["Muon-HLT-McEff"].evaluate(
                variation,
                hlt_path_id_map[self.year],
                muon_abseta if self.nano_version == "9" else muon_eta,
                muon_pt,
            )
            mc_eff = ak.where(in_muons_mask, mc_eff, ak.ones_like(mc_eff))
            mc_eff = ak.unflatten(mc_eff, self.muons_counts)
            mc_eff_leading = ak.firsts(mc_eff)
            mc_eff_subleading = ak.pad_none(mc_eff, target=2)[:, 1]
            full_mc_eff = (
                mc_eff_leading + mc_eff_subleading - mc_eff_leading * mc_eff_subleading
            )
            full_mc_eff = ak.fill_none(full_mc_eff, 1)

            nominal_sf = full_data_eff / full_mc_eff

        return nominal_sf

