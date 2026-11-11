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
            
    def add_trigger_weights(self, id_wp, iso_wp, promptMVA_wp, hlt_paths, dataset):
        """
        add muon iso weights to weights container
        """
        # get nominal scale factors
        print("Getting nominal hlt weights.")
        nominal_weights = self.get_hlt_weights(
            id_wp=id_wp,
            iso_wp=iso_wp,
            promptMVA_wp=promptMVA_wp,
            hlt_paths=hlt_paths,
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

            
    def get_hlt_weights(self, id_wp, iso_wp, promptMVA_wp, hlt_paths, variation):
        """Compute muon HLT weights"""
        print("Muon working points received for applying muon efficiency scale factors:")
        if id_wp: print(f"ID: {id_wp}") 
        if iso_wp: print(f"ISO: {iso_wp}") 
        if promptMVA_wp: print(f"PromptMVA: {promptMVA_wp}") 
            
        flags_with_triggernames = tree()
        flags_with_triggernames = get_config_trigger_flags(hlt_paths) # get trigger flags listed in workflow config
        trigger_flags = load_trigger_flags(Path(__file__)) # load file having for each year the trigger flag and corresponding trigger names
        year_key = get_yearkey_for_trigger_flags(self.year)
        trigger_types = flags_with_triggernames.keys() # types of trigger: single lepton trigger, double lepton trigger ... defined in analysis/corrections/trigger.py get_config_trigger_flags()

        # get for each trigger flag in the workflow config, the proper trigger names 
        # beware that the combination of working points must be available in the correction files
        # for example, only tight ID and ISO are available to be combined with the tthMVA-Run3 promptMVA working point
        for trigger_type in trigger_types:
            
            for trigger_flag in flags_with_triggernames[trigger_type].keys():
                
                trigger_names = trigger_flags[year_key][trigger_flag]
                
                for trigger_name in trigger_names:
                    flags_with_triggernames[trigger_type][trigger_flag][trigger_name]["jsonkey"] = construct_trigger_NUMDEN_name(trigger_flag, trigger_name, self.year, id_wp, iso_wp, promptMVA_wp, hlt_paths)

                    legs = get_trigger_legs(trigger_name)
                    cset_filepaths = get_cset_path(legs, trigger_type, self.year, id_wp, iso_wp, promptMVA_wp)
                    legnumber = 1
                    eff_types = ["mc", "data", "sf"]
                    for leg in legs: 
                        flags_with_triggernames[trigger_type][trigger_flag][trigger_name][f"jsonleg{legnumber}"] = legs[legnumber-1]

                        for eff_type in eff_types:
                            print("cset_filepaths[eff_type]", cset_filepaths[eff_type])
                            flags_with_triggernames[trigger_type][trigger_flag][trigger_name][f"cset_{eff_type}_filepath_leg{legnumber}"] = cset_filepaths[eff_type][legnumber-1]
                        
                        legnumber += 1

                    # hier berekening plaatsen van de scale factors -> we weten nu alle triggers die we moeten nemen
                    # als enkel SingleMu, dan moeten we hier de berekening van voorheen doen
                    # als enkel SingleEle opgegeven in config analoog aan enkel SingleMu
                    # als zowel Single als Double -> combinatie van de trigger effs
        # display final dictionary for checks:
        display(JSON(json.loads(json.dumps(flags_with_triggernames))))

        # check what trigger types we have a trigger for in the workflow config
        apply_triggertype = tree() # store in this dict/tree if we find triggers given in the workflow config: if so, trigger SFs are applied for that type
        for trigger_type in flags_with_triggernames.keys():
            print('trigger_type: ', trigger_type)
            print("flags_with_triggernames[trigger_type].keys(): ", flags_with_triggernames[trigger_type].keys()) 
            if flags_with_triggernames[trigger_type].keys() in [[], None] or flags_with_triggernames[trigger_type] == {}:
                
                print(f"WARNING: {trigger_type} does not have any triggers listed in the workflow config file.")
                apply_triggertype[trigger_type] = False
            else:
                apply_triggertype[trigger_type] = True
                if trigger_type == "Single Lepton Trigger":
                    for trigger_flag in ["SingleMu", "SingleEle"]:
                        if trigger_flag not in flags_with_triggernames[trigger_type].keys():
                            print(f"WARNING: {trigger_flag} does not have any triggers listed in the workflow config file.")
                            apply_triggertype[trigger_flag] = False
                        else:
                            apply_triggertype[trigger_flag] = True
            print('apply_triggertype: ',  apply_triggertype[trigger_type])
        display(JSON(json.loads(json.dumps(apply_triggertype))))

        # make a list of first 4 elements of values in apply_triggertype tree
        apply_triggertype_list = list(apply_triggertype.values())[:4]
        print(apply_triggertype_list)
        match apply_triggertype_list:
            case [True, True, True, False]:
            # both SingleMu and SingleEle, but no double lepton trigger
            # insert function form trigger.py for calculation of scale factor for SingleMu and SingleEle
                print("")
                
            case [True, True, False, False]:
            # only SingleMu and no double lepton trigger
                print("Only SingleMu and no double lepton trigger provided. Applying the SingleMu trigger SF to all events with 1 or more muons.")
                trigger_type = "Single Lepton Trigger"
                trigger_flag = "SingleMu"
                trigger_name = list(flags_with_triggernames[trigger_type][trigger_flag].keys())[0]
                
                for eff_type in eff_types:
                    self.cset_eff[eff_type] = correctionlib.CorrectionSet.from_file(flags_with_triggernames[trigger_type][trigger_flag][trigger_name][f"cset_{eff_type}_filepath_leg1"])
                    jsonkey = flags_with_triggernames[trigger_type][trigger_flag][trigger_name]["jsonkey"]
                    
                if not flags_with_triggernames[trigger_type][trigger_flag][trigger_name][f"cset_sf_filepath_leg1"] == "":
                    has_sf_path = True
                else:
                    has_sf_path = False 
                    
                nominal_sf = calculate_SingleMu_SF(self.events, self.cset_eff, jsonkey, self.year, self.variation, has_sf_path)
                
                # geeft het problemen als je een empty string aan correctionSet.from_file geeft?

            case [True, False, True, False]:
            # only SingleEle and no double lepton trigger
                print("")
            case [False, True]:
            # only a double lepton trigger
                print("Only a double lepton trigger was provided. Are you sure this is correct? At least one single lepton trigger is probably expected.")
            case [True, True, True, True]:
            # calculate scale factor for combination of SingleMu, SingleEle and double lepton
                print('WARNING: Only the single and double lepton trigger scale factor application is implemented, so if you provided also triggers for more than 2 leptons, those efficiencies are not taken into account.')
            case _:
            # combining only SingleMu with a double lepton also not implemented yet, same for SingleEle -> still to be done
                print('WARNING: A combination of triggers was provided for which the trigger scale factors / efficiencies are not available or for which the calculation is not yet implemented. No Trigger Scale Factors were applied.')
                  #, True, True] | [False, False, True, True] | [False, True, False, True] | [False, False, False, True] | [False, True] | [False, True, True] | [False, False, True] | [False, True, False]:
                # the case where there is no Single Lepton Trigger: first element False
                print("aaah")

        # we actually want something with double leptons: so muon or electron. the code below is for applying the SingleMuon trigger SF for the combination of 2 muons 

        return nominal_sf

