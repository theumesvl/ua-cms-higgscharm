import numpy as np
import awkward as ak
import correctionlib
import json
import yaml
from pathlib import Path
from typing import Type
from coffea.analysis_tools import Weights
from analysis.filesets.utils import get_nano_version
from analysis.corrections.utils import correction_files, unflat_sf, tree
from IPython.display import display, JSON

class TriggerWeights:
    """
    Trigger weights class

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
        muons:
            id_wp:
                ID working point {loose, medium, tight, tight_HWW}
            iso_wp:
                Iso working point {loose, medium, tight, tight_HWW}
            promptMVA_wp
                Prompt MVA working point based on the ttH MVA {tight_HWW}, retrained for 2022-2023 (no available in NanoAOD) and promptMVA in 2024 NanoAOD
        electrons
            id_wp:
                ID working point {wp80iso, wp90iso, Medium}
            promptMVA_wp
                Prompt MVA working point based on the ttH MVA {tight_HWW}, retrained for 2022-2023 (no available in NanoAOD) and promptMVA in 2024 NanoAOD

    more info muons: https://cms-analysis-corrections.docs.cern.ch/corrections/MUO/
    more info electrons: https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3#Scale_factors_and_correction_AN1
    """

    def __init__(
        self,
        events: ak.Array,
        weights: Type[Weights],
        year: str,
        variation: str = "nominal",
    ) -> None:
        # general
        self.events = events
        self.weights = weights
        self.year = year
        self.variation = variation
        self.nano_version = get_nano_version(year)

        # muons
        self.use_muon_HWW_files = False
        self.muons = events.selected_muons
        self.flat_muons = ak.flatten(self.muons)
        self.muons_counts = ak.num(self.muons)        

        # electrons
        self.use_electron_HWW_files = False 
        self.electrons = events.selected_electrons
        self.flat_electrons = ak.flatten(self.electrons)
        self.electrons_counts = ak.num(self.electrons)
        
        # get muon id/iso/HLT correction set
        # get muon id/iso/HLT correction set from muon POG
        self.cset = tree()
        self.cset["muon"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon"][year]
        )
        # get muon HLT eff correction set - below the default setting, which is for non-HWW WPs
        self.cset_eff = tree()
        self.cset_eff_hlt = tree()
        self.cset_eff_mc = tree()
        self.cset_eff_data = tree()
        self.cset_eff["muon"]["mc"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon_hlt"][year]
        ) # json file contains only MC and Data trigger efficiencies for tight ID and tight ISO (no SFs)
        self.cset_eff["muon"]["data"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon_hlt"][year]
        ) # json file contains only MC and Data trigger efficiencies for tight ID and tight ISO (no SFs)
        self.cset_eff["muon"]["sf"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon"][year]
        ) # json file contains ID, ISO SFs and also trigger SFs (no efficiencies)

        # get electron HLT eff correction set - below the default setting, which is for Run 3 non-HWW WPs
        self.cset_eff["electron"]["mc"] = correctionlib.CorrectionSet.from_file(
            correction_files["electron_hlt"][year]
        ) # json file contains MC and Data trigger efficiencies as well as trigger SFs for several ID working points for Run 3. For Run 2, these are split in several files and this one only contains the SF
        self.cset_eff["electron"]["data"] = correctionlib.CorrectionSet.from_file(
            correction_files["electron_hlt"][year]
        ) # json file contains MC and Data trigger efficiencies as well as trigger SFs for several ID working points. For Run 2, these are split in several files and this one only contains the SF
        self.cset_eff["electron"]["sf"] = correctionlib.CorrectionSet.from_file(
            correction_files["electron_hlt"][year]
        ) # json file contains MC and Data trigger efficiencies as well as trigger SFs for several ID working points. For Run 2, these are split in several files and this one only contains the SF
        
        self.cset_eff_data["electron"] = self.cset_eff_hlt["Electron-HLT-DataEff"]
        self.cset_eff_mc["electron"] = self.cset_eff_hlt["Electron-HLT-McEff"]
        # run 3
        """
        cset = correctionlib.CorrectionSet.from_file(
            correction_files["electron_hlt"][self.year]
        )
        sf = cset["Electron-HLT-SF"].evaluate(
        data_eff = cset["Electron-HLT-DataEff"].evaluate(
        mc_eff = cset["Electron-HLT-McEff"].evaluate(
        """
        # run 2
        """
        cset = correctionlib.CorrectionSet.from_file(
                correction_files["electron_hlt"][self.year]
            )
        sf = cset[f"HLT_SF_{hlt_paths[self.year]}_MVAiso80ID"].evaluate(
        data_eff = cset_data_eff[
                f"HLT_DataEff_{hlt_paths[self.year]}_MVAiso80ID"
            ].evaluate(
         mc_eff = cset_mc_eff[
                f"HLT_MCEff_{hlt_paths[self.year]}_MVAiso80ID"
            ].evaluate(
        """
    
    def get_config_trigger_flags(self, hlt_paths):
        """
        Extracts the trigger flags for each primary dataset (Data) given in the workflow config file under hlt_paths.
        """
        flags = tree()
        flags["Single Lepton Trigger"] = tree()
        flags["Double Lepton Trigger"] = tree()
        data_samples = ["MuonEG", "Muon", "EGamma"] # data samples used in the analysis
        singleTriggerFlags = ["SingleMu", "SingleEle"]
        doubleTriggerFlags = ["DiMu", "DiEle", "MuEle", "MuEle_HWW"]
        print("Expected data samples to be used are: ", data_samples)
        print("Checking for which data samples trigger flags are provided.")
        print("Trigger flags are provided for: ", list(hlt_paths.keys()))
        for sample in data_samples:
            if sample not in hlt_paths.keys():
                print("WARNING: No trigger flag(s) provided for ", sample)
        for sample in hlt_paths.keys():
            for flag in hlt_paths[sample]:
                if flag in singleTriggerFlags:
                    flags["Single Lepton Trigger"][flag] = tree()
                elif flag in doubleTriggerFlags:
                    flags["Double Lepton Trigger"][flag] = tree()
                else:
                    print(f"WARNING: provided trigger flag {flag} not recognized for applying trigger corrections/weights.")
        return flags
        
    def load_trigger_flags(self, filepath):
        """Load trigger flag definitions from analysis/selections."""
        trigger_flags_path = (
            filepath.resolve().parent.parent
            / "selections"
            / "trigger_flags.yaml"
        )
        
        with trigger_flags_path.open("r") as f:
            return yaml.safe_load(f)
                                        
    def get_yearkey_for_trigger_flags(self):
        """
        Extracts the year key for each campaign
        """
        year_map = {
            "2016preVFP": 2016,
            "2016postVFP": 2016,
            "2017": 2017,
            "2018": 2018,
            "2022preEE": 2022,
            "2022postEE": 2022,
            "2023preBPix": 2023,
            "2023postBPix": 2023,
            "2024": 2024,
        }
        year_key = year_map[self.year]
        return year_key
    
    def get_available_SF_names_DEN(self, particle, leptonWPdict):
        """
        Returns the available scale factor names based on the nano version and the year, as well as what type of working point is being looked at. It raises a Value error if the combination of working points does not correspond to a combination for which trigger Scale Factors are available
        """
        # if no promptMVA working point was given, that is treated as an empty string.
        # dictionary order of levels: nano_version > id_wp > iso_wp > promptMVA_wp: is the trigger name available?: trigger as NUM and whole WP as DEN
        years_with_trigger_SFs  = ["2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
        run2_years = ["2016preVFP", "2016postVFP", "2017", "2018"]
        if self.year not in years_with_trigger_SFs:
            raise ValueError(
                f"There are no muon HLT weights/trigger scale factors for year '{self.year}'"
            )
        if self.year in run2_years:
            raise ValueError(
                f"HLT weights/trigger scale factors are available for year '{self.year}', but not implemented in the framework yet"
            )
        available_trigger_SFs = tree()
        id_wp = ""
        iso_wp = ""
        promptMVA_wp = ""
        
        if particle == "muon":
            id_wp = leptonWPdict["id"]
            iso_wp = leptonWPdict["iso"]
            promptMVA_wp = leptonWPdict["promptMVA"]
            available_trigger_SFs = {
                "9": {
                    "tight_HWW": {
                        "tight_HWW": {
                            "tight_HWW": None,
                            "": "CutBasedIdTight_and_PFIsoTight",
                        },
                    },
                },
                "12": {
                    "tight_HWW": {
                        "tight_HWW": {
                            "tight_HWW": "CutBasedIdTight_and_PFIsoTight",
                            "": "CutBasedIdTight_and_PFIsoTight",
                        },
                    },
                    "tight": {
                        "tight_HWW": {
                            "tight_HWW": "testnametighttightHWWonly", # OK to keep like this: TriggerEff eventually used
                            "": "CutBasedIdTight_and_PFIsoTight",
                        },
                        "tight": {
                            "tight_HWW": "testnametighttighttightHWW", # OK to keep like this: TriggerEff eventually used
                            "": "CutBasedIdTight_and_PFIsoTight",
                        },
                    },
                },
                "15": {
                    "tight_HWW": {
                        "tight_HWW": {
                            "tight_HWW": None,
                            "": "CutBasedIdTight_and_PFIsoTight",
                        },
                    },
                    "tight": {
                        "tight_HWW": {
                            "tight_HWW": None,
                            "": "CutBasedIdTight_and_PFIsoTight",
                        },
                        "tight": {
                            "tight_HWW": None,
                            "": "CutBasedIdTight_and_PFIsoTight",
                        },
                    },
                },
            }
        elif particle == "electron":
            # for iso_wp, the default empty string stays
            id_wp = leptonWPdict["id"]
            promptMVA_wp = leptonWPdict["promptMVA"]
            available_trigger_SFs = {
                "9": {
                    "wp80iso": {
                        "": {
                            "ttHMVA_Run3": None,
                            "": "MVAiso80ID",
                        },
                    },
                    "wp90iso": {
                        "": {
                            "ttHMVA_Run3": None,
                            "": "MVAiso90ID",
                        },
                    },
                },
                "12": {
                    "wp80iso": {
                        "": {
                            "ttHMVA_Run3": None,
                            "": "MVAiso80ID",
                        },
                    },
                    "wp90iso": {
                        "": {
                            "ttHMVA_Run3": None,
                            "": "MVAiso90ID",
                        },
                    },
                    "Medium": {
                        "": {
                            "ttHMVA_Run3":  "testnamemediumttHMVA", # OK to keep like this: TriggerEff eventually
                            "": "", # still to be specified, but not used in our analysis
                        },
                    },
                },
                "15": {
                    "wp80iso": {
                        "": {
                            "ttHMVA_Run3": None,
                            "": "MVAiso80ID",
                        },
                    },
                    "wp90iso": {
                        "": {
                            "ttHMVA_Run3": None,
                            "": "MVAiso90ID",
                        },
                    },
                    "Medium": {
                        "": {
                            "ttHMVA_Run3": "testnamemediumttHMVA", # OK to keep like this: TriggerEff eventually
                            "": "", # still to be specified, but not used in our analysis
                        },
                    },
                },
            }            
        else:
            raise ValueError(
                f"Unknown particle: only muon and electron accepted for determining the trigger DEN name."
            )    
        print(f"Getting {particle} trigger_SF_DEN for nano_version '{self.nano_version}', ID WP '{id_wp}', ISO WP '{iso_wp}', promptMVA WP '{promptMVA_wp}'.")
        trigger_SF_DEN = available_trigger_SFs[self.nano_version][id_wp][iso_wp][promptMVA_wp]
        print(f"trigger_SF_DEN: {trigger_SF_DEN}")
        if trigger_SF_DEN is None:
            raise ValueError(
                f"There are no muon HLT weights for this id wp '{id_wp}', iso wp '{iso_wp}', promptMVA wp '{iso_wp}' combination in year {self.year}."
            )
        else:
            return trigger_SF_DEN
    
    def construct_trigger_NUMDEN_name(self, legs, leg_particles, trigger_flag, trigger_name, triggerSFs_leptonWP_inputs, hlt_paths):
        # the name of the DEN should be the same for each of the triggers as long as they are the same type (for example: SingleMu)
        # inside the HWW trigger efficiency scale factor files, all efficiencies are listed under TriggerEff and the particular name does not even matter - the trigger name or name of the leg is in the json file name wherein the correction factor is stored
        muon_intrigger = False
        electron_intrigger = False
        trigger_muon_DEN_name = ""
        trigger_electron_DEN_name = ""
        muons_for_legs = [particle for particle in leg_particles if particle == "muon"]
        electrons_for_legs = [particle for particle in leg_particles if particle == "electron"]
        nonleptons_for_legs = [particle for particle in leg_particles if not particle in ["muon", "electron"]]
            
        if len(leg_particles) > 2:
            raise ValueError(
                f"Only implemented up to double lepton triggers. For triggers with 3 or more leptons, implementation to be added."
            )
        if len(muons_for_legs) > 0:
            muon_intrigger = True
            print("triggerSFs_leptonWP_inputs['muon']", triggerSFs_leptonWP_inputs['muon'])
            trigger_muon_DEN_name = self.get_available_SF_names_DEN("muon", triggerSFs_leptonWP_inputs["muon"]) 
        if len(electrons_for_legs) > 0:
            electron_intrigger = True
            trigger_electron_DEN_name = self.get_available_SF_names_DEN("electron", triggerSFs_leptonWP_inputs["electron"]) 
        if len(nonleptons_for_legs) > 0:
            raise ValueError(
                f"Unknown particle: only muon and electron accepted for determining the trigger NUMDEN name."
            )
        # denominator defined with muon WP first, then electron WP
        denominator_parts = [
            trigger_muon_DEN_name,
            trigger_electron_DEN_name,
        ]
        
        trigger_DEN_name = "_".join(
            part for part in denominator_parts if part
        )
            
        trigger_NUM_name = trigger_name

        muon_id_wp = triggerSFs_leptonWP_inputs["muon"]["id"]
        muon_iso_wp = triggerSFs_leptonWP_inputs["muon"]["iso"]
        muon_promptMVA_wp = triggerSFs_leptonWP_inputs["muon"]["promptMVA"]

        electron_id_wp = triggerSFs_leptonWP_inputs["electron"]["id"]
        electron_promptMVA_wp = triggerSFs_leptonWP_inputs["electron"]["promptMVA"] 
        if trigger_flag == "SingleMu":
            if muon_promptMVA_wp == "tight_HWW" or muon_id_wp == "tight_HWW":
                # tight_HWW for muon ISO is exactly the same as tight muon ISO, 
                # but for ID we have tight + other selections in the definitions in tight_HWW
                trigger_NUMDEN_name = "TriggerEff"
            else:
                trigger_NUMDEN_name = "NUM_" + trigger_NUM_name + "_DEN_" + trigger_DEN_name
        elif trigger_flag == "SingleEle":
            if electron_promptMVA_wp == "ttHMVA_Run3" or electron_promptMVA_wp == "ttHMVA_HWW":
                trigger_NUMDEN_name = "TriggerEff"
            else:
                trigger_NUM_name = legs[0]
                trigger_NUMDEN_name = "HLT_SF_" + trigger_NUM_name + "_" + trigger_DEN_name
        else:
            if muon_promptMVA_wp == "tight_HWW" or muon_id_wp == "tight_HWW" or electron_promptMVA_wp == "ttHMVA_Run3" or electron_promptMVA_wp == "ttHMVA_HWW":
                # tight_HWW for muon ISO is exactly the same as tight muon ISO, 
                # but for ID we have tight + other selections in the definitions in tight_HWW
                trigger_NUMDEN_name = "TriggerEff"
            else:
                trigger_NUMDEN_name = trigger_NUM_name ##### POSSIBLY this needs to be changed still to also include a DEN_name
        return trigger_NUMDEN_name
    
    def get_trigger_legs(self, trigger_name):
        """
        Extracts the trigger legs from the trigger name: every piece in the trigger name starting with Mu or Ele.
        It also provides what lepton the leg corresponds to: muon or electron
        """
        legs = []
        leg_particles = []
        name_parts = trigger_name.split('_')
        for part in name_parts:
            if part.startswith('Mu') or part.startswith('Ele') or part.startswith('IsoMu'):
                legs.append(part)
        for leg in legs:
            if leg.startswith('Mu') or leg.startswith('IsoMu'):
                leg_particles.append("muon")
            else:
                leg_particles.append("electron")
        return legs, leg_particles
    
    def get_cset_path(self, leg_particles, legs, trigger_type, triggerSFs_leptonWP_inputs):
        """
        Returns the correct correction set path for the trigger efficiency scale factors, based on the year and the lepton working point.
        We take as many centrally derived (POG) efficiencies and scale factors as possible, expecially for the electrons.
        For HWW working points, the muon efficiencies used as the efficiencies derived by the HWW group.
        If only SingleMu triggers are used, the scale factor is the one provided by the POG for 1 muon events (cset["muon"]).
        For 2 muon events, the efficiencies are taken from the files in (cset["muon_hlt"]) and from these, a total efficiency for mc and for data is calculated, to then calculate a scale factor.
        For electrons, this same logic applies (to be checked!!!!)
        For muon-electron events, we always have 1 muon and 1 electron: if only SingleMu and SingleEle are listed in the workflow config file, then only those are used, similarly to 2 muons or 2 electrons events
        When also double lepton triggers are used, we combine the single and double lepton trigger efficiencies into 1 efficiency for data and one for mc and use those to calculate the scale factors.
        """
        # set the correction set basic path based on the working points
        # if ID or promptMVA working point is set to an HWW working point, then HWW files are used
        # for ISO, tight and tight_HWW are exactly the same, so we take the centrally derived files for that if no HWW working points are used 
        muon_id_wp = triggerSFs_leptonWP_inputs["muon"]["id"]
        muon_iso_wp = triggerSFs_leptonWP_inputs["muon"]["iso"]
        muon_promptMVA_wp = triggerSFs_leptonWP_inputs["muon"]["promptMVA"]

        electron_id_wp = triggerSFs_leptonWP_inputs["electron"]["id"]
        electron_promptMVA_wp = triggerSFs_leptonWP_inputs["electron"]["promptMVA"]
        
        eff_types = ["mc", "data", "sf"]

        cset_files = {
            "mc": [],
            "data": [],
            "sf": [],
        }
        year_key = self.get_yearkey_for_trigger_flags() # needed in case other years have other json paths

        if muon_promptMVA_wp == "tight_HWW" or muon_id_wp == "tight_HWW":
            self.use_muon_HWW_files = True
        if electron_promptMVA_wp == "ttHMVA_Run3" or electron_promptMVA_wp == "ttHMVA_HWW":
            self.use_electron_HWW_files = True 
        
        if trigger_type == "Single Lepton Trigger":
            leg = legs[0]
            particle = leg_particles[0]
            jsonfile_paths = {
                particle: {
                    "mc": [],
                    "data": [],
                    "sf": [],
                }
            }
            if particle == "muon":
                if self.use_muon_HWW_files:
                    cset_files["mc"] = correction_files[f"{particle}_hlt_HWW_mc_eff"][self.year]
                    cset_files["data"] = correction_files[f"{particle}_hlt_HWW_data_eff"][self.year]
                    cset_files["sf"] = ""                  
                    for eff_type in eff_types:
                        if eff_type == "sf":
                            # for HWW working points, no scale factors immediately available: need to be calculated with efficiencies
                            jsonfile_name = ""
                        else:
                            jsonfile_name = f"{leg}_pt_eta_efficiency.json"
                            
                        jsonfile_paths[particle][eff_type].append(
                            f"{cset_files[eff_type]}{jsonfile_name}"
                        )
                else:
                    jsonfile_paths[particle]["mc"].append(correction_files[f"{particle}_hlt"][self.year])
                    jsonfile_paths[particle]["data"].append(correction_files[f"{particle}_hlt"][self.year])
                    jsonfile_paths[particle]["sf"].append(correction_files[f"{particle}"][self.year]) # the muon correction set from the POG contains the Single Muon scale factors                 
            elif particle == "electron":
                if self.use_electron_HWW_files:
                    cset_files["mc"] = correction_files[f"{particle}_hlt_HWW_mc_eff"][self.year]
                    cset_files["data"] = correction_files[f"{particle}_hlt_HWW_data_eff"][self.year]
                    cset_files["sf"] = ""
                    for eff_type in eff_types:
                        if eff_type == "sf":
                            # for HWW working points, no scale factors immediately available: need to be calculated with efficiencies
                            jsonfile_name = ""
                        else:
                            jsonfile_name = f"{leg}_pt_eta_efficiency.json"
                        jsonfile_paths[particle][eff_type].append(
                            f"{cset_files[eff_type]}{jsonfile_name}"
                        )
                elif year_key < 2022:
                    jsonfile_paths[particle]["mc"].append(correction_files[f"{particle}_hlt_mc_eff"][self.year])
                    jsonfile_paths[particle]["data"].append(correction_files[f"{particle}_hlt_data_eff"][self.year])
                    jsonfile_paths[particle]["sf"].append(correction_files[f"{particle}_hlt"][self.year])
                else:
                    jsonfile_paths[particle]["mc"].append(correction_files[f"{particle}_hlt"][self.year])
                    jsonfile_paths[particle]["data"].append(correction_files[f"{particle}_hlt"][self.year])
                    jsonfile_paths[particle]["sf"].append(correction_files[f"{particle}_hlt"][self.year]) # the electron correction set from the POG contains the Single Electron scale factors 
            else:
                raise ValueError(
                    f"Unknown leg_particle '{particle}' for determining trigger SF path."
                )                
                    
        elif trigger_type == "Double Lepton Trigger":
            # in this case, always the efficiencies must be used instead of the scale factors
            legnumber = 1
            for leg in legs:
                particle = leg_particles[legnumber-1]
                jsonfile_paths = {
                    particle: {
                        "mc": [],
                        "data": [],
                        "sf": [],
                    }
                }
                if particle == "muon":
                    if self.use_muon_HWW_files:
                        cset_files["mc"] = correction_files[f"{particle}_hlt_HWW_mc_eff"][self.year]
                        cset_files["data"] = correction_files[f"{particle}_hlt_HWW_data_eff"][self.year]
                        cset_files["sf"] = ""                  
                        for eff_type in eff_types:
                            if eff_type == "sf":
                                # for HWW working points, no scale factors immediately available: need to be calculated with efficiencies
                                jsonfile_name_leg1 = ""
                                jsonfile_name_leg2 = ""
                            else:
                                jsonfile_name_leg1 = f"{legs[0]}_{legs[1]}_leg1_pt_eta_efficiency.json"
                                jsonfile_name_leg2 = f"{legs[0]}_{legs[1]}_leg2_pt_eta_efficiency.json"
                            jsonfile_paths[particle][eff_type].append(
                                f"{cset_files[eff_type]}{jsonfile_name_leg1}"
                            )
                            jsonfile_paths[particle][eff_type].append(
                                f"{cset_files[eff_type]}{jsonfile_name_leg2}" 
                            )
                    else:
                        jsonfile_paths[particle]["mc"].append(correction_files[f"{particle}_hlt"][self.year])
                        jsonfile_paths[particle]["data"].append(correction_files[f"{particle}_hlt"][self.year])
                        jsonfile_paths[particle]["sf"].append(correction_files[f"{particle}"][self.year]) # the muon correction set from the POG contains the Single Muon scale factors 
                elif particle == "electron":
                    if self.use_electron_HWW_files:
                        cset_files["mc"] = correction_files[f"{particle}_hlt_HWW_mc_eff"][self.year]
                        cset_files["data"] = correction_files[f"{particle}_hlt_HWW_data_eff"][self.year]
                        cset_files["sf"] = ""
                        for eff_type in eff_types:
                            if eff_type == "sf":
                                # for HWW working points, no scale factors immediately available: need to be calculated with efficiencies
                                jsonfile_name_leg = ""
                            else:
                                jsonfile_name_leg = f"{legs[0]}_{legs[1]}_leg{legnumber}_pt_eta_efficiency.json"
                            jsonfile_paths[particle][eff_type].append(
                                f"{cset_files[eff_type]}{jsonfile_name_leg}"
                            )
                    elif year_key < 2022:
                        jsonfile_paths[particle]["mc"].append(correction_files[f"{particle}_hlt_mc_eff"][self.year])
                        jsonfile_paths[particle]["data"].append(correction_files[f"{particle}_hlt_data_eff"][self.year])
                        jsonfile_paths[particle]["sf"].append(correction_files[f"{particle}_hlt"][self.year])
                    else:
                        jsonfile_paths[particle]["mc"].append(correction_files[f"{particle}_hlt"][self.year])
                        jsonfile_paths[particle]["data"].append(correction_files[f"{particle}_hlt"][self.year])
                        jsonfile_paths[particle]["sf"].append(correction_files[f"{particle}_hlt"][self.year]) # the electron correction set from the POG contains the Single Electron scale factors 
                else:
                    raise ValueError(
                        f"Unknown leg_particle '{particle}' for determining trigger SF path."
                    ) 
                legnumber += 1
        else:
            print(f"WARNING: provided trigger type {trigger_type} not recognized for retrieving correctionset path.")  
        return jsonfile_paths
    
    def get_muons_in_SFbinning_mask(self):
        # gets mask for muons within SF binning
        # later we can change this to be automatically detected by scanning the json file
        if self.use_muon_HWW_files:
            pt_lower_limit = 29.0 if self.year == "2017" else 5.0 # value not exactly known for 2017 or run 2
            pt_upper_limit = 150.0
            abseta_upper_limit = 2.4
        else:
            pt_lower_limit = 29.0 if self.year == "2017" else 26.0
            pt_upper_limit = 199.99
            abseta_upper_limit = 2.399
            
        muon_pt_mask = (self.flat_muons.pt > pt_lower_limit) & (self.flat_muons.pt < pt_upper_limit)
        muon_eta_mask = np.abs(self.flat_muons.eta) < abseta_upper_limit
        
        in_muons_mask = muon_pt_mask & muon_eta_mask
    
        return in_muons_mask

    def get_electrons_in_SFbinning_mask(self):
        # gets mask for muons within SF binning
        # later we can change this to be automatically detected by scanning the json file
        pt_lower_limits = {
            "2016preVFP": 27.0,
            "2016postVFP": 27.0,
            "2017": 35.0,
            "2018": 32.0,
        }
        if self.use_electron_HWW_files:
            pt_lower_limit = 25.0 # not exactly known for run 2. still need to check if this is ok. This is just lower limit in json file, but not checked if also electron WP needs to be taken into account
            pt_upper_limit = 200.0
            abseta_upper_limit = 2.5
            electron_pt_mask = (self.flat_electrons.pt > pt_lower_limit) & (self.flat_electrons.pt < pt_upper_limit)
            electron_eta_mask = np.abs(self.flat_electrons.eta) < abseta_upper_limit
        else:
            pt_upper_limit = 499.9
            abseta_upper_limit = 2.5
            if self.year in pt_lower_limits.keys():
                pt_lower_limit = pt_lower_limits[self.year]
                electron_pt_mask = (self.flat_electrons.pt > pt_lower_limit) & (self.flat_electrons.pt < pt_upper_limit)
                eta_mask = (np.abs(self.flat_electrons.eta + self.flat_electrons.deltaEtaSC) < abseta_upper_limit)
            else:
                # get electrons that pass the id wp, and within SF binning
                pt_lower_limit = 25.0 # still need to check if this is ok. This is just lower limit in json file, but not checked if also electron WP needs to be taken into account
                electron_pt_mask = self.flat_electrons.pt > pt_lower_limit
                electron_eta_mask = np.abs(self.flat_electrons.eta) < abseta_upper_limit
        
        in_electrons_mask = electron_pt_mask & electron_eta_mask
    
        return in_electrons_mask
    
    def get_in_muons_kinematics(self, in_muons_mask):
        # get muons pT and abseta (replace None values with some 'in-limit' value) from the muons in the SF binning
        in_muons = self.flat_muons.mask[in_muons_mask]
        muon_pt = ak.fill_none(in_muons.pt, 30.0)
        muon_eta = ak.fill_none(in_muons.eta, 0)
        muon_abseta = ak.fill_none(np.abs(in_muons.eta), 0)
        return muon_pt, muon_eta, muon_abseta

    def get_in_electrons_kinematics(self, in_electrons_mask):
        # get electrons pT and abseta (replace None values with some 'in-limit' value)
        run2_years  = ["2016preVFP", "2016postVFP", "2017", "2018"]
        in_electrons = self.flat_electrons.mask[in_electrons_mask]
        if self.year in run2_years:
            electron_pt = ak.fill_none(in_electrons.pt, 40)
            electron_eta = ak.fill_none(self.flat_electrons.eta + self.flat_electrons.deltaEtaSC, 0)
        else:
            electron_pt = ak.fill_none(in_electrons.pt, 25)
            electron_eta = ak.fill_none(in_electrons.eta, 0)
        return electron_pt, electron_eta

    def get_electron_variation_key(self, variation, eff_type):
        # keys for the variations in the json files from the EGamma POG, saved locally (electronHlt.json for example)
        # electronHLT.json only contains the nominal value and up and down, so no separation in systematic and statistical uncertainty
        if eff_type == "sf":
            if variation == "nominal":
                variation_key = "sf"
            else:
                variation_key = f"sf{variation}"
        elif eff_type == "eff":
            if variation == "nominal":
                variation_key = "nom"
            else:
                variation_key = variation
        else:
            raise ValueError(f"Unknow efficiency type '{eff_type}' given to get_electron_variation_key.")

        return variation_key
    
    def calculate_SingleMu_SF(self, jsonkey, has_sf_path, variation):
        particle = "muon"
        in_muons_mask = self.get_muons_in_SFbinning_mask()
        muon_pt, muon_eta, muon_abseta = self.get_in_muons_kinematics(in_muons_mask)
        
        kind = "single" if ak.all(self.muons_counts == 1) else "double"
        # for single muon events, take/calculate SF from POG if possible
        # for WPs not available from POG, we use privatly derived or HWW efficiencies to calculate SF
        if self.use_muon_HWW_files:   
            cset_eff_data_eval = self.cset_eff[particle]["data"][jsonkey]
            cset_eff_mc_eval = self.cset_eff[particle]["mc"][jsonkey]
            evaluate_eff_inputlist = [        
                muon_abseta if self.nano_version == "9" else muon_eta,
                muon_pt,
                variation,
            ]
            evaluate_sf_inputlist = evaluate_eff_inputlist
        else:
            cset_eff_data_eval = self.cset_eff[particle]["data"]["Muon-HLT-DataEff"]
            cset_eff_mc_eval = self.cset_eff[particle]["mc"]["Muon-HLT-McEff"]
             
            evaluate_eff_inputlist = [        
                variation,
                jsonkey,
                muon_abseta if self.nano_version == "9" else muon_eta,
                muon_pt,
            ]
            evaluate_sf_inputlist = [        
                muon_abseta if self.nano_version == "9" else muon_eta,
                muon_pt,
                variation,
            ]
            
        if kind == "single":
            if has_sf_path:
                cset_eff_sf_eval = self.cset_eff[particle]["sf"][jsonkey] # For the files with SF, only NUMDEN name needed
                # if SF file is available, take SF straight from there
                sf = cset_eff_sf_eval.evaluate(*evaluate_sf_inputlist) 
                nominal_sf = unflat_sf(sf, in_muons_mask, self.muons_counts)  
            else:
                # if SF file not available, calculate SF from data and mc efficiencies
                data_eff = cset_eff_data_eval.evaluate(*evaluate_eff_inputlist)
                mc_eff = cset_eff_mc_eval.evaluate(*evaluate_eff_inputlist)

                sf = data_eff / mc_eff
                nominal_sf = unflat_sf(sf, in_muons_mask, self.muons_counts) 
        # for double (or more) muon events, always compute SF from data/mc muon hlt efficiencies. Calculation is limited to only 2 muons. For more than 2 muons, same 2 muon calculation is applied at the moment (needs to be changed probably, when actually using triggers with 3 or more legs)
        elif kind == "double":
            data_eff = cset_eff_data_eval.evaluate(*evaluate_eff_inputlist)
            mc_eff = cset_eff_mc_eval.evaluate(*evaluate_eff_inputlist)
                
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

    def calculate_SingleEle_SF(self, jsonkey, has_sf_path, variation):
        particle = "electron"
        in_electrons_mask = self.get_electrons_in_SFbinning_mask() 
        electron_pt, electron_eta = self.get_in_electrons_kinematics(in_electrons_mask)

        kind = "single" if ak.all(self.electrons_counts == 1) else "double"
        # for single electron events, take/calculate SF from POG if possible
        # for WPs not available from POG, we use privatly derived or HWW efficiencies to calculate SF
        if self.year in ["2016postVFP", "2016postVFP", "2017", "2018"]:
            print("WARNING: no implementation of retrieving SingleEle electron trigger SFs or efficiencies yet. Read comments in analysis/corrections/trigger.py in the calculate_SingleEle_SF function.")
        if self.use_electron_HWW_files:   
            cset_eff_data_eval = self.cset_eff[particle]["data"][jsonkey]
            cset_eff_mc_eval = self.cset_eff[particle]["mc"][jsonkey]
            evaluate_eff_inputlist = [        
                electron_eta,
                electron_pt,
                variation,
            ]
            evaluate_sf_inputlist = evaluate_eff_inputlist
        else: 
            cset_eff_data_eval = self.cset_eff[particle]["data"]["Electron-HLT-DataEff"]
            cset_eff_mc_eval = self.cset_eff[particle]["mc"]["Electron-HLT-McEff"]
            year_map = {
                "2016preVFP": "2016preVFP", # for Run 2 years, there is a json file per year, stored locally: no year_map needed #### we need to handle the retrieval of the scale factors / efficiencies different for that: only the NUMDEN name is necessary and it specifies HLT-SF, HLT-DataEff or HLT-McEff at the begining of the NUMDEN name. The locally stored Run 2 files only conatin the nominal Scale Factor value
                "2016postVFP": "2016postVFP",
                "2017": "2017",
                "2018": "2018",
                "2022postEE": "2022Re-recoE+PromptFG",
                "2022preEE": "2022Re-recoBCD",
                "2023preBPix": "2023PromptC",
                "2023postBPix": "2023PromptD",
                "2024": "2024Prompt",
            }
            variation_key_sf = self.get_electron_variation_key(variation, "sf")
            variation_key_eff = self.get_electron_variation_key(variation, "eff")
            evaluate_eff_inputlist = [ 
                year_map[self.year],
                variation_key_eff,
                jsonkey, 
                electron_eta,
                electron_pt,
            ]
            
            evaluate_sf_inputlist = [ 
                year_map[self.year],
                variation_key_sf,
                jsonkey, 
                electron_eta,
                electron_pt,
            ]

        if kind == "single":
            if has_sf_path:
                cset_eff_sf_eval = self.cset_eff[particle]["sf"]["Electron-HLT-SF"]
                # if SF file is available, take SF straight from there
                sf = cset_eff_sf_eval.evaluate(*evaluate_sf_inputlist)
                nominal_sf = unflat_sf(sf, in_electrons_mask, self.electrons_counts)  
            else:
                # if SF file not available, calculate SF from data and mc efficiencies
                data_eff = cset_eff_data_eval.evaluate(*evaluate_eff_inputlist)
                mc_eff = cset_eff_mc_eval.evaluate(*evaluate_eff_inputlist)

                sf = data_eff / mc_eff
                nominal_sf = unflat_sf(sf, in_electrons_mask, self.electrons_counts) 
        # for double (or more) electron events, always compute SF from data/mc electron hlt efficiencies. Calculation is limited to only 2 electrons. For more than 2 electrons, same 2 electron calculation is applied at the moment (needs to be changed probably, when actually using triggers with 3 or more legs)
        elif kind == "double":
            data_eff = cset_eff_data_eval.evaluate(*evaluate_eff_inputlist)
            mc_eff = cset_eff_mc_eval.evaluate(*evaluate_eff_inputlist)
                
            data_eff = ak.where(in_electrons_mask, data_eff, ak.ones_like(data_eff))
            data_eff = ak.unflatten(data_eff, self.electrons_counts)
            data_eff_leading = ak.firsts(data_eff)
            data_eff_subleading = ak.pad_none(data_eff, target=2)[:, 1]
            full_data_eff = (
                data_eff_leading
                + data_eff_subleading
                - data_eff_leading * data_eff_subleading
            )
            full_data_eff = ak.fill_none(full_data_eff, 1)
        
            mc_eff = ak.where(in_electrons_mask, mc_eff, ak.ones_like(mc_eff))
            mc_eff = ak.unflatten(mc_eff, self.electrons_counts)
            mc_eff_leading = ak.firsts(mc_eff)
            mc_eff_subleading = ak.pad_none(mc_eff, target=2)[:, 1]
            full_mc_eff = (
                mc_eff_leading + mc_eff_subleading - mc_eff_leading * mc_eff_subleading
            )
            full_mc_eff = ak.fill_none(full_mc_eff, 1)
        
            nominal_sf = full_data_eff / full_mc_eff
        return nominal_sf
    
    def evaluate_trigger_efficiency(self):
        return 0
        
    def add_trigger_weights(self, triggerSF_particles, triggerSFs_leptonWP_inputs, hlt_paths, dataset):
        """
        add muon and electron trigger weights to weights container
        """
        # get nominal scale factors
        print("Getting nominal hlt weights.")
        nominal_weights = self.get_hlt_weights(
            triggerSF_particles,
            triggerSFs_leptonWP_inputs,
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
        
    def get_hlt_weights(self, triggerSF_particles, triggerSFs_leptonWP_inputs, hlt_paths, variation):
        """Compute muon/electron HLT weights"""
        muon_id_wp = triggerSFs_leptonWP_inputs["muon"]["id"]
        muon_iso_wp = triggerSFs_leptonWP_inputs["muon"]["iso"]
        muon_promptMVA_wp = triggerSFs_leptonWP_inputs["muon"]["promptMVA"]

        electron_id_wp = triggerSFs_leptonWP_inputs["electron"]["id"]
        electron_promptMVA_wp = triggerSFs_leptonWP_inputs["electron"]["promptMVA"]       

        #####
        # need to check once again the bin lower and upper edges for scale factor ANd efficiency json files
        # especially for centrally provided ones: they go sometimes to infinity
        #####
        
        print(f"Muon working points received for applying muon efficiency scale factors:")
        if muon_id_wp: print(f"ID: {muon_id_wp}") 
        if muon_iso_wp: print(f"ISO: {muon_iso_wp}") 
        if muon_promptMVA_wp: print(f"PromptMVA: {muon_promptMVA_wp}") 

        print(f"Electron working points received for applying electron efficiency scale factors:")
        if electron_id_wp: print(f"ID: {electron_id_wp}") 
        if electron_promptMVA_wp: print(f"PromptMVA: {electron_promptMVA_wp}") 

        print(f"List of particles for which to apply trigger SFs: {triggerSF_particles}")
            
        flags_with_triggernames = tree()
        flags_with_triggernames = self.get_config_trigger_flags(hlt_paths) # get trigger flags listed in workflow config
        trigger_flags = self.load_trigger_flags(Path(__file__)) # load file having for each year the trigger flag and corresponding trigger names
        year_key = self.get_yearkey_for_trigger_flags()
        trigger_types = flags_with_triggernames.keys() # types of trigger: single lepton trigger, double lepton trigger ... defined in analysis/corrections/trigger.py get_config_trigger_flags()
        
        # get for each trigger flag in the workflow config, the proper trigger names 
        # beware that the combination of working points must be available in the correction files
        # for example, only tight ID and ISO are available to be combined with the tthMVA-Run3 promptMVA working point
        for trigger_type in trigger_types:
            
            for trigger_flag in flags_with_triggernames[trigger_type].keys():
                
                trigger_names = trigger_flags[year_key][trigger_flag]
                
                for trigger_name in trigger_names:
        
                    legs, leg_particles = self.get_trigger_legs(trigger_name)
                    cset_filepaths = self.get_cset_path(leg_particles, legs, trigger_type, triggerSFs_leptonWP_inputs)
                    
                    legnumber = 1
                    eff_types = ["mc", "data", "sf"]
                    for leg_particle in leg_particles: 
                        if leg_particle in triggerSF_particles:
                            flags_with_triggernames[trigger_type][trigger_flag][trigger_name][f"jsonleg{legnumber}"] = legs[legnumber-1]
                            flags_with_triggernames[trigger_type][trigger_flag][trigger_name][f"jsonkey{legnumber}"] = self.construct_trigger_NUMDEN_name(legs, leg_particles, trigger_flag, trigger_name, triggerSFs_leptonWP_inputs, hlt_paths)
                        else:
                            print(f"WARNING: leg particle {leg_particle} was not intended for trigger SFs to be applied (check your workflow config file), so skipping this leg.")
                            if trigger_type != "Single Lepton Trigger":
                                raise ValueError(
                                    f"A scale factor cannot be computed if one of the legs' efficiencies cannot be used. Make sure the particles corresponding to all trigger legs have trigger set to 'true' in the workflow config file."
                                ) 
                            else:
                                continue
                            
                        # Store the correct path to the efficiency json for each leg 
                        for eff_type in eff_types:
                            print("cset_filepaths[eff_type]", cset_filepaths[leg_particle][eff_type])
                            flags_with_triggernames[trigger_type][trigger_flag][trigger_name][f"cset_{eff_type}_filepath_leg{legnumber}"] = cset_filepaths[leg_particle][eff_type][legnumber-1]
                        
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
                            # check that also the workflow config asked to apply the trigger SF for the corresponding particle
                            if (trigger_flag == "SingleMu") and ("muon" in triggerSF_particles):
                                apply_triggertype[trigger_flag] = True
                            elif (trigger_flag == "SingleEle") and ("electron" in triggerSF_particles):
                                apply_triggertype[trigger_flag] = True
                            else:
                                apply_triggertype[trigger_flag] = False
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
                particle = "muon"
                trigger_type = "Single Lepton Trigger"
                trigger_flag = "SingleMu"
                trigger_name = list(flags_with_triggernames[trigger_type][trigger_flag].keys())[0]
                
                for eff_type in eff_types:
                    filepath = flags_with_triggernames[trigger_type][trigger_flag][trigger_name][
                        f"cset_{eff_type}_filepath_leg1"
                    ]
                
                    if filepath:
                        self.cset_eff[particle][eff_type] = correctionlib.CorrectionSet.from_file(filepath)
                    
                jsonkey = flags_with_triggernames[trigger_type][trigger_flag][trigger_name]["jsonkey1"]
                    #display(JSON(json.loads(json.dumps(self.cset_eff))))
                if not flags_with_triggernames[trigger_type][trigger_flag][trigger_name][f"cset_sf_filepath_leg1"] == "":
                    has_sf_path = True
                else:
                    has_sf_path = False 
                    
                nominal_sf = self.calculate_SingleMu_SF(jsonkey, has_sf_path, variation)
                
                # geeft het problemen als je een empty string aan correctionSet.from_file geeft? Nee
        
            case [True, False, True, False]:
            # only SingleEle and no double lepton trigger
                print("Only SingleEle and no double lepton trigger provided. Applying the SingleEle trigger SF to all events with 1 or more electrons.")
                particle = "electron"
                trigger_type = "Single Lepton Trigger"
                trigger_flag = "SingleEle"
                trigger_name = list(flags_with_triggernames[trigger_type][trigger_flag].keys())[0]
                
                for eff_type in eff_types:
                    filepath = flags_with_triggernames[trigger_type][trigger_flag][trigger_name][
                        f"cset_{eff_type}_filepath_leg1"
                    ]
                
                    if filepath:
                        self.cset_eff[particle][eff_type] = correctionlib.CorrectionSet.from_file(filepath)
                
                jsonkey = flags_with_triggernames[trigger_type][trigger_flag][trigger_name]["jsonkey1"]
                    
                if not flags_with_triggernames[trigger_type][trigger_flag][trigger_name][f"cset_sf_filepath_leg1"] == "":
                    has_sf_path = True
                else:
                    has_sf_path = False 
                    
                nominal_sf = self.calculate_SingleEle_SF(jsonkey, has_sf_path, variation)
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