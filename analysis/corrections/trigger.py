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
        self.muons = events.selected_muons
        self.flat_muons = ak.flatten(self.muons)
        self.muons_counts = ak.num(self.muons)        

        # electrons
        self.electrons = events.selected_electrons
        self.flat_electrons = ak.flatten(self.electrons)
        self.electrons_counts = ak.num(self.electrons)
        
        # get muon id/iso/HLT correction set
        # get muon id/iso/HLT correction set from muon POG
        self.cset = tree()
        self.cset["muon"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon"][year]
        )
        # get muon HLT eff correction set
        self.cset_eff = tree()
        self.cset_eff_hlt = tree()
        self.cset_eff_mc = tree()
        self.cset_eff_data = tree()
        self.cset_eff["muon"]["mc"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon_hlt"][year]
        )
        self.cset_eff["muon"]["data"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon_hlt"][year]
        )
        self.cset_eff["muon"]["sf"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon"][year]
        )
        self.cset_eff_hlt["muon"] = correctionlib.CorrectionSet.from_file(
            correction_files["muon_hlt"][year]
        )
        self.cset_eff_data["muon"] = self.cset_eff_hlt["Muon-HLT-DataEff"]
        self.cset_eff_mc["muon"] = self.cset_eff_hlt["Muon-HLT-McEff"]

        # get electron HLT eff correction set
    
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
    
    def get_available_SF_names_DEN(self, id_wp, iso_wp, promptMVA_wp):
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
                        "tight_HWW": "testnametighttightHWWonly",
                        "": "CutBasedIdTight_and_PFIsoTight",
                    },
                    "tight": {
                        "tight_HWW": "testnametighttighttightHWW",
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
        print(f"Getting trigger_SF_DEN for nano_version '{self.nano_version}', ID WP '{id_wp}', ISO WP '{iso_wp}', promptMVA WP '{promptMVA_wp}'.")
        trigger_SF_DEN = available_trigger_SFs[self.nano_version][id_wp][iso_wp][promptMVA_wp]
        print(f"trigger_SF_DEN: {trigger_SF_DEN}")
        if trigger_SF_DEN is None:
            raise ValueError(
                f"There are no muon HLT weights for this id wp '{id_wp}', iso wp '{iso_wp}', promptMVA wp '{iso_wp}' combination in year {self.year}."
            )
        else:
            return trigger_SF_DEN
    
    def construct_trigger_NUMDEN_name(self, trigger_flag, trigger_name, id_wp, iso_wp, promptMVA_wp, hlt_paths):
        # the name of the DEN should be the same for each of the triggers as long as they are the same type (for example: SingleMu)
        # inside the HWW trigger efficiency scale factor files, all efficiencies are listed under TriggerEff and the particular name does not even matter - the trigger name or name of the leg is in the json file name wherein the correction factor is stored
        trigger_DEN_name = self.get_available_SF_names_DEN(id_wp, iso_wp, promptMVA_wp)
        trigger_NUM_name = trigger_name
        if promptMVA_wp == "tight_HWW" or id_wp == "tight_HWW":
            # tight_HWW for muon ISO is exactly the same as tight muon ISO, 
            # but for ID we have tight + other selections in the definitions in tight_HWW
            trigger_NUMDEN_name = "TriggerEff"
            #trigger_NUMDEN_name = "NUM_" + trigger_NUM_name + "_DEN_" + trigger_DEN_name
        else:
            if trigger_flag in ["SingleMu", "SingleEle"]:
                trigger_NUMDEN_name = "NUM_" + trigger_NUM_name + "_DEN_" + trigger_DEN_name
            else:
                trigger_NUMDEN_name = trigger_NUM_name
        return trigger_NUMDEN_name
    
    def get_trigger_legs(self, trigger_name):
        """
        Extracts the trigger legs from the trigger name: every piece in the trigger name starting with Mu or Ele.
        """
        legs = []
        name_parts = trigger_name.split('_')
        for part in name_parts:
            if part.startswith('Mu') or part.startswith('Ele') or part.startswith('IsoMu'):
                print(part)
                legs.append(part)
        return legs
    
    def get_cset_path(self, particle, legs, trigger_type, id_wp, iso_wp, promptMVA_wp):
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
        eff_types = ["mc", "data", "sf"]
        jsonfile_paths = {
            particle: {
                "mc": [],
                "data": [],
                "sf": [],
            }
        }
        cset_files = {
            "mc": [],
            "data": [],
            "sf": [],
        }
        if promptMVA_wp == "tight_HWW" or id_wp == "tight_HWW":
            cset_files["mc"] = correction_files["muon_hlt_HWW_mc_eff"][self.year]
            cset_files["data"] = correction_files["muon_hlt_HWW_data_eff"][self.year]
            cset_files["sf"] = ""
        
            if trigger_type == "Single Lepton Trigger":
                for eff_type in eff_types:
                    if eff_type == "sf":
                        # for HWW working points, no scale factors immediately available: need to be calculated with efficiencies
                        jsonfile_name = ""
                    else:
                        jsonfile_name = f"{legs[0]}_pt_eta_efficiency.json"
                    jsonfile_paths[particle][eff_type].append(
                        f"{cset_files[eff_type]}{jsonfile_name}"
                    )
                    
            elif trigger_type == "Double Lepton Trigger":
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
                print(f"WARNING: provided trigger type {trigger_type} not recognized for retrieving correctionset path.")
        else:
            for leg in legs:
                # in this case, efficiencies for each leg are obtained from the same file
                jsonfile_paths[particle]["mc"].append(correction_files["muon_hlt"][self.year])
                jsonfile_paths[particle]["data"].append(correction_files["muon_hlt"][self.year])
                jsonfile_paths[particle]["sf"].append(correction_files["muon"][self.year]) # the muon correction set from the POG contains the Single Muon scale factors 
        
        return jsonfile_paths
    
    def get_muons_in_SFbinning_mask(self):
        # gets mask for muons within SF binning
        # later we can change this to be automatically detected by scanning the json file
        pt_lower_limit = 29.0 if self.year == "2017" else 26.0
        muon_pt_mask = (self.flat_muons.pt > pt_lower_limit) & (
            self.flat_muons.pt < 199.99
        )
        muon_eta_mask = np.abs(self.flat_muons.eta) < 2.399
        
        in_muons_mask = muon_pt_mask & muon_eta_mask
    
        return in_muons_mask
    
    def get_in_muons_kinematics(self, in_muons_mask):
        # get muons pT and abseta (replace None values with some 'in-limit' value) from the muons in the SF binning
        in_muons = self.flat_muons.mask[in_muons_mask]
        muon_pt = ak.fill_none(in_muons.pt, 30.0)
        muon_eta = ak.fill_none(in_muons.eta, 0)
        muon_abseta = ak.fill_none(np.abs(in_muons.eta), 0)
        return muon_pt, muon_eta, muon_abseta
    
    def calculate_SingleMu_SF(self, jsonkey, has_sf_path, variation):
        particle = "muon"
        in_muons_mask = self.get_muons_in_SFbinning_mask()
        muon_pt, muon_eta, muon_abseta = self.get_in_muons_kinematics(in_muons_mask)
        
        kind = "single" if ak.all(self.muons_counts == 1) else "double"
        # for single muon events, take/calculate SF from POG if possible
        # for WPs not available from POG, we use privatly derived or HWW efficiencies to calculate SF
        if kind == "single":
            if has_sf_path:
                # if SF file is available, take SF straight from there
                sf = self.cset_eff[particle]["sf"][jsonkey].evaluate(
                    muon_abseta, muon_pt, variation
                )
                nominal_sf = unflat_sf(sf, in_muons_mask, self.muons_counts)  
            else:
                # if SF file not available, calculate SF from data and mc efficiencies
                data_eff = self.cset_eff[particle]["data"]["Muon-HLT-DataEff"].evaluate(
                    variation,
                    jsonkey,
                    muon_abseta if self.nano_version == "9" else muon_eta,
                    muon_pt,
                )
                mc_eff = self.cset_eff[particle]["mc"]["Muon-HLT-McEff"].evaluate(
                    variation,
                    jsonkey,
                    muon_abseta if self.nano_version == "9" else muon_eta,
                    muon_pt,
                )               
                sf = data_eff / mc_eff
                nominal_sf = unflat_sf(sf, in_muons_mask, self.muons_counts) 
        # for double (or more) muon events, always compute SF from data/mc muon hlt efficiencies. Calculation is limited to only 2 muons. For more than 2 muons, same 2 muon calculation is applied
        elif kind == "double":
            data_eff = self.cset_eff[particle]["data"]["Muon-HLT-DataEff"].evaluate(
                variation,
                jsonkey,
                muon_abseta if self.nano_version == "9" else muon_eta,
                muon_pt,
            )
            mc_eff = self.cset_eff[particle]["mc"]["Muon-HLT-McEff"].evaluate(
                variation,
                jsonkey,
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
    
    def evaluate_trigger_efficiency(self):
        return 0
    
    def set_trigger_correction_file(self):
        return 0
    
    
    # copied from muon.py
    def add_trigger_weights(self, id_wp, iso_wp, promptMVA_wp, hlt_paths, dataset):
        """
        add muon iso weights to weights container
        """
        # get nominal scale factors
        print("Getting nominal hlt weights.")
        nominal_weights = self.get_hlt_weights(
            particle="muon", ### probably need to remove particle again and rather move to individual id_wp ... for each particle, as in get_hlt_weights, we want to calculate the complete scale factor
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
        
    def get_hlt_weights(self, particle, id_wp, iso_wp, promptMVA_wp, hlt_paths, variation):
        """Compute muon/electron HLT weights"""
        print(f"Particle: {particle} - working points received for applying muon efficiency scale factors:")
        if id_wp: print(f"ID: {id_wp}") 
        if iso_wp: print(f"ISO: {iso_wp}") 
        if promptMVA_wp: print(f"PromptMVA: {promptMVA_wp}") 
            
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
                    flags_with_triggernames[trigger_type][trigger_flag][trigger_name]["jsonkey"] = self.construct_trigger_NUMDEN_name(trigger_flag, trigger_name, id_wp, iso_wp, promptMVA_wp, hlt_paths)
        
                    legs = self.get_trigger_legs(trigger_name)
                    cset_filepaths = self.get_cset_path(particle, legs, trigger_type, id_wp, iso_wp, promptMVA_wp)
                    legnumber = 1
                    eff_types = ["mc", "data", "sf"]
                    for leg in legs: 
                        flags_with_triggernames[trigger_type][trigger_flag][trigger_name][f"jsonleg{legnumber}"] = legs[legnumber-1]
        
                        for eff_type in eff_types:
                            print("cset_filepaths[eff_type]", cset_filepaths[particle][eff_type])
                            flags_with_triggernames[trigger_type][trigger_flag][trigger_name][f"cset_{eff_type}_filepath_leg{legnumber}"] = cset_filepaths[particle][eff_type][legnumber-1]
                        
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
                    
                nominal_sf = self.calculate_SingleMu_SF(jsonkey, has_sf_path, variation)
                
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