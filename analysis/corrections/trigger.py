import numpy as np
import awkward as ak
import yaml
from collections import defaultdict
from analysis.filesets.utils import get_nano_version
from analysis.corrections.utils import correction_files, unflat_sf

def tree():
    return defaultdict(tree)

def get_config_trigger_flags(hlt_paths):
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
    
def load_trigger_flags(filepath):
    """Load trigger flag definitions from analysis/selections."""
    trigger_flags_path = (
        filepath.resolve().parent.parent
        / "selections"
        / "trigger_flags.yaml"
    )

    with trigger_flags_path.open("r") as f:
        return yaml.safe_load(f)


def get_yearkey_for_trigger_flags(year):
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
    year_key = year_map[year]
    return year_key

def get_available_SF_names_DEN(year, id_wp, iso_wp, promptMVA_wp):
    """
    Returns the available scale factor names based on the nano version and the year, as well as what type of working point is being looked at. It raises a Value error if the combination of working points does not correspond to a combination for which trigger Scale Factors are available
    """
    # if no promptMVA working point was given, that is treated as an empty string.
    # dictionary order of levels: nanoversion > id_wp > iso_wp > promptMVA_wp: is the trigger name available?: trigger as NUM and whole WP as DEN
    years_with_trigger_SFs  = ["2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
    run2_years = ["2016preVFP", "2016postVFP", "2017", "2018"]
    if year not in years_with_trigger_SFs:
        raise ValueError(
            f"There are no muon HLT weights/trigger scale factors for year '{year}'"
        )
    if year in run2_years:
        raise ValueError(
            f"HLT weights/trigger scale factors are available for year '{year}', but not implemented in the framework yet"
        )
    nanoversion = get_nano_version(year)
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
    print(f"Getting trigger_SF_DEN for nanoversion '{nanoversion}', ID WP '{id_wp}', ISO WP '{iso_wp}', promptMVA WP '{promptMVA_wp}'.")
    trigger_SF_DEN = available_trigger_SFs[nanoversion][id_wp][iso_wp][promptMVA_wp]
    print(f"trigger_SF_DEN: {trigger_SF_DEN}")
    if trigger_SF_DEN is None:
        raise ValueError(
            f"There are no muon HLT weights for this id wp '{id_wp}', iso wp '{iso_wp}', promptMVA wp '{iso_wp}' combination in year {year}."
        )
    else:
        return trigger_SF_DEN

    
def construct_trigger_NUMDEN_name(trigger_flag, trigger_name, year, id_wp, iso_wp, promptMVA_wp, hlt_paths):
    # the name of the DEN should be the same for each of the triggers as long as they are the same type (for example: SingleMu)
    # inside the HWW trigger efficiency scale factor files, all efficiencies are listed under TriggerEff and the particular name does not even matter - the trigger name or name of the leg is in the json file name wherein the correction factor is stored
    trigger_DEN_name = get_available_SF_names_DEN(year, id_wp, iso_wp, promptMVA_wp)
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

def get_trigger_legs(trigger_name):
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

def get_cset_path(legs, trigger_type, year, id_wp, iso_wp, promptMVA_wp):
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
        "mc": [],
        "data": [],
        "sf": [],
    }
    cset_files = {
        "mc": [],
        "data": [],
        "sf": [],
    }
    if promptMVA_wp == "tight_HWW" or id_wp == "tight_HWW":
        cset_files["mc"] = correction_files["muon_hlt_HWW_mc_eff"][year]
        cset_files["data"] = correction_files["muon_hlt_HWW_data_eff"][year]
        cset_files["sf"] = ""

        if trigger_type == "Single Lepton Trigger":
            for eff_type in eff_types:
                if eff_type == "sf":
                    # for HWW working points, no scale factors immediately available: need to be calculated with efficiencies
                    jsonfile_name = ""
                else:
                    jsonfile_name = f"{legs[0]}_pt_eta_efficiency.json"
                jsonfile_paths[eff_type].append(
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
                jsonfile_paths[eff_type].append(
                    f"{cset_files[eff_type]}{jsonfile_name_leg1}"
                )
                jsonfile_paths[eff_type].append(
                    f"{cset_files[eff_type]}{jsonfile_name_leg2}"
                )
        else:
            print(f"WARNING: provided trigger type {trigger_type} not recognized for retrieving correctionset path.")
    else:
        for leg in legs:
            # in this case, efficiencies for each leg are obtained from the same file
            jsonfile_paths["mc"].append(correction_files["muon_hlt"][year])
            jsonfile_paths["data"].append(correction_files["muon_hlt"][year])
            jsonfile_paths["sf"].append(correction_files["muon"][year]) # the muon correction set from the POG contains the Single Muon scale factors 

    return jsonfile_paths

def get_muons_in_SFbinning_mask(year, flat_muons):
    # gets mask for muons within SF binning
    # later we can change this to be automatically detected by scanning the json file
    pt_lower_limit = 29.0 if year == "2017" else 26.0
    muon_pt_mask = (flat_muons.pt > pt_lower_limit) & (
        flat_muons.pt < 199.99
    )
    muon_eta_mask = np.abs(flat_muons.eta) < 2.399

    in_muons_mask = muon_pt_mask & muon_eta_mask
    
    return in_muons_mask

def get_in_muons_kinematics(flat_muons, in_muons_mask):
    # get muons pT and abseta (replace None values with some 'in-limit' value) from the muons in the SF binning
    in_muons = flat_muons.mask[in_muons_mask]
    muon_pt = ak.fill_none(in_muons.pt, 30.0)
    muon_eta = ak.fill_none(in_muons.eta, 0)
    muon_abseta = ak.fill_none(np.abs(in_muons.eta), 0)
    return muon_pt, muon_eta, muon_abseta
    
def calculate_SingleMu_SF(events, cset_eff, jsonkey, year, variation, has_sf_path):
    muons = events.selected_muons
    flat_muons = ak.flatten(muons)
    muons_counts = ak.num(muons)
    in_muons_mask = get_muons_in_SFbinning_mask(year, flat_muons)
    muon_pt, muon_eta, muon_abseta = get_in_muons_kinematics(flat_muons, in_muons_mask)
    
    kind = "single" if ak.all(ak.num(muons) == 1) else "double"
    # for single muon events, take/calculate SF from POG if possible
    # for WPs not available from POG, we use privatly derived or HWW efficiencies to calculate SF
    if kind == "single":
        if has_sf_path:
            # if SF file is available, take SF straight from there
            sf = cset_eff["sf"][jsonkey].evaluate(
                muon_abseta, muon_pt, variation
            )
            nominal_sf = unflat_sf(sf, in_muons_mask, muons_counts)  
        else:
            # if SF file not available, calculate SF from data and mc efficiencies
            data_eff = cset_eff["data"]["Muon-HLT-DataEff"].evaluate(
                variation,
                jsonkey,
                muon_abseta if self.nano_version == "9" else muon_eta,
                muon_pt,
            )
            mc_eff = cset_eff["mc"]["Muon-HLT-McEff"].evaluate(
                variation,
                jsonkey,
                muon_abseta if self.nano_version == "9" else muon_eta,
                muon_pt,
            )               
            sf = data_eff / mc_eff
            nominal_sf = unflat_sf(sf, in_muons_mask, muons_counts) 
    # for double (or more) muon events, always compute SF from data/mc muon hlt efficiencies. Calculation is limited to only 2 muons. For more than 2 muons, same 2 muon calculation is applied
    elif kind == "double":
        data_eff = cset_eff["data"]["Muon-HLT-DataEff"].evaluate(
            variation,
            jsonkey,
            muon_abseta if self.nano_version == "9" else muon_eta,
            muon_pt,
        )
        mc_eff = self.cset_eff["mc"]["Muon-HLT-McEff"].evaluate(
            variation,
            jsonkey,
            muon_abseta if self.nano_version == "9" else muon_eta,
            muon_pt,
        )
            
        data_eff = ak.where(in_muons_mask, data_eff, ak.ones_like(data_eff))
        data_eff = ak.unflatten(data_eff, muons_counts)
        data_eff_leading = ak.firsts(data_eff)
        data_eff_subleading = ak.pad_none(data_eff, target=2)[:, 1]
        full_data_eff = (
            data_eff_leading
            + data_eff_subleading
            - data_eff_leading * data_eff_subleading
        )
        full_data_eff = ak.fill_none(full_data_eff, 1)

        mc_eff = ak.where(in_muons_mask, mc_eff, ak.ones_like(mc_eff))
        mc_eff = ak.unflatten(mc_eff, muons_counts)
        mc_eff_leading = ak.firsts(mc_eff)
        mc_eff_subleading = ak.pad_none(mc_eff, target=2)[:, 1]
        full_mc_eff = (
            mc_eff_leading + mc_eff_subleading - mc_eff_leading * mc_eff_subleading
        )
        full_mc_eff = ak.fill_none(full_mc_eff, 1)

        nominal_sf = full_data_eff / full_mc_eff
    return nominal_sf

def evaluate_trigger_efficiency():
    return 0

def set_trigger_correction_file():
    return 0


