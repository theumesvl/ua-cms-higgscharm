# Corrections

### Trigger Scale Factors / HLT weigths
For now, mainly two files in `analysis/corrections/` are used for this: `muon.py` and `trigger.py`. In `trigger.py` all the helper functions for applying the trigger scale factors are stored and in `muon.py` the actual scale factors are applied.
Later, when combining with electron legs and electron trigger scale factors, it might be better to have one dedicated file for the triggers and move the content for the triggers fully in one file for electrons as well as muons. For now this is done in `muon.py`, so there we actually also look at electron scale factors. `electron.py` also still contains a part on applying Single Electron trigger scale factors, which eventually should be fully removed.

In `muon.py`, we start with listing the ID, ISO and promptMVA working points for the muons. ID and ISO working points are required to be given in the workflow config file. `analysis/corrections/correction_manager.py` will otherwise give an error in case the trigger correction option is set to `true` in the workflow config file. Is no promptMVA WP is specified, the `correction_manager.py` will set this `promptMVA_wp` to an empty string and just continue functioning normally.
Further in `muon.py`, a whole dictionalry/tree `flags_with_triggernames` is constructed with all the trigger info: the type of trigger, the flags, the corresponding trigger names, the names of each trigger leg and the key used in the json file to specify the efficiency or scale factor value.
You will end up with a structure like this for `flags_with_triggernames`, depending on which lepton working points you specified:
```yaml
Single Lepton Trigger:
  SingleMu: 
    IsoMu24:
      jsonkey: "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"
      jsonleg1: "IsoMu24"
      cset_mc_filepath_leg1: "/home/cms-jovyan/ua-cms-higgscharm/analysis/data/muon_HWW_hlt/2022preEE/MC_eff/IsoMu24_pt_eta_efficiency.json"
      cset_data_filepath_leg1: "/home/cms-jovyan/ua-cms-higgscharm/analysis/data/muon_HWW_hlt/2022preEE/Data_eff/IsoMu24_pt_eta_efficiency.json"
      cset_sf_filepath_leg1: ""
  SingleEle:
    Ele30_WPTight_Gsf:
      jsonkey: "TriggerEff"
      jsonleg1: "Ele30"
Double Lepton Trigger:
  MuEle_HWW:
    Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ:  
      jsonkey: "TriggerEff"
      jsonleg1: "Mu12"
      cset_mc_filepath_leg1: "/home/cms-jovyan/ua-cms-higgscharm/analysis/data/muon_HWW_hlt/2022preEE/MC_eff/Mu23_Ele12_leg1_pt_eta_efficiency.json"
      cset_data_filepath_leg1: "/home/cms-jovyan/ua-cms-higgscharm/analysis/data/muon_HWW_hlt/2022preEE/Data_eff/Mu23_Ele12_leg1_pt_eta_efficiency.json"
      cset_sf_filepath_leg1: ""
      jsonleg2: "Ele23"
      cset_mc_filepath_leg2: "/home/cms-jovyan/ua-cms-higgscharm/analysis/data/muon_HWW_hlt/2022preEE/MC_eff/Mu23_Ele12_leg2_pt_eta_efficiency.json"
      cset_data_filepath_leg2: "/home/cms-jovyan/ua-cms-higgscharm/analysis/data/muon_HWW_hlt/2022preEE/Data_eff/Mu23_Ele12_leg2_pt_eta_efficiency.json"
      cset_sf_filepath_leg2: ""
    Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL:  
      jsonkey: "TriggerEff"
      jsonleg1: "Mu23"
      cset_mc_filepath_leg1: "/home/cms-jovyan/ua-cms-higgscharm/analysis/data/muon_HWW_hlt/2022preEE/MC_eff/Mu12_Ele23_leg1_pt_eta_efficiency.json"
      cset_data_filepath_leg1: "/home/cms-jovyan/ua-cms-higgscharm/analysis/data/muon_HWW_hlt/2022preEE/Data_eff/Mu12_Ele23_leg1_pt_eta_efficiency.json"
      cset_sf_filepath_leg1: ""
      jsonleg2: "Ele12"
      cset_mc_filepath_leg2: "/home/cms-jovyan/ua-cms-higgscharm/analysis/data/muon_HWW_hlt/2022preEE/MC_eff/Mu12_Ele23_leg2_pt_eta_efficiency.json"
      cset_data_filepath_leg2: "/home/cms-jovyan/ua-cms-higgscharm/analysis/data/muon_HWW_hlt/2022preEE/Data_eff/Mu12_Ele23_leg2_pt_eta_efficiency.json"
      cset_sf_filepath_leg2: ""
```

We distinguish two types for the triggers: `Single Lepton Trigger` and `Double Lepton Trigger`. 
The trigger flags are obtained from the `hlt_paths` in the workflow config file using `get_config_trigger_flags()` in `trigger.py`. In there is specified which trigger flags correspond to which trigger type. For trigger flags not falling under one of the two types, a warning is printed. These will not be added to the dictionary, so if necessary, they need to be added still. The same function will also check if a trigger flag was specified for each of the Data samples listed in that function and generate a warning otherwise (it will still continue normally). The Data samples are set by default as `MuonEG`, `Muon` and `EGamma`.  

For each trigger flag, the corresponding trigger names are added. These are obtained from `analysis/selections/trigger_flags.yaml` for each year (campaigns 2022preEE and 2022postEE are both taken as year 2022 and 2023preBPix and 2023postBPix are both taken as year 2023).

To get the trigger legs (`jsonlegs`), we split the name of the trigger at each underscore `_` and look for parts starting with 'Mu', 'Ele' or 'IsoMu'. These names are later used to get the correct efficiency json file (depending on the lepton working points). The `jsonkey` is just listed in `get_available_SF_names_DEN()` in `trigger.py` for each NanoAOD version and working point (for now only for muons). This is the key you need in the efficiency json file to evaluate the effciency or scale factor value. For the efficienies derived by the HWW group, this is always "TriggerEff". In other files, this is usually a name specifying the numerator (NUM) and denominator (DEN) used when deriving the efficiencies with the tag and probe method.

Beware that the working points for which the trigger scale factors are applied still need to be specified separately from the definition of the muon and electron object in the workflow config file.

Now the trick is to get the correct the correct efficiency file attached to each leg and combine the legs, double and single lepton trigger efficiencies into one trigger efficiency data and one for mc and then one scale factor.