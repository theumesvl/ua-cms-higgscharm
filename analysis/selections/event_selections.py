import yaml
import json
import math
from pathlib import Path
import numpy as np
import awkward as ak
import importlib.resources
from coffea.lumi_tools import LumiMask
from analysis.selections.trigger import trigger_mask, trigger_match_mask, zzto4l_trigger
from analysis.tthMVA.cache import TTHMVACache

def get_muon_tthMVA_cached(events, year, dataset, filename):
    """
    Get retrained Muon tthMVA scores using the persistent cache.

    For Run 3 years for which the retrained MVA is used, scores are
    calculated only for muons that are not already present in the cache.
    """

    if year not in [
        "2022preEE",
        "2022postEE",
        "2023preBPix",
        "2023postBPix",
    ]:
        raise ValueError(
            f"Input year is '{year}', but cached retrained tthMVA "
            "scores are only defined for 2022preEE, 2022postEE, "
            "2023preBPix and 2023postBPix."
        )

    weightfile = (
        Path.cwd()
        / "analysis/data/tthMVA_2022-2023_retrained"
        / "Muon-mvaTTH.2022EE.weights.json"
    )

    cache = TTHMVACache(weightfile)

    return cache.get_scores(
        events,
        dataset,
        filename,
    )

def get_lumi_mask(events, year):
    year_map = {
        "2016": "analysis/data/lumimask/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
        "2017": "analysis/data/lumimask/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
        "2018": "analysis/data/lumimask/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
        "2022": "analysis/data/lumimask/Cert_Collisions2022_355100_362760_Golden.txt",
        "2023": "analysis/data/lumimask/Cert_Collisions2023_366442_370790_Golden.txt",
        "2024": "analysis/data/lumimask/Cert_Collisions2024_378981_386951_Golden.txt",
    }
    for key in year_map:
        if year.startswith(key):
            goldenjson = year_map[key]
            break
    else:
        raise ValueError(f"Unrecognized year format: '{year}'")

    if hasattr(events, "genWeight"):
        lumi_mask = np.ones(len(events), dtype="bool")
    else:
        lumi_info = LumiMask(goldenjson)
        lumi_mask = lumi_info(events.run, events.luminosityBlock)
    return lumi_mask == 1


def get_zzto4l_trigger_mask(events, hlt_paths, dataset_key, year):
    return zzto4l_trigger(events, hlt_paths, dataset_key, year)


def get_trigger_mask(events, hlt_paths, dataset_key, year):
    return trigger_mask(events, hlt_paths, dataset_key, year)


def get_trigger_match_mask(events, hlt_paths, year, leptons):
    mask = trigger_match_mask(events, hlt_paths, year, leptons)
    return ak.sum(mask, axis=-1) > 0


def get_metfilters_mask(events, year):
    with importlib.resources.path("analysis.data.met", "metfilters.json") as path:
        with open(path, "r") as handle:
            metfilters = json.load(handle)[year]
    metfilters_mask = np.ones(len(events), dtype="bool")
    metfilterkey = "mc" if hasattr(events, "genWeight") else "data"
    for mf in metfilters[metfilterkey]:
        if mf in events.Flag.fields:
            metfilters_mask = metfilters_mask & events.Flag[mf]
    return metfilters_mask


def get_stitching_mask(events, dataset, dataset_key, ht_value):
    stitching_mask = np.ones(len(events), dtype="bool")
    if dataset.startswith(dataset_key):
        stitching_mask = events.LHE.HT < ht_value
    return stitching_mask
