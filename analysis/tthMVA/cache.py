from pathlib import Path
import hashlib

import awkward as ak
import numpy as np
import pandas as pd
import os
import uuid
from filelock import FileLock

from analysis.tthMVA.evaluator import TMVAGradBDT


cwd = str(Path.cwd())
if cwd.startswith('/eos'):
    userpath = cwd.split('user')[-1] # part of the path with first letter and username
    userpath = userpath.split('/', 3)[3]
    userpath = userpath.split('/', 2)[0] + '/' + userpath.split('/', 2)[1]
elif cwd.startswith('/afs'):
    userpath = cwd.split('user')[-1].strip('/') # part of the path with first letter and username
    userpath = userpath.split('/', 2)[0] + '/' + userpath.split('/', 2)[1]
    print(userpath)
elif cwd.startswith('/home'):
    print("Beware that when running from coffea-casa, you have to manually put in the target directory's 'userpath' for writing the tth MVA scores in analysis/tthMVA/cache.py")
    userpath = 't/tvanlaer'
else:
    raise ValueError(
        "You ran from an unexpected directory. Please run from SWAN (/eos/...), coffea-casa (/home/...) or lxplus (/.../user/...)."
        "Beware that when running from coffea-casa, you have to manually put in the target directory's 'userpath' for writing the tth MVA scores in analysis/tthMVA/cache.py"
    )
cache_base = Path(f"/eos/user/{userpath}/higgscharm/analysis/data/tthMVA_2022-2023_retrained/tthMVA_cache")
#cache_base = Path("/home/cms-jovyan/ua-cms-higgscharm/analysis/data/tthMVA_2022-2023_retrained/tthMVA_cache")
print(f"Your current path to write the tthMVA score cache to is '{cache_base}'")


class TTHMVACache:

    def __init__(self, model_file):

        self.model_file = str(model_file)

        self.model_name = Path(model_file).stem

        self.evaluator = TMVAGradBDT(
            model_file
        )

    # ------------------------------------------------------------
    # Cache filename
    # ------------------------------------------------------------

    def _cache_file(
        self,
        dataset,
        filename,
    ):

        print("dataset:", dataset)
        print("filename:", filename)

        filename = str(filename)

        file_id = Path(filename).name

        cache_base_sample = (
            cache_base
            / dataset
        )

        cache_base_sample.mkdir(parents=True, exist_ok=True)

        print(f"parquet cache file path: {cache_base_sample}/{file_id}.parquet")

        return (
            cache_base_sample
            / f"{file_id}.parquet"
        )

    # ------------------------------------------------------------
    # Calculate scores
    # ------------------------------------------------------------

    def _calculate_scores(
        self,
        events,
        indices,
        particle,
    ):

        values = []

        if particle == "muons":
            for event_idx, mu_idx in indices:

                mu = events.Muon[event_idx, mu_idx]

                # ----------------------------------------------------
                # Variables directly from NanoAOD
                # ----------------------------------------------------

                pt = float(mu.pt)
                eta = float(mu.eta)
                pfRelIso03_all = float(mu.pfRelIso03_all)
                miniPFRelIso_chg = float(mu.miniPFRelIso_chg)
                jetNDauCharged = float(mu.jetNDauCharged)
                jetPtRelv2 = float(mu.jetPtRelv2)
                sip3d = float(mu.sip3d)
                dxy = float(mu.dxy)
                dz = float(mu.dz)
                segmentComp = float(mu.segmentComp)

                # ----------------------------------------------------
                # Derived variables
                # ----------------------------------------------------
                miniRelIsoNeutral = (float(mu.miniPFRelIso_all)- miniPFRelIso_chg)

                jet_idx = int(mu.jetIdx)

                if jet_idx >= 0:
                    jetBTagDeepFlavB = float(events.Jet[event_idx,jet_idx,].btagDeepFlavB)
                else:
                    jetBTagDeepFlavB = 0.0

                jetPtRatio = min(1.0 / (1.0 + float(mu.jetRelIso)),1.5)

                log_dxy = np.log(abs(dxy))
                log_dz = np.log(abs(dz))

                # ----------------------------------------------------
                # MVA input vector
                # ----------------------------------------------------

                x = [
                    pt,
                    eta,
                    pfRelIso03_all,
                    miniPFRelIso_chg,
                    miniRelIsoNeutral,
                    jetNDauCharged,
                    jetPtRelv2,
                    jetBTagDeepFlavB,
                    jetPtRatio,
                    sip3d,
                    log_dxy,
                    log_dz,
                    segmentComp,
                ]

                values.append(x)

        elif particle == "electrons":
            for event_idx, ele_idx in indices:

                ele = events.Electron[event_idx, ele_idx]

                # ----------------------------------------------------
                # Variables directly from NanoAOD
                # ----------------------------------------------------

                pt = float(ele.pt)
                eta = float(ele.eta)
                pfRelIso03_all = float(ele.pfRelIso03_all)
                miniPFRelIso_chg = float(ele.miniPFRelIso_chg)
                jetNDauCharged = float(ele.jetNDauCharged)
                jetPtRelv2 = float(ele.jetPtRelv2)
                sip3d = float(ele.sip3d)
                dxy = float(ele.dxy)
                dz = float(ele.dz)
                mvaIso = float(ele.mvaIso)

                # ----------------------------------------------------
                # Derived variables
                # ----------------------------------------------------
                miniRelIsoNeutral = (float(ele.miniPFRelIso_all)- miniPFRelIso_chg)

                jet_idx = int(ele.jetIdx)

                if jet_idx >= 0:
                    jetBTagDeepFlavB = float(events.Jet[event_idx,jet_idx,].btagDeepFlavB)
                else:
                    jetBTagDeepFlavB = 0.0

                jetPtRatio = min(1.0 / (1.0 + float(ele.jetRelIso)),1.5)

                log_dxy = np.log(abs(dxy))
                log_dz = np.log(abs(dz))

                # ----------------------------------------------------
                # MVA input vector
                # ----------------------------------------------------

                x = [
                    pt,
                    eta,
                    pfRelIso03_all,
                    miniPFRelIso_chg,
                    miniRelIsoNeutral,
                    jetNDauCharged,
                    jetPtRelv2,
                    jetBTagDeepFlavB,
                    jetPtRatio,
                    sip3d,
                    log_dxy,
                    log_dz,
                    mvaIso,
                ]

                values.append(x)
        else:
            raise ValueError(
                f"Input partcile is '{particle}', but the only valid options are 'muons' or 'electrons'"
                "Check the object_selections.py and the hww config file to debug."
            )

        if not values:
            return np.empty(
                0,
                dtype=np.float32,
            )

        values = np.asarray(
            values,
            dtype=np.float64,
        )

        return self.evaluator.evaluate_array(
            values
        )


    # ------------------------------------------------------------
    # Get scores
    # ------------------------------------------------------------

    def get_scores(
        self,
        events,
        dataset,
        filename,
        particle,
    ):

        cache_file = self._cache_file(
            dataset,
            filename,
        )

        # --------------------------------------------------------
        # Event identifiers
        # --------------------------------------------------------

        run = np.asarray(
            ak.to_numpy(events.run),
            dtype=np.int64,
        )

        lumi = np.asarray(
            ak.to_numpy(events.luminosityBlock),
            dtype=np.int64,
        )

        event = np.asarray(
            ak.to_numpy(events.event),
            dtype=np.int64,
        )

        # --------------------------------------------------------
        # Read cache
        # --------------------------------------------------------

        if cache_file.exists():

            cache = pd.read_parquet(
                cache_file
            )

        else:

            cache = pd.DataFrame(
                columns=[
                    "run",
                    "lumi",
                    "event",
                    "muon",
                    "muon_tthMVA",
                    "electron",
                    "electron_tthMVA",
                ]
            )

        # --------------------------------------------------------
        # Build cache lookup
        # --------------------------------------------------------
        if particle == "muons": 
            muon_cache = cache[cache["muon"].notna()]
            cached_keys = set(
                zip(
                    muon_cache["run"].astype(np.int64),
                    muon_cache["lumi"].astype(np.int64),
                    muon_cache["event"].astype(np.int64),
                    muon_cache["muon"].astype(np.int64),
                )
            )
        elif particle == "electrons":  
            electron_cache = cache[cache["electron"].notna()]
            
            cached_keys = set(
                zip(
                    electron_cache["run"].astype(np.int64),
                    electron_cache["lumi"].astype(np.int64),
                    electron_cache["event"].astype(np.int64),
                    electron_cache["electron"].astype(np.int64),
                )
            )
        else:
            raise ValueError(
                f"Input particle is '{particle}', but the only valid options are 'muons' or 'electrons'"
                "Check the object_selections.py and the hww config file to debug."
            )

        # --------------------------------------------------------
        # Find missing muons and electrons
        # --------------------------------------------------------

        missing_muon_keys = []
        missing_muon_indices = []
        missing_electron_keys = []
        missing_electron_indices = []

        for ievt in range(len(events)):
            
            if particle == "muons":
                for imu in range(
                    len(events.Muon[ievt])
                ):

                    key = (
                        int(run[ievt]),
                        int(lumi[ievt]),
                        int(event[ievt]),
                        int(imu),
                    )

                    if key not in cached_keys:

                        missing_muon_keys.append(key)

                        missing_muon_indices.append(
                            (
                                ievt,
                                imu,
                            )
                        )
            
            elif particle == "electrons":
                for iele in range(
                    len(events.Electron[ievt])
                ):

                    key = (
                        int(run[ievt]),
                        int(lumi[ievt]),
                        int(event[ievt]),
                        int(iele),
                    )

                    if key not in cached_keys:

                        missing_electron_keys.append(key)

                        missing_electron_indices.append(
                            (
                                ievt,
                                iele,
                            )
                        )
            
            else:
                raise ValueError(
                    f"Input particle is '{particle}', but the only valid options are 'muons' or 'electrons'"
                    "Check the object_selections.py and the hww config file to debug."
                )

        # --------------------------------------------------------
        # Calculate missing scores
        # -------------------------------------------------------- 
        new_rows = []
        
        if missing_muon_indices:
            print("Calculating missing muon MVA scores.")
            scores = self._calculate_scores(
                events,
                missing_muon_indices,
                "muons",
            )

            muon_cache = pd.DataFrame(
                {
                    "run": [
                        key[0]
                        for key in missing_muon_keys
                    ],
                    "lumi": [
                        key[1]
                        for key in missing_muon_keys
                    ],
                    "event": [
                        key[2]
                        for key in missing_muon_keys
                    ],
                    "muon": [
                        key[3]
                        for key in missing_muon_keys
                    ],
                    "muon_tthMVA": scores,
                    "electron": [
                        None
                        for key in missing_muon_keys
                    ],
                    "electron_tthMVA": [ 
                        None
                        for key in missing_muon_keys
                    ],
                }
            )
            new_rows.append(muon_cache)

        if missing_electron_indices:
            print("Calculating missing electron MVA scores.")
            scores = self._calculate_scores(
                events,
                missing_electron_indices,
                "electrons",
            )

            electron_cache = pd.DataFrame(
                {
                    "run": [
                        key[0]
                        for key in missing_electron_keys
                    ],
                    "lumi": [
                        key[1]
                        for key in missing_electron_keys
                    ],
                    "event": [
                        key[2]
                        for key in missing_electron_keys
                    ],
                    "muon": [
                        None
                        for key in missing_electron_keys
                    ],
                    "muon_tthMVA": [
                        None
                        for key in missing_electron_keys
                    ],
                    "electron": [
                        key[3]
                        for key in missing_electron_keys
                    ],
                    "electron_tthMVA": scores,
                }
            )
            new_rows.append(electron_cache)
        
        # --------------------------------------------------------
        # Add newly calculated scores to cache
        # --------------------------------------------------------
        if new_rows:

            # Combine all scores calculated by this worker
            new_cache = pd.concat(
                new_rows,
                ignore_index=True,
            )

            # ----------------------------------------------------
            # Write worker-specific temporary cache
            # ----------------------------------------------------
            tmp_file = cache_file.with_name(
                f"{cache_file.stem}.{os.getpid()}.{uuid.uuid4().hex}.tmp.parquet"
            )

            print(  
                f"Writing {len(new_cache)} new MVA scores "
                f"to temporary cache: {tmp_file}"
            )

            new_cache.to_parquet(
                tmp_file,
                index=False,
            )

            # ----------------------------------------------------
            # Merge temporary cache into shared cache
            # ----------------------------------------------------
            lock_file = Path(str(cache_file) + ".lock")

            print(f"Waiting for cache lock: {lock_file}")

            with FileLock(str(lock_file)):

                print(f"Acquired cache lock: {lock_file}")

                # Another worker may have updated the cache while
                # this worker was calculating its MVA scores.
                if cache_file.exists():

                    cache = pd.read_parquet(
                        cache_file
                    )

                else:

                    cache = pd.DataFrame(
                        columns=[
                            "run",
                            "lumi",
                            "event",
                            "muon",
                            "muon_tthMVA",
                            "electron",
                            "electron_tthMVA",
                        ]
                    )

                # Read this worker's newly calculated scores
                tmp = pd.read_parquet(
                    tmp_file
                )

                print(
                    f"Merging {len(tmp)} new rows "
                    f"into cache containing {len(cache)} rows."
                )

                if cache.empty:
                    cache = tmp.copy()
                else:
                    cache = pd.concat(
                        [cache, tmp],
                        ignore_index=True,
                    )

                # ------------------------------------------------
                # Remove duplicate lepton entries
                #
                # The combination of run/lumi/event/muon/electron
                # uniquely identifies a lepton.
                # ------------------------------------------------
                cache = cache.drop_duplicates(
                    subset=[
                        "run",
                        "lumi",
                        "event",
                        "muon",
                        "electron",
                    ],
                    keep="last",
                )

                print(
                    f"Saving merged cache with {len(cache)} rows."
                )

                cache.to_parquet(
                    cache_file,
                    index=False,
                )

                # Temporary file is no longer needed
                tmp_file.unlink()

                print(
                    f"Released cache lock: {lock_file}"
                )


        # --------------------------------------------------------
        # Build lookup
        # --------------------------------------------------------
        
        if particle == "muons":
            lookup = {
                (
                    int(row.run),
                    int(row.lumi),
                    int(row.event),
                    int(row.muon),
                ): np.float32(row.muon_tthMVA)
                for row in cache.itertuples(index=False)
                if pd.notna(row.muon)
            }

        elif particle == "electrons":
            lookup = {
                (
                    int(row.run),
                    int(row.lumi),
                    int(row.event),
                    int(row.electron),
                ): np.float32(row.electron_tthMVA)
                for row in cache.itertuples(index=False)
                if pd.notna(row.electron)
            }

        else:
            raise ValueError(
                f"Input partcile is '{particle}', but the only valid options are 'muons' or 'electrons'"
                "Check the object_selections.py and the hww config file to debug."
            )

        # --------------------------------------------------------
        # Construct Awkward array
        # --------------------------------------------------------

        output = []

        for ievt in range(len(events)):

            scores = []
            if particle == "muons":
                for imu in range(
                    len(events.Muon[ievt])
                ):

                    key = (
                        int(run[ievt]),
                        int(lumi[ievt]),
                        int(event[ievt]),
                        int(imu),
                    )

                    scores.append(
                        lookup[key]
                    )
            elif particle == "electrons":
                for iele in range(
                    len(events.Electron[ievt])
                ):

                    key = (
                        int(run[ievt]),
                        int(lumi[ievt]),
                        int(event[ievt]),
                        int(iele),
                    )

                    scores.append(
                        lookup[key]
                    )
            else:
                raise ValueError(
                    f"Input partcile is '{particle}', but the only valid options are 'muons' or 'electrons'"
                    "Check the object_selections.py and the hww config file to debug."
                )

            output.append(scores)

        return ak.Array(output)
