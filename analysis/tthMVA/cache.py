from pathlib import Path
import hashlib

import awkward as ak
import numpy as np
import pandas as pd

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

        identifier = (
            f"{dataset}::{filename}"
        )

        file_hash = hashlib.sha1(
            identifier.encode()
        ).hexdigest()

        cache_dir = (
            cache_base
            / self.model_name
        )

        cache_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        return (
            cache_dir
            / f"{file_hash}.parquet"
        )

    # ------------------------------------------------------------
    # Calculate scores
    # ------------------------------------------------------------

    def _calculate_scores(
        self,
        events,
        muon_indices,
    ):

        values = []

        for event_idx, mu_idx in muon_indices:

            mu = events.Muon[
                event_idx,
                mu_idx,
            ]

            jet_idx = int(mu.jetIdx)

            if jet_idx >= 0:

                jet_btag = float(
                    events.Jet[
                        event_idx,
                        jet_idx,
                    ].btagDeepFlavB
                )

            else:

                jet_btag = 0.0

            x = [
                float(mu.pt),
                float(mu.eta),
                float(mu.pfRelIso03_all),
                float(mu.miniPFRelIso_chg),

                float(
                    mu.miniPFRelIso_all
                    - mu.miniPFRelIso_chg
                ),

                float(mu.jetNDauCharged),
                float(mu.jetPtRelv2),

                jet_btag,

                min(
                    1.0 / (
                        1.0
                        + float(mu.jetRelIso)
                    ),
                    1.5,
                ),

                float(mu.sip3d),

                np.log(
                    abs(float(mu.dxy))
                ),

                np.log(
                    abs(float(mu.dz))
                ),

                float(mu.segmentComp),
            ]

            values.append(x)

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
                    "tthMVA",
                ]
            )

        # --------------------------------------------------------
        # Build cache lookup
        # --------------------------------------------------------

        cached_keys = set(
            zip(
                cache["run"].astype(
                    np.int64
                ),
                cache["lumi"].astype(
                    np.int64
                ),
                cache["event"].astype(
                    np.int64
                ),
                cache["muon"].astype(
                    np.int64
                ),
            )
        )

        # --------------------------------------------------------
        # Find missing muons
        # --------------------------------------------------------

        missing_keys = []
        missing_indices = []

        for ievt in range(len(events)):

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

                    missing_keys.append(key)

                    missing_indices.append(
                        (
                            ievt,
                            imu,
                        )
                    )

        # --------------------------------------------------------
        # Calculate missing scores
        # --------------------------------------------------------
        
        if missing_indices:

            scores = self._calculate_scores(
                events,
                missing_indices,
            )

            new_cache = pd.DataFrame(
                {
                    "run": [
                        key[0]
                        for key in missing_keys
                    ],
                    "lumi": [
                        key[1]
                        for key in missing_keys
                    ],
                    "event": [
                        key[2]
                        for key in missing_keys
                    ],
                    "muon": [
                        key[3]
                        for key in missing_keys
                    ],
                    "tthMVA": scores,
                }
            )

            cache = pd.concat(
                [
                    cache,
                    new_cache,
                ],
                ignore_index=True,
            )

            cache.to_parquet(
                cache_file,
                index=False,
            )

        # --------------------------------------------------------
        # Build lookup
        # --------------------------------------------------------

        lookup = {
            (
                int(row.run),
                int(row.lumi),
                int(row.event),
                int(row.muon),
            ): np.float32(row.tthMVA)

            for row in cache.itertuples(
                index=False
            )
        }

        # --------------------------------------------------------
        # Construct Awkward array
        # --------------------------------------------------------

        output = []

        for ievt in range(len(events)):

            scores = []

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

            output.append(scores)

        return ak.Array(output)
