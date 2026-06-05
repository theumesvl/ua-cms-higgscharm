import copy
import yaml
import glob
import logging
import numpy as np
import pandas as pd
import dask.dataframe as dd
from pathlib import Path
from coffea.util import load, save
from coffea.processor import accumulate
from analysis.filesets.utils import get_dataset_config
from analysis.histograms import HistBuilder, fill_histogram
from analysis.postprocess.utils import (
    print_header,
    get_variations_keys,
    find_kin_and_axis,
    get_lumi_weight,
    accumulate_histograms,
    accumulate_metadata,
    get_process_dict,
    save_cutflows,
    accumulate_and_save_cutflows,
    add_mva_labels,
)


def fill_histograms_from_parquets(
    year, sample, categories, workflow_config, output_dir, skipmerging
):
    """Build and fill histograms from parquet files for a given sample"""
    dataset_config = get_dataset_config(year)
    histogram_config = workflow_config.histogram_config
    variables = list(histogram_config.axes.keys())
    histograms = HistBuilder(workflow_config).build_histogram()
    process_dict = get_process_dict(output_dir, year, categories)

    for category in categories:
        logging.info(f"Filling {sample} histograms from parquets for category {category}.")

        # merge sample parquets
        sample_df_folder = output_dir / category
        sample_df_folder.mkdir(parents=True, exist_ok=True)
            
        sample_df_file = sample_df_folder / f"{sample}.parquet"
        if sample_df_file.exists() and skipmerging:
            logging.info(f"{sample}.parquet already exists and skipmerging flag provided, so using that file for the filling.\nNo more parquet file merging to be done for this sample.")
            sample_df = pd.read_parquet(sample_df_file)       
        else:
            if sample_df_file.exists():
                logging.info(f"{sample}.parquet already exists, but no skipmerging flag provided. \nSo overwriting existing file and merging all parquet files in parquets_{sample}/{category} again.")
            else:
                logging.info(f"{sample}.parquet does not exist yet. Creating it by merging all parquet files in parquets_{sample}/{category}.")
            sample_parquets = glob.glob(
                f"{output_dir}/parquets_{sample}/{category}/*.parquet"
            )
            
            # make sure there are no empty parquet files (could be from low statistics). Skip empty ones
            valid_parquets = []
            sumlengths = 0 # sum of lengths of dataframes to calculate averages
            averagelength = 0

            for parquet_file in sample_parquets:
                try:
                    df = pd.read_parquet(parquet_file)
            
                    # skip empty parquet files
                    if len(df) == 0:
                        logging.info(f"Skipping empty parquet file: {parquet_file}. Probably empty due to low statistics in this selection. \nIf the average length (will be printed below) is high, check for other possible issues.")
                        continue
                    sumlengths += len(df)
                    valid_parquets.append(parquet_file)
            
                except Exception as e:
                    logging.warning(f"Could not read {parquet_file}: {e}")
           
            if len(valid_parquets) == 0:
                logging.warning(
                    f"No valid parquet files found for sample {sample} in category {category}. They are all empty for each partition. \nPlease check that the parquet outputs were actually produced for this sample and category."
                )
                sample_df = pd.read_parquet(sample_parquets[0])
                logging.info(f"Because all parquets were empty, just copied the first partition parquet file to {sample}.parquet.")
            else:
                averagelength = sumlengths / len(valid_parquets)
                logging.info(f"The average number of events in each partition for sample {sample} is {averagelength}.")
                sample_df = dd.read_parquet(
                    valid_parquets, engine="pyarrow", calculate_divisions=False
                ).compute()
            sample_df = sample_df.replace({None: np.nan})
            sample_df.to_parquet(
                f"{output_dir}/{category}/{sample}.parquet", engine="pyarrow", index=False
            )
            logging.info(f"Saved {sample}.parquet in directory '{category}'.")

        # build variables map
        variables_map = {}
        variables_mask_map = {}
        for variable in variables:
            if variable in sample_df.columns:
                variable_array = sample_df[variable].values
            else:
                logging.info(f"Could not find variable {variable} for sample {sample}")
            if variable_array.dtype.type is np.object_:
                variable_array = np.array(
                    [x if x is not None else np.nan for x in variable_array], dtype=bool
                )
            variables_map[variable] = variable_array

        # compute nominal weights
        partial_weights = list(
            set(
                [
                    w.replace("Up", "").replace("Down", "")
                    for w in sample_df.columns
                    if w.startswith("weight") and "nominal" not in w
                ]
            )
        )
        nominal_weights = sample_df[partial_weights].prod(axis=1).values
        if len(partial_weights) > 0:
            logging.info(
                f"weights: {[w.replace('weight_','') for w in partial_weights]}"
            )

        # fill nominal histograms
        sample_histograms = copy.deepcopy(histograms)
        fill_args = {
            "histograms": sample_histograms,
            "histogram_config": histogram_config,
            "variables_map": variables_map,
            "category": category,
            "flow": True,
            "weights": nominal_weights,
            "variation": "nominal",
        }
        fill_histogram(**fill_args)

        # fill syst variation histograms
        if dataset_config[sample]["era"] in ["mc", "signal"]:
            for syst in partial_weights:
                for variation in ["Up", "Down"]:
                    syst_name = f"{syst}{variation}"
                    if syst_name in sample_df.columns:
                        fill_args["weights"] = sample_df[syst_name].values
                        fill_args["variation"] = syst_name.replace("weight_", "")
                        fill_histogram(**fill_args)

    return sample_histograms


def save_histograms_by_sample(
    grouped_outputs,
    sample,
    year,
    output_dir,
    categories,
    workflow_config,
    nocutflow,
    output_format,
    skipmerging,
):
    """Accumulate, scale, and save histograms for a single sample"""
    print_header(f"Processing {sample} outputs (SAMPLE)")

    # get histograms
    if output_format == "coffea":
        histograms = accumulate_histograms(grouped_outputs, sample)
    elif output_format == "parquet":
        histograms = fill_histograms_from_parquets(
            year, sample, categories, workflow_config, output_dir, skipmerging
        )
    else:
        raise ValueError(f"Unsupported output_format: {output_format}")

    # accumulate metadata and compute lumi weight
    metadata = accumulate_metadata(grouped_outputs, sample)
    weight = get_lumi_weight(year, sample, metadata)

    # scale histograms by lumi-xsec weight
    scaled_histograms = {
        variable: histograms[variable] * weight for variable in histograms
    }
    save(scaled_histograms, Path(output_dir) / f"{sample}.coffea")

    # save cutflows if requested
    if not nocutflow:
        save_cutflows(metadata, categories, sample, weight, output_dir)


def save_histograms_by_process(
    process: str,
    output_dir: str,
    process_samples_map: dict,
    categories: list,
    nocutflow: bool,
    output_format: str,
    add_mva_labels_flag: bool = False,
):
    """Accumulate and save all outputs for a given physics process.

    Parameters
    ----------
    process : str
        Process name (e.g., "tt", "DY+Jets")
    output_dir : str
        Output directory path
    process_samples_map : dict
        Mapping of process names to sample lists
    categories : list
        List of category names
    nocutflow : bool
        If True, skip saving cutflow tables
    output_format : str
        Output format ("coffea" or "parquet")
    add_mva_labels_flag : bool
        If True, add MVA training labels to parquet files (hww workflow only)
    """
    print_header(f"Processing {process} outputs (PROCESS)")

    # accumulate and save all histograms into a single dictionary
    coffea_files = []
    for sample in process_samples_map[process]:
        coffea_files += glob.glob(f"{output_dir}/{sample}.coffea", recursive=True)

    logging.info(f"Accumulating histograms for process {process}")
    hist_to_accumulate = [load(f) for f in coffea_files]
    output_histograms = {process: accumulate(hist_to_accumulate)}
    save(output_histograms, Path(output_dir) / f"combined_{process}.coffea")

    # accumulate and save all parquets into a single parquet file per process per category
    if output_format == "parquet":
        for category in categories:
            logging.info(f"Accumulating sample parquets for process {process} and category {category}") 
            parquet_files = []
            for sample in process_samples_map[process]:
                parquet_files.append(f"{output_dir}/{category}/{sample}.parquet")

            # failsafe in case some parquet files are completely empty and do not even have the column headers with all variables
            valid_parquets = []
            sumlengths = 0 # sum of lengths of dataframes to calculate averages
            averagelength = 0
            for parquet_file in parquet_files:
                try:
                    df = pd.read_parquet(parquet_file)
            
                    # skip empty parquet files
                    if len(df) == 0:
                        logging.info(f"Skipping empty parquet file: {parquet_file}. Probably empty due to low statistics in this selection. \nThis issue has been raised already in the merging step (if skipmerging flag was not given).")
                        continue
                    sumlengths += len(df)
                    valid_parquets.append(parquet_file)
            
                except Exception as e:
                    logging.warning(f"Could not read {parquet_file}: {e}")
           
            if len(valid_parquets) == 0:
                logging.warning(
                    f"No valid parquet files found for process {process} in category {category}. They are all empty for each sample. \nPlease check that the parquet outputs were actually produced. It is uncommon that all samples (if there are multiple) for a process have no events that pass the selection."
                )
                process_df = pd.read_parquet(parquet_files[0])
                logging.info(f"Because all parquets were empty, just copied the first sample parquet file to combined_{process}.parquet.")
            else:
                averagelength = sumlengths / len(valid_parquets)
                logging.info(f"The average number of events in each sample for process {process} is {averagelength}.")
                process_df = dd.read_parquet(
                    valid_parquets, engine="pyarrow", calculate_divisions=False
                ).compute()
            process_df = process_df.replace({None: np.nan})
    
            print("List of parquet files merged: ", parquet_files)
    
            # Add MVA training labels if requested (hww workflow)
            if add_mva_labels_flag:
                logging.info(f"Adding MVA labels for process {process}")
                process_df = add_mva_labels(process_df, process)
     
            process_df.to_parquet(Path(output_dir) / category / f"combined_{process}.parquet") 
            logging.info(f"Saved parquet file for process {process} in {category}/combined_{process}.parquet") 

    # accumulate and save cutflows if requested
    if not nocutflow:
        accumulate_and_save_cutflows(
            process, process_samples_map, output_dir, categories
        )
