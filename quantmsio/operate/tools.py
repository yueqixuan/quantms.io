import os
import re
from collections import defaultdict
from pathlib import Path
from typing import Generator, List, Optional, Union
import logging

import ahocorasick
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from Bio import SeqIO
from pyopenms import FASTAFile

from quantmsio.core.common import FEATURE_SCHEMA, IBAQ_SCHEMA, IBAQ_USECOLS, PSM_SCHEMA
from quantmsio.core.openms import OpenMSHandler
from quantmsio.core.sdrf import SDRFHandler
from quantmsio.operate.query import Query, map_spectrum_mz
from quantmsio.utils.file_utils import (
    close_file,
    load_de_or_ae,
    save_file,
    save_slice_file,
)
from quantmsio.utils.pride_utils import get_unanimous_name


def init_save_info(parquet_path: str) -> tuple[dict, None, str]:
    pqwriters: dict = {}
    pqwriter_no_part = None
    filename = os.path.basename(parquet_path)
    return pqwriters, pqwriter_no_part, filename


def generate_psms_of_spectrum(
    parquet_path: str,
    mzml_directory: str,
    output_folder: str,
    file_num: int,
    partitions: Optional[List[str]] = None,
):
    """
    parquet_path: parquet file path
    mzml_directory: mzml file directory path
    """
    pqwriters, pqwriter_no_part, filename = init_save_info(parquet_path)
    p = Query(parquet_path)
    for refs, table in p.iter_file(file_num=file_num):
        mzml_handlers = {ref: OpenMSHandler() for ref in refs}
        table[["number_peaks", "mz_array", "intensity_array"]] = table[
            ["reference_file_name", "scan"]
        ].apply(
            lambda x: map_spectrum_mz(
                x["reference_file_name"],
                x["scan"],
                mzml_handlers,
                mzml_directory,
            ),
            axis=1,
            result_type="expand",
        )
        # Initialize new fragment annotation arrays as None for backward compatibility
        table["charge_array"] = None
        table["ion_type_array"] = None
        table["ion_mobility_array"] = None
        pqwriters, pqwriter_no_part = save_parquet_file(
            partitions,
            table,
            output_folder,
            filename,
            pqwriters,
            pqwriter_no_part,
            PSM_SCHEMA,
        )
    close_file(pqwriters=pqwriters, pqwriter=pqwriter_no_part)


def save_parquet_file(
    partitions: list,
    table: pd.DataFrame,
    output_folder: str,
    filename: str,
    pqwriters: dict | None = None,
    pqwriter_no_part: pq.ParquetWriter | None = None,
    schema: pa.Schema = FEATURE_SCHEMA,
):

    if pqwriters is None:
        pqwriters: dict = {}
    if partitions and len(partitions) > 0:
        for key, df in table.groupby(partitions):
            parquet_table = pa.Table.from_pandas(df, schema=schema)
            pqwriters = save_slice_file(
                parquet_table, pqwriters, output_folder, key, filename
            )
        return pqwriters, pqwriter_no_part
    else:
        parquet_table = pa.Table.from_pandas(table, schema=schema)
        pqwriter_no_part = save_file(
            parquet_table=parquet_table,
            pqwriter=pqwriter_no_part,
            output_folder=output_folder,
            filename=filename,
        )
        return pqwriters, pqwriter_no_part


def generate_feature_of_gene(
    parquet_path: str,
    fasta: str,
    output_folder: str,
    file_num: int,
    partitions: list | None = None,
    species: str = "human",
):
    pqwriters, pqwriter_no_part, filename = init_save_info(parquet_path)
    p = Query(parquet_path)
    map_gene_names = p.get_protein_to_gene_map(fasta)
    for _, table in p.iter_file(file_num=file_num):
        table = p.inject_gene_msg(table, map_gene_names, species)
        pqwriters, pqwriter_no_part = save_parquet_file(
            partitions, table, output_folder, filename, pqwriters, pqwriter_no_part
        )
    close_file(pqwriters, pqwriter_no_part)


def map_protein_for_tsv(
    path: str, fasta: str, output_path: str, map_parameter: str
) -> None:
    """
    according fasta database to map the proteins accessions to uniprot names.
    :param path: de_path or ae_path
    :param fasta: Reference fasta database
    :param output_path: output file path
    :param map_parameter: map_protein_name or map_protein_accession
    retrun: None
    """
    map_protein_names = defaultdict(set)
    if map_parameter == "map_protein_name":
        for seq_record in SeqIO.parse(fasta, "fasta"):
            accession, name = seq_record.id.split("|")[1:]
            map_protein_names[seq_record.id].add(name)
            map_protein_names[accession].add(name)
            map_protein_names[name].add(name)
    else:
        for seq_record in SeqIO.parse(fasta, "fasta"):
            accession, name = seq_record.id.split("|")[1:]
            map_protein_names[seq_record.id].add(accession)
            map_protein_names[accession].add(accession)
            map_protein_names[name].add(accession)
    df, content = load_de_or_ae(path)
    df["protein"] = df["protein"].apply(
        lambda x: get_unanimous_name(x, map_protein_names)
    )
    content += df.columns.str.cat(sep="\t") + "\n"
    for _, row in df.iterrows():
        content += "\t".join(map(str, row)).strip() + "\n"
    with open(output_path, "w", encoding="utf8") as f:
        f.write(content)


# Before Code:
def get_ahocorasick(mods_dict: dict) -> ahocorasick.Automaton:
    automaton = ahocorasick.Automaton()
    for key in mods_dict.keys():
        key = "(" + key + ")"
        automaton.add_word(key, key)
    automaton.make_automaton()
    return automaton


def get_modification_details(
    seq: str,
    mods_dict: dict,
    automaton: ahocorasick.Automaton,
    select_mods: list | None = None,
) -> tuple[str, list]:
    if "(" not in seq:
        return seq, []
    total: int = 0
    modifications: list = []
    modification_details: list = []
    peptidoform: str = ""
    pre: int = 0
    for item in automaton.iter(seq):
        total += len(item[1])
        modification = item[1][1:-1]
        position = item[0] - total + 1
        name = re.search(r"([^ ]+)\s?", modification)
        name = name.group(1)
        if position == 0:
            peptidoform = f"[{name}]-{peptidoform}"
        elif item[0] + 1 == len(seq):
            peptidoform += seq[pre : item[0] - len(item[1]) + 1]
            peptidoform += f"-[{name}]"
        else:
            peptidoform += seq[pre : item[0] - len(item[1]) + 1]
            peptidoform += f"[{name}]"
        pre = item[0] + 1
        if modification in modifications:
            index = modifications.index(modification)
            modification_details[index]["fields"].append(
                {"position": position, "localization_probability": 1.0}
            )
        elif modification in select_mods:
            modifications.append(modification)
            modification_details.append(
                {
                    "modification_name": mods_dict[modification][0],
                    "fields": [{"position": position, "localization_probability": 1.0}],
                }
            )
    peptidoform += seq[pre:]
    return (peptidoform, modification_details)


def get_peptide_map(unique_peptides: list, fasta: str) -> defaultdict:
    peptide_map: defaultdict = defaultdict(list)
    automaton = ahocorasick.Automaton()
    for sequence in unique_peptides:
        automaton.add_word(sequence, sequence)
    automaton.make_automaton()

    fasta_proteins: list = list()
    FASTAFile().load(fasta, fasta_proteins)

    for entry in fasta_proteins:
        accession = entry.identifier.split("|")[1]
        for match in automaton.iter(entry.sequence):
            peptide = match[1]
            if accession not in peptide_map[peptide]:
                peptide_map[peptide].append(accession)
    return peptide_map


def map_peptide_to_protein(
    parquet_file: str,
    fasta: str,
    output_folder: str,
    filename: str,
    label: str = "feature",
) -> None:
    p = Query(parquet_file)
    unique_peptides = p.get_unique_peptides()
    peptide_map = get_peptide_map(unique_peptides, fasta)
    pqwriter = None
    for table in p.iter_chunk(batch_size=2000000):
        table["pg_accessions"] = table["sequence"].map(peptide_map)
        table = table[table["pg_accessions"].apply(lambda x: len(x) > 0)]
        table.loc[:, "unique"] = (
            table["pg_accessions"]
            .apply(lambda x: 0 if len(x) > 1 else 1)
            .astype(np.int32)
        )
        if label == "feature":
            parquet_table = pa.Table.from_pandas(table, schema=FEATURE_SCHEMA)
            pqwriter = save_file(parquet_table, pqwriter, output_folder, filename)
        else:
            parquet_table = pa.Table.from_pandas(table, schema=IBAQ_SCHEMA)
            pqwriter = save_file(parquet_table, pqwriter, output_folder, filename)
    close_file(None, pqwriter)


def get_field_schema(parquet_path: str) -> pa.Schema:
    schema = pq.read_schema(parquet_path)
    return schema


PROTEIN_ACCESSION = r"\|([^|]*)\|"


def get_protein_accession(proteins: Optional[str] = None) -> list:
    proteins = str(proteins)
    if "|" in proteins:
        return re.findall(PROTEIN_ACCESSION, proteins)
    else:
        return [re.split(r"[;,]", proteins)[0]]


def transform_ibaq(df: pd.DataFrame) -> pd.DataFrame:
    def transform(row: pd.Series) -> tuple:
        map_dict = row["intensities"]
        return map_dict["sample_accession"], map_dict["channel"], map_dict["intensity"]

    df = df.explode("intensities")

    # Check for NA in the "intensities" column
    if df["intensities"].isna().any():
        logging.warning(
            "[transform_ibaq]: The 'intensities' column contains NaN values."
        )
        df.dropna(subset=["intensities"], inplace=True)

    df.reset_index(drop=True, inplace=True)
    df[["sample_accession", "channel", "intensity"]] = df[["intensities"]].apply(
        transform, axis=1, result_type="expand"
    )
    df.drop(["intensities"], axis=1, inplace=True)
    return df


def genereate_ibaq_feature(
    sdrf_path: Union[Path, str], parquet_path: Union[Path, str]
) -> Generator[pa.Table, None, None]:
    sdrf_parser = SDRFHandler(sdrf_path)
    sdrf_df = sdrf_parser.transform_sdrf()
    sdrf_df.drop(columns=["sample_accession"], inplace=True)
    experiment_type = sdrf_parser.get_experiment_type_from_sdrf()
    p = Query(parquet_path)
    for _, df in p.iter_file(file_num=10, columns=IBAQ_USECOLS):
        df = transform_ibaq(df)
        if experiment_type != "LFQ":
            df = pd.merge(
                df,
                sdrf_df,
                left_on=["reference_file_name", "channel"],
                right_on=["reference_file_name", "channel"],
                how="left",
            )
        else:
            df = pd.merge(
                df,
                sdrf_df,
                left_on=["reference_file_name"],
                right_on=["reference_file_name"],
                how="left",
            )
            df["channel"] = "LFQ"
        df["fraction"] = df["fraction"].astype(str)
        feature = pa.Table.from_pandas(df, schema=IBAQ_SCHEMA)
        yield feature


def write_ibaq_feature(
    sdrf_path: Union[Path, str],
    parquet_path: Union[Path, str],
    output_path: Union[Path, str],
) -> None:

    logger = logging.getLogger("transform.ibaq")

    # Log input and output paths
    logger.info(f"Input SDRF file: {sdrf_path}")
    logger.info(f"Input feature file: {parquet_path}")
    logger.info(f"Output path: {output_path}")

    pqwriter = None
    for feature in genereate_ibaq_feature(sdrf_path, parquet_path):
        if not pqwriter:
            pqwriter = pq.ParquetWriter(output_path, feature.schema)
        pqwriter.write_table(feature)
    if pqwriter:
        pqwriter.close()

    logger.info("The iBAQ conversion has been completed.")
