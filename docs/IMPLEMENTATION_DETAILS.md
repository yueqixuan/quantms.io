# Implementation Details

This document provides detailed technical information about the internal workings of the quantms.io library, focusing on data processing pipelines and transformation logic.

## Optimized Protein Groups Processing

### Overview

The `MzTabProteinGroups` class provides high-performance protein quantification using DuckDB SQL aggregation, achieving 4-6x performance improvements over traditional pandas-based processing.

### Key Features

#### 1. Context Manager Resource Management

All resources are automatically cleaned up using context managers:

```python
from quantmsio.core.quantms.pg import MzTabProteinGroups

# Automatic cleanup of temporary files and database connections
with MzTabProteinGroups('data.mzTab.gz') as pg:
    result = pg.quantify_from_msstats_optimized(
        msstats_path='msstats.csv',
        sdrf_path='experiment.sdrf.tsv'
    )
```

#### 2. DuckDB SQL Optimization

The optimized pipeline uses SQL aggregation instead of pandas operations:

```python
def quantify_from_msstats_optimized(self, msstats_path, sdrf_path, ...):
    # Step 1: Create protein groups lookup table
    protein_groups_info = self._create_protein_groups_table_optimized()
    
    # Step 2: Initialize MsstatsIN with context manager
    with MsstatsIN(msstats_path, sdrf_path, ...) as msstats_in:
        # Step 3: Create optimized SQL views
        self._create_msstats_protein_join_optimized(msstats_in, protein_groups_info)
        
        # Step 4: Process in batches using SQL aggregation
        for file_batch in self._get_file_batches_optimized(msstats_in, file_num):
            sql = self._get_protein_aggregation_sql(experiment_type, file_batch)
            batch_results = msstats_in._duckdb.execute(sql).df()
```

#### 3. SQL Aggregation Strategy

Instead of expensive pandas operations, the system uses optimized SQL:

```sql
SELECT 
    anchor_protein,
    pg_accessions,
    reference_file_name,
    channel,
    SUM(intensity) as total_intensity,
    COUNT(*) as feature_count,
    COUNT(DISTINCT peptidoform) as unique_peptide_count,
    MAX(intensity) as max_intensity,
    AVG(intensity) as avg_intensity
FROM processed_msstats_with_pg 
WHERE reference_file_name = ANY({file_batch})
AND anchor_protein IS NOT NULL
AND intensity > 0
GROUP BY anchor_protein, pg_accessions, reference_file_name, channel
```

#### 4. Performance Improvements

- **Before:** 7+ minutes for large datasets
- **After:** 10-45 seconds for the same datasets
- **Memory usage:** Reduced by eliminating pandas aggregation
- **Database files:** Automatically cleaned up after processing

### Resource Management

#### Temporary File Cleanup

The system tracks and cleans up temporary files:

```python
def cleanup(self):
    """Clean up any temporary files and resources."""
    # Close file handles
    for file_handle in self._file_handles:
        if not file_handle.closed:
            file_handle.close()
    
    # Remove temporary files
    for temp_file in self._temp_files:
        if os.path.exists(temp_file):
            os.unlink(temp_file)
```

#### DuckDB Database Cleanup

DuckDB databases are automatically cleaned up:

```python
class MsstatsIN(DuckDB):
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        self.destroy_duckdb_database()
    
    def __del__(self):
        """Cleanup database views and tables."""
        try:
            if hasattr(self, "_duckdb") and self._duckdb:
                self._duckdb.execute("DROP VIEW IF EXISTS processed_msstats")
                self._duckdb.execute("DROP TABLE IF EXISTS protein_groups")
                self.destroy_duckdb_database()
        except:
            pass
```

### Temporary File Patterns

The system creates temporary files with specific patterns that are automatically ignored by Git:

- **mzTab temp files:** `mztab_temp*.mzTab*`
- **DuckDB databases:** `*duckdb*.db*`
- **Test outputs:** `test_output*/`, `output*/`

## MsstatsIN, SDRF, and Channel Mapping Integration

### Overview

The quantms.io library processes mass spectrometry data from different sources (MaxQuant, FragPipe, DIA-NN) and transforms it into a standardized format. A critical component of this process is the proper handling of multiplexed experiments (TMT/iTRAQ) where multiple samples are analyzed simultaneously using isobaric tags.

### Architecture Components

#### 1. MsstatsIN Class (`quantmsio.core.quantms.msstats_in.py`)

The `MsstatsIN` class is responsible for:
- Reading and processing msstats CSV files (including gzipped files via DuckDB)
- Transforming raw numeric channel data into standardized channel names
- Aggregating channel intensities per unique peptide feature
- Integrating with SDRF metadata for experimental design information

#### 2. SDRF Handler (`quantmsio.core.sdrf.py`)

The SDRF (Sample and Data Relationship Format) handler:
- Parses experimental metadata from SDRF files
- Detects experiment type (LFQ, TMT6, TMT10, TMT11, TMT16, iTRAQ4, iTRAQ8)
- Creates sample-to-channel mappings
- Provides experimental design context

#### 3. Channel Constants (`quantmsio.utils.constants.py`)

Defines the standardized channel naming conventions for different multiplexing methods.

### Channel Mapping Logic

#### Step 1: Experiment Type Detection

The SDRF file is analyzed to determine the experiment type based on the `comment[label]` column:

```python
def get_experiment_type_from_sdrf(self):
    labeling_values = get_complex_value_sdrf_column(self.sdrf_table, self.LABELING)
    labeling_values = [i.upper() for i in labeling_values]
    
    if len([i for i in labeling_values if "TMT" in i]) > 0:
        if len(labeling_values) == 10:
            return "TMT10"
        elif len(labeling_values) == 11:
            return "TMT11"
        elif len(labeling_values) == 16:
            return "TMT16"
        elif len(labeling_values) == 6:
            return "TMT6"
```

**Example SDRF TMT11 labels:**
```
comment[label]
TMT126
TMT127N
TMT127C
TMT128N
TMT128C
TMT129N
TMT129C
TMT130N
TMT130C
TMT131N
TMT131C
```

**Result:** `experiment_type = "TMT11"`

#### Step 2: Raw Data Structure

In msstats CSV files, channels are represented as **numeric values**:

```csv
Channel,Reference,ProteinName,PeptideSequence,Charge,Intensity
1,a05058.mzML_...,sp|P55011|S12A2_HUMAN,.(TMT6plex)PEPTIDE...,3,188745.40
2,a05058.mzML_...,sp|P55011|S12A2_HUMAN,.(TMT6plex)PEPTIDE...,3,183611.90
3,a05058.mzML_...,sp|P55011|S12A2_HUMAN,.(TMT6plex)PEPTIDE...,3,194053.50
...
11,a05058.mzML_...,sp|P55011|S12A2_HUMAN,.(TMT6plex)PEPTIDE...,3,199489.80
```

#### Step 3: Channel Number to Name Mapping

The transformation uses array-based lookup with zero-based indexing:

```python
if "TMT" in self.experiment_type:
    msstats["channel"] = msstats["channel"].apply(
        lambda row: TMT_CHANNELS[self.experiment_type][row - 1]
    )
```

**TMT11 Channel Mapping Table:**

| Channel Number | Array Index | TMT Name | Mass (Da) |
|---------------|-------------|----------|-----------|
| 1             | 0           | TMT126   | 126.127726|
| 2             | 1           | TMT127N  | 127.124761|
| 3             | 2           | TMT127C  | 127.131081|
| 4             | 3           | TMT128N  | 128.128116|
| 5             | 4           | TMT128C  | 128.134436|
| 6             | 5           | TMT129N  | 129.131471|
| 7             | 6           | TMT129C  | 129.137790|
| 8             | 7           | TMT130N  | 130.134825|
| 9             | 8           | TMT130C  | 130.141145|
| 10            | 9           | TMT131N  | 131.138180|
| 11            | 10          | TMT131C  | 131.144500|

**Mapping Formula:** `TMT_CHANNELS["TMT11"][channel_number - 1]`

#### Step 4: Channel Aggregation

For each unique peptide feature (defined by `reference_file_name + peptidoform + precursor_charge`), all channel intensities are aggregated:

```python
def transform_experiment(self, msstats):
    intensities_map = {}
    
    def get_intensities_map(rows):
        key = rows["map"]  # unique peptide identifier
        sample_key = rows["reference_file_name"] + "-" + rows["channel"]
        if key not in intensities_map:
            intensities_map[key] = []
        intensities_map[key].append({
            "sample_accession": self._sample_map[sample_key],
            "channel": rows["channel"],
            "intensity": rows["intensity"],
        })
    
    # Aggregate all channels for each unique peptide
    msstats[select_cols].apply(get_intensities_map, axis=1)
    
    # Remove duplicate rows (keeping one row per peptide)
    msstats.drop_duplicates(
        subset=["reference_file_name", "peptidoform", "precursor_charge"],
        inplace=True
    )
    
    # Map aggregated intensities back to each peptide
    msstats.loc[:, "intensities"] = msstats["map"].map(intensities_map)
```

### Data Flow Example

#### Input (Raw msstats data):
```
Channel | Reference | ProteinName | PeptideSequence | Intensity
1       | a05058    | P55011      | PEPTIDEK        | 188745.40
2       | a05058    | P55011      | PEPTIDEK        | 183611.90
3       | a05058    | P55011      | PEPTIDEK        | 194053.50
...
11      | a05058    | P55011      | PEPTIDEK        | 199489.80
```

#### After Channel Mapping:
```
channel | reference_file_name | pg_accessions | peptidoform | intensity
TMT126  | a05058             | P55011        | PEPTIDEK    | 188745.40
TMT127N | a05058             | P55011        | PEPTIDEK    | 183611.90
TMT127C | a05058             | P55011        | PEPTIDEK    | 194053.50
...
TMT131C | a05058             | P55011        | PEPTIDEK    | 199489.80
```

#### After Aggregation (Final output):
```
reference_file_name | pg_accessions | peptidoform | channel | intensities
a05058             | P55011        | PEPTIDEK    | TMT126  | [
                                                           |   {"channel": "TMT126", "intensity": 188745.40, "sample_accession": "Sample1"},
                                                           |   {"channel": "TMT127N", "intensity": 183611.90, "sample_accession": "Sample2"},
                                                           |   {"channel": "TMT127C", "intensity": 194053.50, "sample_accession": "Sample3"},
                                                           |   ...
                                                           |   {"channel": "TMT131C", "intensity": 199489.80, "sample_accession": "Sample11"}
                                                           | ]
```

### SDRF Integration

#### Sample Mapping

The SDRF handler creates a mapping between file-channel combinations and sample accessions:

```python
def get_sample_map_run(self):
    sdrf = self.sdrf_table[["source name", "comment[data file]", "comment[label]"]].copy()
    sdrf["comment[data file]"] = sdrf["comment[data file]"].str.split(".").str[0]
    
    if self.get_experiment_type_from_sdrf() != "LFQ":
        sdrf.loc[:, "map_sample"] = (
            sdrf["comment[data file]"] + "-" + sdrf["comment[label]"]
        )
    else:
        sdrf.loc[:, "map_sample"] = sdrf["comment[data file]"] + "-LFQ"
    
    sample_map = sdrf.to_dict()["source name"]
    return sample_map
```

**Example Sample Map:**
```python
{
    "a05058-TMT126": "PXD007683-Sample-1",
    "a05058-TMT127N": "PXD007683-Sample-2",
    "a05058-TMT127C": "PXD007683-Sample-3",
    ...
    "a05058-TMT131C": "PXD007683-Sample-11"
}
```

### Supported Multiplexing Methods

#### TMT (Tandem Mass Tags)
- **TMT6**: 6 channels (TMT126, TMT127, TMT128, TMT129, TMT130, TMT131)
- **TMT10**: 10 channels (includes N and C variants)
- **TMT11**: 11 channels (TMT126 through TMT131C)
- **TMT16**: 16 channels (extends to TMT134N)

#### iTRAQ (Isobaric Tags for Relative and Absolute Quantitation)
- **iTRAQ4**: 4 channels (iTRAQ114, iTRAQ115, iTRAQ116, iTRAQ117)
- **iTRAQ8**: 8 channels (iTRAQ113 through iTRAQ121)

### Performance Considerations

#### DuckDB Integration
- Native support for gzipped CSV files (no temporary file extraction needed)
- Efficient batch processing for large datasets
- Configurable memory limits and thread usage
- **SQL-based aggregation** instead of pandas operations for 4-6x speed improvement
- **Automatic resource cleanup** prevents memory leaks and file accumulation

#### Batch Processing
- Data is processed in configurable batches by reference files
- Memory-efficient handling of large datasets
- Parallel processing capabilities
- **Optimized SQL joins** with indexed protein group lookups
- **Eliminated expensive LIKE operations** in favor of exact matching

#### Resource Management
- **Context managers** ensure automatic cleanup of database connections
- **Temporary file tracking** prevents accumulation of intermediate files
- **Multiple safety nets** (context manager, destructor, parent cleanup)
- **Git ignore patterns** for temporary files (`mztab_temp*.mzTab*`, `*duckdb*.db*`)

#### Memory Optimization
- **Eliminated ARRAY_AGG operations** that caused memory spikes
- **Simplified SQL queries** using basic aggregation functions
- **Streaming processing** for large datasets
- **Immediate cleanup** of intermediate results

### Error Handling

#### Validation Checks
- Experiment type detection validates supported TMT/iTRAQ configurations
- Channel number bounds checking during mapping
- Sample map uniqueness validation

#### Common Issues
- **Missing SDRF labels**: Experiment type cannot be determined
- **Channel out of bounds**: Raw channel numbers exceed expected range
- **Duplicate sample mappings**: SDRF contains conflicting sample assignments

### Usage Example

```python
from quantmsio.core.quantms.msstats_in import MsstatsIN

# Initialize with gzipped msstats file and SDRF
msstats_in = MsstatsIN(
    report_path="data/msstats_in.csv.gz",
    sdrf_path="metadata/experiment.sdrf.tsv"
)

# Process data in batches
for batch in msstats_in.generate_msstats_in(file_num=5):
    # Each batch contains aggregated channel data
    print(f"Batch size: {len(batch)} features")
    print(f"Experiment type: {msstats_in.experiment_type}")
    
    # Access aggregated intensities
    for _, row in batch.iterrows():
        intensities = row['intensities']  # List of channel intensities
        for intensity_info in intensities:
            channel = intensity_info['channel']      # e.g., "TMT126"
            intensity = intensity_info['intensity']  # Raw intensity value
            sample = intensity_info['sample_accession']  # Sample identifier
```

This implementation ensures proper handling of multiplexed data while maintaining compatibility with different experimental designs and mass spectrometry platforms. 