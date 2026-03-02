# Quality Control Examples

Generate comprehensive QC reports with statistics and visualizations.

## Comprehensive QC Report

Generate a complete QC report using the Python API.

```python
import qpx
import matplotlib.pyplot as plt
from pathlib import Path

OUTPUT_DIR = Path("output")
QC_DIR = Path("qc_report")
QC_DIR.mkdir(parents=True, exist_ok=True)

# Open dataset
ds = qpx.open(str(OUTPUT_DIR))

# --- Statistics ---
report_lines = ["QPX Quality Control Report", "=" * 50, ""]

# PSM statistics
if hasattr(ds, "psm") and ds.psm.count() > 0:
    psm_df = ds.psm.to_df()
    report_lines.append("PSM Statistics:")
    report_lines.append(f"  Total PSMs: {len(psm_df):,}")
    report_lines.append(f"  Unique peptides: {psm_df['sequence'].nunique():,}")
    report_lines.append(f"  Unique proteins: {psm_df['protein_accessions'].nunique():,}")
    report_lines.append(f"  Runs: {psm_df['run_file_name'].nunique()}")
    report_lines.append("")

# Feature statistics
if hasattr(ds, "feature") and ds.feature.count() > 0:
    report_lines.append("Feature Statistics:")
    report_lines.append(f"  Total features: {ds.feature.count():,}")
    report_lines.append("")

# Protein group statistics
if hasattr(ds, "pg") and ds.pg.count() > 0:
    report_lines.append("Protein Group Statistics:")
    report_lines.append(f"  Total protein groups: {ds.pg.count():,}")
    report_lines.append("")

report_text = "\n".join(report_lines)
(QC_DIR / "report.txt").write_text(report_text)
print(report_text)

# --- QC Visualizations ---
# See the Visualization Guide for more plot types:
#   https://bigbio.github.io/qpx/guide/visualize/

if hasattr(ds, "feature") and ds.feature.count() > 0:
    feature_df = ds.feature.to_df()

    # Intensity box plot
    fig, ax = plt.subplots(figsize=(12, 6))
    intensity_cols = [c for c in feature_df.columns if c.startswith("sample_")]
    feature_df[intensity_cols].apply(lambda x: x[x > 0].apply("log10")).boxplot(ax=ax)
    ax.set_ylabel("log10(Intensity)")
    ax.set_title("Intensity Distribution")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    fig.savefig(QC_DIR / "intensity_boxplot.svg", format="svg")
    plt.close(fig)

print(f"QC report generated: {QC_DIR}")
```

## Related Documentation

- **[CLI Reference](../guide/index.md)** - Complete command documentation
- **[Format Specification](../spec/index.md)** - Data format details
- **[GitHub Repository](https://github.com/bigbio/qpx)** - Source code and more examples

---

**Need more examples?** Check the [`tests/examples/`](https://github.com/bigbio/qpx/tree/main/tests/examples) directory in the repository for real data files you can use for testing.

[← Back to Examples Overview](index.md)
