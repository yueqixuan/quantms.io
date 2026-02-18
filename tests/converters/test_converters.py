"""Converter adapter tests with synthetic data."""

import pyarrow.parquet as pq


class TestFragPipeFeatureAdapter:
    def _write_ion_tsv(self, tmp_path):
        """Write a minimal combined_ion.tsv."""
        tsv = tmp_path / "combined_ion.tsv"
        lines = [
            "Peptide Sequence\tModified Sequence\tCharge\tM/Z\tProtein\tProtein ID\tGene\tAssigned Modifications\texperiment_1 Intensity\texperiment_2 Intensity",
            "PEPTIDEK\tPEPTIDEK\t2\t450.25\tsp|P12345|PROT_HUMAN\tP12345\tBRCA1\t\t1000.0\t2000.0",
            "ANOTHERK\tANOTHERK\t3\t500.30\tsp|P67890|PROT2_HUMAN\tP67890\tTP53\t7M(15.9949)\t3000.0\t0.0",
        ]
        tsv.write_text("\n".join(lines) + "\n")
        return tsv

    def _write_peptide_tsv(self, tmp_path):
        """Write a minimal combined_peptide.tsv."""
        tsv = tmp_path / "combined_peptide.tsv"
        lines = [
            "Peptide Sequence\tModified Peptide\tCharges\tProtein\tProtein ID\tGene\texperiment_1 Intensity\texperiment_2 Intensity",
            "PEPTIDEK\tPEPTIDEK\t2,3\tsp|P12345|PROT_HUMAN\tP12345\tBRCA1\t1000.0\t2000.0",
        ]
        tsv.write_text("\n".join(lines) + "\n")
        return tsv

    def test_convert_combined_ion(self, tmp_path):
        from qpx.converters.fragpipe.feature_adapter import FragPipeFeatureAdapter

        tsv = self._write_ion_tsv(tmp_path)
        output = tmp_path / "test.feature.parquet"
        with FragPipeFeatureAdapter() as adapter:
            adapter.convert(feature_path=str(tsv), output_path=str(output))
        assert output.exists()
        table = pq.read_table(output)
        # row1 -> 2 experiments (1000 + 2000), row2 -> 1 experiment (3000, experiment_2 = 0)
        assert table.num_rows == 3
        assert "sequence" in table.schema.names
        assert "anchor_protein" in table.schema.names

    def test_convert_combined_peptide(self, tmp_path):
        from qpx.converters.fragpipe.feature_adapter import FragPipeFeatureAdapter

        tsv = self._write_peptide_tsv(tmp_path)
        output = tmp_path / "test.feature.parquet"
        with FragPipeFeatureAdapter() as adapter:
            adapter.convert(feature_path=str(tsv), output_path=str(output))
        assert output.exists()
        table = pq.read_table(output)
        # 1 peptide x 2 experiments x 2 charges = 4 rows
        assert table.num_rows == 4
        assert "sequence" in table.schema.names

    def test_ion_modifications_parsed(self, tmp_path):
        from qpx.converters.fragpipe.feature_adapter import FragPipeFeatureAdapter

        tsv = self._write_ion_tsv(tmp_path)
        output = tmp_path / "test.feature.parquet"
        with FragPipeFeatureAdapter() as adapter:
            adapter.convert(feature_path=str(tsv), output_path=str(output))
        table = pq.read_table(output)
        # The second peptide (ANOTHERK) has mods, appears in experiment_1 only
        df = table.to_pandas()
        mod_rows = df[df["sequence"] == "ANOTHERK"]
        assert len(mod_rows) == 1
        mods = mod_rows.iloc[0]["modifications"]
        assert mods is not None
        assert len(mods) == 1


class TestFragPipePsmAdapter:
    def test_convert_basic(self, tmp_path):
        from qpx.converters.fragpipe.psm_adapter import FragPipePsmAdapter

        tsv = tmp_path / "psm.tsv"
        lines = [
            "Peptide\tCharge\tSpectrum\tProtein ID\tModified Peptide\tPeptideProphet Probability\tCalculated M/Z\tCalibrated Observed M/Z\tRetention\tAssigned Modifications",
            "PEPTIDEK\t2\tfile1.0001.0001.2\tP12345\tPEPTIDEK\t0.99\t450.25\t450.26\t120.5\t",
        ]
        tsv.write_text("\n".join(lines) + "\n")

        output = tmp_path / "psm.parquet"
        with FragPipePsmAdapter() as adapter:
            adapter.convert(psm_path=str(tsv), output_path=str(output))
        assert output.exists()
        table = pq.read_table(output)
        assert table.num_rows == 1
        assert table.column("sequence")[0].as_py() == "PEPTIDEK"


class TestModificationParsing:
    def test_parse_modifications_from_peptidoform_unimod(self):
        from qpx.converters.quantms.psm_adapter import _parse_modifications_from_peptidoform

        mods_meta = {
            "UNIMOD:35": ("Oxidation", ["M"], ["Anywhere"]),
            "UNIMOD:4": ("Carbamidomethyl", ["C"], ["Anywhere"]),
        }
        result = _parse_modifications_from_peptidoform(
            peptidoform="M[UNIMOD:35]PEPTIDEC[UNIMOD:4]K",
            sequence="MPEPTIDECK",
            modifications_meta=mods_meta,
        )
        assert result is not None
        assert len(result) == 2
        ox = next(m for m in result if m["name"] == "Oxidation")
        assert ox["accession"] == "UNIMOD:35"
        assert ox["positions"][0]["position"] == 1
        assert ox["positions"][0]["amino_acid"] == "M"

    def test_parse_modifications_no_mods(self):
        from qpx.converters.quantms.psm_adapter import _parse_modifications_from_peptidoform

        result = _parse_modifications_from_peptidoform("PEPTIDEK", "PEPTIDEK", {})
        assert result is None

    def test_parse_modifications_mass_shift(self):
        from qpx.converters.quantms.psm_adapter import _parse_modifications_from_peptidoform

        mods_meta = {"UNIMOD:35": ("Oxidation", ["M"], ["Anywhere"])}
        result = _parse_modifications_from_peptidoform(
            peptidoform="M[+15.9949]PEPTIDEK",
            sequence="MPEPTIDEK",
            modifications_meta=mods_meta,
        )
        assert result is not None
        assert len(result) == 1
        assert result[0]["name"] == "Oxidation"
        assert result[0]["accession"] == "UNIMOD:35"

    def test_parse_modifications_nterm(self):
        from qpx.converters.quantms.psm_adapter import _parse_modifications_from_peptidoform

        mods_meta = {"UNIMOD:1": ("Acetyl", ["X"], ["N-term"])}
        result = _parse_modifications_from_peptidoform(
            peptidoform="[UNIMOD:1]-PEPTIDEK",
            sequence="PEPTIDEK",
            modifications_meta=mods_meta,
        )
        assert result is not None
        assert len(result) == 1
        assert result[0]["positions"][0]["position"] == 0  # N-terminal


class TestFragPipePgAdapter:
    def test_convert_basic(self, tmp_path):
        from qpx.converters.fragpipe.pg_adapter import FragPipePgAdapter

        tsv = tmp_path / "combined_protein.tsv"
        lines = [
            "Protein\tProtein ID\tGene\tDescription\tCombined Total Peptides\tCombined Unique Peptides\tCombined Spectral Count\texperiment_1 Total Intensity\texperiment_1 Spectral Count",
            "sp|P12345|PROT_HUMAN\tP12345\tBRCA1\tProtein desc\t10\t5\t20\t50000.0\t20",
        ]
        tsv.write_text("\n".join(lines) + "\n")

        output = tmp_path / "pg.parquet"
        with FragPipePgAdapter() as adapter:
            adapter.convert(protein_path=str(tsv), output_path=str(output))
        assert output.exists()
        table = pq.read_table(output)
        assert table.num_rows == 1
