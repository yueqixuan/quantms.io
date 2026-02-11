# Avro Schema Reference

All QPX views have formal schema definitions in Apache Avro format. These schemas define the structure and types of each field and serve as the canonical, language-independent schema definition for the QPX format.

!!! info "Avro vs Parquet"
    Avro schemas (`.avsc`) define the **logical data model** -- the field names, types, nullability, and documentation for each view. The actual data is serialized in **Apache Parquet** format for efficient columnar storage and querying. The Avro schemas are the single source of truth for the QPX data model, while Parquet handles the physical storage layer.

## How to use these schemas

- **Validation**: Use these schemas to validate that a QPX file conforms to the expected structure.
- **Code generation**: Avro schemas can be used to generate data classes in Python, Java, and other languages.
- **Documentation**: Each field includes a `doc` attribute describing its purpose and semantics.

## Schemas

=== "PSM"

    ```json title="psm.avsc"
    --8<-- "docs/psm.avsc"
    ```

=== "Feature"

    ```json title="feature.avsc"
    --8<-- "docs/feature.avsc"
    ```

=== "Peptide"

    ```json title="peptide.avsc"
    --8<-- "docs/peptide.avsc"
    ```

=== "Protein Group"

    ```json title="pg.avsc"
    --8<-- "docs/pg.avsc"
    ```

=== "Protein"

    ```json title="protein.avsc"
    --8<-- "docs/protein.avsc"
    ```

=== "Mass Spectra"

    ```json title="mz.avsc"
    --8<-- "docs/mz.avsc"
    ```

=== "Absolute Expression"

    ```json title="absolute.avsc"
    --8<-- "docs/absolute.avsc"
    ```

=== "Differential Expression"

    ```json title="differential.avsc"
    --8<-- "docs/differential.avsc"
    ```

!!! tip "Viewing schemas programmatically"
    You can load and inspect any Avro schema in Python using the `fastavro` library:

    ```python
    import json

    with open("psm.avsc") as f:
        schema = json.load(f)

    for field in schema["fields"]:
        print(f"{field['name']}: {field['type']}")
    ```
