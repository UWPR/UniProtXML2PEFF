# UniProtXML2PEFF

UniProtXML2PEFF is a command-line tool designed to convert UniProt XML files into [PEFF (PSI extended FASTA format)](https://www.psidev.info/peff), enabling better compatibility with proteomics tools such as [Comet](https://github.com/UWPR/Comet). 

This tool processes **sequence variants** and **modifications** in UniProt XML files, accurately encoding them into PEFF format with `\VariantSimple`, `\VariantComplex`, and `\ModResUnimod` annotations. The result is optimized for use in proteomics database searching.

---

## Why Use UniProtXML2PEFF?

PEFF files directly retrieved from UniProt via their API (e.g., `curl "https://www.ebi.ac.uk/proteins/api/proteins?query=organism_id:9606&format=peff" -o human.peff`) may not always be suitable for Comet or other proteomics tools. These API-provided entries often use placeholders like `SGRP` (Sequence GRouP) for `\VariantSimple` values. For example:

```
\VariantSimple=(8|L|SGRP)(137|Q|SGRP)(179|R|SGRP)
```

In this context, `SGRP` does not provide specific residue substitutions and is therefore not useful for Comet searching. Similarly, modifications such as phosphorylations or acetylations might lack the specificity required for meaningful proteomics searches.

UniProtXML2PEFF addresses these issues by:
- Extracting specific sequence variant annotations directly from UniProt XML.
- Mapping modifications (e.g., phosphorylations, methylations) to **Unimod identifiers** for encoding in `\ModResUnimod`.
- Skipping or reporting entries that lack sufficient information for proper PEFF annotations.

---

## Features

- Converts UniProt XML files to PEFF.
- Processes:
  - **Simple substitutions (e.g., A → V)** into `\VariantSimple`.
  - **Complex variants** such as deletions, insertions, and multi-residue substitutions into `\VariantComplex`.
  - **Post-translational modifications (PTMs)** such as phosphorylations and acetylations are encoded into `\ModResUnimod` with Unimod identifiers.
- Logs skipped or unsupported variants and modifications for audit and troubleshooting.

---

## Installation

To build this tool, you'll need a C++ compiler that supports C++11 and the TinyXML-2 library.

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/UniProtXML2PEFF.git
   cd UniProtXML2PEFF
   ```

2. Build the executable:
   ```bash
   make
   ```

   This will generate the `UniProtXML2PEFF.exe` executable.

---

## Usage

The tool requires an input UniProt XML file and an output PEFF file path.

### Command Syntax
```
./UniProtXML2PEFF.exe <input.xml> <output.peff> [--strict]
```

- **`<input.xml>`**: The path to the input UniProt XML file.
- **`<output.peff>`**: The desired path for the PEFF output file.
- **`--strict`** (optional): Enforces strict handling of annotations. If unsupported variants or modifications are encountered, the program will exit with an error.

### Example
```bash
./UniProtXML2PEFF.exe test.xml test.peff
```

---

## Input Format

### Supported XML Structure
The input file should adhere to the standard UniProt XML format. The tool expects `<feature>` elements under `<entry>` with `type` values such as:

- `sequence variant`
- `mutagenesis site`
- `modified residue`

Each `<feature>` element can include:
- **Substitutions**:
  ```xml
  <feature type="sequence variant" description="Substitution A to V.">
      <original>A</original>
      <variation>V</variation>
      <location>
          <position position="15"/>
      </location>
  </feature>
  ```

- **Deletions**:
  ```xml
  <feature type="sequence variant" description="Deletion of residues.">
      <location>
          <begin position="20"/>
          <end position="22"/>
      </location>
  </feature>
  ```

- **Modifications**:
  Marked as `modified residue` in UniProt XML, these modifications are mapped to Unimod identifiers in the PEFF output.

  Example XML:
  ```xml
  <feature type="modified residue" description="Phosphothreonine.">
      <location>
          <position position="10"/>
      </location>
  </feature>
  ```

---

## Output Format

### PEFF Format Annotations
The generated PEFF files include:
- **`VariantSimple`**: Encodes single residue substitutions.
  Example:
  ```
  \VariantSimple=(15|A|V)
  ```

- **`VariantComplex`**: Encodes multi-residue changes, deletions, or insertions.
  Examples:
  ```
  \VariantComplex=(20|22|)
  \VariantComplex=(6|6|GP)
  ```

- **`ModResUnimod`**: Encodes PTMs with Unimod accession numbers.
  - The tool maps modification descriptions (e.g., `Phosphothreonine`) to Unimod identifiers (`UNIMOD:21`) using a predefined mapping.
  Example:
  ```
  \ModResUnimod=(10|UNIMOD:21)(20|UNIMOD:1)
  ```

### Example Output
```
# PEFF 1.0
# Database=UniProt
# VariantSimple=true
# VariantComplex=true
# ModResUnimod=true
>tr|TEST123| \VariantSimple=(15|A|V) \VariantComplex=(20|22|) \ModResUnimod=(10|UNIMOD:21)
AAAAAGGGGG
```

---

## ModResUnimod Mapping

The tool uses a predefined map of modification descriptions to Unimod accessions. For example:
```
"Phosphothreonine" → UNIMOD:21
"Phosphoserine"    → UNIMOD:21
"Acetylation"      → UNIMOD:1
```

Modifications in the XML that do not match this map will be skipped unless the `--strict` option is used.

---

## Audit Logs

The tool generates two audit files to track skipped or processed annotations:
1. **`variant_skipped.csv`**:
   Logs the reasons for skipped variants and modifications, including:
   - Unsupported feature types.
   - Modifications or variants with missing location data.

2. **`variant_complex.csv`**:
   Summarizes the types and counts of `\VariantComplex` entries generated.

---

## Troubleshooting

### Common Issues
1. **Missing Modifications in Output**:
   - Check the `variant_skipped.csv` file to confirm if modifications were skipped due to missing Unimod mappings.

2. **Strict Mode Fails**:
   - Use the `--strict` option for debugging. If the tool fails, inspect the logs for skipped features.

3. **SGRP in Retrieved PEFF Files**:
   - If you retrieve PEFF files directly from UniProt and encounter entries like `\VariantSimple=(8|L|SGRP)`, use this tool to regenerate the PEFF with specific substitutions or variants.

### Debug Output
Run with standard error redirection to debug processing:
```bash
./UniProtXML2PEFF.exe input.xml output.peff 2> debug.log
```

---

## License

This project is open-source and distributed under the MIT License.
