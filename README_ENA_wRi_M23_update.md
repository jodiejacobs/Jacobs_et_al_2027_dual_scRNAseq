# README: Updating ENA Assembly ERZ28976196 (PRJEB107302) with New Annotation

This documents the working procedure used to submit an updated, annotated
version of the wRi Merrill 23 *Wolbachia* assembly to ENA, replacing the
original `wRi_Merrill23_polished` submission. Confirmed successful — this is
the real sequence of commands, including the fixes needed along the way.

## Known values (reuse these for future updates)

| Field | Value |
|---|---|
| Study accession | `PRJEB107302` |
| Secondary study accession | `ERP188321` |
| Sample accession | `ERS28645251` |
| Original analysis (ERZ) | `ERZ28976196` |
| Original assembly name | `wRi_Merrill23_polished` (retired — do not reuse) |
| Locus tag prefix | `M23` (registered on the project — reuse, don't re-register) |
| Assembly type | `clone or isolate` |
| Coverage | `57` |
| Assembly program | `Flye` |
| Sequencing platform | `OXFORD_NANOPORE` |
| Min gap length | `10` |
| Molecule type | `genomic DNA` |
| Chromosome name | `wRi_chromosome` |
| Chromosome type/topology | `Circular-Chromosome` |
| Bakta database | `/private/groups/russelllab/dbs/bakta_db` |
| Webin-CLI jar | `/private/home/jomojaco/ena_submission/webin-cli.jar` |
| Webin username | `Webin-74134` |

Note on chromosome naming: `wRi_chromosome` contains the substring
"chromosome", which current ENA naming rules disallow for *new* chromosome
names. It passed here because it's a reused name from the original
submission — updates carry it forward rather than re-validating it as new.

## Working directory layout

Everything (manifest, flat file, chromosome list) lives together in one
submission directory, e.g.:
```
/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/ena_submission/genome/16S_Assembly_polished/
```
Webin-CLI resolves manifest file paths relative to this directory — files
referenced in `manifest.txt` must physically be here, not just in the Bakta
output folder.

---

## Step 1: Annotate the new assembly with Bakta

```bash
mamba activate bakta

bakta --db /private/groups/russelllab/dbs/bakta_db \
  --locus-tag M23 \
  --genus Wolbachia --species "sp. wRi" --strain "Merrill 23" \
  --complete \
  --output bakta_out \
  --prefix wRi_M23_v3 \
  /private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/v3_assembly_fix_16S/manual_insertion/20260402_wRi_M23_manual_insertion_pilon.fasta
```

Output of interest: `bakta_out/wRi_M23_v3.embl`

Bakta named the single circular chromosome **`contig_1`** (check with
`grep -n "^ID\|^AC" wRi_M23_v3.embl` — this name varies run to run, always
re-check it).

Copy the flat file into the submission directory:
```bash
cp bakta_out/wRi_M23_v3.embl /private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/ena_submission/genome/16S_Assembly_polished/
```

---

## Step 2: Fix the flat file — add the ENA `AC *` line

Bakta produces standard EMBL format, which does **not** include the special
`AC * _{entry_name}` line ENA's Webin-CLI requires to identify the sequence
before it has a real accession. This has to be inserted manually every time.

```bash
cd /private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/ena_submission/genome/16S_Assembly_polished/

sed -i '3a XX\nAC * _contig_1' wRi_M23_v3.embl
```

(Adjust `_contig_1` to match whatever sequence name Bakta actually gave you
this run.)

Confirm:
```bash
head -8 wRi_M23_v3.embl
```
Should read:
```
ID   contig_1; ; circular; DNA; ; UNC; 1447516 BP.
XX
AC   contig_1;
XX
AC * _contig_1
XX
DE   Wolbachia sp. wri Merrill 23 contig_1, whole genome shotgun sequence
...
```

Gzip it:
```bash
gzip wRi_M23_v3.embl
```

---

## Step 3: Build the chromosome list file

The `OBJECT_NAME` (column 1) must exactly match the sequence name from the
`AC *` line — **without** the leading underscore. It does *not* need to
match the original submission's old sequence ID.

```bash
echo -e "contig_1\twRi_chromosome\tCircular-Chromosome" > chromosome_list.txt
gzip -f chromosome_list.txt
```

Error this catches if wrong: `ERROR: Sequenceless chromosomes are not
allowed in assembly : <name>` — means column 1 doesn't match any sequence
actually in the flat file.

---

## Step 4: Write the Webin-CLI manifest

`manifest.txt`:
```
STUDY           PRJEB107302
SAMPLE          ERS28645251
ASSEMBLYNAME    wRi_Merrill23_polished_v2
ASSEMBLY_TYPE   clone or isolate
COVERAGE        57
PROGRAM         Flye
PLATFORM        OXFORD_NANOPORE
MINGAPLENGTH    10
MOLECULETYPE    genomic DNA
FLATFILE        wRi_M23_v3.embl.gz
CHROMOSOME_LIST chromosome_list.txt.gz
```

`ASSEMBLYNAME` must be unique each time — increment it for future updates
(e.g. `_v3`, `_v4`).

Error this catches if a data file is missing/misplaced: `ERROR: Could not
read data file: "<filename>"` — check the file is physically in this
directory (see Working directory layout, above).

---

## Step 5: Validate

```bash
java -jar /private/home/jomojaco/ena_submission/webin-cli.jar \
  -username Webin-74134 -password <your_webin_password> \
  -context genome -manifest manifest.txt -validate
```

On failure, check the report(s) under:
```
genome/<ASSEMBLYNAME>/validate/
```

Errors hit and fixed during this submission:
- `Could not read data file` → flat file wasn't copied into the submission directory yet (Step 1, last command).
- `Sequenceless chromosomes are not allowed in assembly` → chromosome list still referenced the old/wrong sequence name (Step 3).

---

## Step 6: Submit

Once validation is clean:
```bash
java -jar /private/home/jomojaco/ena_submission/webin-cli.jar \
  -username Webin-74134 -password <your_webin_password> \
  -context genome -manifest manifest.txt -submit
```

Returns a new ERZ accession. Assembly version increments; GCA accession
stays the same. Confirmation and final accessions also arrive by email;
allow up to a week for the public record to update.

---

## Key constraints to remember for next time

- Same study (`PRJEB107302`) + sample (`ERS28645251`) pair — required for this to register as an update, not a new assembly.
- Chromosome name (`wRi_chromosome`) must be preserved exactly.
- `ASSEMBLYNAME` must be new/unique every submission.
- Locus tag prefix `M23` is already registered — reuse, never re-register.
- Bakta's EMBL output always needs the `AC *` line added by hand (Step 2) — this isn't a one-off bug, it's a gap between standard EMBL and ENA's submission format.
- The chromosome list's sequence name must match whatever Bakta happens to call the contig on that run (check every time — it isn't guaranteed to be `contig_1`).
- Don't commit your Webin password to any script or shared file — pass it at the command line each time.
