# Updating ENA Assembly ERZ28976196 (PRJEB107302) with New Annotation

Target assembly to submit:
`/private/groups/russelllab/jodie/Jacobs_et_al_2026_de_novo_wRi_merrill_23_assembly/v3_assembly_fix_16S/manual_insertion/20260402_wRi_M23_manual_insertion_pilon.fasta`

This is an **update** to an existing annotated assembly (not a new submission). ENA rules require reusing the same study, sample, and chromosome name(s) as the original.

## What we already know

| Field | Value |
|---|---|
| Study accession | `PRJEB107302` |
| Secondary study accession | `ERP188321` |
| Sample accession | `ERS28645251` |
| Original analysis (ERZ) | `ERZ28976196` |
| Original assembly name | `wRi_Merrill23_polished` (taken — must use a new name) |
| Locus tag prefix | `M23` (already registered — reuse as-is) |
| Assembly type | `clone or isolate` |
| Original coverage | `57` |
| Assembly program | `Flye` |
| Sequencing platform | `OXFORD_NANOPORE` |
| Min gap length | `10` |
| Molecule type | `genomic DNA` |
| Bakta database path | `/private/groups/russelllab/dbs/bakta_db` |
| Original chromosome name | `wRi_chromosome` |
| Original chromosome type/topology | `Circular-Chromosome` |

Confirmed from the original chromosome list:
```
wRi_Merrill23_polished    wRi_chromosome    Circular-Chromosome
```
(column 1 = old sequence ID, column 2 = chromosome name, column 3 = topology-chromosome_type)

**Still needed:** the sequence ID (header) used in your *new* FASTA/Bakta EMBL output — this replaces column 1 above and won't match the old `wRi_Merrill23_polished` id since it's a new assembly run. Check it with:
```bash
grep "^>" 20260402_wRi_M23_manual_insertion_pilon.fasta
```
(or, once Bakta finishes, the `AC * ` line in `bakta_out/wRi_M23_v3.embl`)

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

Confirm the `/locus_tag` values in this file use the `M23` prefix and are unique — this is the most common validation failure.

---

## Step 2: Build the chromosome list file

Create `chromosome_list.txt`, reusing the confirmed chromosome name/type and swapping in your new sequence ID (from the `grep "^>"` check above):

```
<your_new_seq_id>	wRi_chromosome	Circular-Chromosome
```

Gzip it:
```bash
gzip chromosome_list.txt
```

---

## Step 3: Gzip the flat file

```bash
gzip -k bakta_out/wRi_M23_v3.embl
```

---

## Step 4: Write the Webin-CLI manifest

Create `manifest.txt`:

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

Notes:
- `ASSEMBLYNAME` must be new/unique — `wRi_Merrill23_polished` is already taken by the original.
- Update `COVERAGE` if the v3 assembly's actual coverage differs from the original.
- Put `manifest.txt`, the gzipped flat file, and the gzipped chromosome list in the same directory.

---

## Step 5: Download Webin-CLI

```bash
wget https://github.com/enasequence/webin-cli/releases/latest/download/webin-cli.jar
```

(Or check the [releases page](https://github.com/enasequence/webin-cli/releases) for the current filename/version.)

---

## Step 6: Validate

```bash
java -jar webin-cli.jar \
  -username Webin-XXXXX -password YYYYYYY \
  -context genome -manifest manifest.txt -validate
```

If this fails, check the report file under `genome/<ASSEMBLYNAME>/validate/` for specific errors (commonly: locus_tag issues, chromosome name mismatch, or partiality/translation errors in CDS features). Fix and re-run `-validate` until clean.

---

## Step 7: Submit

```bash
java -jar webin-cli.jar \
  -username Webin-XXXXX -password YYYYYYY \
  -context genome -manifest manifest.txt -submit
```

This returns a new ERZ accession. The assembly version increments; the GCA accession stays the same.

---

## Step 8: Verify

- Check the [Webin Portal](https://www.ebi.ac.uk/ena/submit/webin/) for submission status.
- Confirmation and accession numbers also arrive by email.
- Allow up to a week for the public record to reflect the update.

---

## Key constraints to remember

- Same study + sample pair as the original — non-negotiable, this is what marks it as an "update" rather than a new assembly.
- Original chromosome name must be preserved exactly (Step 1).
- `ASSEMBLYNAME` must be unique — never reuse a previous one.
- Locus tag prefix `M23` is already registered — do not re-register or change it.
- Submit at least 24 hours after any prior submission's processing completed, to avoid pipeline conflicts.
- Current ENA naming rules disallow "chromosome" as a substring in a chromosome name, but `wRi_chromosome` was already accepted in the original submission. Since this is an update reusing that exact name, it should carry over — but if Webin-CLI flags it during validation, this is a known possible snag to raise with ENA helpdesk rather than something to fix by renaming (renaming would break update continuity).
