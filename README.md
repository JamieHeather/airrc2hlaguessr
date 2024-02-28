# airrc2hlaguessr

### 2024, JH @ MGH

This is a simple script to convert TCR repertoire files in the [AIRR Community standard schema](https://docs.airr-community.org) (as proposed in [Vander Heiden *et al.*, 2018](https://doi.org/10.3389/fimmu.2018.02206)) into that required for input into the [HLAGuessr](https://github.com/statbiophys/HLAGuessr) tool (currently in [pre-print](https://doi.org/10.1101/2024.01.25.577228)).

#### Requisites

Python >= 3.6 ish

No non-standard packages

#### Usage

Download the script and place it in an appropriate file with or near your AIRR-C-format TCR data. Then either pipe a list of files to the script, or explicitly provide names or paths via the `-in` flag.

```bash
# Simply pipe multiple files in ...
ls *-all-my-files.tsv | python airrc2hlaguessr.py 

# ... or use the specific input flag
python airrc2hlaguessr.py -in some-specific-file.tsv.gz
python airrc2hlaguessr.py -in file1.tsv,file2.tsv
```

The script will then iterate through the list of input files, converting them to the necessary format:

* A three column tsv file consisting of `cdr3aa`, `v_family`, and `Patient` 
  * `cdr3aa` is simply the CDR3 junction amino acid sequence (conserved C to F residue inclusive)
  * `v_family` is the *family* name of the V gene used, e.g. 'TRAV1' not 'TRAV1-2'
  * `Patient` is a unique identifier for the individual
* The rows then become the *unique* and **productive** entries found in the input file (based off the `productive` field of the AIRR-C input)
* Output rows are also sorted on `cdr3aa` sequence

#### Additional commands

The following flags can be used for additional functionality:

| Flag  | Field               | Description                                                                                                                                                                                                           |
|-------|---------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-ia`  | `--ignore_ambiguous` | If a CDR3 was provided with an ambiguous comma-delimited V gene list (e.g. 'TRBV7,TRBV11') the script ordinarily write an entry out for *both* possibilities. This flag overrides that behaviour and just skips them. |
| `-z`  | `--compress`        | Automatically gzip compress the output.                                                                                                                                                                               |
| `-kd` | `--keep_duplicates` | Keep duplicate copies of TCRs if present in the input file. Doesn't seem to make a difference to prediction.                                                                                                          |
| `-ts` | `--truncate_str`    | Useful for truncating long filenames to get better Patient IDs. E.g. use `-ts _` for files named SomeID_XYZ.tsv.                                                                                                      |
| `-c`  | `--chain_filter`    | Include only one chain's rearrangements. `-c A` to keep alpha, `-c B` to keep beta.                                                                                                                                   |
| `-n` | `--name_override` | Explicitly provide a string to use in the `Patient` ID field, instead of inferring from the input file name.                                                                                                          |

#### Disclaimery stuff

I'm not affiliated with the authors of `HLAGuessr`, and nor am I responsible in any way for the AIRR-C format, so questions about either should be directed to their respective authors. This is just a tool to convert files from the most widely used (non-commercial) TCR rearrangement format to the one required for a cool new tool.