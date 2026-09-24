![PyPI - Version](https://img.shields.io/pypi/v/fastapy?label=version&color=blue)
![tests workflow](https://github.com/aziele/fastapy/actions/workflows/run-tests.yml/badge.svg)
![OS](https://img.shields.io/badge/OS-Linux%20MacOS%20Windows-7373e3)
[![Python versions](https://img.shields.io/pypi/pyversions/fastapy)](https://pypi.org/project/fastapy/)
[![License](https://img.shields.io/github/license/aziele/fastapy)](https://github.com/aziele/fastapy/blob/main/LICENSE)

# fastapy

A lightweight Python package to read and write sequence records in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format), with support for compressed files and no external dependencies.

## Key features

- **No dependencies:** Uses only the Python standard library.
- **Streaming:** Reads and writes records one by one without loading entire files into memory.
- **Compression:** Reads gzip, bzip2, and ZIP files; writes gzip and bzip2 files.
- **Path support:** Accepts string paths and `pathlib.Path` objects.

## Installation

Requires Python 3.8 or later. Install from [PyPI](https://pypi.org/project/fastapy/):

```bash
pip install fastapy
```

Or install directly from GitHub:

```bash
pip install "git+https://github.com/aziele/fastapy.git"
```

Because `fastapy` has no third-party dependencies, you can also use it directly by downloading or cloning the repository:

```bash
git clone https://github.com/aziele/fastapy.git
cd fastapy
python -c "import fastapy"
```

## Quick start

Read a FASTA file and save sequences of at least 100 residues to a compressed file:

```python
import fastapy

records = (
    record
    for record in fastapy.parse("input.fasta")
    if len(record) >= 100
)
count = fastapy.write(records, "filtered.fasta.gz")
print(f"Wrote {count} records")
```

Records are processed one by one without loading the entire file into memory.

## Usage

| Object or function | Purpose |
|---|---|
| `Record(id, seq, desc=None)` | Create a record from an identifier, sequence, and optional description |
| `parse(filename)` | Iterate over records in a plain or compressed FASTA file |
| `parse_handle(handle)` | Iterate over FASTA records from an open text handle |
| `read(filename)` | Return the first record from a FASTA file |
| `to_dict(records)` | Create a dictionary mapping sequence identifiers to records |
| `write(records, destination)` | Write records to a file |

### Reading records

Use `parse()` to iterate over FASTA records:

```python
import fastapy

for record in fastapy.parse("tests/test.fasta"):
    print(record.id, len(record), record.seq[:10], record.desc)
```

Output:

```
NP_002433.1  362   METDAPQPGL   RNA-binding protein Musashi homolog 1 [Homo sapiens]
ENO94161.1    79   MKLLISGLGP   RRM domain-containing RNA-binding protein
sequence     292   MKLSKIALMM
```

Input can be a filename or `pathlib.Path`. Compression is selected from the extension: `.gz` or `.gzip` (gzip), `.bz2` (bzip2), or `.zip` (ZIP). For ZIP archives, the first entry is read and must be a FASTA text file.

```python
for record in fastapy.parse("tests/test.fasta.gz"):
    print(record.id)

# NP_002433.1
# ENO94161.1
# sequence
```


For an open text handle, use `parse_handle()`:

```python
with open("input.fasta") as handle:
    for record in fastapy.parse_handle(handle):
        print(record.id)
```

### Reading the first record

Use `read()` to return the first record. Additional records are ignored; a file containing no records raises `ValueError`.

```python
record = fastapy.read("tests/test.fasta")
print(record.id)      # NP_002433.1
```

### Working with records

A `Record` stores a sequence identifier (`id`), sequence (`seq`), and optional description (`desc`). The `description` property returns the complete FASTA header, including `>`.

```python
record = fastapy.Record(
    id='NP_950171.2', 
    seq='MEEEAETEEQQRFSYQQRLKAAVHYTVGCLCEEVALDKEMQFSKQTIAAISELTFRQCENFAKDLEMFASICRKRQE',
    desc='APITD1-CORT protein isoform 2 [Homo sapiens]'
)

print(record.id)            # NP_950171.2
print(record.desc)          # APITD1-CORT protein isoform 2 [Homo sapiens]
print(record.seq)           # MEEEAE..
print(record.description)   # >NP_950171.2 G APITD1-CORT protein isoform 2 [Homo sapiens]
print(len(record))          # 77
print('EEEA' in record)     # True
```

Iterating over a record yields its sequence characters. 

```python
print(list(record)[:10])
# ['M', 'E', 'E', 'E', 'A', 'E', 'T', 'E', 'E', 'Q']
```

By default, sequences are wrapped at 70 characters per line. Use `wrap` to set a different line width, or `wrap=None` (or `wrap=0`) for no wrapping.

```python
print(record)
# >NP_950171.2 APITD1-CORT protein isoform 2 [Homo sapiens]
# MEEEAETEEQQRFSYQQRLKAAVHYTVGCLCEEVALDKEMQFSKQTIAAISELTFRQCENFAKDLEMFAS
# ICRKRQE

print(record.format(wrap=30), end="")
# >NP_950171.2 APITD1-CORT protein isoform 2 [Homo sapiens]
# MEEEAETEEQQRFSYQQRLKAAVHYTVGCL
# CEEVALDKEMQFSKQTIAAISELTFRQCEN
# FAKDLEMFASICRKRQE

print(record.format(wrap=None), end="")
# >NP_950171.2 APITD1-CORT protein isoform 2 [Homo sapiens]
# MEEEAETEEQQRFSYQQRLKAAVHYTVGCLCEEVALDKEMQFSKQTIAAISELTFRQCENFAKDLEMFASICRKRQE
```

The default line width is 70 characters. Use `wrap=0` or `wrap=None` for no wrapping.

### Storing records in a list or dictionary

Use a list for repeated access or positional indexing:

```python
records = list(fastapy.parse("input.fasta"))
print(records[0].id)   # First record id
print(records[-1].id)  # Last record id
```

Use `to_dict()` to index records by identifier. Duplicate identifiers raise `ValueError`.

```python
records = fastapy.to_dict(fastapy.parse("tests/test.fasta"))
print(records["NP_002433.1"].description)
# NP_002433.1 RNA-binding protein Musashi homolog 1 [Homo sapiens]
```

Both approaches store all records in memory.

### Writing records

Use `write()` to save a single record or an iterable of records to a filename, `pathlib.Path`, or open text handle. The function returns the number of records written.

```python
# Write a single record
record = fastapy.Record(id="seq1", seq="ATCG")
fastapy.write(record, "sequence.fasta")

# Write a list of records
records = [record, fastapy.Record(id="seq2", seq="GGTA")]
count = fastapy.write(records, "sequences.fasta", wrap=80)
print(count)  # 2

# Write records stored in a dictionary
records_by_id = fastapy.to_dict(records)
fastapy.write(records_by_id.values(), "sequences.fasta.gz")

# Write directly from a generator
fastapy.write(fastapy.parse("input.fasta"), "output.fasta.bz2", wrap=0)
```

Compression and options:

- Output compression is selected from the extension: `.gz` or `.gzip` for gzip, and `.bz2` for bzip2.
- ZIP output is not supported; a `.zip` destination raises `ValueError`.
- Sequences are wrapped to 70 characters by default. Use `wrap=0` or `wrap=None` for no wrapping.
- Existing output files are overwritten.
- When you pass an open text handle, `write()` does not close it.

```python
with open("sequences.fasta", "w") as handle:
    fastapy.write(records, handle)
```

## Testing

From the repository root, run:

```bash
pip install pytest
pytest
```

## License

Distributed under the terms of the [GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html).