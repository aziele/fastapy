![PyPI - Version](https://img.shields.io/pypi/v/fastapy?label=version&color=blue)
![tests workflow](https://github.com/aziele/fastapy/actions/workflows/run-tests.yml/badge.svg)
![OS](https://img.shields.io/badge/OS-Linux%20MacOS%20Windows-7373e3)

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
python
```

```python
>>> import fastapy
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
| `write(records, destination, wrap=70)` | Write records to a file or text handle and return the number written |

### Reading records

Use `parse()` to iterate over FASTA records:

```python
import fastapy

for record in fastapy.parse("input.fasta"):
    print(record.id, len(record), record.desc)
```

Input can be a filename or `pathlib.Path`. Compression is selected from the extension: `.gz` or `.gzip` (gzip), `.bz2` (bzip2), or `.zip` (ZIP). For ZIP archives, the first entry is read and must be a FASTA text file.

For an open text handle, use `parse_handle()`:

```python
with open("input.fasta") as handle:
    for record in fastapy.parse_handle(handle):
        print(record.id)
```

### Reading the first record

Use `read()` to return the first record. Additional records are ignored; a file containing no records raises `ValueError`.

```python
record = fastapy.read("input.fasta")
print(record.id)
```

### Working with records

A `Record` stores a sequence identifier (`id`), sequence (`seq`), and optional description (`desc`). The `description` property returns the complete FASTA header, including `>`.

```python
record = fastapy.Record(id="seq1", seq="MRELEAKAT", desc="Example protein")

print(record.id)           # seq1
print(record.seq)          # MRELEAKAT
print(record.desc)         # Example protein
print(record.description)  # >seq1 Example protein
print(len(record))         # 9
print("LEA" in record)      # True
```

Iterating over a record yields its sequence characters. Use `print(record)` to display FASTA, or `record.format()` to obtain a FASTA string:

```python
print(record.format(wrap=3), end="")
# >seq1 Example protein
# MRE
# LEA
# KAT
```

The default line width is 70 characters. Use `wrap=0` or `wrap=None` for no wrapping.

### Storing records in a list or dictionary

Use a list for repeated access or positional indexing:

```python
records = list(fastapy.parse("input.fasta"))
print(records[0].id)   # First record
print(records[-1].id)  # Last record
```

Use `to_dict()` to index records by identifier. Duplicate identifiers raise `ValueError`.

```python
records = fastapy.to_dict(fastapy.parse("input.fasta"))
print(records["seq1"])  # Replace with an identifier from your file
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