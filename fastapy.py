"""A lightweight Python module to read and write FASTA sequence records"""

import bz2
import gzip
import io
import pathlib
import typing
import zipfile

__version__ = "1.1.0"

DEFAULT_WRAP = 70


class Record:
    """Object representing a FASTA (aka Pearson) record.

    Attributes:
        id  (str)         : Sequence identifier
        seq (str)         : Sequence
        description (str) : Description line (defline)
    """

    def __init__(self, id: str, seq: str, desc: typing.Optional[str] = None):
        """Creates a Record.

        Example:
        >>> record = Record(id="NP_055309.2", 
        ...                 seq="MRELEAKAT",
        ...                 desc="TNRC6A")
        >>> print(record)
        >NP_055309.2 TNRC6A
        MRELEAKAT
        """
        self.id = id
        self.seq = seq
        self.desc = desc

    @property
    def description(self) -> str:
        """Returns a description line (defline) of FASTA record.

        Example:
        >>> record = Record(id="NP_055309.2", seq="MRELEAKAT", desc="TNRC6A")
        >>> print(record.description)
        >NP_055309.2 TNRC6A

        >>> record = Record(id="seqid", seq="ATCGA")
        >>> print(record.description)
        >seqid
        """
        lst = [f">{self.id}"]
        if self.desc:
            lst.append(f"{self.desc}")
        return " ".join(lst)

    def __iter__(self):
        """Iterates over the characters in the sequence.

        Example:
        >>> record = Record(id="NP_055309.2", seq="MRELEAKAT", desc="TNRC6A")
        >>> for amino_acid in record:
        ...     print(amino_acid)
        M
        R
        E
        L
        E

        This is equivalent to iterating over the sequence directly:
        >>> for amino_acid in record.seq:
        ...     print(amino_acid)
        M
        R
        E
        L
        E
        """
        return iter(self.seq)

    def __contains__(self, char):
        """Implements the 'in' keyword to search the sequence.

        Example:
        >>> record = Record(id="NP_055309.2", seq="MRELEAKAT", desc="TNRC6A")
        >>> print("M" in record)
        True
        """
        return char in self.seq

    def __str__(self):
        """Returns the record as a string in the FASTA format.

        Example:
        >>> record = Record(id="NP_055309.2", seq="MRELEAKAT", desc="TNRC6A")
        >>> print(record)
        >NP_055309.2 TNRC6A
        MRELEAKAT
        """
        return self.format().rstrip()

    def __len__(self):
        """Returns the length of the sequence.

        Example:
        >>> record = Record(id="NP_055309.2", seq="MRELEAKAT")
        >>> len(record)
        9
        """
        return len(self.seq)

    def format(self, wrap:int = DEFAULT_WRAP) -> str:
        """Returns a formatted FASTA record.

        Args:
        wrap (int): line length to wrap sequence lines.
          Default: 70 characters, use zero (or None) for no wrapping.

        Example:
        >>> record = Record(id="NP_055309.2", seq="MRELEAKAT", desc="TNRC6A")

        >>> print(record.format())
        >NP_055309.2 TNRC6A
        MRELEAKAT

        >>> print(record.format(wrap=3))
        >NP_055309.2 TNRC6A
        MRE
        LEA
        KAT
        """
        lst = [self.description, "\n"]
        if wrap:
            for i in range(0, len(self.seq), wrap):
                lst.append(f'{self.seq[i:i + wrap]}\n')
        else:
            lst.append(self.seq)
            lst.append('\n')
        return "".join(lst)


def parse_handle(handle) -> Record:
    """Iterates over `Record` objects from a FASTA file handle.

    Args:
        handle: a handle (file-like objects) that contains FASTA sequences

    Returns:
        A generator that yields `Record` objects.

    Raises:
        This function does not raise any exception.

    Example:
        with open("test/test.fasta") as handle:
            for record in parse_handle(handle):
                print(record.id, record.seq, record.description)
    """
    seqid = None
    desc = None
    seq = []
    for line in handle:
        if line.startswith(">"):
            if seq:
                yield Record(seqid, "".join(seq), desc)
                seq.clear()
            seqid = line.split()[0][1:]
            desc = line[len(seqid) + 1:].strip()
        else:
            seq.append("".join(line.split()))
    if seq:
        yield Record(seqid, "".join(seq), desc)


def parse(filename: typing.Union[str, pathlib.Path]) -> Record:
    """Determines the compression type of a file and yields FASTA records.

    Args:
        filename: a name or pathlib.Path of a file containing FASTA sequences

    Returns:
        A generator that yields `Record` objects.

    Raises:
        FileNotFoundError: If the input file cannot be found.
        OSError: If an operating system error occurs while opening the file.
    """
    compression_type = get_compression_type(filename)
    if compression_type == "bz2":
        with bz2.open(filename, "rt") as fh:
            yield from parse_handle(fh)
    elif compression_type == "gz":
        with gzip.open(filename, "rt") as fh:
            yield from parse_handle(fh)
    elif compression_type == "zip":
        with zipfile.ZipFile(filename) as z:
            # Assuming the first file in the archive is the one we want
            inner_filename = z.namelist()[0]
            with io.TextIOWrapper(z.open(inner_filename)) as fh:
                yield from parse_handle(fh)
    else:
        with open(filename, "rt") as fh:
            yield from parse_handle(fh)


def read(filename: typing.Union[str, pathlib.Path]) -> Record:
    """Reads a single `Record` object from a FASTA file.
    
    Args:
        filename: a name or pathlib.Path of a file containing FASTA sequences

    Returns:
        A single `Record` object.

    Raises:
        ValueError: If no FASTA records are found in the file.
    """
    iterator = parse(filename)
    try:
        record = next(iterator)
    except StopIteration:
        raise ValueError(f"No records found in file: '{filename}'") from None
    return record


def write(
    records: typing.Union[Record, typing.Iterable[Record]],
    destination: typing.Union[str, pathlib.Path, typing.TextIO],
    wrap: typing.Optional[int] = DEFAULT_WRAP,
) -> int:
    """Writes FASTA records to a file or text handle.

    Args:
        records: a single `Record` or an iterable of records.
        destination: a filename, pathlib.Path, or writable text handle.
            Supports plain, gzip, and bzip2 files. Files are overwritten;
            supplied handles remain open.
        wrap (int): line length, default 70; use zero or None for no wrapping.

    Returns:
        The number of records written.

    Example:
    >>> record = Record(id="seq1", seq="ATCG")
    >>> write(record, "sequence.fasta")
    1
    >>> records = [record, Record(id="seq2", seq="GGTA")]
    >>> write(records, "sequences.fasta")
    2
    >>> records = to_dict(records)
    >>> write(records.values(), "sequences.fasta.gz")
    2
    """
    if wrap is not None:
        if isinstance(wrap, bool) or not isinstance(wrap, int):
            raise TypeError("wrap must be an integer or None")
        if wrap < 0:
            raise ValueError("wrap must be non-negative")

    if isinstance(records, Record):
        records = (records,)

    def write_records(handle):
        count = 0
        for record in records:
            if not isinstance(record, Record):
                raise TypeError("Expected Record objects in the iterable")
            handle.write(record.format(wrap=wrap))
            count += 1
        return count

    if isinstance(destination, (str, pathlib.Path)):
        compression = get_compression_type(destination)
        if compression == "zip":
            raise ValueError("ZIP output is not supported")
        opener = {"gz": gzip.open, "bz2": bz2.open}.get(compression, open)
        with opener(destination, "wt", encoding="utf-8", newline="\n") as handle:
            return write_records(handle)

    return write_records(destination)


def to_dict(records) -> dict:
    """Turns a generator or list of `Record` objects into a dictionary.

    Args:
        records: a generator that yields `Record` objects, or a list of 
          `Record` objects

    Returns:
        A dict mapping sequence id (key) to `Record` object (value).

    Raises:
        ValuError: If duplicate record ids are found.

    Example:
    >>> import fasta
    >>> record_dict = fasta.to_dict(fasta.parse("test.fasta"))
    >>> print(sorted(record_dict.keys()))
    ["ENO94161.1", "NP_002433.1", "sequence"]
    >>> print(record_dict["ENO94161.1"].description)
    RRM domain-containing RNA-binding protein
    >>> len(pdict)
    3
    """
    d = {}
    for record in records:
        if record.id in d:
            raise ValueError(f"Duplicated key: '{record.id}'")
        d[record.id] = record
    return d


COMPRESSION_EXTENSIONS = {
    ".gz": "gz",
    ".gzip": "gz",
    ".bz2": "bz2",
    ".zip": "zip",
}


def get_compression_type(
        filename: typing.Union[str, pathlib.Path]
    ) -> typing.Union[str, None]:
    """Returns the compression type of a file based on its extension.

    Args:
        filename: a name or pathlib.Path of a file containing FASTA sequences

    Returns:
        A string representing the compression type of the file, or None 
        if the compression type could not be determined.
    """
    ext = pathlib.Path(filename).suffix.lower()
    return COMPRESSION_EXTENSIONS.get(ext)