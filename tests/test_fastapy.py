import bz2
import gzip
import io
import pathlib

import pytest

import fastapy as fp

TEST_DIR = pathlib.Path(__file__).resolve().parent
TEST_FASTA = TEST_DIR / "test.fasta"


def test_parse_sequence_with_internal_whitespace():
    fasta_content = ">seq1 example\nAT CG\tAT\nGG\r\n TA\n"
    handle = io.StringIO(fasta_content)
    records = list(fp.parse_handle(handle))
    assert len(records) == 1
    assert records[0].seq == "ATCGATGGTA"


def test_record_id():
    lst = ["NP_002433.1", "ENO94161.1", "sequence"]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert record.id == lst[i]


def test_record_len():
    lst = [362, 79, 292]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert len(record) == lst[i]


def test_record_desc():
    lst = [
        "RNA-binding protein Musashi homolog 1 [Homo sapiens]",
        "RRM domain-containing RNA-binding protein",
        ""
    ]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert record.desc == lst[i]


def test_record_description():
    lst = [
        ">NP_002433.1 RNA-binding protein Musashi homolog 1 [Homo sapiens]",
        ">ENO94161.1 RRM domain-containing RNA-binding protein",
        ">sequence"
    ]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert record.description == lst[i]


def test_record_iter():
    lst = [list("METDA"), list("MKLLI"), list("MKLSK")]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert list(record)[:5] == lst[i]


def test_record_in():
    lst = ["METDA", "MKLLI", "MKLSK"]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert lst[i] in record


def test_record_format():
    lst = [
      # First sequence record
      (">NP_002433.1 RNA-binding protein Musashi homolog 1 [Homo sapiens]\n",
       "METDAPQPGLASPDSPHDPCKMFIGGLSWQTTQEGLREYFGQFGEVKECLVMRDPLTKRS\n",
       "RGFGFVTFMDQAGVDKVLAQSRHELDSKTIDPKVAFPRRAQPKMVTRTKKIFVGGLSVNT\n",
       "TVEDVKQYFEQFGKVDDAMLMFDKTTNRHRGFGFVTFESEDIVEKVCEIHFHEINNKMVE\n",
       "CKKAQPKEVMSPTGSARGRSRVMPYGMDAFMLGIGMLGYPGFQATTYASRSYTGLAPGYT\n",
       "YQFPEFRVERTPLPSAPVLPELTAIPLTAYGPMAAAAAAAAVVRGTGSHPWTMAPPPGST\n",
       "PSRTGGFLGTTSPGPMAELYGAANQDSGVSSYISAASPAPSTGFGHSLGGPLIATAFTNG\n",
       "YH\n"),
      # Second sequence record
      (">ENO94161.1 RRM domain-containing RNA-binding protein\n",
       "MKLLISGLGPDTDLDTLRERMRHFGPVLDILVVREGDPERPWFIIDMDITPDVATEVARR\n",
       "IDGIYFHGSFVHARVMLHD\n"),
      # Third sequence record
      (">sequence\n",
       "MKLSKIALMMATLAASSAAWSHGYIEVPESRAYKCKLGSNTDCGRAQWEPQSVEQVSGFP\n",
       "GGATPLDGQLASGGVNGFESLDRQGVNVWALNTMKPGPQTFTWYHTAKHKTNNWRYYITK\n",
       "QDWDVNKPLSREAFEKEPFCEIDGHAKPPKDREVHQCVVPERTGYQVIYGVWEDASQPLT\n",
       "MALSNVSEGHHMLKVIASNDNGQSIQPDIENFNLEAESTGGGGGNGDYNFVFPNALSKYT\n",
       "AGTTVLQPKDGKVYQCKPFPYSGYCMQWNSGATHFEPGVGSNWQDAWILKK*\n")
    ]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert record.format(wrap=60) == "".join(lst[i])


def test_get_compression_type_plain():
    assert fp.get_compression_type(TEST_FASTA) is None


def test_get_compression_type_gz():
    file_type = fp.get_compression_type(TEST_DIR / "test.fasta.gz")
    assert file_type == "gz"
    file_type_upper = fp.get_compression_type(TEST_DIR / "test.fasta.GZ")
    assert file_type_upper == "gz"


def test_get_compression_type_bz2():
    file_type = fp.get_compression_type(TEST_DIR / "test.fasta.bz2")
    assert file_type == "bz2"
    file_type_upper = fp.get_compression_type(TEST_DIR / "test.fasta.BZ2")
    assert file_type_upper == "bz2"


def test_get_compression_type_zip():
    file_type = fp.get_compression_type(TEST_DIR / "test.fasta.zip")
    assert file_type == "zip"


def test_parse_fasta_file():
    lst = [r.id for r in fp.parse(TEST_DIR / "test.fasta.gz")]
    assert len(lst) == 3
    assert set(lst) == {"NP_002433.1", "ENO94161.1", "sequence"}


def test_parse_empty_file():
    lst = [rec for rec in fp.parse(TEST_DIR / "empty_file.fasta")]
    assert len(lst) == 0


def test_parse_missing_file():
    with pytest.raises(FileNotFoundError):
        list(fp.parse("non_existent_file.fasta"))


def test_parse_gz_file():
    record = list(fp.parse(TEST_DIR / "test.fasta.gz"))[0]
    assert record.id == "NP_002433.1"
    assert len(record) == 362


def test_parse_bz2_file():
    record = list(fp.parse(TEST_DIR / "test.fasta.bz2"))[0]
    assert record.id == "NP_002433.1"
    assert len(record) == 362


def test_parse_zip_file():
    record = list(fp.parse(TEST_DIR / "test.fasta.zip"))[0]
    assert record.id == "NP_002433.1"
    assert len(record) == 362


def test_read():
    record = fp.read(TEST_FASTA)
    assert record.id == "NP_002433.1"
    assert len(record) == 362


def test_read_gz_file():
    record = fp.read(TEST_DIR / "test.fasta.gz")
    assert record.id == "NP_002433.1"
    assert len(record) == 362


def test_read_zip_file():
    record = fp.read(TEST_DIR / "test.fasta.zip")
    assert record.id == "NP_002433.1"
    assert len(record) == 362


def test_read_bz2_file():
    record = fp.read(TEST_DIR / "test.fasta.bz2")
    assert record.id == "NP_002433.1"
    assert len(record) == 362


def test_to_dict():
    d = fp.to_dict(fp.parse(TEST_FASTA))
    assert len(d) == 3
    assert len(d["ENO94161.1"]) == 79


def test_to_dict_duplicate_records():
    records = [
        fp.Record(id="id1", seq="ATGC"),
        fp.Record(id="id2", seq="CGTA"),
        fp.Record(id="id1", seq="ATGC"),
    ]
    with pytest.raises(ValueError):
        fp.to_dict(records)


def test_to_dict_empty_records():
    assert fp.to_dict([]) == {}


@pytest.mark.parametrize("as_string", [False, True])
def test_write_single_record(tmp_path, as_string):
    path = tmp_path / "output.fasta"
    destination = str(path) if as_string else path
    record = fp.Record(id="seq1", seq="ATCG", desc="example")

    assert fp.write(record, destination) == 1
    assert path.read_bytes() == b">seq1 example\nATCG\n"


@pytest.mark.parametrize("container", ["list", "tuple", "generator", "dict_values"])
def test_write_records(tmp_path, container):
    records = [fp.Record("seq1", "ATCG"), fp.Record("seq2", "GGTA")]

    if container == "tuple":
        records = tuple(records)
    elif container == "generator":
        records = (record for record in records)
    elif container == "dict_values":
        records = fp.to_dict(records).values()

    path = tmp_path / "output.fasta"
    assert fp.write(records, path) == 2
    assert path.read_text() == ">seq1\nATCG\n>seq2\nGGTA\n"


@pytest.mark.parametrize("wrap, expected", [
    (2, ">seq1\nAT\nCG\nA\n"),
    (0, ">seq1\nATCGA\n"),
    (None, ">seq1\nATCGA\n"),
])
def test_write_wrap(tmp_path, wrap, expected):
    path = tmp_path / "output.fasta"

    assert fp.write(fp.Record("seq1", "ATCGA"), path, wrap=wrap) == 1
    assert path.read_text() == expected


def test_write_default_wrap(tmp_path):
    path = tmp_path / "output.fasta"

    assert fp.write(fp.Record("seq1", "A" * (fp.DEFAULT_WRAP + 1)), path) == 1
    assert path.read_text() == ">seq1\n" + "A" * fp.DEFAULT_WRAP + "\nA\n"


@pytest.mark.parametrize("suffix, opener", [
    (".gz", gzip.open),
    (".gzip", gzip.open),
    (".GZ", gzip.open),
    (".bz2", bz2.open),
    (".BZ2", bz2.open),
])
def test_write_compressed(tmp_path, suffix, opener):
    path = tmp_path / ("output.fasta" + suffix)
    records = [fp.Record("seq1", "ATCG"), fp.Record("seq2", "GGTA")]

    assert fp.write(records, path) == 2
    with opener(path, "rt", encoding="utf-8") as handle:
        assert handle.read() == ">seq1\nATCG\n>seq2\nGGTA\n"


def test_write_empty_records(tmp_path):
    path = tmp_path / "output.fasta"

    assert fp.write([], path) == 0
    assert path.read_bytes() == b""


def test_write_overwrites_file(tmp_path):
    path = tmp_path / "output.fasta"
    path.write_text("Previous contents that should be removed.\n")

    assert fp.write(fp.Record("seq1", "AT"), path) == 1
    assert path.read_text() == ">seq1\nAT\n"


@pytest.mark.parametrize("wrap, error", [
    (-1, ValueError),
    (1.5, TypeError),
    ("70", TypeError),
    (True, TypeError),
])
def test_write_invalid_wrap(tmp_path, wrap, error):
    path = tmp_path / "output.fasta"

    with pytest.raises(error):
        fp.write(fp.Record("seq1", "ATCG"), path, wrap=wrap)
    assert not path.exists()


def test_write_zip_not_supported(tmp_path):
    path = tmp_path / "output.fasta.zip"

    with pytest.raises(ValueError, match="ZIP"):
        fp.write(fp.Record("seq1", "ATCG"), path)
    assert not path.exists()


@pytest.mark.parametrize("records", [
    ["ATCG"],
    [None],
    {"seq1": fp.Record("seq1", "ATCG")},
])
def test_write_invalid_records(records):
    handle = io.StringIO()

    with pytest.raises(TypeError, match="Expected Record objects in the iterable"):
        fp.write(records, handle)
    assert not handle.closed