import pathlib
import pytest

import fastapy as fp

TEST_DIR = pathlib.Path(__file__).resolve().parent
TEST_FASTA = TEST_DIR / 'test.fasta'


def test_record_id():
    lst = ['NP_002433.1', 'ENO94161.1', 'sequence']
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert record.id == lst[i]


def test_record_len():
    lst = [362, 79, 292]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert len(record) == lst[i]


def test_record_desc():
    lst = [
        'RNA-binding protein Musashi homolog 1 [Homo sapiens]',
        'RRM domain-containing RNA-binding protein',
        ''
    ]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert record.desc == lst[i]


def test_record_description():
    lst = [
        '>NP_002433.1 RNA-binding protein Musashi homolog 1 [Homo sapiens]',
        '>ENO94161.1 RRM domain-containing RNA-binding protein',
        '>sequence'
    ]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert record.description == lst[i]


def test_record_iter():
    lst = [list('METDA'), list('MKLLI'), list('MKLSK')]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert list(record)[:5] == lst[i]


def test_record_in():
    lst = ['METDA', 'MKLLI', 'MKLSK']
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert lst[i] in record


def test_record_format():
    lst = [
      # First sequence record
      ('>NP_002433.1 RNA-binding protein Musashi homolog 1 [Homo sapiens]\n',
       'METDAPQPGLASPDSPHDPCKMFIGGLSWQTTQEGLREYFGQFGEVKECLVMRDPLTKRS\n',
       'RGFGFVTFMDQAGVDKVLAQSRHELDSKTIDPKVAFPRRAQPKMVTRTKKIFVGGLSVNT\n',
       'TVEDVKQYFEQFGKVDDAMLMFDKTTNRHRGFGFVTFESEDIVEKVCEIHFHEINNKMVE\n',
       'CKKAQPKEVMSPTGSARGRSRVMPYGMDAFMLGIGMLGYPGFQATTYASRSYTGLAPGYT\n',
       'YQFPEFRVERTPLPSAPVLPELTAIPLTAYGPMAAAAAAAAVVRGTGSHPWTMAPPPGST\n',
       'PSRTGGFLGTTSPGPMAELYGAANQDSGVSSYISAASPAPSTGFGHSLGGPLIATAFTNG\n',
       'YH\n'),
      # Second sequence record
      ('>ENO94161.1 RRM domain-containing RNA-binding protein\n',
       'MKLLISGLGPDTDLDTLRERMRHFGPVLDILVVREGDPERPWFIIDMDITPDVATEVARR\n',
       'IDGIYFHGSFVHARVMLHD\n'),
      # Third sequence record
      ('>sequence\n',
       'MKLSKIALMMATLAASSAAWSHGYIEVPESRAYKCKLGSNTDCGRAQWEPQSVEQVSGFP\n',
       'GGATPLDGQLASGGVNGFESLDRQGVNVWALNTMKPGPQTFTWYHTAKHKTNNWRYYITK\n',
       'QDWDVNKPLSREAFEKEPFCEIDGHAKPPKDREVHQCVVPERTGYQVIYGVWEDASQPLT\n',
       'MALSNVSEGHHMLKVIASNDNGQSIQPDIENFNLEAESTGGGGGNGDYNFVFPNALSKYT\n',
       'AGTTVLQPKDGKVYQCKPFPYSGYCMQWNSGATHFEPGVGSNWQDAWILKK*\n')
    ]
    for i, record in enumerate(fp.parse(TEST_FASTA)):
        assert record.format(wrap=60) == "".join(lst[i])


def test_get_compression_type_plain():
    assert fp.get_compression_type(TEST_FASTA) is None


def test_get_compression_type_gz():
    file_type = fp.get_compression_type(TEST_DIR / 'test.fasta.gz')
    assert file_type == 'gz'


def test_get_compression_type_bz2():
    file_type = fp.get_compression_type(TEST_DIR / 'test.fasta.bz2')
    assert file_type == 'bz2'


def test_get_compression_type_zip():
    file_type = fp.get_compression_type(TEST_DIR / 'test.fasta.zip')
    assert file_type == 'zip'


def test_parse_fasta_file():
    lst = [r.id for r in fp.parse(TEST_DIR / 'test.fasta.gz')]
    assert len(lst) == 3
    assert set(lst) == {'NP_002433.1', 'ENO94161.1', 'sequence'}


def test_parse_empty_file():
    lst = [rec for rec in fp.parse(TEST_DIR / 'empty_file.fasta')]
    assert len(lst) == 0


def test_parse_missing_file():
    with pytest.raises(FileNotFoundError):
        list(fp.parse("non_existent_file.fasta"))


def test_parse_gz_file():
    record = list(fp.parse(TEST_DIR / 'test.fasta.gz'))[0]
    assert record.id == 'NP_002433.1'
    assert len(record) == 362


def test_parse_bz2_file():
    record = list(fp.parse(TEST_DIR / 'test.fasta.bz2'))[0]
    assert record.id == 'NP_002433.1'
    assert len(record) == 362


def test_parse_zip_file():
    record = list(fp.parse(TEST_DIR / 'test.fasta.zip'))[0]
    assert record.id == 'NP_002433.1'
    assert len(record) == 362


def test_read():
    record = fp.read(TEST_FASTA)
    assert record.id == 'NP_002433.1'
    assert len(record) == 362


def test_read_gz_file():
    record = fp.read(TEST_DIR / 'test.fasta.gz')
    assert record.id == 'NP_002433.1'
    assert len(record) == 362


def test_read_zip_file():
    record = fp.read(TEST_DIR / 'test.fasta.zip')
    assert record.id == 'NP_002433.1'
    assert len(record) == 362


def test_read_bz2_file():
    record = fp.read(TEST_DIR / 'test.fasta.bz2')
    assert record.id == 'NP_002433.1'
    assert len(record) == 362


def test_to_dict():
    d = fp.to_dict(fp.parse(TEST_FASTA))
    assert len(d) == 3
    assert len(d['ENO94161.1']) == 79


def test_to_dict_duplicate_records():
    records = [
        fp.Record(id='id1', seq='ATGC'),
        fp.Record(id='id2', seq='CGTA'),
        fp.Record(id='id1', seq='ATGC'),
    ]
    with pytest.raises(ValueError):
        fp.to_dict(records)


def test_to_dict_empty_records():
    assert fp.to_dict([]) == {}
