#!/usr/bin/env python3

import pathlib
import unittest

import fastapy as fp

TEST_DIR = pathlib.Path(__file__).resolve().parent

class TestFasta(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.filename = TEST_DIR / 'test.fasta'

    def test_record_id(self):
        lst = ['NP_002433.1', 'ENO94161.1', 'sequence']
        for i, record in enumerate(fp.parse(self.filename)):
            self.assertEqual(record.id, lst[i])

    def test_record_len(self):
        lst = [362, 79, 292]
        for i, record in enumerate(fp.parse(self.filename)):
            self.assertEqual(len(record), lst[i])

    def test_record_desc(self):
        lst = [
            'RNA-binding protein Musashi homolog 1 [Homo sapiens]',
            'RRM domain-containing RNA-binding protein',
            ''
        ]
        for i, record in enumerate(fp.parse(self.filename)):
            self.assertEqual(record.desc, lst[i])

    def test_record_description(self):
        lst = [
            '>NP_002433.1 RNA-binding protein Musashi homolog 1 [Homo sapiens]',
            '>ENO94161.1 RRM domain-containing RNA-binding protein',
            '>sequence'
        ]
        for i, record in enumerate(fp.parse(self.filename)):
            self.assertEqual(record.description, lst[i])

    def test_record_iter(self):
        lst = [list('METDA'), list('MKLLI'), list('MKLSK')]
        for i, record in enumerate(fp.parse(self.filename)):
            self.assertEqual(list(record)[:5], lst[i])

    def test_record_in(self):
        lst = ['METDA', 'MKLLI', 'MKLSK']
        for i, record in enumerate(fp.parse(self.filename)):
            self.assertTrue(lst[i] in record)   

    def test_record_format(self):
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
          # Thir sequence record
          ('>sequence\n',
           'MKLSKIALMMATLAASSAAWSHGYIEVPESRAYKCKLGSNTDCGRAQWEPQSVEQVSGFP\n',
           'GGATPLDGQLASGGVNGFESLDRQGVNVWALNTMKPGPQTFTWYHTAKHKTNNWRYYITK\n',
           'QDWDVNKPLSREAFEKEPFCEIDGHAKPPKDREVHQCVVPERTGYQVIYGVWEDASQPLT\n',
           'MALSNVSEGHHMLKVIASNDNGQSIQPDIENFNLEAESTGGGGGNGDYNFVFPNALSKYT\n',
           'AGTTVLQPKDGKVYQCKPFPYSGYCMQWNSGATHFEPGVGSNWQDAWILKK*\n')
        ]
        for i, record in enumerate(fp.parse(self.filename)):
            self.assertEqual(record.format(wrap=60), "".join(lst[i]))

    def test_get_compression_type_plain(self):
        self.assertIsNone(fp.get_compression_type(self.filename))

    def test_get_compression_type_gz(self):
        file_type = fp.get_compression_type(TEST_DIR / 'test.fasta.gz')
        self.assertEqual(file_type, 'gz')

    def test_get_compression_type_bz2(self):
        file_type = fp.get_compression_type(TEST_DIR / 'test.fasta.bz2')
        self.assertEqual(file_type, 'bz2')

    def test_get_compression_type_zip(self):
        file_type = fp.get_compression_type(TEST_DIR / 'test.fasta.zip')
        self.assertEqual(file_type, 'zip')

    def test_parse_fasta_file(self):
        lst = [r.id for r in fp.parse(TEST_DIR / 'test.fasta.gz')]
        self.assertEqual(len(lst), 3)
        self.assertEqual(set(lst), {'NP_002433.1', 'ENO94161.1', 'sequence'})

    def test_parse_empty_file(self):
        lst = [rec for rec in fp.parse(TEST_DIR / 'empty_file.fasta')]
        self.assertEqual(len(lst), 0)

    def test_parse_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            list(fp.parse("non_existent_file.fasta"))

    def test_parse_gz_file(self):
        record = list(fp.parse(TEST_DIR / 'test.fasta.gz'))[0]
        self.assertEqual(record.id, 'NP_002433.1')
        self.assertEqual(len(record), 362)

    def test_parse_bz2_file(self):
        record = list(fp.parse(TEST_DIR / 'test.fasta.bz2'))[0]
        self.assertEqual(record.id, 'NP_002433.1')
        self.assertEqual(len(record), 362)

    def test_parse_zip_file(self):
        record = list(fp.parse(TEST_DIR / 'test.fasta.zip'))[0]
        self.assertEqual(record.id, 'NP_002433.1')
        self.assertEqual(len(record), 362)

    def test_read(self):
        record = fp.read(TEST_DIR / 'test.fasta.zip')
        self.assertEqual(record.id, 'NP_002433.1')
        self.assertEqual(len(record), 362)        

    def test_to_dict(self):
        d = fp.to_dict(fp.parse(self.filename))
        self.assertEqual(len(d), 3)
        self.assertEqual(len(d['ENO94161.1']), 79)

    def test_to_dict_duplicate_records(self):
        records = [
          fp.Record(id='id1', seq='ATGC'),
          fp.Record(id='id2', seq='CGTA'),
          fp.Record(id='id1', seq='ATGC'),
        ]
        with self.assertRaises(ValueError):
            fp.to_dict(records)

    def test_to_dict_empty_records(self):
        self.assertEqual(fp.to_dict([]), {})


unittest.main()