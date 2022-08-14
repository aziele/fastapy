#!/usr/bin/env python3

import pathlib
import unittest

import fasta


class TestFasta(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_dir = pathlib.Path('test')
        cls.filename = cls.test_dir / 'test.fasta'

    def test_record_id(self):
        lst = ['NP_002433.1', 'ENO94161.1', 'sequence']
        for i, record in enumerate(fasta.parse(self.filename)):
            self.assertEqual(record.id, lst[i])

    def test_record_len(self):
        lst = [362, 79, 292]
        for i, record in enumerate(fasta.parse(self.filename)):
            self.assertEqual(len(record), lst[i])

    def test_record_desc(self):
        lst = [
            'RNA-binding protein Musashi homolog 1 [Homo sapiens]',
            'RRM domain-containing RNA-binding protein',
            ''
        ]
        for i, record in enumerate(fasta.parse(self.filename)):
            self.assertEqual(record.desc, lst[i])

    def test_record_description(self):
        lst = [
            '>NP_002433.1 RNA-binding protein Musashi homolog 1 [Homo sapiens]',
            '>ENO94161.1 RRM domain-containing RNA-binding protein',
            '>sequence'
        ]
        for i, record in enumerate(fasta.parse(self.filename)):
            self.assertEqual(record.description, lst[i])

    def test_record_iter(self):
        lst = [list('METDA'), list('MKLLI'), list('MKLSK')]
        for i, record in enumerate(fasta.parse(self.filename)):
            self.assertEqual(list(record)[:5], lst[i])

    def test_record_in(self):
        lst = ['METDA', 'MKLLI', 'MKLSK']
        for i, record in enumerate(fasta.parse(self.filename)):
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
        for i, record in enumerate(fasta.parse(self.filename)):
            self.assertEqual(record.format(wrap=60), "".join(lst[i]))

    def test_parse_gz_file(self):
        record = list(fasta.parse(self.test_dir / 'test.fasta.gz'))[0]
        self.assertEqual(record.id, 'NP_002433.1')
        self.assertEqual(len(record), 362)

    def test_read_bz2_file(self):
        record = list(fasta.parse(self.test_dir / 'test.fasta.bz2'))[0]
        self.assertEqual(record.id, 'NP_002433.1')
        self.assertEqual(len(record), 362)

    def test_to_dict1(self):
        d = fasta.to_dict(fasta.parse(self.filename))
        self.assertEqual(len(d), 3)
        self.assertEqual(len(d['ENO94161.1']), 79)

unittest.main()