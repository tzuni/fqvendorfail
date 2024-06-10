#!/usr/bin/env python3

import unittest

from fqvendorfail.vfail import (
    Block,
    get_outfile,
    get_outfile_dict,
    iter_linescan,
    vendorfail,
)

from fqvendorfail import __main__ as MOD


class ThisTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.runner = CliRunner()

    def test_pass(self):
        result = self.runner.invoke(MOD.main)
        self.assertEqual(result.exit_code, 0)


    assert result == expected


def test_full_fracfail():
    fastq_list = [fracfail_R1, fracfail_R2]
    outprefix = 'testout'
    outfiles = [f'{outprefix}_R1.fq.gz', f'{outprefix}_R2.fq.gz']

    ex_R1 = b'''@K00270:69:HN2HGBBXX:1:1101:1722:1244 1:N:0:CGTACG
NCTTTTTATTAGGAACCAGGGGAATGAGCTGCTTATCCCTCTATAACAGTCTAGAGCAGGTCATCAGGCCCAGGA
+
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJ7FAJJJFJJJJJJJJJJJJJFJJJJJJJFJJJJJJJJ7FJJJJJJA
@K00270:69:HN2HGBBXX:1:1101:1722:1244 1:N:0:CGTACG
NCTTTTTATTAGGAACCAGGGGAATGAGCTGCTTATCCCTCTATAACAGTCTAGAGCAGGTCATCAGGCCCAGGA
+
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJ7FAJJJFJJJJJJJJJJJJJFJJJJJJJFJJJJJJJJ7FJJJJJJA
@K00270:69:HN2HGBBXX:1:1101:1722:1244 1:N:0:CGTACG
NCTTTTTATTAGGAACCAGGGGAATGAGCTGCTTATCCCTCTATAACAGTCTAGAGCAGGTCATCAGGCCCAGGA
+
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJ7FAJJJFJJJJJJJJJJJJJFJJJJJJJFJJJJJJJJ7FJJJJJJA
@K00270:69:HN2HGBBXX:1:1101:1722:1244 1:N:0:CGTACG
NCTTTTTATTAGGAACCAGGGGAATGAGCTGCTTATCCCTCTATAACAGTCTAGAGCAGGTCATCAGGCCCAGGA
+
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJ7FAJJJFJJJJJJJJJJJJJFJJJJJJJFJJJJJJJJ7FJJJJJJA
@K00270:69:HN2HGBBXX:1:1101:1722:1244 1:N:0:CGTACG
NCTTTTTATTAGGAACCAGGGGAATGAGCTGCTTATCCCTCTATAACAGTCTAGAGCAGGTCATCAGGCCCAGGA
+
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJ7FAJJJFJJJJJJJJJJJJJFJJJJJJJFJJJJJJJJ7FJJJJJJA
'''
    ex_R2 = b'''@K00270:69:HN2HGBBXX:1:1101:1722:1244 2:N:0:CGTACG
GTTTTGACTAACACCAGTTCCTGCCAACCTCTGTTGCCACCACCTTTCCTTCCAGGCCCTAAGCACGTGCAGCAA
+
AAFFFJFJJJJJJJJJJJJJJJJJJ7JJJJJJJJJJJJJJJJJJFFJJFJJFJJJ<FAJJJJJFFJAAJJFAJA-
@K00270:69:HN2HGBBXX:1:1101:1722:1244 2:N:0:CGTACG
GTTTTGACTAACACCAGTTCCTGCCAACCTCTGTTGCCACCACCTTTCCTTCCAGGCCCTAAGCACGTGCAGCAA
+
AAFFFJFJJJJJJJJJJJJJJJJJJ7JJJJJJJJJJJJJJJJJJFFJJFJJFJJJ<FAJJJJJFFJAAJJFAJA-
@K00270:69:HN2HGBBXX:1:1101:1722:1244 2:N:0:CGTACG
GTTTTGACTAACACCAGTTCCTGCCAACCTCTGTTGCCACCACCTTTCCTTCCAGGCCCTAAGCACGTGCAGCAA
+
AAFFFJFJJJJJJJJJJJJJJJJJJ7JJJJJJJJJJJJJJJJJJFFJJFJJFJJJ<FAJJJJJFFJAAJJFAJA-
@K00270:69:HN2HGBBXX:1:1101:1722:1244 2:N:0:CGTACG
GTTTTGACTAACACCAGTTCCTGCCAACCTCTGTTGCCACCACCTTTCCTTCCAGGCCCTAAGCACGTGCAGCAA
+
AAFFFJFJJJJJJJJJJJJJJJJJJ7JJJJJJJJJJJJJJJJJJFFJJFJJFJJJ<FAJJJJJFFJAAJJFAJA-
@K00270:69:HN2HGBBXX:1:1101:1722:1244 2:N:0:CGTACG
GTTTTGACTAACACCAGTTCCTGCCAACCTCTGTTGCCACCACCTTTCCTTCCAGGCCCTAAGCACGTGCAGCAA
+
AAFFFJFJJJJJJJJJJJJJJJJJJ7JJJJJJJJJJJJJJJJJJFFJJFJJFJJJ<FAJJJJJFFJAAJJFAJA-
'''
    expected = [ex_R1.splitlines(keepends=True), ex_R2.splitlines(keepends=True)]

    with mock.patch.object(logger, 'debug') as mock_debug:
        vendorfail(fastq_list, outprefix, mock_debug)

    results = []
    for fq in outfiles:
        with gzip.open(fq) as fh:
            results.append(fh.readlines())
        os.remove(fq)
    assert results == expected


def test_get_outfile():
    prefix = 'testprefix'
    ordinal = 3
    result = get_outfile(prefix, ordinal)
    expected = 'testprefix_R3.fq.gz'
    assert result == expected


def test_get_outfile_dict():
    fq_list = ['somename_R1.fq.gz', 'somename_R2.fq.gz']
    out_prefix = 'testprefix'
    pair_ordinals = ['1', '2']

    expected = {
        fq_list[0]: f'{out_prefix}_R1.fq.gz',
        fq_list[1]: f'{out_prefix}_R2.fq.gz',
    }

    result = get_outfile_dict(fq_list, out_prefix, pair_ordinals)

    assert result == expected
