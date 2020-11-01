import unittest
from karma.kmer import KmerClustering

class TestKmer:

    def test_is_palindrome(self):
        assert KmerClustering.is_palindrome("ACGT") is False
        assert KmerClustering.is_palindrome("AAAA") is True
