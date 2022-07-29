import os
import py2bit

class Test():
    fname = os.path.dirname(py2bit.__file__) + "/py2bitTest/foo.2bit"

    def testOpenClose(self):
        tb = py2bit.open(self.fname, True)
        assert(tb is not None)
        tb.close()

    def testChroms(self):
        tb = py2bit.open(self.fname, True)
        chroms = tb.chroms()
        correct = {'chr1': 150, 'chr2': 100}
        for k, v in chroms.items():
            assert(correct[k] == v)
        assert(tb.chroms("chr1") == 150)
        assert(tb.chroms("c") is None)
        tb.close()

    def testInfo(self):
        tb = py2bit.open(self.fname, True)
        correct = {'file size': 161, 'nChroms': 2, 'sequence length': 250, 'hard-masked length': 150, 'soft-masked length': 8}
        check = tb.info()
        assert(len(correct) == len(check))
        for k, v in check.items():
            assert(correct[k] == v)
        tb.close()

    def testSequence(self):
        tb = py2bit.open(self.fname, True)
        assert(tb.sequence("chr1") == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATCGATCGTAGCTAGCTAGCTAGCTGATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
        assert(tb.sequence("chr1", 1, 2) == "N")
        assert(tb.sequence("chr1", 1, 3) == "NN")
        assert(tb.sequence("chr1", 0, 1000) == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATCGATCGTAGCTAGCTAGCTAGCTGATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
        assert(tb.sequence("chr1", 24, 74) == "NNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGTagctagctGATC")
        tb.close()

    def testBases(self):
        tb = py2bit.open(self.fname, True)
        assert(tb.bases("chr1") == {'A': 0.08, 'C': 0.08, 'T': 0.08666666666666667, 'G': 0.08666666666666667})
        assert(tb.bases("chr1", 24, 74) == {'A': 0.12, 'C': 0.12, 'T': 0.12, 'G': 0.12})
        assert(tb.bases("chr1", 24, 74, False) == {'A': 6, 'C': 6, 'T': 6, 'G': 6})
        assert(tb.bases("chr2", 10, 20) == {'A': 0.2, 'C': 0.2, 'T': 0.3, 'G': 0.3})
        assert(tb.bases("chr2", 10, 20, False) == {'A': 2, 'C': 2, 'T': 3, 'G': 3})
        tb.close()

    def testHardMaskedBlocks(self):
        tb = py2bit.open(self.fname, True)
        assert(tb.hardMaskedBlocks("chr1") == [(0, 50), (100, 150)])
        assert(tb.hardMaskedBlocks("chr1", 25, 75) == [(0, 50)])
        assert(tb.hardMaskedBlocks("chr1", 75, 100) == [])
        assert(tb.hardMaskedBlocks("chr1", 75, 101) == [(100, 150)])
        assert(tb.hardMaskedBlocks("chr2") == [(50, 100)])
        tb.close()

    def testSoftMaskedBlocks(self):
        tb = py2bit.open(self.fname, storeMasked=True)
        assert(tb.softMaskedBlocks("chr1") == [(62, 70)])
        assert(tb.softMaskedBlocks("chr1", 0, 50) == [])
        tb.close()
