module Test where

import GenUtils
import GenStructs
import DNAStruct
import HigherOrderStructs
import RepUtils
import Replication
import Control.Applicative
import qualified Data.List as L

geneA :: DNA
geneA = DNA {strands = createDNA firstStrand secondStrand, separatedStrands = ([], []), linkingNumber = ln0 bForm , basePairsPerTurn = Just 5.7, 
			dnaATP = False, dnaADP = True,  helicaseBound = False, replicationOrigin = (Nothing,Nothing), polymeraseBounded = False}
		where 
			firstStrand  = "ACAAGATGCCATTGTACCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCGCTGCCCTGCCCCTGGAGGGTGGCCCCACCGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAATAAGGAAAAGCAGCCTCCTGACTTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGAGAGGAAGCTCGGGAGGTGGCCAGGCGGCAGGAAGGCGCACCCCCCCAGCAATCCGCGCGCCGGGACAGAATGCCCTGCAGGAACTTCTTCTGGAAGACCTTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAGTTTAATTACAGACCTGAA"
			secondStrand = complementarySequence firstStrand

geneB :: DNA
geneB = DNA {strands = createDNA firstStrand secondStrand, separatedStrands = ([], []), linkingNumber = ln0 bForm , basePairsPerTurn = Just 6.2, 
			dnaATP = False, dnaADP = True,  helicaseBound = False, replicationOrigin = (Nothing,Nothing), polymeraseBounded = False}
		where 
			firstStrand  = "GACTGATCGATCGATGGCTAGATCGATCGATTTAGAGATCGCTAGCTAGGGGGGATAGCTAGCTAGCTAGCTAGCTAGCCCGCTAGCTAGATCGATCGACTACGACTACGACTGGGGGGGATAGTCGGCGCGCGAGCGACTACGTGACTACTGCATCGATCGACTAGCTAGCTACGACTGACGATCGATCGATCGATCGACGGCGCATCAATATAGCTAGCTAGGGATATAAAAAAAAAAAAAAAAATAGATAAGGAATTAAAGAATCGCGCCGCGCGCGCGCGATATAGCTACAGCTACG"
			secondStrand = complementarySequence firstStrand

geneC :: DNA
geneC = DNA {strands = createDNA firstStrand secondStrand, separatedStrands = ([], []), linkingNumber = ln0 bForm , basePairsPerTurn = Just 1.5, 
			dnaATP = False, dnaADP = True,  helicaseBound = False, replicationOrigin = (Nothing,Nothing), polymeraseBounded = False}
		where 
			firstStrand  = "TAGATCGCTAGCTAGCTAGCTACGATCAGCATACGACTCGATCGACTAGACTCAGCGCGCGCGCGCATCGACACATCAGCATCAGCCTCGGGGGAGGGAAAAGAAGAAATCGACAGCACTACGAGCTCCGCGCGCGCTATATATATAGCGCATCAGCTACGATCGATCGACGCGGCATCGATCGACTACGCTCCGATCGACATCGCCTGCGCGCGCGCGATATATAGCATCGACTAGC"
			secondStrand = complementarySequence firstStrand

dnaSq :: [[[[[[DNA]]]]]]
dnaSq = [[[[[[geneA,geneB,geneC,bForm,bForm,bForm], [bForm,bForm, geneA,geneB]]]]], [[[[[geneA, bForm], [geneB, geneC]]]]]]

mitosis :: ChromosomeSet -> [ChromosomeSet]
mitosis chromosomeSet = return chromosomeSet >>= replicateCell >>= replicateCell >>= replicateCell

meiosisI :: ChromosomeSet -> [ChromosomeSet]
meiosisI chromosomeSet = return chromosomeSet >>= crossoverCell >>= replicateCell 

meiosisII :: ChromosomeSet -> [ChromosomeSet]
meiosisII chromosomeSet = return chromosomeSet >>= crossoverCell >>= replicateCell >>= crossoverCell >>= replicateCell