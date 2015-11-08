module DNAStruct where

import GenUtils
import GenStructs

import Data.List
import Control.Applicative 

createDNA :: [Nucleotide] -> [Nucleotide] -> [(Nucleotide, Nucleotide)] -- creates a DNA molecule from two (complementary) nucleotide sequences. 
createDNA firstStrand secondStrand 
		| check firstStrand secondStrand = zipWith (,) firstStrand secondStrand 
		| otherwise = error "Strands are not aligned properly."
			
ln0 :: DNA -> Maybe Double -- linking number in bForm DNA is calculated by dividing the length of the dna by the number of base pairs per turn. 
ln0 dna = (/) <$> (fromIntegral <$> lengthOfBasePairs) <*> (basePairsPerTurn dna)
			where 
				lengthOfBasePairs = Just $ (length $ strands dna) 

superHelicalDensity :: DNA -> Maybe Double -- superhelical density of DNA. See README for details. 
superHelicalDensity a = (/) <$> ((-) <$> linkno <*> relaxLinkNo) <*> (relaxLinkNo) 
			where
				linkno = linkingNumber a
				relaxLinkNo = ln0 bForm

typeOfSuperCoil :: DNA -> Maybe Double -- returns (-1) if DNA is a negative supercoil, and 1 if DNA is a positive superCoil. 
typeOfSuperCoil dna = signum <$> superHelicalDensity dna

bpPerTurn :: DNA -> Maybe Double  -- the number of base pairs per helical turn of DNA, calculated by dividing the length of the dna by the linking number.
bpPerTurn dna = (/) <$> (fromIntegral <$> lengthOfBasePairs ) <*> (linkingNumber dna)
			where 
				lengthOfBasePairs = Just $ (length $ strands dna) 

wind :: Double -> DNA -> DNA -- wind the dna in the positive or negative direction by increasing or decreasing the linking number, respectively. 
wind lnChange dna = 
	let 
		d = dna {linkingNumber = (+ lnChange) <$> (linkingNumber dna)}
		g = d {basePairsPerTurn = bpPerTurn d} 
	in 
		if (basePairsPerTurn g) >= (Just 1) then g else error "Too many twists -- You've broken your DNA!"

bForm :: DNA -- DNA in bform, with a linking number calculated by ln0, and 10.5 base pairs per helical turn. 
bForm = DNA {strands = createDNA firstStrand secondStrand, separatedStrands = ([], []), linkingNumber = ln0 bForm ,
			 basePairsPerTurn = Just 10.5, dnaATP = False, dnaADP = True,  helicaseBound = False, 
			 replicationOrigin = (Nothing,Nothing), polymeraseBounded = False}
		where 
			firstStrand  = "ACAAGATGCCATTGTACCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCGCTGCCCTGCCCCTGGAGGGTGGCCCCACCGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAATAAGGAAAAGCAGCCTCCTGACTTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGAGAGGAAGCTCGGGAGGTGGCCAGGCGGCAGGAAGGCGCACCCCCCCAGCAATCCGCGCGCCGGGACAGAATGCCCTGCAGGAACTTCTTCTGGAAGACCTTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAGTTTAATTACAGACCTGAA"
			secondStrand = complementarySequence firstStrand