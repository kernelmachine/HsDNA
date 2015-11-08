module RepUtils where

import GenStructs
import GenUtils
import DNAStruct
import HigherOrderStructs

import Control.Applicative
import Data.List 
import qualified Data.List.Utils as LU

adptoATP :: DNA -> DNA -- Activation 
adptoATP dna = dna {dnaATP = True}

atptoADP :: DNA -> DNA -- Dectivation 
atptoADP dna = dna {dnaATP = False}

bindHelicase :: DNA -> DNA 
bindHelicase dna = dna {helicaseBound = True}

unBindHelicase :: DNA -> DNA
unBindHelicase dna = dna {helicaseBound = False}

split :: [(Nucleotide, Nucleotide)] -> ([Nucleotide], [Nucleotide]) -- Change representation of DNA strands from list of tuples to tuple of lists
split [] = ([], [])
split (x:xs) = (map fst (x:xs), map snd (x:xs))

helicase :: DNA -> DNA -- input helicase-unzipped dna into separated Strands field for easy access.
helicase dna = dna {separatedStrands = split (strands dna)}

getReplicationOrigin :: DNA -> DNA -- replication origin begins at repeated sequence of "A=T" bindings, and ends 50 nucleotides later (just a guess)
getReplicationOrigin dna = dna { replicationOrigin = (x, (+50) <$> x) }
				where x = head $ findATATAT (compressPairs $ strands dna)

bindPolymerase :: DNA -> DNA 
bindPolymerase dna = dna {polymeraseBounded = True}

unBindPolymerase :: DNA -> DNA 
unBindPolymerase dna = dna {polymeraseBounded = False}

duplicateDNA :: DNA -> [DNA] -- semi-conservative replication. 
duplicateDNA dna = [dna{strands = createDNA (fst $ separatedStrands dna) (complementarySequence $ fst $ separatedStrands dna), 
						separatedStrands = ([],[])}, 
						dna{strands = createDNA (snd $ separatedStrands dna) (complementarySequence $ snd $ separatedStrands dna), 
						separatedStrands = ([],[])} ] 

topisomeraseII :: DNA -> DNA -- topisomeraseII positively winds the DNA after replication to decrease strain. 
topisomeraseII dna = wind (2) dna