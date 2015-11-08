module HigherOrderStructs where

import GenUtils
import GenStructs
import DNAStruct

import Data.List
import Control.Applicative 

createNucleosome :: DNA -> Maybe Int ->  Nucleosome
createNucleosome k pos = Nucleosome { startPos = pos,
									endPos = lengthOfNucleosome
									}
						where 
							lengthOfNucleosome = if d < wrapMaybe (length $ strands k) then d else Just $ (length $ strands k)
							d = (+ 146) <$> pos
							wrapMaybe a = Just a

-- Every nucleosomeAssembly consists of nucleosomes mapped to repeated sequences of "A=T" bindings. 
-- Each nucleosomeassembly binding decreases the linking number by 1.	

createNucleosomeAssembly :: DNA -> NucleosomeAssembly 
createNucleosomeAssembly k = NucleosomeAssembly {  dna = wind (-1) k,  nucleosomes = map (\x -> createNucleosome k x) positions}
	where
		positions = findATATAT (compressPairs $ strands k)

createNucleosomeAssemblies :: [DNA] -> [NucleosomeAssembly]
createNucleosomeAssemblies = map createNucleosomeAssembly 

-- Sequential assembly of higher order DNA structures

createThirtyNMLoop :: [NucleosomeAssembly] -> ThirtyNMLoop
createThirtyNMLoop ns = ThirtyNMLoop { nucleosomeAssemblies = ns} 

createRosette :: [ThirtyNMLoop] -> Rosette
createRosette ts = Rosette { thirtyNMLoops = ts}

createCoil :: [Rosette] -> Coil
createCoil rs = Coil { rosettes = rs}

createChromatid :: [Coil] -> Chromatid
createChromatid cl = Chromatid { coils =  cl }

createChromosome :: [Chromatid] -> Chromosome
createChromosome cd = Chromosome { chromatids = cd }

createChromosomeSet :: [Chromosome]  -> ChromosomeSet
createChromosomeSet cs = ChromosomeSet { chromosomes = cs}

-- create a chromosome from a sequence of genes.

chromosome    ::    [[[[[DNA]]]]] -> Chromosome
chromosome dnaSq = 
	let 
		ns dnaSq                =  (map . map . map . map . map) createNucleosomeAssembly dnaSq
		loops nucleoassemblies  =  (map . map . map . map) createThirtyNMLoop nucleoassemblies
		roses ts                =  (map . map . map) createRosette ts
		coilings rs             =  (map . map) createCoil rs
		ctids cs                =  (map) createChromatid cs
	in 
		createChromosome $ ctids . coilings . roses . loops . ns $ dnaSq 

-- create a chromosome set from a sequence of sequences of genes.

chromosomeSet :: [[[[[[DNA]]]]]] -> ChromosomeSet
chromosomeSet dnaSq = createChromosomeSet $ map chromosome dnaSq 

-- Access the DNA of a chromosome. 

getChromosomeDNA :: Chromosome -> [[[[[DNA]]]]]
getChromosomeDNA chromosome = 
	let 
		coilings cds     = (map) coils cds 
		roses cs         = (map . map) rosettes cs
		loops rs 		 = (map . map . map) thirtyNMLoops rs
		ns ls 		     = (map . map . map . map) nucleosomeAssemblies ls
		ds n             = (map . map . map . map . map) dna n
 	in 
 		ds $ ns . loops . roses . coilings .  chromatids $ chromosome

-- Access the DNA of a chromosomeSet 

getAllDNA :: ChromosomeSet -> [[[[[[DNA]]]]]]
getAllDNA chrset = map getChromosomeDNA (chromosomes chrset)