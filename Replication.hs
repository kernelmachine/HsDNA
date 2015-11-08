module Replication where

import GenUtils
import GenStructs
import DNAStruct
import HigherOrderStructs
import RepUtils 

import Control.Applicative
import Data.List 
import qualified Data.List.Utils as LU

-- Initiation phase of Replication. See README for details. 
-- Remember, underwinding executed for enzyme access to particular base pairs has a linking difference of about -0.05 to -0.07,
-- meaning about 5 - 7 percent of helical turns are removed. 

initiation :: Monad m => DNA -> m DNA
initiation dna = do
	dna' <- return $ adptoATP dna
	dnaPrime2 <- return $ wind (-0.05 * (25)) dna' 
	dnaPrime3 <-  return $ bindHelicase dnaPrime2 
	dnaPrime4 <-  return $ helicase dnaPrime3 
	dnaPrime5 <- return $ topisomeraseII dnaPrime4 
	dnaPrime6 <-  return $ atptoADP dnaPrime5
	dnaPrime7 <-  return $ unBindHelicase dnaPrime6
	return dnaPrime7

-- Elongation phase of Replication. See README for details. 

elongation :: Monad m => DNA -> m [DNA]
elongation dna = do
	dna'  <-  return $ bindPolymerase dna  
	dna'' <-  return $ duplicateDNA dna'
	dnaPrime3 <- return $ map (unBindPolymerase) dna''
	return dnaPrime3

-- Chromosomal crossover. Takes 50 pieces of substructure and exchanges it another chromosome. 

crossoverChromosomes:: [Chromosome]-> [Chromosome]
crossoverChromosomes [chr, chr1] = [chromosome $ replaceMap chr, chromosome $ replaceMap1 chr1] 
					where 
						p x = map (take 50) (getChromosomeDNA x)
						replaceMap x =LU.replace (p chr) (p chr1) (getChromosomeDNA x)
						replaceMap1 y = LU.replace (p chr1) (p chr) (getChromosomeDNA y)

-- Chromosomal crossover mapped over entire chromosome set. crossovers only occur between homologues, 
-- which are determined here through similarity function (called 'residue' here). See README for details.

crossoverCell :: ChromosomeSet -> [ChromosomeSet]
crossoverCell cs = map createChromosomeSet $ map crossoverChromosomes (homologues cs)
				where
					homologues cs = filter (\x -> (length x) == 2) $ groupBy (\x y -> residue x y < 0.50) (chromosomes cs)
					residue x y = 0.2

--Sequential replication of DNA and higher order structures. 

replicateGene :: DNA -> [DNA]
replicateGene dna = [head strands, strands !! 1]
				where strands = concat $ return dna >>= initiation >>= elongation

replicateNucleosomeAssembly :: NucleosomeAssembly -> [NucleosomeAssembly]
replicateNucleosomeAssembly ns = map createNucleosomeAssembly $ replicateGene (dna ns) 

replicateThirtyNMLoop :: ThirtyNMLoop -> [ThirtyNMLoop]
replicateThirtyNMLoop tl =  map createThirtyNMLoop $ [map head k, map (\x -> x !! 1) k]
					where k = map replicateNucleosomeAssembly (nucleosomeAssemblies tl)

replicateRosette :: Rosette -> [Rosette]
replicateRosette rs =  map createRosette $ [map head k, map (\x -> x !! 1) k]
					where k =  map replicateThirtyNMLoop (thirtyNMLoops rs)

replicateCoil :: Coil -> [Coil]
replicateCoil cl = map createCoil $ [map head k, map (\x -> x !! 1) k]
					where k =  map replicateRosette (rosettes cl)

replicateChromatid :: Chromatid -> [Chromatid]
replicateChromatid cd = map createChromatid $ [map head k, map (\x -> x !! 1) k]
					where k = map replicateCoil (coils cd) 

replicateChromosomes :: Chromosome -> [Chromosome]
replicateChromosomes cs = map createChromosome $ [map head k, map (\x -> x !! 1) k]
					where k = map replicateChromatid (chromatids cs) 

replicateCell:: ChromosomeSet -> [ChromosomeSet]
replicateCell ct = map createChromosomeSet $ [map head k, map (\x -> x !! 1) k]
					where k = map replicateChromosomes (chromosomes ct) 