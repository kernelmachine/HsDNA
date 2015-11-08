module GenStructs where

type Nucleotide = Char

data DNA = DNA {
	strands ::  [(Nucleotide, Nucleotide)],
	separatedStrands :: ([Nucleotide],[Nucleotide]), -- space in struct for separated strands
	linkingNumber :: Maybe Double,
	basePairsPerTurn :: Maybe Double,
	helicaseBound :: Bool,
	dnaATP :: Bool,  -- is ATP bound to DNA-A protein?
	dnaADP :: Bool,  -- is ADP bound to DNA-A protein?
	replicationOrigin :: (Maybe Int, Maybe Int), -- position of replication origin
	polymeraseBounded :: Bool
	}
	deriving (Show, Eq)

data Nucleosome = Nucleosome {
	startPos :: Maybe Int, -- start position of nucleosome on dna seq
	endPos :: Maybe Int   -- end position of nucleosome on dna seq
}
	deriving (Show, Eq)

data NucleosomeAssembly = NucleosomeAssembly { 
			dna :: DNA, 
			nucleosomes :: [Nucleosome]
			}
			deriving (Show, Eq)

newtype ThirtyNMLoop = ThirtyNMLoop {nucleosomeAssemblies :: [NucleosomeAssembly]} deriving (Show, Eq)

newtype Rosette = Rosette {thirtyNMLoops :: [ThirtyNMLoop]} deriving (Show, Eq)

newtype Coil = Coil {rosettes :: [Rosette]} deriving (Show, Eq)

newtype Chromatid = Chromatid {coils :: [Coil]} deriving (Show, Eq)

newtype Chromosome = Chromosome {chromatids :: [Chromatid]} deriving (Show, Eq)

newtype ChromosomeSet = ChromosomeSet {chromosomes :: [Chromosome]} deriving (Show, Eq)