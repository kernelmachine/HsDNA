module GenUtils where

import GenStructs

import Control.Applicative
import Prelude hiding (length,null)
import Data.List (map, nub,filter,length)
import Data.ByteString (ByteString,null,isPrefixOf,pack,length)
import Data.ByteString.Internal (c2w)
import Data.ByteString.Unsafe 


complementarySequence :: [Nucleotide] -> [Nucleotide]
complementarySequence [] = []
complementarySequence (x:xs)= figure x : complementarySequence xs
				where 
					figure x 
						| x == 'A' = 'T'
						| x == 'T' = 'A'
						| x == 'C' = 'G'
						| x == 'G' = 'C'
						| otherwise = error (charToString x  ++ " is not a nucleotide")

charToString :: Char -> String
charToString c = [c]

check :: [Nucleotide] -> [Nucleotide] -> Bool
check fs ss = Data.List.length fs == Data.List.length ss  

findNucleotidePatterns :: ByteString -> ByteString -> [Maybe Int] -- Used to find positions of particular nucleotide patterns in DNA sequence.   
findNucleotidePatterns pat str
    | null pat         = map wrapMaybe [ 0 .. Data.ByteString.length str]
    | otherwise        = search 0 str
  where
    search n s
        | null s             = []
        | pat `isPrefixOf` s = Just n : search (n+1) (unsafeTail s)
        | otherwise          = Nothing :search (n+1) (unsafeTail s)
    wrapMaybe a = Just a

-- Nucleosomes tend to bind to sequences of repeated "A=T" bindings. 'findATATAT' uses 'findNucleotidePatterns' to return positions of such sequences,
-- filtering out 'Nothing' and dividing by 2 because the inputted string contains nucleotides from both strands.

findATATAT :: [Nucleotide] -> [Maybe Int]
findATATAT dna = map (fmap ( `div` 2)) $ nub $ filter (\x -> x /= Nothing) $ findNucleotidePatterns (at) (pack $ fmap c2w dna) 
				++ findNucleotidePatterns (ta) (pack $ fmap c2w dna)
		where 
			at = pack (fmap c2w "ATATAT")
			ta = pack (fmap c2w "TATATA")

compressPairs :: [(Nucleotide, Nucleotide)] -> [Nucleotide] -- changes the representation of DNA from a list of tuples to a string.
compressPairs [] = []
compressPairs (x:xs) = ((charToString $ fst x) ++ (charToString $ snd x)) ++ compressPairs xs