library(AhoCorasickTrie)

if (suppressPackageStartupMessages(require(microbenchmark)))
{
  set.seed(0)
  # generate some random protein sequences
  bigProteins = replicate(10, paste(sample(strsplit("ACDEFGHIKLMNPQRSTUVWY", "")[[1]], 100, replace=T), collapse=""))

  set.seed(0)
  # generate random tryptic peptdides from the above proteins
  bigPeptides = c(replicate(100, sapply(bigProteins,
                                        function(protein)
                                        {
                                          sites=sort(sample(c(0, gregexpr("[KR]", protein, perl=T)[[1]], nchar(protein)), 2, replace=F))
                                          if (sites[2]-sites[1] < 30 & sites[2]-sites[1] > 5) substr(protein, sites[1]+1, sites[2])
                                        })),
                  recursive=T)

  print(microbenchmark(gregexpr = lapply(bigPeptides, gregexpr, bigProteins),
                       AhoCorasickASCII = AhoCorasickSearch(bigPeptides, bigProteins),
                       AhoCorasickAminoAcid = AhoCorasickSearch(bigPeptides, bigProteins, alphabet="aminoacid"),
                       AhoCorasickGroupByKeyword = AhoCorasickSearch(bigPeptides, bigProteins, groupByKeyword=T),
                       times=3))

  if (suppressPackageStartupMessages(require(Biostrings)))
  {
    set.seed(0)
    # generate some random DNA sequences
    bigDNA = replicate(10, paste(sample(strsplit("AGTC", "")[[1]], 100, replace=T), collapse=""))

    set.seed(0)
    #generate random 10mers from the above DNA sequences
    bigProbes = c(replicate(100, sapply(bigDNA,
                                        function(dna)
                                        {
                                          site=sort(sample(seq(1, nchar(dna)-10, 10), 1, replace=F))
                                          substr(dna, site, site+10)
                                        })),
                  recursive=T)

    dict = DNAStringSet(bigProbes)
    dna = DNAString(paste(bigDNA, collapse="."))
    pdict = PDict(dict)
    print(microbenchmark(createDNAStringSet = DNAStringSet(bigProbes),
                         createDNAString = DNAString(paste(bigDNA, collapse=".")),
                         createPDict = PDict(dict),
                         matchPDict = matchPDict(pdict, dna),
                         AhoCorasickASCII = AhoCorasickSearch(bigProbes, bigDNA),
                         AhoCorasickNucleicAcid = AhoCorasickSearch(bigProbes, bigDNA, alphabet="nucleicacid"),
                         times=3))
  }
}
