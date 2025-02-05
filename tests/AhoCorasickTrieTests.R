library(AhoCorasickTrie)
listEquals = function(a, b) { is.null(unlist(a)) && is.null(unlist(b)) || !is.null(a) && !is.null(b) && all(unlist(a) == unlist(b)) }

# simple search of multiple keywords in a single text
keywords = c("Abra", "cadabra", "is", "the", "Magic", "Word")
oneSearch = AhoCorasickSearch(keywords, "Is Abracadabra the Magic Word?")
stopifnot(listEquals(oneSearch[[1]][[1]], list(keyword="Abra", offset=4)))
stopifnot(listEquals(oneSearch[[1]][[2]], list(keyword="cadabra", offset=8)))
stopifnot(listEquals(oneSearch[[1]][[3]], list(keyword="the", offset=16)))
stopifnot(listEquals(oneSearch[[1]][[4]], list(keyword="Magic", offset=20)))
stopifnot(listEquals(oneSearch[[1]][[5]], list(keyword="Word", offset=26)))

# search a list of lists
# * sublists are accessed by index
# * texts are accessed by index
# * non-matched texts are kept (to preserve index order)
listSearch = AhoCorasickSearchList(keywords, list(c("What in", "the world"), c("is"), "secret about", "the Magic Word?"))
stopifnot(listEquals(listSearch[[1]][[1]], list()))
stopifnot(listEquals(listSearch[[1]][[2]][[1]], list(keyword="the", offset=1)))
stopifnot(listEquals(listSearch[[2]][[1]][[1]], list(keyword="is", offset=1)))
stopifnot(listEquals(listSearch[[3]], list()))
stopifnot(listEquals(listSearch[[4]][[1]][[1]], list(keyword="the", offset=1)))
stopifnot(listEquals(listSearch[[4]][[1]][[2]], list(keyword="Magic", offset=5)))
stopifnot(listEquals(listSearch[[4]][[1]][[3]], list(keyword="Word", offset=11)))

# named search of a list of lists
# * sublists are accessed by name
# * matched texts are accessed by name
# * non-matched texts are dropped
namedSearch = AhoCorasickSearchList(keywords, list(subject=c(phrase1="What in", phrase2="the world"),
                                                   verb=c(phrase1="is"),
                                                   predicate1=c(phrase1="secret about"),
                                                   predicate2=c(phrase1="the Magic Word?")))
stopifnot(listEquals(namedSearch$subject$phrase2[[1]], list(keyword="the", offset=1)))
stopifnot(listEquals(namedSearch$verb$phrase1[[1]], list(keyword="is", offset=1)))
stopifnot(listEquals(namedSearch$predicate1, list()))
stopifnot(listEquals(namedSearch$predicate2$phrase1[[1]], list(keyword="the", offset=1)))
stopifnot(listEquals(namedSearch$predicate2$phrase1[[2]], list(keyword="Magic", offset=5)))
stopifnot(listEquals(namedSearch$predicate2$phrase1[[3]], list(keyword="Word", offset=11)))

# named search of multiple texts in a single list with keyword grouping and aminoacid alphabet
# * all matches to a keyword are accessed by name
# * non-matched keywords are dropped
proteins = c(protein1="PEPTIDEPEPTIDEDADADARARARARAKEKEKEKEPEPTIDE",
             protein2="DERPADERPAPEWPEWPEEPEERAWRAWWARRAGTAGPEPTIDEKESEQUENCE")
peptides = c("PEPTIDE", "DERPA", "SEQUENCE", "KEKE", "PEPPIE")
peptideSearch = AhoCorasickSearch(peptides, proteins, alphabet="aminoacid", groupByKeyword=T)
stopifnot(listEquals(peptideSearch$PEPTIDE, list(list(keyword="protein1", offset=1),
                                                 list(keyword="protein1", offset=8),
                                                 list(keyword="protein1", offset=37),
                                                 list(keyword="protein2", offset=38))))
stopifnot(listEquals(peptideSearch$DERPA, list(list(keyword="protein2", offset=1),
                                               list(keyword="protein2", offset=6))))
stopifnot(listEquals(peptideSearch$SEQUENCE, list(list(keyword="protein2", offset=47))))
stopifnot(listEquals(peptideSearch$KEKE, list(list(keyword="protein1", offset=29),
                                              list(keyword="protein1", offset=31),
                                              list(keyword="protein1", offset=33))))
stopifnot(listEquals(peptideSearch$PEPPIE, NULL))

# grouping by keyword without text names: offsets are given without reference to the text
names(proteins) = NULL
peptideSearch = AhoCorasickSearch(peptides, proteins, groupByKeyword=T)
stopifnot(listEquals(peptideSearch$PEPTIDE, list(1, 8, 37, 38)))
stopifnot(listEquals(peptideSearch$DERPA, list(1, 6)))
stopifnot(listEquals(peptideSearch$SEQUENCE, list(47)))
stopifnot(listEquals(peptideSearch$KEKE, list(29, 31, 33)))

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

  if (0) #suppressPackageStartupMessages(require(Biostrings)))
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

    #dict = DNAStringSet(bigProbes)
    #dna = DNAString(paste(bigDNA, collapse="."))
    #pdict = PDict(dict)
    print(microbenchmark(#createDNAStringSet = DNAStringSet(bigProbes),
                         #createDNAString = DNAString(paste(bigDNA, collapse=".")),
                         #createPDict = PDict(dict),
                         #matchPDict = matchPDict(pdict, dna),
                         AhoCorasickASCII = AhoCorasickSearch(bigProbes, bigDNA),
                         AhoCorasickNucleicAcid = AhoCorasickSearch(bigProbes, bigDNA, alphabet="nucleicacid"),
                         times=3))
  }
}
