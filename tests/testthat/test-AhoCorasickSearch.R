library(testthat)
library(AhoCorasickTrie)

context("AhoCorasickSearch")

# simple search of multiple keywords in a single text
test_that("simple search of multiple keywords in a single text", {
  keywords = c("Abra", "cadabra", "is", "the", "Magic", "Word")
  oneSearch = AhoCorasickSearch(keywords, "Is Abracadabra the Magic Word?")
  expect_equal(oneSearch[[1]][[1]], list(Keyword="Abra", Offset=4))
  expect_equal(oneSearch[[1]][[2]], list(Keyword="cadabra", Offset=8))
  expect_equal(oneSearch[[1]][[3]], list(Keyword="the", Offset=16))
  expect_equal(oneSearch[[1]][[4]], list(Keyword="Magic", Offset=20))
  expect_equal(oneSearch[[1]][[5]], list(Keyword="Word", Offset=26))
})


# named search of multiple texts in a single list with keyword grouping and aminoacid alphabet
# * all matches to a keyword are accessed by name
# * non-matched keywords are dropped
proteins = c(protein1="PEPTIDEPEPTIDEDADADARARARARAKEKEKEKEPEPTIDE",
             protein2="DERPADERPAPEWPEWPEEPEERAWRAWWARRAGTAGPEPTIDEKESEQUENCE")
peptides = c("PEPTIDE", "DERPA", "SEQUENCE", "KEKE", "PEPPIE")

test_that("named search of multiple texts in a single list with keyword grouping and aminoacid alphabet", {
  peptideSearch = AhoCorasickSearch(peptides, proteins, alphabet="aminoacid", groupByKeyword=T)
  expect_equal(peptideSearch$PEPTIDE, list(list(Text="protein1", Offset=1),
                                           list(Text="protein1", Offset=8),
                                           list(Text="protein1", Offset=37),
                                           list(Text="protein2", Offset=38)))
  expect_equal(peptideSearch$DERPA, list(list(Text="protein2", Offset=1),
                                         list(Text="protein2", Offset=6)))
  expect_equal(peptideSearch$SEQUENCE, list(list(Text="protein2", Offset=47)))
  expect_equal(peptideSearch$KEKE, list(list(Text="protein1", Offset=29),
                                        list(Text="protein1", Offset=31),
                                        list(Text="protein1", Offset=33)))
  expect_equal(peptideSearch$PEPPIE, NULL)
})


# grouping by keyword without text names: offsets are given without reference to the text
test_that("grouping by keyword without text names", {
  names(proteins) = NULL
  peptideSearch = AhoCorasickSearch(peptides, proteins, groupByKeyword=T)
  expect_equal(peptideSearch$PEPTIDE, list(1, 8, 37, 38))
  expect_equal(peptideSearch$DERPA, list(1, 6))
  expect_equal(peptideSearch$SEQUENCE, list(47))
  expect_equal(peptideSearch$KEKE, list(29, 31, 33))
})

