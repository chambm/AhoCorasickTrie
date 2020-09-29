library(testthat)
library(AhoCorasickTrie)

context("AhoCorasickSearchList")

keywords = c("Abra", "cadabra", "is", "the", "Magic", "Word")

# search a list of lists
# * sublists are accessed by index
# * texts are accessed by index
# * non-matched texts are kept (to preserve index order)
test_that("search a list of lists", {
  listSearch = AhoCorasickSearchList(keywords, list(c("What in", "the world"), c("is"), "secret about", "the Magic Word?"))
  expect_equal(listSearch[[1]][[1]], list())
  expect_equal(listSearch[[1]][[2]][[1]], list(Keyword="the", Offset=1))
  expect_equal(listSearch[[2]][[1]][[1]], list(Keyword="is", Offset=1))
  expect_equal(listSearch[[3]][[1]], list())
  expect_equal(listSearch[[4]][[1]][[1]], list(Keyword="the", Offset=1))
  expect_equal(listSearch[[4]][[1]][[2]], list(Keyword="Magic", Offset=5))
  expect_equal(listSearch[[4]][[1]][[3]], list(Keyword="Word", Offset=11))
})



# named search of a list of lists
# * sublists are accessed by name
# * matched texts are accessed by name
# * non-matched texts are dropped
# namedSearch = AhoCorasickSearchList(keywords, list(subject=c(phrase1="What in", phrase2="the world"),
#                                                    verb=c(phrase1="is"),
#                                                    predicate1=c(phrase1="secret about"),
#                                                    predicate2=c(phrase1="the Magic Word?")))
# stopifnot(listEquals(namedSearch$subject$phrase2[[1]], list(keyword="the", offset=1)))
# stopifnot(listEquals(namedSearch$verb$phrase1[[1]], list(keyword="is", offset=1)))
# stopifnot(listEquals(namedSearch$predicate1, list()))
# stopifnot(listEquals(namedSearch$predicate2$phrase1[[1]], list(keyword="the", offset=1)))
# stopifnot(listEquals(namedSearch$predicate2$phrase1[[2]], list(keyword="Magic", offset=5)))
# stopifnot(listEquals(namedSearch$predicate2$phrase1[[3]], list(keyword="Word", offset=11)))
test_that("named search of a list of lists", {
  namedSearch = AhoCorasickSearchList(keywords, list(subject=c(phrase1="What in", phrase2="the world"),
                                                     verb=c(phrase1="is"),
                                                     predicate1=c(phrase1="secret about"),
                                                     predicate2=c(phrase1="the Magic Word?")))
  expect_equal(namedSearch$subject$phrase2[[1]], list(Keyword="the", Offset=1))
  expect_equal(namedSearch$verb$phrase1[[1]], list(Keyword="is", Offset=1))
  expect_equal(unlist(namedSearch$predicate1), NULL)
  expect_equal(namedSearch$predicate2$phrase1[[1]], list(Keyword="the", Offset=1))
  expect_equal(namedSearch$predicate2$phrase1[[2]], list(Keyword="Magic", Offset=5))
  expect_equal(namedSearch$predicate2$phrase1[[3]], list(Keyword="Word", Offset=11))
})
