#' AhoCorasickTrie: fast searching for multiple keywords in multiple texts
#'
#' @docType package
#' @name AhoCorasickTrie
#' @importFrom Rcpp evalCpp
#' @useDynLib AhoCorasickTrie
#' @description Builds an Aho-Corasick trie from one or more keywords and uses it to
#'   search one or more texts. For a large number of keywords, Aho-Corasick is much faster
#'   than a naive approach (such as \code{lapply(keywords, gregexpr, text)}).
NULL
