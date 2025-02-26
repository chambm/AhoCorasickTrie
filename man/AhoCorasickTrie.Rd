% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AhoCorasickTrie.R
\docType{package}
\name{AhoCorasickTrie}
\alias{AhoCorasickTrie-package}
\alias{AhoCorasickTrie}
\title{AhoCorasickTrie: fast searching for multiple keywords in multiple texts}
\description{
Builds an Aho-Corasick trie from one or more keywords and uses it to
  search one or more texts. For a large number of keywords, Aho-Corasick is much faster
  than a naive approach (such as \code{lapply(keywords, gregexpr, text)}).
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/chambm/AhoCorasickTrie}
  \item Report bugs at \url{https://github.com/chambm/AhoCorasickTrie/issues}
}

}
\author{
\strong{Maintainer}: Matt Chambers \email{matt.chambers42@gmail.com}

Authors:
\itemize{
  \item Tomas Petricek [copyright holder]
}

Other contributors:
\itemize{
  \item Vanderbilt University [copyright holder]
}

}
