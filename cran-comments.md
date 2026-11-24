## Submission summary

This is a bug-fix release. A keyword containing a character outside the
selected alphabet used to index past the end of an array while the trie was
built, which crashed R. It now gives an R error.

## R CMD check results

0 errors | 0 warnings | 0 notes

Checked with `R CMD check --as-cran` on GitHub Actions:

* ubuntu-latest: R-devel, R-release, R-oldrel-1
* macos-latest: R-release
* windows-latest: R-release
* rocker/r-devel-san (R-devel with ASAN and UBSAN): no sanitizer errors

## Reverse dependencies

There is 1 reverse dependency (prozor). Its `R CMD check` against this
version shows no new problems.
