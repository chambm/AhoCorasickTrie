# AhoCorasickTrie 0.1.4

* A keyword containing a character outside the selected alphabet (e.g. `[` or
  a lowercase letter with `alphabet = "aminoacid"`) now gives an R error instead
  of crashing R. A text character one past the end of the alphabet is also now
  caught (#3).

# AhoCorasickTrie 0.1.3

* Fixed CRAN check issues (UBSAN, C++ `SystemRequirements`).
* Removed the Biostrings suggestion.
