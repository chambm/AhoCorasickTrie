#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

//
// Original author: Tomas Petricek (http://www.codeproject.com/KB/recipes/ahocorasick.aspx)
// Copyright 2005 Tomas Petricek
//
// Adapted to C++ by: Matt Chambers <matt.chambers .@. vanderbilt.edu>
// Copyright 2016 Vanderbilt University
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include <string>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <exception>
#include <algorithm>

namespace freicore {

struct ascii_translator
{
  static int size() { return 128; }
  static char translate(int index) { return static_cast<char>(index); }
  static int translate(char symbol) { return static_cast<int>(symbol); }
};

struct AminoAcidTranslator
{
  static int size() { return 26; }
  static int translate(char aa) { return aa - 'A'; };
  static char translate(int index) { return static_cast<char>(index)+'A'; }
};

struct NucleicAcidTranslator
{
  static constexpr const char* ALPHABET = "ACGTUMRWSYKVHDBN-+.";
  static int size() { return 19; }
  static char translate(int index) { return ALPHABET[index]; }

  static int translate(char na)
  {
    switch(na)
    {
      default:
      case 'A': return 0;
      case 'C': return 1;
      case 'G': return 2;
      case 'T': return 3;
      case 'U': return 4;
      case 'M': return 5;
      case 'R': return 6;
      case 'W': return 7;
      case 'S': return 8;
      case 'Y': return 9;
      case 'K': return 10;
      case 'V': return 11;
      case 'H': return 12;
      case 'D': return 13;
      case 'B': return 14;
      case 'N': return 15;
      case '-': return 16;
      case '+': return 17;
      case '.': return 18;
    }
  }
};


class AhoCorasickTrie
{
  public:
  typedef std::shared_ptr<std::string> shared_keytype;

  struct SearchResult
  {
    size_t offset() const { return _offset; }
    const shared_keytype& keyword() const { return _keyword; }

    SearchResult(size_t offset, const shared_keytype& keyword)
      : _offset(offset), _keyword(keyword)
    {}

    protected:
    size_t _offset;
    shared_keytype _keyword;
    friend class AhoCorasickTrie;
  };

  virtual ~AhoCorasickTrie() {};

  /// inserts a single shared_string and rebuilds the trie
  virtual void insert(const shared_keytype& keyword) = 0;

  /// returns the first instance of a keyword in the text
  virtual SearchResult find_first(const std::string& text) = 0;

  /// returns all instances of all keywords in the text
  virtual std::vector<SearchResult> find_all(const std::string& text) = 0;

  virtual size_t size() const = 0;

  virtual bool empty() const = 0;

  virtual void clear() = 0;
};


template <typename SymbolTranslator = ascii_translator>
class AhoCorasickTrieImpl : public AhoCorasickTrie
{
public:

  typedef std::shared_ptr<std::string> shared_keytype;

  /// default constructor
  AhoCorasickTrieImpl() : _emptyHash(0), _root(0) {}

  template <typename FwdIterator>
  AhoCorasickTrieImpl(FwdIterator begin, FwdIterator end) : _emptyHash(0), _root(0)
  {
    insert(begin, end);
  }

  virtual ~AhoCorasickTrieImpl()
  {
    if (_root)
      delete _root;
    if (_emptyHash)
      delete[] _emptyHash;
  }

  /// inserts a range of shared_string and rebuilds the trie
  template <typename FwdIterator>
  void insert(FwdIterator begin, FwdIterator end)
  {
    for (; begin != end; ++begin) _insert(*begin);
    _build();
  }

  /// inserts a single shared_string and rebuilds the trie
  void insert(const shared_keytype& keyword)
  {
    _insert(keyword);
    _build();
  }

  /// returns the first instance of a keyword in the text
  SearchResult find_first(const std::string& text)
  {
    if (!_root)
      return SearchResult(text.length(), shared_keytype());

    typename Node::Ptr ptr = _root;
    size_t offset = 0;

    while (offset < text.length())
    {
      typename Node::Ptr transition = 0;
      while (!transition)
      {
        int index = SymbolTranslator::translate(text[offset]);
        if (index < 0 || index > SymbolTranslator::size())
          throw std::out_of_range(std::string("[AhoCorasickTrie::find_first] character '") + text[offset] + "' is not in the trie's alphabet");

        transition = ptr->transitionHash[index];
        if (ptr == _root) break;
        if (!transition) ptr = ptr->failure;
      }
      if (transition) ptr = transition;

      if (!ptr->results->empty())
      {
        const shared_keytype& result = *ptr->results->begin();
        return SearchResult(offset - result->length() + 1, result);
      }
      ++offset;
    }
    return SearchResult(text.length(), shared_keytype());
  }

  /// returns all instances of all keywords in the text
  std::vector<SearchResult> find_all(const std::string& text)
  {
    std::vector<SearchResult> results;

    if (!_root)
      return results;

    typename Node::Ptr ptr = _root;
    size_t offset = 0;

    size_t textLength = text.length();
    while (offset < textLength)
    {
      typename Node::Ptr transition = 0;
      while (!transition)
      {
        int index = SymbolTranslator::translate(text[offset]);
        if (index < 0 || index > SymbolTranslator::size())
          throw std::out_of_range(std::string("[AhoCorasickTrie::find_all] character '") + text[offset] + "' is not in the trie's alphabet (" + text.substr(std::max(int(offset), 5) - 5, 10) + ")");

        transition = ptr->transitionHash[index];
        if (ptr == _root) break;
        if (!transition) ptr = ptr->failure;
      }
      if (transition) ptr = transition;

      auto resultSet = *ptr->results;
      for (auto& result : resultSet)
        results.push_back(SearchResult(offset - static_cast<const std::string&>(*result).length() + 1, result));
      ++offset;
    }
    return results;
  }

  size_t size() const
  {
    return _keywords.size();
  }

  bool empty() const
  {
    return size() == 0;
  }

  void clear()
  {
    _keywords.clear();

    if (_root) { delete _root; _root = 0; }
    if (_emptyHash) { delete[] _emptyHash; _emptyHash = 0; }
  }

private:

  struct SharedKeyTypeFastLessThan
  {
    bool operator() (const shared_keytype& lhs, const shared_keytype& rhs) const
    {
      const std::string& lhsStr = static_cast<const std::string&>(*lhs);
      const std::string& rhsStr = static_cast<const std::string&>(*rhs);
      if (lhsStr.length() == rhsStr.length())
        return lhsStr < rhsStr;
      return lhsStr.length() < rhsStr.length();
    }
  };

  typedef std::set<shared_keytype, SharedKeyTypeFastLessThan> SharedKeyTypeSet;

  struct Node
  {
    typedef Node* Ptr;

    Node(SharedKeyTypeSet& emptyResults, Ptr*& emptyHash, Ptr parent, char c)
      : value(c),
        parent(parent),
        failure(0),
        results(&emptyResults),
        transitionHash(emptyHash),
        emptyHash(emptyHash)
    {}

    ~Node()
    {
      if (!results->empty())
        delete results;
      if (transitionHash != emptyHash)
      {
        for (int i = 0; i < SymbolTranslator::size(); ++i)
          if (transitionHash[i])
            delete transitionHash[i];
          delete[] transitionHash;
      }
    }

    void addResult(const shared_keytype& result)
    {
      if (results->empty())
        results = new SharedKeyTypeSet;
      results->insert(result);
    }

    void addTransition(char c, Ptr& node)
    {
      if (transitionHash == emptyHash)
      {
        transitionHash = new Ptr[SymbolTranslator::size()];
        std::fill(transitionHash, transitionHash + SymbolTranslator::size(), Ptr());
      }
      transitionHash[SymbolTranslator::translate(c)] = node;
    }

    char value;
    Ptr parent;
    Ptr failure;
    SharedKeyTypeSet* results;
    Ptr* transitionHash;
    Ptr*& emptyHash;
  };

  void _insert(const shared_keytype& keyword)
  {
    _isDirty = true;
    _keywords.insert(keyword);
  }

  SharedKeyTypeSet _emptyResults;
  typename Node::Ptr* _emptyHash;

  void _build()
  {
    if (!_isDirty)
      return;

    if (!_emptyHash)
    {
      _emptyHash = new typename Node::Ptr[SymbolTranslator::size()];
      std::fill(_emptyHash, _emptyHash + SymbolTranslator::size(), typename Node::Ptr());
    }

    // Build keyword tree and transition function
    _root = new Node(_emptyResults, _emptyHash, typename Node::Ptr(), SymbolTranslator::translate(0));
    for (auto& keyword : _keywords)
    {
      // add pattern to tree
      typename Node::Ptr node = _root;
      for (auto& c : static_cast<const std::string&>(*keyword))
      {
        typename Node::Ptr newNode = node->transitionHash[SymbolTranslator::translate(c)];

        if (!newNode)
        {
          newNode = new Node(_emptyResults, _emptyHash, node, c);
          node->addTransition(c, newNode);
        }
        node = newNode;
      }
      node->addResult(keyword);
    }

    // Find failure functions
    std::vector<typename Node::Ptr> nodes;

    // level 1 nodes - fail to root node
    for (int i = 0; i < SymbolTranslator::size(); ++i)
    {
      const typename Node::Ptr& depth1Node = _root->transitionHash[i];
      if (!depth1Node)
        continue;

      depth1Node->failure = _root;
      for (int j = 0; j < SymbolTranslator::size(); ++j)
      {
        const typename Node::Ptr& depth2Node = depth1Node->transitionHash[j];
        if (depth2Node)
          nodes.push_back(depth2Node);
      }
    }

    // other nodes - using BFS
    while (!nodes.empty())
    {
      std::vector<typename Node::Ptr> nextLevelNodes;
      for (auto& node : nodes)
      {
        typename Node::Ptr r = node->parent->failure;
        char c = node->value;

        while (r && !r->transitionHash[SymbolTranslator::translate(c)]) r = r->failure;
        if (!r)
          node->failure = _root;
        else
        {
          node->failure = r->transitionHash[SymbolTranslator::translate(c)];
          for (auto& result : *node->failure->results)
            node->addResult(result);
        }

        // add child nodes to BFS list
        for (int i = 0; i < SymbolTranslator::size(); ++i)
        {
          const typename Node::Ptr& transition = node->transitionHash[i];
          if (transition)
            nextLevelNodes.push_back(transition);
        }
      }
      nodes = nextLevelNodes;
    }
    _root->failure = _root;

    _isDirty = false;
  }

  bool _isDirty;
  SharedKeyTypeSet _keywords;
  typename Node::Ptr _root;
};

enum class Alphabet
{
  ASCII, AminoAcid, NucleicAcid
};

template <typename FwdIterator>
std::unique_ptr<AhoCorasickTrie> CreateAhoCorasickTrie(FwdIterator begin, FwdIterator end, Alphabet alphabet = freicore::Alphabet::ASCII)
{
  switch(alphabet)
  {
    default:
    case Alphabet::ASCII: return std::unique_ptr<AhoCorasickTrie>(new AhoCorasickTrieImpl<>(begin, end));
    case Alphabet::AminoAcid: return std::unique_ptr<AhoCorasickTrie>(new AhoCorasickTrieImpl<AminoAcidTranslator>(begin, end));
    case Alphabet::NucleicAcid: return std::unique_ptr<AhoCorasickTrie>(new AhoCorasickTrieImpl<NucleicAcidTranslator>(begin, end));
  }
  return std::unique_ptr<AhoCorasickTrie>();
}

} // namespace freicore


//' Fast searching for one or more keywords in a list of texts
//'
//' @param keywords Character vector of one or more keywords
//' @param textList List of lists, each sublist with one or more texts to search
//' @param alphabet Alphabet to use; one of \code{ascii}, \code{aminoacid}, or \code{nucleicacid}
//' @param groupByKeyword If true, matches are grouped by keyword (instead of by text)
//' @param iterationFeedback When set to a positive integer \code{i}, console output will indicate when searching every \code{i}th text
//' @return List of lists of matches, grouped by either text or by keyword (each list of texts gets its own list of matches)
//' @description Builds an Aho-Corasick trie from one or more keywords and uses it to search a list of
//'   one or more texts. For a large number of keywords, Aho-Corasick is much faster
//'   than a naive approach (such as \code{lapply(keywords, gregexpr, text)}).
//'
//'   Use \code{\link{AhoCorasickSearchList}} instead of \code{\link{AhoCorasickSearch}} when you want to keep the matches
//'   of each input sublist separate. If the sublists of the input list have names, the resulting list of lists
//'   will use those names, but sublists with no matches will still be in the resulting list.
//'   If the texts of the sublists have names, the resulting sublists of matches will use
//'   those names, and the texts with no matches will be dropped. If the input texts do
//'   not have names, then the resulting sublists of matches will be in the same order as the
//'   input texts, and non-matched texts will be kept to preserve that order. Thus, it is more
//'   efficient to use named input texts (so non-matched texts can be dropped).
//'
//'   The default alphabet allows all 128 ASCII characters in the keywords and the texts.
//'   Characters outside this range will cause an error. A more efficient trie is possible
//'   if the alphabet size can be reduced. For example, DNA sequences use at most 19 distinct
//'   characters and usually only 4; protein sequences use at most 26 distinct characters and
//'   usually only 20. Set the \code{alphabet} parameter if a reduced alphabet is appropriate.
//'
//'   UTF-8 (Unicode) matching is not currently supported.
//' @seealso
//' \itemize{
//' \item \href{https://www.codeproject.com/Articles/12383/Aho-Corasick-string-matching-in-C}{Aho-Corasick string matching in C#} for the article this package is based on
//' \item \code{\link[Biostrings]{matchPDict}} for a more memory efficient, but DNA-only, implementation of the algorithm
//' }
//' @examples
//' listEquals = function(a, b) { is.null(unlist(a)) && is.null(unlist(b)) ||
//'                               !is.null(a) && !is.null(b) && all(unlist(a) == unlist(b)) }
//' keywords = c("Abra", "cadabra", "is", "the", "Magic", "Word")
//'
//' # 1. Search a list of lists without names
//' # * sublists are accessed by index
//' # * texts are accessed by index
//' # * non-matched texts are kept (input index order is preserved)
//' listSearch = AhoCorasickSearchList(keywords,
//'                                    list(c("What in", "the world"),
//'                                         c("is"),
//'                                         "secret about",
//'                                         "the Magic Word?"))
//' stopifnot(listEquals(listSearch[[1]][[1]], list()))
//' stopifnot(listEquals(listSearch[[1]][[2]][[1]], list(keyword="the", offset=1)))
//' stopifnot(listEquals(listSearch[[2]][[1]][[1]], list(keyword="is", offset=1)))
//' stopifnot(listEquals(listSearch[[3]], list()))
//' stopifnot(listEquals(listSearch[[4]][[1]][[1]], list(keyword="the", offset=1)))
//' stopifnot(listEquals(listSearch[[4]][[1]][[2]], list(keyword="Magic", offset=5)))
//' stopifnot(listEquals(listSearch[[4]][[1]][[3]], list(keyword="Word", offset=11)))
//'
//' # 2. Search a named list of named lists
//' # * sublists are accessed by name
//' # * matched texts are accessed by name
//' # * non-matched texts are dropped
//' namedSearch = AhoCorasickSearchList(keywords,
//'                                     list(subject=c(phrase1="What in", phrase2="the world"),
//'                                          verb=c(phrase1="is"),
//'                                          predicate1=c(phrase1="secret about"),
//'                                          predicate2=c(phrase1="the Magic Word?")))
//' stopifnot(listEquals(namedSearch$subject$phrase2[[1]], list(keyword="the", offset=1)))
//' stopifnot(listEquals(namedSearch$verb$phrase1[[1]], list(keyword="is", offset=1)))
//' stopifnot(listEquals(namedSearch$predicate1, list()))
//' stopifnot(listEquals(namedSearch$predicate2$phrase1[[1]], list(keyword="the", offset=1)))
//' stopifnot(listEquals(namedSearch$predicate2$phrase1[[2]], list(keyword="Magic", offset=5)))
//' stopifnot(listEquals(namedSearch$predicate2$phrase1[[3]], list(keyword="Word", offset=11)))
//' @export
// [[Rcpp::export]]
Rcpp::List AhoCorasickSearchList(Rcpp::StringVector keywords,
                                 Rcpp::List textList,
                                 std::string alphabet = "ascii",
                                 bool groupByKeyword = false,
                                 int iterationFeedback = 0)
{
  using freicore::AhoCorasickTrie;
  using freicore::CreateAhoCorasickTrie;
  using freicore::Alphabet;

  Alphabet alphabetEnum = alphabet == "ascii" ? Alphabet::ASCII :
                          alphabet == "aminoacid" ? Alphabet::AminoAcid :
                          alphabet == "nucleicacid" ? Alphabet::NucleicAcid :
    throw std::runtime_error("alphabet must be one of {ascii, aminoacid, nucleicacid}");

  std::vector<std::shared_ptr<std::string> > keywordsStd(keywords.size());
  for (int i=0; i < keywords.size(); ++i)
    keywordsStd[i] = std::make_shared<std::string>(keywords(i));

  std::unique_ptr<AhoCorasickTrie> trie;
  trie = CreateAhoCorasickTrie(keywordsStd.begin(), keywordsStd.end(), alphabetEnum);
  Rcpp::List resultPerList;
  auto listNames = textList.names();

  for (int i=0, textListSize=textList.size(); i < textListSize; ++i)
  {
    if (groupByKeyword)
    {
      Rcpp::List resultPerKeyword;
      if (Rf_isNull(textList(i)))
      {
        resultPerList.push_back(resultPerKeyword);
        continue;
      }

      std::map<std::shared_ptr<std::string>, Rcpp::List> groupByKeywordMap;

      Rcpp::StringVector texts = Rcpp::as<Rcpp::StringVector>(textList(i));
      bool hasNames = texts.hasAttribute("names");
      Rcpp::StringVector textNames;
      if (hasNames)
        textNames = texts.names();
      Rcpp::StringVector resultNames;

      for (int k=0, textsSize=texts.size(); k < textsSize; ++k)
      {
        if (iterationFeedback > 0 && (k==0 || k+1==textsSize || ((k+1) % iterationFeedback)==0))
        {
          if (textListSize > 1)
            Rcpp::Rcout << "Searching TextList " << (i+1) << "/" << textListSize << ", Text " << (k+1) << "/" << textsSize << std::endl;
          else
            Rcpp::Rcout << "Searching Text " << (k+1) << "/" << textsSize << std::endl;
          R_FlushConsole();
          R_ProcessEvents();
          R_CheckUserInterrupt();
        }

        std::vector<AhoCorasickTrie::SearchResult> results = trie->find_all(std::string(texts(k)));

        if (results.empty())
          continue;

        for (size_t j=0; j < results.size(); ++j)
        {
          Rcpp::List& keywordResult = groupByKeywordMap[results[j].keyword()];
          if (hasNames)
            keywordResult.push_back(Rcpp::List::create(Rcpp::_["Text"] = Rcpp::String(textNames(k)),
                                                       Rcpp::_["Offset"] = results[j].offset()+1)); // one-based offset
          else
            keywordResult.push_back(results[j].offset()+1); // one-based offset
        }
      }

      for (auto keywordPair : groupByKeywordMap)
      {
        resultPerKeyword.push_back(keywordPair.second);
        resultNames.push_back(*keywordPair.first);
      }

      resultPerKeyword.names() = resultNames;
      resultPerList.push_back(resultPerKeyword);
    }
    else
    {
      Rcpp::List resultPerText;
      if (Rf_isNull(textList(i)))
      {
        resultPerList.push_back(resultPerText);
        continue;
      }

      Rcpp::StringVector texts = Rcpp::as<Rcpp::StringVector>(textList(i));
      bool hasNames = texts.hasAttribute("names");
      Rcpp::StringVector textNames;
      if (hasNames)
        textNames = texts.names();
      Rcpp::StringVector resultNames;

      for (int k=0, textsSize=texts.size(); k < textsSize; ++k)
      {
        if (iterationFeedback > 0 && (k==0 || k+1==textsSize || ((k+1) % iterationFeedback)==0))
        {
          if (textListSize > 1)
            Rcpp::Rcout << "Searching TextList " << (i+1) << "/" << textListSize << ", Text " << (k+1) << "/" << textsSize << std::endl;
          else
            Rcpp::Rcout << "Searching Text " << (k+1) << "/" << textsSize << std::endl;
          R_FlushConsole();
          R_ProcessEvents();
          R_CheckUserInterrupt();
        }

        std::vector<AhoCorasickTrie::SearchResult> results = trie->find_all(std::string(texts(k)));

        if (hasNames && results.empty())
          continue;

        Rcpp::List textResult(results.size());
        for (size_t j=0; j < results.size(); ++j)
          textResult[j] = Rcpp::List::create(Rcpp::_["Keyword"] = *results[j].keyword(),
                                             Rcpp::_["Offset"] = results[j].offset()+1); // one-based offset

        resultPerText.push_back(textResult);
        if (hasNames)
          resultNames.push_back(textNames(k));
      }

      if (hasNames)
        resultPerText.names() = resultNames;
      resultPerList.push_back(resultPerText);
    }
  }
  resultPerList.names() = listNames;
  return resultPerList;
}

//' Fast searching for one or more keywords in one or more texts
//'
//' @param text Character vector of one or more texts to search
//' @inheritParams AhoCorasickSearchList
//' @return List of matches, grouped by either text or by keyword
//' @description Builds an Aho-Corasick trie from one or more keywords and uses it to
//'   search one or more texts. For a large number of keywords, Aho-Corasick is much faster
//'   than a naive approach (such as \code{lapply(keywords, gregexpr, text)}).
//'
//'   Use \code{\link{AhoCorasickSearchList}} instead of \code{\link{AhoCorasickSearch}} when you want to keep the matches
//'   of each input text separate. If the input texts have names, the resulting list of matches will include those
//'   names and non-matched texts will be excluded from the results. If the input texts do
//'   not have names, then the resulting list of matches will be in the same order as the
//'   input texts, and non-matched texts will be kept to preserve that order. Thus, it is more
//'   efficient to use named input texts (so non-matched texts can be dropped).
//'
//'   The default alphabet allows all 128 ASCII characters in the keywords and the texts.
//'   Characters outside this range will cause an error. A more efficient trie is possible
//'   if the alphabet size can be reduced. For example, DNA sequences use at most 19 distinct
//'   characters and usually only 4; protein sequences use at most 26 distinct characters and
//'   usually only 20. Set the \code{alphabet} parameter if a reduced alphabet is appropriate.
//'
//'   UTF-8 (Unicode) matching is not currently supported.
//' @seealso
//' \itemize{
//' \item \href{https://www.codeproject.com/Articles/12383/Aho-Corasick-string-matching-in-C}{Aho-Corasick string matching in C#} for the article this package is based on
//' \item \code{\link[Biostrings]{matchPDict}} for a more memory efficient, but DNA-only, implementation of the algorithm
//' }
//' @examples
//' listEquals = function(a, b) { is.null(unlist(a)) && is.null(unlist(b)) ||
//'                               !is.null(a) && !is.null(b) && all(unlist(a) == unlist(b)) }
//'
//' # 1. Search for multiple keywords in a single text
//' keywords = c("Abra", "cadabra", "is", "the", "Magic", "Word")
//' oneSearch = AhoCorasickSearch(keywords, "Is Abracadabra the Magic Word?")
//' stopifnot(listEquals(oneSearch[[1]][[1]], list(keyword="Abra", offset=4)))
//' stopifnot(listEquals(oneSearch[[1]][[2]], list(keyword="cadabra", offset=8)))
//' stopifnot(listEquals(oneSearch[[1]][[3]], list(keyword="the", offset=16)))
//' stopifnot(listEquals(oneSearch[[1]][[4]], list(keyword="Magic", offset=20)))
//' stopifnot(listEquals(oneSearch[[1]][[5]], list(keyword="Word", offset=26)))
//'
//' # 2. Search multiple named texts in a named list with keyword grouping and aminoacid alphabet
//' # * all matches to a keyword are accessed by name
//' # * non-matched keywords are dropped
//' proteins = c(protein1="PEPTIDEPEPTIDEDADADARARARARAKEKEKEKEPEPTIDE",
//'              protein2="DERPADERPAPEWPEWPEEPEERAWRAWWARRAGTAGPEPTIDEKESEQUENCE")
//' peptides = c("PEPTIDE", "DERPA", "SEQUENCE", "KEKE", "PEPPIE")
//'
//' peptideSearch = AhoCorasickSearch(peptides, proteins, alphabet="aminoacid", groupByKeyword=TRUE)
//' stopifnot(listEquals(peptideSearch$PEPTIDE, list(list(keyword="protein1", offset=1),
//'                                                  list(keyword="protein1", offset=8),
//'                                                  list(keyword="protein1", offset=37),
//'                                                  list(keyword="protein2", offset=38))))
//' stopifnot(listEquals(peptideSearch$DERPA, list(list(keyword="protein2", offset=1),
//'                                                list(keyword="protein2", offset=6))))
//' stopifnot(listEquals(peptideSearch$SEQUENCE, list(list(keyword="protein2", offset=47))))
//' stopifnot(listEquals(peptideSearch$KEKE, list(list(keyword="protein1", offset=29),
//'                                               list(keyword="protein1", offset=31),
//'                                               list(keyword="protein1", offset=33))))
//' stopifnot(listEquals(peptideSearch$PEPPIE, NULL))
//'
//' # 3. Grouping by keyword without text names: offsets are given without reference to the text
//' names(proteins) = NULL
//' peptideSearch = AhoCorasickSearch(peptides, proteins, groupByKeyword=TRUE)
//' stopifnot(listEquals(peptideSearch$PEPTIDE, list(1, 8, 37, 38)))
//' stopifnot(listEquals(peptideSearch$DERPA, list(1, 6)))
//' stopifnot(listEquals(peptideSearch$SEQUENCE, list(47)))
//' stopifnot(listEquals(peptideSearch$KEKE, list(29, 31, 33)))
//' @export
// [[Rcpp::export]]
Rcpp::List AhoCorasickSearch(Rcpp::StringVector keywords,
                             Rcpp::StringVector text,
                             std::string alphabet = "ascii",
                             bool groupByKeyword = false,
                             int iterationFeedback = 0)
{
  return AhoCorasickSearchList(keywords, Rcpp::List::create(text), alphabet, groupByKeyword, iterationFeedback)(0);
}

/*** R
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

  biggerDNA = rep(bigDNA, times=2000)
  foo=AhoCorasickSearch(bigProbes, biggerDNA, alphabet="nucleicacid", groupByKeyword=T, iterationFeedback=5000)
}

*/
