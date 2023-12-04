#include "Huffman.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>

HuffmanNode::HuffmanNode() {
    left = right = nullptr;
    id = '$';
    freq = 0;
}

HuffmanNode::HuffmanNode(char data, unsigned freq) {
    left = right = nullptr;
    this->id = data;
    this->freq = freq;
}

// For comparison of
// two heap nodes (needed in min heap)
struct compare {
    bool operator()(HuffmanNode* l, HuffmanNode* r) {
        return (l->freq > r->freq);
    }
};

// prototypes
int binary_to_decimal(const std::string& in);
std::string decimal_to_binary(int in);
uint32_t getMinimalAlphabetElementSize(uint32_t count);
uint32_t getMinimalAlphabetElementBitSize(uint32_t count);

const uint32_t maxAlphabetSize = 256;

Huffman::Huffman(const std::string& in, const std::string& out) {
    in_file_name = in;
    out_file_name = out;
}

void Huffman::encode() {
    const auto begin = std::chrono::high_resolution_clock::now();
    std::priority_queue<HuffmanNode*, std::vector<HuffmanNode*>, compare> minHeap;
    node_ptr root = nullptr;
    const auto freq = calcFreq(); // get frequencies of all letters in file
    for (const auto& it : freq) {
        minHeap.push(new HuffmanNode(it.first, it.second)); // fill initial state for priority queue
    }
    // create the Huffman tree with highest frequecy characher being leaf from bottom to top
    while (minHeap.size() > 1) {
        std::cout << '1' << std::endl;
        root = new HuffmanNode;
        root->freq = 0;
        root->left = minHeap.top();
        root->freq += minHeap.top()->freq;
        minHeap.pop();
        root->right = minHeap.top();
        root->freq += minHeap.top()->freq;
        minHeap.pop();
        minHeap.push(root);
    }
    buildAlphabet(root, ""); // collect letters and their Huffman bit representation
    saveEncoded(); // save to output file
    const auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time ellapsed (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
        << std::endl;
    showStatistics();
    //showAlphabet();
}

void Huffman::decode() {
    const auto begin = std::chrono::high_resolution_clock::now();
    std::fstream in_file;
    in_file.open(in_file_name, std::ios::in | std::ios::binary);
    unsigned char size;
    in_file.read(reinterpret_cast<char*>(&size), 1); // get number of nodes in Huffman tree and restore it
    // increase by 1 to work properly with 256 alphabet
    uint32_t fixedSize = size + 1;
    unsigned char* h_code_c = new unsigned char[getMinimalAlphabetElementSize(fixedSize)];
    for (uint32_t i = 0; i < fixedSize; i++) {
        std::cout << '1' << std::endl;
        char a_code;
        in_file.read(&a_code, 1); // get ascii char representation
        in_file.read(reinterpret_cast<char*>(h_code_c),
            getMinimalAlphabetElementSize(fixedSize)); // obtain the binary code
        std::string h_code_s;
        for (uint32_t k = 0; k < getMinimalAlphabetElementSize(fixedSize); k++) { // obtain the oringinal 128-bit binary string
            h_code_s += decimal_to_binary(h_code_c[k]);
        }
        int j = 0;
        while (h_code_s[j] == '0') { // delete the added '000����1' to get the real Huffman code
            j++;
        }
        h_code_s = h_code_s.substr(j + 1); // remove leading '1' which signals Huffman code start
        reverseAlphabet[h_code_s] = a_code;
    }
    delete[] h_code_c;
    std::streampos encodedStart = in_file.tellg(); // save current position - encoded text starts here
    // jump to the last one byte to get the number of '0' append to the string at last
    in_file.seekg(-1, std::ios::end);
    char count0;
    in_file.read(&count0, 1);
    in_file.seekg(encodedStart, std::ios::beg);  // jump to the position where text starts
    std::vector<unsigned char> text;
    unsigned char textseg;
    in_file.read(reinterpret_cast<char*>(&textseg), 1);
    while (!in_file.eof()) {  // get encoded text byte by byte
        text.push_back(textseg);
        in_file.read(reinterpret_cast<char*>(&textseg), 1);
    }
    text.pop_back();  // remove last byte with zero bit count
    in_file.close();
    saveDecoded(text, count0);
    const auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time ellapsed (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
        << std::endl;
    showStatistics();
    //showReverseAlphabet();
}

void Huffman::saveEncoded() {
    // save alphabet
    // alphabet size in 1 byte
    // elements array
    // element contains char represented and byte with Huffman code in format "(optional zeros)1(Huffman code)"
    if (alphabet.size() == 0) {
        std::cout << "Nothing to save, exiting..." << std::endl;
        return;
    }
    // init with alphabet size (256 is max for 1 byte char, but byte contains 255 at max, so decrease by 1)
    std::string result;
    result.push_back((char)(alphabet.size() - 1));
    for (const auto& a : alphabet) {
        result.push_back(a.first);
        std::string s(getMinimalAlphabetElementBitSize(alphabet.size()) - 1 - a.second.size(),
            '0');  // set the codes with a fixed maxAlphabetSize-bit string form[000����1 + real code]
        s.push_back('1');    //'1' indicates the start of Huffman code
        s.append(a.second);
        result.push_back((char)binary_to_decimal(s.substr(0, 8)));
        for (uint32_t i = 0; i < getMinimalAlphabetElementSize(alphabet.size()) - 1;
            i++) {  // cut into 8-bit binary codes that can convert into saving char needed for binary file
            s = s.substr(8);
            result.push_back((char)binary_to_decimal(s.substr(0, 8)));
        }
    }
    // save encoded text
    std::ifstream input(in_file_name, std::ios::binary);
    std::vector<char> buffer(std::istreambuf_iterator<char>(input), {});
    std::string s;
    for (const auto& c : buffer) {
        s += alphabet[c];
        while (s.size() > 8) {  // cut into 8-bit binary codes that can convert into saving char needed for
                                // binary file
            result.push_back((char)binary_to_decimal(s.substr(0, 8)));
            s = s.substr(8);
        }
    }
    // after last encoded symbol put count of empty bits in last byte
    int count = 8 - s.size();
    if (s.size() < 8) {  // append number of 'count' '0' to the last few codes to create the last byte of text
        s.append(count, '0');
    }
    result.push_back((char)binary_to_decimal(s));  // save number of 'count' at last
    result.push_back(count);
    std::ofstream output(out_file_name, std::ios::out | std::ios::binary);
    output << result;
}

void Huffman::saveDecoded(const std::vector<unsigned char>& text, char count0) {
    std::string processed;
    std::fstream out_file;
    out_file.open(out_file_name, std::ios::out | std::ios::binary);
    for (size_t i = 0; i < text.size(); i++) {    // translate the Huffman code
        processed += decimal_to_binary(text[i]);  // add to processed text new byte
        if (i == text.size() - 1) {
            processed = processed.substr(0, processed.size() - count0);
        }
        std::string search;
        bool searchAgain = true;
        while (searchAgain) {
            searchAgain = false;
            for (size_t j = 0; j < processed.size(); j++) {
                search += processed[j];
                const auto it = reverseAlphabet.find(search);
                if (it != reverseAlphabet.end()) {
                    out_file.put(reverseAlphabet[search]);
                    processed = processed.substr(search.size());
                    search.clear();
                    searchAgain = !processed.empty();
                    break;
                }
            }
        }
    }
    std::cout << std::endl;  // print new line after progress bar
    out_file.close();
}

// recursively collect codes for every leaf (symbol) in tree
void Huffman::buildAlphabet(node_ptr node, std::string code) {
    if (node->left == nullptr && node->right == nullptr) {
        alphabet[node->id] = code;
    }
    else {
        buildAlphabet(node->left, code + '0');
        buildAlphabet(node->right, code + '1');
    }
}

// helper function for letter's frequencies discovery
std::map<char, uint32_t> Huffman::calcFreq() {
    std::ifstream input(in_file_name, std::ios::binary);
    std::vector<char> buffer(std::istreambuf_iterator<char>(input), {});
    std::map<char, uint32_t> freq;
    for (const auto& c : buffer) {
        freq[c]++;
    }
    return freq;
}

// console output with additional info >>>>
void Huffman::showAlphabet() {
    std::cout << "Alphabet (" << alphabet.size() << " elements)" << std::endl;
    for (const auto& a : alphabet) {
        std::cout << "\'" << a.first << "\'\t" << uint32_t((unsigned char)(a.first)) << '\t' << a.second << std::endl;
    }
}

// sort unordered_map by value >>>>>
template <typename A, typename B>
using wrapped_pair = std::pair<std::reference_wrapper<A>, std::reference_wrapper<B>>;

struct cmp_wrapped_pair {
    template <typename A, typename B>
    bool operator()(const A& a, const B& b) const {
        return a.first.get() < b.first.get() || (!(b.first.get() < a.first.get()) && a.second.get() < b.second.get());
    }
};

template <typename MAP_TYPE>
auto sorted_snapshot_of(const MAP_TYPE& map) { // auto: function return type deduction is C++14
    using KEY = typename MAP_TYPE::key_type;
    using VALUE = typename MAP_TYPE::mapped_type;
    std::multiset<wrapped_pair<const VALUE, const KEY>, cmp_wrapped_pair> snapshot;

    for (auto& pair : map)
        snapshot.emplace(pair.second, pair.first);
    return snapshot;
}
// sort unordered_map by value <<<<<

void Huffman::showReverseAlphabet() {
    std::cout << "ReverseAlphabet (" << reverseAlphabet.size() << " elements)" << std::endl;
    for (const auto& a : sorted_snapshot_of(reverseAlphabet)) {
        std::cout << "\'" << a.first.get() << "\'\t" << uint32_t((unsigned char)(a.first.get())) << '\t' << a.second.get()
            << std::endl;
    }
}

void Huffman::showStatistics() {
    std::filesystem::path in = std::filesystem::current_path() / in_file_name;
    std::cout << "Input file size = " << std::filesystem::file_size(in) << std::endl;
    std::filesystem::path out = std::filesystem::current_path() / out_file_name;
    std::cout << "Output file size = " << std::filesystem::file_size(out) << std::endl;
    std::cout << "Output to input ratio is "
        << std::filesystem::file_size(out) / (double)(std::filesystem::file_size(in)) * 100 << "%" << std::endl;
}
// console output with additional info <<<<<

// converters >>>>>
// converts string of '1' and '0' bits to int value
int binary_to_decimal(const std::string& in) {
    int result = 0;
    for (size_t i = 0; i < in.size(); i++)
        result = result * 2 + in[i] - '0';
    return result;
}

// converts int value to bit string
std::string decimal_to_binary(int in) {
    std::string temp = "";
    std::string result = "";
    while (in) {
        temp += ('0' + in % 2);
        in /= 2;
    }
    result.append(8 - temp.size(), '0');  // append '0' ahead to let the result become fixed length of 8
    for (int i = temp.size() - 1; i >= 0; i--) {
        result += temp[i];
    }
    return result;
}
// converters <<<<<

// how many bytes needed to represent symbol using Huffman code with _count_ elements in alphabet
uint32_t getMinimalAlphabetElementSize(uint32_t count) {
    return count / 8 + 1;
}

// get bites for operation above
uint32_t getMinimalAlphabetElementBitSize(uint32_t count) {
    return getMinimalAlphabetElementSize(count) * 8;
}
