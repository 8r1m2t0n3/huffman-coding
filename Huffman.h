#pragma once
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

// A Huffman tree node
struct HuffmanNode {
	// One of the input characters
	char id;
	// Frequency of the character
	unsigned freq;
	// Left and right child
	HuffmanNode* left, * right;

	HuffmanNode();
	HuffmanNode(char data, unsigned freq);
};

typedef HuffmanNode* node_ptr;

class Huffman {
public:
	Huffman(const std::string& in, const std::string& out);
	void encode();
	void decode();

private:
	void saveEncoded();
	void saveDecoded(const std::vector<unsigned char>& text, char count0);

	// recursively collect codes for every leaf (symbol) in tree
	void buildAlphabet(node_ptr node, std::string code);

	// helper function for letter's frequencies discovery
	std::map<char, uint32_t> calcFreq();

	// console output with additional info >>>>
	void showAlphabet();
	void showReverseAlphabet();
	void showStatistics();
	// console output with additional info <<<<<

	std::string in_file_name;
	std::string out_file_name;
	std::map<char, std::string> alphabet;
	// std::unordered_map<std::string, char> reverseAlphabet;
	std::unordered_map<std::string, char> reverseAlphabet;
};
