# cancer_mutation_calling
use vardict, strelka, and mutect, keep mutations by consensus voting, force call downstream mutations, final output in maflite format for Firehose use

To compile forcecall.cpp: g++ -std=c++11 forcecall.cpp -o forcecall
