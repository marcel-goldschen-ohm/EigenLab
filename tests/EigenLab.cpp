// --*-  Mode: C++; c-basic-offset:4; indent-tabs-mode:t; tab-width:4 -*--
// EigenLab
// Version: 1.0.0
// Author: Michael Tesch
// Email:  tesch1@gmail.com
// Copyright (c) 2018 by Michael Tesch
// Licence: MIT
//----------------------------------------
#include <regex>
#include <getopt.h>

int eigenlab_debug = 0;
long int eigenlab_maxmatrix = 128*1024*1024; // 128 MB
#define EIGENLAB_DEBUG eigenlab_debug
#define EIGENLAB_MAXMATRIX eigenlab_maxmatrix
#include "EigenLab.h"

//
// Simple command-line program to evaluate expressions using EigenLab.h
//

typedef EigenLab::ParserXd ParserType;


// Evaluate a single expression
bool eval(ParserType & parser, const std::string & expr)
{
	try {
		size_t pos = 0;
		while (pos < expr.size()) {
			size_t end = expr.find_first_of("|", pos);
			if (end == std::string::npos)
				end = expr.size();
			std::cout << parser.eval(expr.substr(pos, end-pos)).matrix() << "\n";
			pos = end + 1;
		}
	}
	catch (std::exception & ex) {
		std::cerr << "err:" << ex.what() << "\n";
	}
	return true;
}

void usage()
{
	std::cerr << "usage:\n"
		" -e,--expr <expression>    | parse the given expression\n"
		" -f,--file <filename>      | parse each line of filename as expression\n"
		" -l,--limit <bytes>        | limit matrices to this many bytes\n"
		" -h,--help                 | print this help message and exit\n"
		" -v,--verbose              | be more verbose\n";
	exit(-1);
}

int main(int argc, char *argv[])
{
	ParserType parser;
	
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"expr", 	required_argument, 0,  'e' },
			{"file", 	required_argument, 0,  'f' },
			{"limit", 	required_argument, 0,  'l' },
			{"help",    no_argument,       0,  'h' },
			{"verbose", no_argument,       0,  'v' },
			{0,         0,                 0,  0 }
		};

		int c = getopt_long(argc, argv, "e:f:hv", long_options, &option_index);
		if (c == -1)
			break;

		switch (c) {
		case 0:
			printf("option %s", long_options[option_index].name);
			if (optarg)
				printf(" with arg %s", optarg);
			printf("\n");
			break;

		case 'f': {
			printf("expressions from file '%s':\n", optarg);
			if (!strcmp(optarg, "-"))
				strcpy(optarg, "/dev/stdin");
			if (FILE * fp = fopen(optarg, "r")) {
				char * line = NULL;
				size_t len = 0;
				ssize_t nread;
				while ((nread = getline(&line, &len, fp)) != -1) {
					printf(">> %s", line);
					eval(parser, line);
				}

				free(line);
				fclose(fp);
			}
			else
				perror("fopen");
		}
			break;

		case 'e':
			printf("expression '%s'\n", optarg);
			eval(parser, optarg);
			break;

		case 'l':
			eigenlab_maxmatrix = atol(optarg);
			printf("max matrix bytes '%s' = %ld\n", optarg, eigenlab_maxmatrix);
			break;

		case 'v':
			eigenlab_debug++;
			break;

		case 'h':
		case '?':
			usage();
			break;

		default:
			printf("?? getopt returned character code 0%o ??\n", c);
		}
	}
	
	if (optind < argc) {
		printf("Error, unrecognized command line arguments: ");
		while (optind < argc)
			printf("%s ", argv[optind++]);
		printf("\n");
		return -1;
	}

	return 0;
}
