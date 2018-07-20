// --*-  Mode: C++; c-basic-offset:4; indent-tabs-mode:t; tab-width:4 -*--
// EigenLab
// Version: 1.0.0
// Author: Michael Tesch
// Email:  tesch1@gmail.com
// Copyright (c) 2018 by Michael Tesch
// Licence: MIT
//----------------------------------------
#include <getopt.h>
#include "EigenLab.h"

//
// Simple command-line program to evaluate expressions using EigenLab.h
//

typedef EigenLab::ParserXd ParserType;


// Evaluate a single expression
bool eval(ParserType & parser, const std::string & expr)
{
	try {
		std::cout << parser.eval(expr).matrix() << "\n";
	}
	catch (std::exception & ex) {
		std::cerr << "err:" << ex.what() << "\n";
	}
	return true;
}



int main(int argc, char *argv[])
{
	ParserType parser;
	
	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"file", 	required_argument, 0,  'f' },
			{"expr", 	required_argument, 0,  'e' },
			{"verbose", no_argument,       0,  'v' },
			{0,         0,                 0,  0 }
		};

		int c = getopt_long(argc, argv, "f:e:v", long_options, &option_index);
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
					printf("%s:\n", line);
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

		case '?':
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
