//----------------------------------------
// EigenLab Speed Test
// Author: Dr. Marcel Paz Goldschen-Ohm (based on scripts by Ilja Honkonen <ilja.j.honkonen@nasa.gov>)
// Email:  marcel.goldschen@gmail.com
// Copyright (c) 2014 by Dr. Marcel Paz Goldschen-Ohm.
// Licence: MIT
//----------------------------------------

#include "EigenLab.h"
#include "mpParser.h"
#include <array>
#include <chrono>
#include <iostream>
#include <vector>

#ifndef EIGEN_NDEBUG
#	define EIGEN_NDEBUG
#endif

int main(int argc, const char * argv[])
{
	// Test muParser.
	// Evaluate "i+1" one million times.
	{
		mup::ParserX parser;//(mup::pckCOMMON | mup::pckNON_COMPLEX | mup::pckMATRIX | mup::pckUNIT);
		mup::Value value1{0};
		mup::Variable variable1{& value1};
		parser.DefineVar("v1", variable1);
		parser.SetExpr("v1 + 1");
		const auto timeStart = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < 1000000; i++) {
			value1 = i;
			if(parser.Eval().GetInteger() != i + 1) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Wrong value for variable: " << parser.Eval().GetInteger()
				<< ", should be: " << i + 1
				<< std::endl;
				abort();
			}
		}
		const auto timeEnd = std::chrono::high_resolution_clock::now();
		const auto total = std::chrono::duration_cast<std::chrono::duration<double> >(timeEnd - timeStart).count();
		std::cout << "1e6 evaluations of `i+1` took muParserX (s): " << total << std::endl;
	}
	
	// Test EigenLab.
	// Evaluate "i+1" one million times with expression caching off.
	{
		EigenLab::ParserXd parser;
		double variable1;
		parser.var("v1").setShared(& variable1);
		EigenLab::ValueXd result;
		parser.setCacheExpressions(false);
		const auto timeStart = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < 1000000; i++) {
			variable1 = i;
			result = parser.eval("v1 + 1");
			if(result.matrix()(0, 0) != i + 1) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Wrong value for variable: " << result.matrix()(0, 0)
				<< ", should be: " << i + 1
				<< std::endl;
				abort();
			}
		}
		const auto timeEnd = std::chrono::high_resolution_clock::now();
		const auto total = std::chrono::duration_cast<std::chrono::duration<double> >(timeEnd - timeStart).count();
		std::cout << "1e6 evaluations of `i+1` (expression caching off) took EigenLab (s): " << total << std::endl;
	}
	
	// Test EigenLab.
	// Evaluate "i+1" one million times with expression caching on.
	{
		EigenLab::ParserXd parser;
		double variable1;
		parser.var("v1").setShared(& variable1);
		EigenLab::ValueXd result;
		parser.setCacheExpressions(true);
		const auto timeStart = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < 1000000; i++) {
			variable1 = i;
			result = parser.eval("v1 + 1");
			if(result.matrix()(0, 0) != i + 1) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Wrong value for variable: " << result.matrix()(0, 0)
				<< ", should be: " << i + 1
				<< std::endl;
				abort();
			}
		}
		const auto timeEnd = std::chrono::high_resolution_clock::now();
		const auto total = std::chrono::duration_cast<std::chrono::duration<double> >(timeEnd - timeStart).count();
		std::cout << "1e6 evaluations of `i+1` (expression caching on) took EigenLab (s): " << total << std::endl;
	}
	
	// Test EigenLab.
	// Evaluate "i+1" on an array of one million values.
	{
		EigenLab::ParserXd parser;
		const auto timeStart = std::chrono::high_resolution_clock::now();
		Eigen::MatrixXd variable1(1, 1000000);
		parser.var("v1").setShared(variable1);
		for(int i = 0; i < 1000000; i++) {
			variable1(0, i) = i;
		}
		parser.eval("v1 = v1 + 1");
		for(int i = 0; i < 1000000; i++) {
			if(variable1(i) != i + 1) {
				std::cerr <<  __FILE__ << "(" << __LINE__<< ") "
				<< "Wrong value for variable: " << variable1(i)
				<< ", should be: " << i + 1
				<< std::endl;
				abort();
			}
		}
		const auto timeEnd = std::chrono::high_resolution_clock::now();
		const auto total = std::chrono::duration_cast<std::chrono::duration<double> >(timeEnd - timeStart).count();
		std::cout << "Creation of array i=[0:1e6-1] and evaluation of `i+1` took EigenLab (s): " << total << std::endl;
	}
	
    return 0;
}
