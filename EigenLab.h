//----------------------------------------
// EigenLab
// Version: 0.9.0
// Author: Dr. Marcel Paz Goldschen-Ohm
// Email:  marcel.goldschen@gmail.com
// Copyright (c) 2014 by Dr. Marcel Paz Goldschen-Ohm.
// Licence: MIT
//----------------------------------------

#ifndef EigenLab_H
#define EigenLab_H

#include <Eigen/Dense>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <iostream>
// Define both DEBUG and EIGENLAB_DEBUG for step-by-step equation parsing printouts.
#ifdef DEBUG
//#define EIGENLAB_DEBUG
#ifdef EIGENLAB_DEBUG
#include <iostream>
#endif
#endif

namespace EigenLab
{
	//----------------------------------------
	// A wrapper for a matrix whose data is either stored locally or shared.
	//
	// Template typename Derived can be any dynamically sized matrix type supported by Eigen.
	//
	// !!! matrix() promises to ALWAYS return a map to the matrix data whether it's
	// stored locally or externally in some shared memory.
	//
	// !!! local() provides direct access to the local data, but this data is
	// ONLY valid when isLocal() is true. In most cases, you're best off
	// accessing the matrix data via matrix() instead.
	//----------------------------------------
	template <typename Derived = Eigen::MatrixXd>
	class Value
	{
	private:
		// Local matrix data.
		Derived mLocal;
		
		// Map to shared matrix data (map points to mLocal if the data is local).
		// !!! This map promises to ALWAYS point to the matrix data whether it's
		// stored locally in mLocal or externally in some shared memory.
		Eigen::Map<Derived> mShared;
		
		// Flag indicating whether the local data is being used.
		bool mIsLocal;
		
	public:
		// Access mapped data (whether its local or shared).
		// !!! matrix() promises to ALWAYS return a map to the matrix data whether it's
		// stored locally in mLocal or externally in some shared memory.
		inline Eigen::Map<Derived> & matrix() { return mShared; }
		inline const Eigen::Map<Derived> & matrix() const { return mShared; }
		
		// Access local data.
		// !!! WARNING! This data is ONLY valid if isLocal() is true.
		// !!! WARNING! If you change the local data via this method, you MUST call mapLocal() immediately afterwards.
		// In most cases, you're best off accessing the matrix data via matrix() instead.
		inline Derived & local() { return mLocal; }
		inline const Derived & local() const { return mLocal; }
		
		// Is mapped data local?
		inline bool isLocal() const { return mIsLocal; }
		
		// Set mapped data to point to local data.
		inline void mapLocal() { new (& mShared) Eigen::Map<Derived>(mLocal.data(), mLocal.rows(), mLocal.cols()); mIsLocal = true; }
		
		//  Copy shared data to local data (if needed).
		inline void copySharedToLocal() { if(!isLocal()) { mLocal = mShared; mapLocal(); } }
		
		// Set local data.
		Value() : mLocal(1, 1), mShared(mLocal.data(), mLocal.rows(), mLocal.cols()), mIsLocal(true) {}
		Value(const typename Derived::Scalar s) : mLocal(Derived::Constant(1, 1, s)), mShared(mLocal.data(), mLocal.rows(), mLocal.cols()), mIsLocal(true) {}
		Value(const Derived & mat) : mLocal(mat), mShared(mLocal.data(), mLocal.rows(), mLocal.cols()), mIsLocal(true) {}
		inline void setLocal(const typename Derived::Scalar s) { mLocal = Derived::Constant(1, 1, s); mapLocal(); }
		inline void setLocal(const Eigen::MatrixBase<Derived> & mat) { mLocal = mat; mapLocal(); }
		inline void setLocal(const Value & val) { mLocal = val.matrix(); mapLocal(); }
		inline void setLocal(const typename Derived::Scalar * data, size_t rows = 1, size_t cols = 1) { setShared(data, rows, cols); copySharedToLocal(); }
		inline Value & operator = (const typename Derived::Scalar s) { setLocal(s); return (* this); }
		inline Value & operator = (const Derived & mat) { setLocal(mat); return (* this); }
		
		// Set shared data.
		Value(const typename Derived::Scalar * data, size_t rows = 1, size_t cols = 1) : mShared(const_cast<typename Derived::Scalar *>(data), rows, cols), mIsLocal(false) {}
		inline void setShared(const typename Derived::Scalar * data, size_t rows = 1, size_t cols = 1) { new (& mShared) Eigen::Map<Derived>(const_cast<typename Derived::Scalar *>(data), rows, cols); mIsLocal = false; }
		inline void setShared(const Derived & mat) { setShared(mat.data(), mat.rows(), mat.cols()); }
		inline void setShared(const Value & val) { setShared(val.matrix().data(), val.matrix().rows(), val.matrix().cols()); }
		
		// Set to local or shared data dependent on whether val maps its own local data or some other shared data.
		Value(const Value & val) : mLocal(1, 1), mShared(mLocal.data(), mLocal.rows(), mLocal.cols()) { (* this) = val; }
		inline Value & operator = (const Value & val) { if(val.isLocal()) { setLocal(val); } else { setShared(val); } return (* this); }
	};
	typedef Value<Eigen::MatrixXd> ValueXd;
	typedef Value<Eigen::MatrixXf> ValueXf;
	typedef Value<Eigen::MatrixXi> ValueXi;
	
	//----------------------------------------
	// Equation parser.
	//
	// Template typename Derived can be any dynamically sized matrix type supported by Eigen.
	//----------------------------------------
	template <typename Derived = Eigen::MatrixXd>
	class Parser
	{
	public:
		// A map to hold named values.
		typedef std::map<std::string, Value<Derived> > ValueMap;
		
	private:
		// Variables are stored in a map of named values.
		ValueMap mVariables;
		
		// Operator symbols and function names used by the parser.
		std::string mOperators1, mOperators2;
		std::vector<std::string> mFunctions;
		
		// Expressions are parsed by first splitting them into chunks.
		struct Chunk {
			std::string field;
			int type;
			Value<Derived> value;
			int row0, col0, rows, cols;
			Chunk(const std::string & str = "", int t = -1, const Value<Derived> & val = Value<Derived>()) : field(str), type(t), value(val), row0(-1), col0(-1), rows(-1), cols(-1) {}
		};
		enum ChunkType { VALUE = 0, VARIABLE, OPERATOR, FUNCTION };
		typedef std::vector<Chunk> ChunkArray;
		
	public:
		// Access to named variables.
		// !!! WARNING! var(name) will create the variable name if it does not already exist.
		inline ValueMap & vars() { return mVariables; }
		inline Value<Derived> & var(const std::string & name) { return mVariables[name]; }
		
		// Check if a variable exists.
		inline bool hasVar(const std::string & name) { return isVariable(name); }
		
		// Delete a variable.
		inline void clearVar(const std::string & name) { typename ValueMap::iterator it = mVariables.find(name); if(it != mVariables.end()) mVariables.erase(it); }
		
		Parser();
		
		// Evaluate an expression and return the result in a value wrapper.
		Value<Derived> eval(const std::string & expression);
		
	private:
		void splitEquationIntoChunks(const std::string & expression, ChunkArray & chunks, std::string & code);
		std::string::const_iterator findClosingBracket(const std::string & str, const std::string::const_iterator openingBracket, const char closingBracket) const;
		std::vector<std::string> splitArguments(const std::string & str, const char delimeter) const;
		void evalIndexRange(const std::string & str, int * first, int * last);
		void evalMatrixExpression(const std::string & str, Value<Derived> & mat);
		void evalFunction(const std::string & name, std::vector<std::string> & args, Value<Derived> & result);
		void evalNumericRange(const std::string & str, Value<Derived> & mat);
		inline bool isVariable(const std::string & name) const { return mVariables.count(name) > 0; }
		inline bool isOperator(const char c) const { return (std::find(mOperators1.begin(), mOperators1.end(), c) != mOperators1.end()); }
		bool isOperator(const std::string & str) const;
		inline bool isFunction(const std::string & str) const { return (std::find(mFunctions.begin(), mFunctions.end(), str) != mFunctions.end()); }
		void evalIndices(ChunkArray & chunks);
		void evalNegations(ChunkArray & chunks);
		void evalPowers(ChunkArray & chunks);
		void evalMultiplication(ChunkArray & chunks);
		void evalAddition(ChunkArray & chunks);
		void evalAssignment(ChunkArray & chunks);
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		void printChunks(ChunkArray & chunks, size_t maxRows = 2, size_t maxCols = 2, int precision = 0);
		void printVars(size_t maxRows = 2, size_t maxCols = 2, int precision = 0);
		std::string textRepresentation(Value<Derived> & val, size_t maxRows = 2, size_t maxCols = 2, int precision = 0);
#endif
#endif
		
	public:
		static std::string trim(const std::string & str);
		static std::vector<std::string> split(const std::string & str, const char delimeter);
		template <typename T> static bool isNumber(const std::string & str, T * num = 0);
		template <typename T> static T stringToNumber(const std::string & str);
		template <typename T> static std::string numberToString(T num, int precision = 0);
	};
	typedef Parser<Eigen::MatrixXd> ParserXd;
	typedef Parser<Eigen::MatrixXf> ParserXf;
	typedef Parser<Eigen::MatrixXi> ParserXi;
	
	//----------------------------------------
	// Function definitions.
	//----------------------------------------
	template <typename Derived>
	Parser<Derived>::Parser() : mOperators1("+-*/^()[]="), mOperators2(".+.-.*./.^")
	{
		// Coefficient-wise operations.
		mFunctions.push_back("abs");
		mFunctions.push_back("sqrt");
		mFunctions.push_back("exp");
		mFunctions.push_back("log");
		mFunctions.push_back("log10");
		mFunctions.push_back("sin");
		mFunctions.push_back("cos");
		mFunctions.push_back("tan");
		mFunctions.push_back("asin");
		mFunctions.push_back("acos");

		// Matrix reduction operations.
		mFunctions.push_back("trace");
		mFunctions.push_back("norm");
		mFunctions.push_back("size");
		mFunctions.push_back("min");
		mFunctions.push_back("max");
		mFunctions.push_back("mean");
		mFunctions.push_back("sum");
		mFunctions.push_back("prod");

		// Matrix operations.
		mFunctions.push_back("transpose");
		mFunctions.push_back("conjugate");
		mFunctions.push_back("adjoint");
	}
	
	template <typename Derived>
	Value<Derived> Parser<Derived>::eval(const std::string & expression)
	{
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		std::cout << "---" << std::endl;
		std::cout << "EXPRESSION: " << expression << std::endl;
#endif
#endif
		ChunkArray chunks;
		std::string code;
		splitEquationIntoChunks(trim(expression), chunks, code);
		evalIndices(chunks);
		evalNegations(chunks);
		evalPowers(chunks);
		evalMultiplication(chunks);
		evalAddition(chunks);
		evalAssignment(chunks);
		if(chunks.size() != 1)
			throw std::runtime_error("Failed to reduce expression '" + expression + "' to a single value.");
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		std::cout << "---" << std::endl;
#endif
#endif
		if(chunks[0].type == VARIABLE) {
			if(!isVariable(chunks[0].field))
				throw std::runtime_error("Unknown variable '" + chunks[0].field + "'.");
			return mVariables[chunks[0].field];
		}
		return chunks[0].value;
	}
	
	template <typename Derived>
	void Parser<Derived>::splitEquationIntoChunks(const std::string & expression, ChunkArray & chunks, std::string & code)
	{
		for(std::string::const_iterator it = expression.begin(); it != expression.end();)
		{
			int prevType = (chunks.size() ? chunks.back().type : -1);
			char ci = (* it);
			if(isspace(ci)) {
				// Ignore whitespace.
				it++;
			} else if(ci == '(' && (prevType == VALUE || prevType == VARIABLE)) {
				// Index group.
				std::string::const_iterator jt = findClosingBracket(expression, it, ')');
				if(jt == expression.end())
					throw std::runtime_error("Missing closing bracket for '" + std::string(it, jt) + "'.");
				std::string field = std::string(it + 1, jt); // Outer brackets stripped.
				int first, last;
				std::vector<std::string> args = splitArguments(field, ',');
				if(args.size() == 1) {
					if(chunks.back().value.matrix().cols() == 1) {
						evalIndexRange(args[0], & first, & last);
						chunks.back().row0 = first;
						chunks.back().col0 = 0;
						chunks.back().rows = (last == -1 ? int(chunks.back().value.matrix().rows()) : last + 1) - first;
						chunks.back().cols = 1;
					} else if(chunks.back().value.matrix().rows() == 1) {
						evalIndexRange(args[0], & first, & last);
						chunks.back().row0 = 0;
						chunks.back().col0 = first;
						chunks.back().rows = 1;
						chunks.back().cols = (last == -1 ? int(chunks.back().value.matrix().cols()) : last + 1) - first;
					} else {
						throw std::runtime_error("Missing row or column indices for '(" + chunks.back().field + "(" + field + ")'.");
					}
				} else if(args.size() == 2) {
					evalIndexRange(args[0], & first, & last);
					chunks.back().row0 = first;
					chunks.back().rows = (last == -1 ? int(chunks.back().value.matrix().rows()) : last + 1) - first;
					evalIndexRange(args[1], & first, & last);
					chunks.back().col0 = first;
					chunks.back().cols = (last == -1 ? int(chunks.back().value.matrix().cols()) : last + 1) - first;
				} else {
					throw std::runtime_error("Invalid index expression '" + chunks.back().field + "(" + field + ")'.");
				}
				code += "i";
				it = jt + 1;
			} else if(ci == '(' && prevType == FUNCTION) {
				// Function argument group.
				std::string::const_iterator jt = findClosingBracket(expression, it, ')');
				if(jt == expression.end())
					throw std::runtime_error("Missing closing bracket for '" + std::string(it, jt) + "'.");
				std::string field = std::string(it + 1, jt); // Outer brackets stripped.
				std::vector<std::string> args = splitArguments(field, ',');
				evalFunction(chunks.back().field, args, chunks.back().value);
				chunks.back().field += "(" + field + ")";
				chunks.back().type = VALUE;
				code += (chunks.back().value.matrix().size() == 1 ? "n" : "M");
				it = jt + 1;
			} else if(ci == '(') {
				// Recursively evaluate group to a single value.
				std::string::const_iterator jt = findClosingBracket(expression, it, ')');
				if(jt == expression.end())
					throw std::runtime_error("Missing closing bracket for '" + std::string(it, jt) + "'.");
				std::string field = std::string(it + 1, jt); // Outer brackets stripped.
				chunks.push_back(Chunk(field, prevType = VALUE, eval(field)));
				code += (chunks.back().value.matrix().size() == 1 ? "n" : "M");
				it = jt + 1;
			} else if(ci == '[') {
				// Evaluate matrix.
				if(prevType == VALUE || prevType == VARIABLE)
					throw std::runtime_error("Invalid operation '" + chunks.back().field + std::string(1, ci) + "'.");
				std::string::const_iterator jt = findClosingBracket(expression, it, ']');
				if(jt == expression.end())
					throw std::runtime_error("Missing closing bracket for '" + std::string(it, jt) + "'.");
				std::string field = std::string(it + 1, jt); // Outer brackets stripped.
				chunks.push_back(Chunk("[" + field + "]", prevType = VALUE));
				evalMatrixExpression(field, chunks.back().value);
				code += (chunks.back().value.matrix().size() == 1 ? "n" : "M");
				it = jt + 1;
			} else if(it + 1 != expression.end() && isOperator(std::string(it, it + 2))) {
				// Double character operator.
				std::string field = std::string(it, it + 2);
				chunks.push_back(Chunk(field, prevType = OPERATOR));
				code += field;
				it += 2;
			} else if(isOperator(ci)) {
				// Single character operator.
				std::string field = std::string(1, ci);
				chunks.push_back(Chunk(field, prevType = OPERATOR));
				code += field;
				it++;
			} else {
				// Non-operator: value range, number, function, or variable name.
				std::string::const_iterator jt = it + 1;
				for(; jt != expression.end(); jt++) {
					if(isOperator(*jt) || (jt + 1 != expression.end() && isOperator(std::string(jt, jt + 2))))
						break;
				}
				std::string field = trim(std::string(it, jt));
				if(prevType == VALUE || prevType == VARIABLE)
					throw std::runtime_error("Invalid operation '" + chunks.back().field + field + "'.");
				double num;
				if(field.find(":") != std::string::npos) {
					// Numeric range.
					chunks.push_back(Chunk(field, prevType = VALUE));
					evalNumericRange(field, chunks.back().value);
					code += (chunks.back().value.matrix().size() == 1 ? "n" : "M");
				} else if(isNumber<double>(field, & num)) {
					// Number.
					chunks.push_back(Chunk(field, prevType = VALUE, Value<Derived>(num)));
					code += "n";
				} else if(isVariable(field)) {
					// Local variable.
					chunks.push_back(Chunk(field, prevType = VARIABLE));
					code += (mVariables[field].matrix().size() == 1 ? "vn" : "vM");
				} else if(isFunction(field)) {
					// Function.
					chunks.push_back(Chunk(field, prevType = FUNCTION));
				} else {
					// New undefined variable.
					chunks.push_back(Chunk(field, prevType = VARIABLE));
					code += "vU";
				}
				it = jt;
			}
		} // it
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		std::cout << "CHUNKS: "; printChunks(chunks); std::cout << std::endl;
		std::cout << "CODE: " << code << std::endl;
#endif
#endif
	}
	
	template <typename Derived>
	std::string::const_iterator Parser<Derived>::findClosingBracket(const std::string & str, const std::string::const_iterator openingBracket, const char closingBracket) const
	{
		int depth = 1;
		std::string::const_iterator it = openingBracket + 1;
		for(; it != str.end(); it++) {
			if((*it) == (*openingBracket)) depth++;
			else if((*it) == closingBracket) depth--;
			if(depth == 0) return it;
		}
		return str.end();
	}

	template <typename Derived>
	std::vector<std::string> Parser<Derived>::splitArguments(const std::string & str, const char delimeter) const
	{
		std::vector<std::string> args;
		std::string::const_iterator i0 = str.begin();
		for(std::string::const_iterator it = str.begin(); it != str.end(); it++) {
			if((*it) == '(') it = findClosingBracket(str, it, ')');
			else if((*it) == '[') it = findClosingBracket(str, it, ']');
			else if((*it) == delimeter) {
				args.push_back(trim(std::string(i0, it)));
				i0 = it + 1;
			}
		}
		args.push_back(std::string(i0, str.end()));
		return args;
	}

	template <typename Derived>
	void Parser<Derived>::evalIndexRange(const std::string & str, int * first, int * last)
	{
		if(str.empty())
			throw std::runtime_error("Empty index range.");
		ValueXi valuei;
		ParserXi parseri;
		for(std::string::const_iterator it = str.begin(); it != str.end(); it++) {
			if((*it) == ':') {
				std::string firstStr = trim(std::string(str.begin(), it));
				std::string lastStr = trim(std::string(it + 1, str.end()));
				if(firstStr.empty() && lastStr.empty()) {
					(* first) = 0;
					(* last) = -1;
					return;
				}
				if(firstStr.empty() || lastStr.empty())
					throw std::runtime_error("Missing indices for '" + str + "'.");
				valuei = parseri.eval(firstStr);
				if(valuei.matrix().size() != 1)
					throw std::runtime_error("Invalid indices '" + str + "'.");
				(* first) = valuei.matrix()(0, 0);
				if(lastStr == "end")
					(* last) = -1;
				else {
					valuei = parseri.eval(lastStr);
					if(valuei.matrix().size() != 1)
						throw std::runtime_error("Invalid indices '" + str + "'.");
					(* last) = valuei.matrix()(0, 0);
				}
				return;
			}
		}
		valuei = parseri.eval(str);
		if(valuei.matrix().size() != 1)
			throw std::runtime_error("Invalid index '" + str + "'.");
		(* first) = valuei.matrix()(0, 0);
		(* last) = (* first);
	}

	template <typename Derived>
	void Parser<Derived>::evalMatrixExpression(const std::string & str, Value<Derived> & mat)
	{
		// !!! Expression may NOT include outer brackets, although brackets for individual rows are OK.
		std::vector<std::string> rows = splitArguments(str, ';');
		std::vector<std::vector<typename Derived::Scalar> > temp;
		Value<Derived> submatrix;
		size_t row0 = 0, col0 = 0, nrows = 0, ncols = 0;
		for(size_t i = 0; i < rows.size(); i++) {
			// Strip row brackets if they exist.
			if(rows[i][0] == '[' && rows[i].back() == ']') rows[i] = rows[i].substr(1, int(rows[i].size()) - 2);
			std::vector<std::string> cols = splitArguments(rows[i], ',');
			col0 = 0;
			ncols = 0;
			for(size_t j = 0; j < cols.size(); j++) {
				submatrix = eval(cols[j]);
				if(j > 0 && submatrix.matrix().cols() != nrows)
					throw std::runtime_error("Invalid matrix definition '[" + str + "]'. Successive column entries '" + cols[int(j) - 1] + "' and '" + cols[j] + "' do not have the same number of rows.");
				nrows = submatrix.matrix().rows();
				ncols += submatrix.matrix().cols();
				temp.resize(row0 + submatrix.matrix().rows());
				for(size_t row = 0; row < submatrix.matrix().rows(); row++) {
					temp[row0 + row].resize(col0 + submatrix.matrix().cols());
					for(size_t col = 0; col < submatrix.matrix().cols(); col++)
						temp[row0 + row][col0 + col] = submatrix.matrix()(row, col);
				}
				col0 += submatrix.matrix().cols();
			}
			if(row0 > 0 && ncols != temp[int(row0) - 1].size())
				throw std::runtime_error("Invalid matrix definition '[" + str + "]'. Successive row entries '" + rows[int(i) - 1] + "' and '" + rows[i] + "' do not have the same number of columns.");
			row0 += nrows;
		}
		nrows = temp.size();
		if(nrows == 0) return;
		ncols = temp[0].size();
		mat.setLocal(Derived(nrows, ncols));
		for(size_t row = 0; row < nrows; row++) {
			for(size_t col = 0; col < ncols; col++)
				mat.matrix()(row, col) = temp[row][col];
		}
		mat.mapLocal();
	}
	
	template <typename Derived>
	void Parser<Derived>::evalFunction(const std::string & name, std::vector<std::string> & args, Value<Derived> & result)
	{
		if(args.size() == 1) {
			Value<Derived> arg = eval(args[0]);
			if(name == "abs") {
				result.local() = arg.matrix().array().abs();
				result.mapLocal();
				return;
			} else if(name == "sqrt") {
				result.local() = arg.matrix().array().sqrt();
				result.mapLocal();
				return;
			} else if(name == "exp") {
				result.local() = arg.matrix().array().exp();
				result.mapLocal();
				return;
			} else if(name == "log") {
				result.local() = arg.matrix().array().log();
				result.mapLocal();
				return;
			} else if(name == "log10") {
				result.local() = arg.matrix().array().log() * (1.0 / log(10));
				result.mapLocal();
				return;
			} else if(name == "sin") {
				result.local() = arg.matrix().array().sin();
				result.mapLocal();
				return;
			} else if(name == "cos") {
				result.local() = arg.matrix().array().cos();
				result.mapLocal();
				return;
			} else if(name == "tan") {
				result.local() = arg.matrix().array().tan();
				result.mapLocal();
				return;
			} else if(name == "acos") {
				result.local() = arg.matrix().array().acos();
				result.mapLocal();
				return;
			} else if(name == "asin") {
				result.local() = arg.matrix().array().asin();
				result.mapLocal();
				return;
			} else if(name == "trace") {
				result.setLocal(arg.matrix().trace());
				return;
			} else if(name == "norm") {
				result.setLocal(arg.matrix().norm());
				return;
			} else if(name == "min") {
				result.setLocal(arg.matrix().minCoeff());
				return;
			} else if(name == "max") {
				result.setLocal(arg.matrix().maxCoeff());
				return;
			} else if(name == "mean") {
				result.setLocal(arg.matrix().mean());
				return;
			} else if(name == "sum") {
				result.setLocal(arg.matrix().sum());
				return;
			} else if(name == "prod") {
				result.setLocal(arg.matrix().prod());
				return;
			} else if(name == "transpose") {
				result.local() = arg.matrix().transpose();
				result.mapLocal();
				return;
			} else if(name == "conjugate") {
				result.local() = arg.matrix().conjugate();
				result.mapLocal();
				return;
			} else if(name == "adjoint") {
				result.local() = arg.matrix().adjoint();
				result.mapLocal();
				return;
			}
			else {
				throw std::runtime_error("Invalid function '" + name + "(" + args[0] + ")'.");
			}
		} else if(args.size() == 2) {
			Value<Derived> arg0 = eval(args[0]);
			ParserXi parseri;
			ValueXi arg1 = parseri.eval(args[1]);
			if(arg1.matrix().size() != 1)
				throw std::runtime_error("Invalid dimension argument for function '" + name + "(" + args[0] + "," + args[1] + ")'.");
			int dim = arg1.matrix()(0, 0);
			if(dim != 0 && dim != 1)
				throw std::runtime_error("Invalid dimension argument for function '" + name + "(" + args[0] + "," + args[1] + ")'.");
			if(name == "size") {
				if(dim == 0) {
					result.setLocal((typename Derived::Scalar) arg0.matrix().rows());
					return;
				} else if(dim == 1) {
					result.setLocal((typename Derived::Scalar) arg0.matrix().cols());
					return;
				}
			} else if(name == "min") {
				if(dim == 0) {
					result.local() = arg0.matrix().colwise().minCoeff();
					result.mapLocal();
					return;
				} else if(dim == 1) {
					result.local() = arg0.matrix().rowwise().minCoeff();
					result.mapLocal();
					return;
				}
			} else if(name == "max") {
				if(dim == 0) {
					result.local() = arg0.matrix().colwise().maxCoeff();
					result.mapLocal();
					return;
				} else if(dim == 1) {
					result.local() = arg0.matrix().rowwise().maxCoeff();
					result.mapLocal();
					return;
				}
			} else if(name == "mean") {
				if(dim == 0) {
					result.local() = arg0.matrix().colwise().mean();
					result.mapLocal();
					return;
				} else if(dim == 1) {
					result.local() = arg0.matrix().rowwise().mean();
					result.mapLocal();
					return;
				}
			} else if(name == "sum") {
				if(dim == 0) {
					result.local() = arg0.matrix().colwise().sum();
					result.mapLocal();
					return;
				} else if(dim == 1) {
					result.local() = arg0.matrix().rowwise().sum();
					result.mapLocal();
					return;
				}
			} else if(name == "prod") {
				if(dim == 0) {
					result.local() = arg0.matrix().colwise().prod();
					result.mapLocal();
					return;
				} else if(dim == 1) {
					result.local() = arg0.matrix().rowwise().prod();
					result.mapLocal();
					return;
				}
			} else {
				throw std::runtime_error("Invalid function '" + name + "(" + args[0] + "," + args[1] + ")'.");
			}
		}
		std::string argsStr = "(";
		for(size_t i = 0; i < args.size(); i++) {
			if(i > 0) argsStr += ",";
			argsStr += args[i];
		}
		argsStr += ")";
		throw std::runtime_error("Invalid function/arguments for '" + name + argsStr + "'.");
	}
	
	template <typename Derived>
	void Parser<Derived>::evalNumericRange(const std::string & str, Value<Derived> & mat)
	{
		size_t pos = str.find(":");
		if(pos == std::string::npos)
			throw std::runtime_error("Invalid numeric range '" + str + "'.");
		size_t pos2 = str.substr(pos + 1).find(":");
		if(pos2 == std::string::npos) {
			// first:last
			std::string firstStr = str.substr(0, pos);
			std::string lastStr = str.substr(pos + 1);
			Value<Derived> first = eval(firstStr);
			Value<Derived> last = eval(lastStr);
			if(first.matrix().size() != 1 || last.matrix().size() != 1 || first.matrix()(0, 0) > last.matrix()(0, 0))
				throw std::runtime_error("Invalid numeric range '" + str + "'.");
			int n = 1 + floor(last.matrix()(0, 0) - first.matrix()(0, 0));
			mat.local().resize(1, n);
			for(int i = 0; i < n; i++)
				mat.local()(0, i) = first.matrix()(0, 0) + i;
			mat.mapLocal();
		} else {
			// first:step:last
			pos2 += pos + 1;
			std::string firstStr = str.substr(0, pos);
			std::string stepStr = str.substr(pos + 1, pos2 - pos - 1);
			std::string lastStr = str.substr(pos2 + 1);
			Value<Derived> first = eval(firstStr);
			Value<Derived> step = eval(stepStr);
			Value<Derived> last = eval(lastStr);
			if(first.matrix().size() != 1 || step.matrix().size() != 1 || last.matrix().size() != 1)
				throw std::runtime_error("Invalid numeric range '" + str + "'.");
			if(first.matrix()(0, 0) == last.matrix()(0, 0)) {
				mat = first.matrix()(0, 0);
			} else if(first.matrix()(0, 0) < last.matrix()(0, 0) && step.matrix()(0, 0) > 0) {
				int n = 1 + floor((last.matrix()(0, 0) - first.matrix()(0, 0)) / step.matrix()(0, 0));
				mat.local().resize(1, n);
				for(int i = 0; i < n; i++)
					mat.local()(0, i) = first.matrix()(0, 0) + i * step.matrix()(0, 0);
				mat.mapLocal();
			} else if(first.matrix()(0, 0) > last.matrix()(0, 0) && step.matrix()(0, 0) < 0) {
				int n = 1 + floor((last.matrix()(0, 0) - first.matrix()(0, 0)) / step.matrix()(0, 0));
				mat.local().resize(1, n);
				for(int i = 0; i < n; i++)
					mat.local()(0, i) = first.matrix()(0, 0) + i * step.matrix()(0, 0);
				mat.mapLocal();
			} else {
				throw std::runtime_error("Invalid numeric range '" + str + "'.");
			}
		}
	}
	
	template <typename Derived>
	bool Parser<Derived>::isOperator(const std::string & str) const
	{
		if(str.size() == 1)
			return isOperator(str[0]);
		else if(str.size() == 2) {
			size_t pos = mOperators2.find(str);
			return (pos != std::string::npos && pos % 2 == 0);
		}
		return false;
	}
	
	template <typename Derived>
	void Parser<Derived>::evalIndices(ChunkArray & chunks)
	{
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		bool operationPerformed = false;
#endif
#endif
		for(typename ChunkArray::iterator it = chunks.begin(); it != chunks.end(); it++) {
			if(it->row0 != -1 && (it->type == VALUE || (it->type == VARIABLE && it->row0 != -1 && (it + 1 == chunks.end() || (it + 1)->type != OPERATOR || (it + 1)->field != "=")))) {
				if(it->type == VALUE) {
					Derived temp = it->value.local().block(it->row0, it->col0, it->rows, it->cols);
					it->value.local() = temp;
					it->value.mapLocal();
				} else { //if(it->type == VARIABLE) {
					if(!isVariable(it->field))
						throw std::runtime_error("Attempted indexing into uninitialized variable '" + it->field + "'.");
					it->value.local() = mVariables[it->field].matrix().block(it->row0, it->col0, it->rows, it->cols);
					it->value.mapLocal();
					it->type = VALUE;
				}
				it->row0 = -1;
				it->col0 = -1;
				it->rows = -1;
				it->cols = -1;
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
				operationPerformed = true;
#endif
#endif
			}
		}
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		if(operationPerformed) { std::cout << "i: "; printChunks(chunks); std::cout << std::endl; }
#endif
#endif
	}
	
	template <typename Derived>
	void Parser<Derived>::evalNegations(ChunkArray & chunks)
	{
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		bool operationPerformed = false;
#endif
#endif
		if(chunks.size() < 2) return;
		typename ChunkArray::iterator lhs = chunks.begin(), op = chunks.begin(), rhs = op + 1;
		for(; lhs != chunks.end() && op != chunks.end() && rhs != chunks.end();) {
			if(op->type == OPERATOR && op->field == "-" && (op == chunks.begin() || (lhs->type != VALUE && lhs->type != VARIABLE)) && (rhs->type == VALUE || rhs->type == VARIABLE)) {
				if(rhs->type == VALUE)
					rhs->value.matrix().array() *= -1;
				else if(rhs->type == VARIABLE) {
					if(!isVariable(rhs->field))
						throw std::runtime_error("Attempted operation '" + op->field + rhs->field + "' on uninitialized variable '" + rhs->field + "'.");
					rhs->value.local() = mVariables[rhs->field].matrix().array() * -1;
					rhs->value.mapLocal();
					rhs->type = VALUE;
				}
				chunks.erase(op);
				op = rhs + 1;
				lhs = op - 1;
				rhs = op + 1;
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
				operationPerformed = true;
#endif
#endif
			} else {
				lhs = op;
				op = rhs;
				rhs++;
			}
		}
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		if(operationPerformed) { std::cout << "-: "; printChunks(chunks); std::cout << std::endl; }
#endif
#endif
	}

	template <typename Derived>
	void Parser<Derived>::evalPowers(ChunkArray & chunks)
	{
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		bool operationPerformed = false;
#endif
#endif
		if(chunks.size() < 3) return;
		typename ChunkArray::iterator lhs = chunks.begin(), op = lhs + 1, rhs = op + 1;
		for(; lhs != chunks.end() && op != chunks.end() && rhs != chunks.end();) {
			if((op->type == OPERATOR) && (op->field == "^" || op->field == ".^")) {
				if(lhs->type == VARIABLE) {
					if(!isVariable(lhs->field))
						throw std::runtime_error("Attempted operation '" + lhs->field + op->field + rhs->field + "' on uninitialized variable '" + lhs->field + "'.");
					lhs->value.setShared(mVariables[lhs->field]);
				}
				if(rhs->type == VARIABLE) {
					if(!isVariable(rhs->field))
						throw std::runtime_error("Attempted operation '" + lhs->field + op->field + rhs->field + "' on uninitialized variable '" + rhs->field + "'.");
					rhs->value.setShared(mVariables[rhs->field]);
				}
				if(rhs->value.matrix().size() == 1) {
					lhs->value.local() = lhs->value.matrix().array().pow(rhs->value.matrix()(0, 0));
					lhs->value.mapLocal();
					lhs->type = VALUE;
				} else if(lhs->value.matrix().size() == 1) {
					typename Derived::Scalar temp = lhs->value.matrix()(0, 0);
					lhs->value.local().resize(rhs->value.matrix().rows(), rhs->value.matrix().cols());
					for(size_t row = 0; row < rhs->value.matrix().rows(); row++) {
						for(size_t col = 0; col < rhs->value.matrix().cols(); col++)
							lhs->value.local()(row, col) = pow(temp, rhs->value.matrix()(row, col));
					}
					lhs->value.mapLocal();
					lhs->type = VALUE;
				} else if(op->field == ".^" && lhs->value.matrix().rows() == rhs->value.matrix().rows() && lhs->value.matrix().cols() == rhs->value.matrix().cols()) {
					lhs->value.local().resize(rhs->value.matrix().rows(), rhs->value.matrix().cols());
					for(size_t row = 0; row < rhs->value.matrix().rows(); row++) {
						for(size_t col = 0; col < rhs->value.matrix().cols(); col++)
							lhs->value.local()(row, col) = pow(lhs->value.matrix()(row, col), rhs->value.matrix()(row, col));
					}
					lhs->value.mapLocal();
					lhs->type = VALUE;
				} else {
					throw std::runtime_error("Invalid operand dimensions for operation '" + lhs->field + op->field + rhs->field + "'.");
				}
				chunks.erase(op, rhs + 1);
				op = lhs + 1;
				rhs = op + 1;
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
				operationPerformed = true;
#endif
#endif
			} else {
				lhs = op;
				op = rhs;
				rhs++;
			}
		}
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		if(operationPerformed) { std::cout << "^: "; printChunks(chunks); std::cout << std::endl; }
#endif
#endif
	}

	template <typename Derived>
	void Parser<Derived>::evalMultiplication(ChunkArray & chunks)
	{
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		bool operationPerformed = false;
#endif
#endif
		if(chunks.size() < 3) return;
		typename ChunkArray::iterator lhs = chunks.begin(), op = lhs + 1, rhs = op + 1;
		for(; lhs != chunks.end() && op != chunks.end() && rhs != chunks.end();) {
			if((op->type == OPERATOR) && (op->field == "*" || op->field == "/" || op->field == ".*" || op->field == "./")) {
				if(lhs->type == VARIABLE) {
					if(!isVariable(lhs->field))
						throw std::runtime_error("Attempted operation '" + lhs->field + op->field + rhs->field + "' on uninitialized variable '" + lhs->field + "'.");
					lhs->value.setShared(mVariables[lhs->field]);
				}
				if(rhs->type == VARIABLE) {
					if(!isVariable(rhs->field))
						throw std::runtime_error("Attempted operation '" + lhs->field + op->field + rhs->field + "' on uninitialized variable '" + rhs->field + "'.");
					rhs->value.setShared(mVariables[rhs->field]);
				}
				if(rhs->value.matrix().size() == 1) {
					if(lhs->value.isLocal()) {
						if(op->field == "*" || op->field == ".*")
							lhs->value.local().array() *= rhs->value.matrix()(0, 0);
						else // if(op->field == "/" || op->field == "./")
							lhs->value.local().array() /= rhs->value.matrix()(0, 0);
					} else {
						if(op->field == "*" || op->field == ".*")
							lhs->value.local() = lhs->value.matrix().array() * rhs->value.matrix()(0, 0);
						else // if(op->field == "/" || op->field == "./")
							lhs->value.local() = lhs->value.matrix().array() / rhs->value.matrix()(0, 0);
						lhs->value.mapLocal();
						lhs->type = VALUE;
					}
				} else if(lhs->value.matrix().size() == 1) {
					typename Derived::Scalar temp = lhs->value.matrix()(0, 0);
					if(op->field == "*" || op->field == ".*")
						lhs->value.local() = rhs->value.matrix().array() * temp;
					else // if(op->field == "/" || op->field == "./")
						lhs->value.local() = Derived::Constant(rhs->value.matrix().rows(), rhs->value.matrix().cols(), temp).array() / rhs->value.matrix().array();
					lhs->value.mapLocal();
					lhs->type = VALUE;
				} else if((op->field == ".*" || op->field == "./") && lhs->value.matrix().rows() == rhs->value.matrix().rows() && lhs->value.matrix().cols() == rhs->value.matrix().cols()) {
					if(lhs->value.isLocal()) {
						if(op->field == ".*")
							lhs->value.local().array() *= rhs->value.matrix().array();
						else // if(op->field == "./")
							lhs->value.local().array() /= rhs->value.matrix().array();
					} else {
						if(op->field == ".*")
							lhs->value.local() = lhs->value.matrix().array() * rhs->value.matrix().array();
						else // if(op->field == "./")
							lhs->value.local() = lhs->value.matrix().array() / rhs->value.matrix().array();
						lhs->value.mapLocal();
						lhs->type = VALUE;
					}
				} else if(op->field == "*" && lhs->value.matrix().cols() == rhs->value.matrix().rows()) {
					if(lhs->value.isLocal()) {
						lhs->value.local() *= rhs->value.matrix();
						lhs->value.mapLocal();
					} else {
						lhs->value.local() = lhs->value.matrix() * rhs->value.matrix();
						lhs->value.mapLocal();
						lhs->type = VALUE;
					}
				} else {
					throw std::runtime_error("Invalid operand dimensions for operation '" + lhs->field + op->field + rhs->field + "'.");
				}
				chunks.erase(op, rhs + 1);
				op = lhs + 1;
				rhs = op + 1;
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
				operationPerformed = true;
#endif
#endif
			} else {
				lhs = op;
				op = rhs;
				rhs++;
			}
		}
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		if(operationPerformed) { std::cout << "*: "; printChunks(chunks); std::cout << std::endl; }
#endif
#endif
	}

	template <typename Derived>
	void Parser<Derived>::evalAddition(ChunkArray & chunks)
	{
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		bool operationPerformed = false;
#endif
#endif
		if(chunks.size() < 3) return;
		typename ChunkArray::iterator lhs = chunks.begin(), op = lhs + 1, rhs = op + 1;
		for(; lhs != chunks.end() && op != chunks.end() && rhs != chunks.end();) {
			if((op->type == OPERATOR) && (op->field == "+" || op->field == "-" || op->field == ".+" || op->field == ".-")) {
				if(lhs->type == VARIABLE) {
					if(!isVariable(lhs->field))
						throw std::runtime_error("Attempted operation '" + lhs->field + op->field + rhs->field + "' on uninitialized variable '" + lhs->field + "'.");
					lhs->value.setShared(mVariables[lhs->field]);
				}
				if(rhs->type == VARIABLE) {
					if(!isVariable(rhs->field))
						throw std::runtime_error("Attempted operation '" + lhs->field + op->field + rhs->field + "' on uninitialized variable '" + rhs->field + "'.");
					rhs->value.setShared(mVariables[rhs->field]);
				}
				if(rhs->value.matrix().size() == 1) {
					if(lhs->value.isLocal()) {
						if(op->field == "+" || op->field == ".+")
							lhs->value.local().array() += rhs->value.matrix()(0, 0);
						else // if(op->field == "-" || op->field == ".-")
							lhs->value.local().array() -= rhs->value.matrix()(0, 0);
					} else {
						if(op->field == "+" || op->field == ".+")
							lhs->value.local() = lhs->value.matrix().array() + rhs->value.matrix()(0, 0);
						else // if(op->field == "-" || op->field == ".-")
							lhs->value.local() = lhs->value.matrix().array() - rhs->value.matrix()(0, 0);
						lhs->value.mapLocal();
						lhs->type = VALUE;
					}
				} else if(lhs->value.matrix().size() == 1) {
					typename Derived::Scalar temp = lhs->value.matrix()(0, 0);
					if(op->field == "+" || op->field == ".+")
						lhs->value.local() = rhs->value.matrix().array() + temp;
					else // if(op->field == "-" || op->field == ".-")
						lhs->value.local() = Derived::Constant(rhs->value.matrix().rows(), rhs->value.matrix().cols(), temp).array() - rhs->value.matrix().array();
					lhs->value.mapLocal();
					lhs->type = VALUE;
				} else if(lhs->value.matrix().rows() == rhs->value.matrix().rows() && lhs->value.matrix().cols() == rhs->value.matrix().cols()) {
					if(lhs->value.isLocal()) {
						if(op->field == "+" || op->field == ".+")
							lhs->value.local().array() += rhs->value.matrix().array();
						else // if(op->field == "-" || op->field == ".-")
							lhs->value.local().array() -= rhs->value.matrix().array();
					} else {
						if(op->field == "+" || op->field == ".+")
							lhs->value.local() = lhs->value.matrix().array() + rhs->value.matrix().array();
						else // if(op->field == "-" || op->field == ".-")
							lhs->value.local() = lhs->value.matrix().array() - rhs->value.matrix().array();
						lhs->value.mapLocal();
						lhs->type = VALUE;
					}
				} else {
					throw std::runtime_error("Invalid operand dimensions for operation '" + lhs->field + op->field + rhs->field + "'.");
				}
				chunks.erase(op, rhs + 1);
				op = lhs + 1;
				rhs = op + 1;
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
				operationPerformed = true;
#endif
#endif
			} else {
				lhs = op;
				op = rhs;
				rhs++;
			}
		}
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		if(operationPerformed) { std::cout << "+: "; printChunks(chunks); std::cout << std::endl; }
#endif
#endif
	}

	template <typename Derived>
	void Parser<Derived>::evalAssignment(ChunkArray & chunks)
	{
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		bool operationPerformed = false;
#endif
#endif
		if(chunks.size() < 3) return;
		typename ChunkArray::iterator rhs = chunks.end() - 1, op = rhs - 1, lhs = op - 1;
		for(; op != chunks.begin() && rhs != chunks.begin();) {
			if(op->type == OPERATOR && op->field == "=" && (lhs->type == VALUE || lhs->type == VARIABLE) && (rhs->type == VALUE || rhs->type == VARIABLE)) {
				if(rhs->type == VARIABLE) {
					if(!isVariable(rhs->field))
						throw std::runtime_error("Attempted operation '" + lhs->field + op->field + rhs->field + "' on uninitialized variable '" + rhs->field + "'.");
					rhs->value.setShared(mVariables[rhs->field]);
				}
				if(lhs->type == VALUE) {
					lhs->value.local() = rhs->value.matrix();
					lhs->value.mapLocal();
				} else { //if(lhs->type == VARIABLE) {
					if(isVariable(lhs->field)) {
						lhs->value.setShared(mVariables[lhs->field]);
						if(lhs->row0 == -1) {
							if(lhs->value.matrix().rows() == rhs->value.matrix().rows() && lhs->value.matrix().cols() == rhs->value.matrix().cols()) {
								lhs->value.matrix() = rhs->value.matrix();
							} else {
								mVariables[lhs->field].local() = rhs->value.matrix();
								mVariables[lhs->field].mapLocal();
							}
						} else { //if(lhs->row0 != -1) {
							lhs->value.matrix().block(lhs->row0, lhs->col0, lhs->rows, lhs->cols) = rhs->value.matrix();
						}
					} else {
						mVariables[lhs->field].local() = rhs->value.matrix();
						mVariables[lhs->field].mapLocal();
					}
				}
				chunks.erase(op, rhs + 1);
				rhs = lhs;
				op = rhs - 1;
				lhs = op - 1;
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
				operationPerformed = true;
#endif
#endif
			} else {
				rhs = op;
				op = lhs;
				lhs--;
			}
		}
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
		if(operationPerformed) { std::cout << "=: "; printChunks(chunks); std::cout << std::endl; }
#endif
#endif
	}
	
#ifdef DEBUG
#ifdef EIGENLAB_DEBUG
	template <typename Derived>
	void Parser<Derived>::printChunks(ChunkArray & chunks, size_t maxRows, size_t maxCols, int precision)
	{
		std::cout << "__";
		for(typename ChunkArray::iterator it = chunks.begin(); it != chunks.end(); it++) {
			switch(it->type) {
				case VALUE:
					std::cout << textRepresentation(it->value, maxRows, maxCols, precision);
					if(it->row0 != -1)
						std::cout << "(" << it->row0 << ":" << it->row0 + it->rows - 1 << "," << it->col0 << ":" << it->col0 + it->cols - 1 << ")";
					break;
				case VARIABLE:
					std::cout << it->field;// << "=" << textRepresentation(mVariables[it->field], maxRows, maxCols, precision);
					if(it->row0 != -1)
						std::cout << "(" << it->row0 << ":" << it->row0 + it->rows - 1 << "," << it->col0 << ":" << it->col0 + it->cols - 1 << ")";
					break;
				case OPERATOR:
					std::cout << it->field;
					break;
				case FUNCTION:
					std::cout << "f()=" << it->field;
					break;
			}
			std::cout << "__";
		}
	}
	
	template <typename Derived>
	void Parser<Derived>::printVars(size_t maxRows, size_t maxCols, int precision)
	{
		for(typename ValueMap::iterator it = mVariables.begin(); it != mVariables.end(); it++)
			std::cout << it->first << " (" << it->second.matrix().rows() << "x" << it->second.matrix().cols() << ") = " << textRepresentation(it->second, maxRows, maxCols, precision) << std::endl;
	}
	
	template <typename Derived>
	std::string Parser<Derived>::textRepresentation(Value<Derived> & val, size_t maxRows, size_t maxCols, int precision)
	{
		if(val.matrix().size() == 1)
			return numberToString<typename Derived::Scalar>(val.matrix()(0, 0), precision);
		else {
			std::string str = "[";
			for(size_t row = 0; row < val.matrix().rows() && row < maxRows; row++) {
				str += (row > 0 ? ";[" : "[");
				for(size_t col = 0; col < val.matrix().cols() && col < maxCols; col++) {
					if(col > 0) str += ",";
					str += numberToString<typename Derived::Scalar>(val.matrix()(row, col), precision);
				}
				str += (val.matrix().cols() > maxCols ? "...]" : "]");
			}
			str += (val.matrix().rows() > maxRows ? "...]" : "]");
			return str;
		}
	}
#endif // #ifdef EIGENLAB_DEBUG
#endif // #ifdef DEBUG
	
	template <typename Derived>
	std::string Parser<Derived>::trim(const std::string & str)
	{
		if(str.size() == 0) return str;
		std::string::const_iterator first = str.begin();
		std::string::const_iterator last = str.end() - 1;
		while(first < last && isspace(*first)) first++;
		while(first < last && isspace(*last)) last--;
		return std::string(first, last + 1);
	}
	
	template <typename Derived>
	std::vector<std::string> Parser<Derived>::split(const std::string & str, const char delimeter)
	{
		std::vector<std::string> args;
		std::string::const_iterator i0 = str.begin();
		for(std::string::const_iterator it = str.begin(); it != str.end(); it++) {
			if((*it) == delimeter) {
				args.push_back(trim(std::string(i0, it)));
				i0 = it + 1;
			}
		}
		args.push_back(std::string(i0, str.end()));
		return args;
	}
	
	template <typename Derived>
	template <typename T>
	bool Parser<Derived>::isNumber(const std::string & str, T * num)
	{
		std::istringstream iss(str);
		if(num)
			iss >> (* num);
		else {
			T number;
			iss >> number;
		}
		return (!iss.fail() && !iss.bad() && iss.eof());
	}
	
	template <typename Derived>
	template<typename T>
	T Parser<Derived>::stringToNumber(const std::string & str)
	{
		std::istringstream iss(str);
		T number;
		iss >> number;
		if(iss.fail() || iss.bad() || !iss.eof())
			throw std::runtime_error("Failed to convert " + str + " to a number.");
		return number;
	}
	
	template <typename Derived>
	template<typename T>
	std::string Parser<Derived>::numberToString(T num, int precision)
	{
		std::ostringstream oss;
		if(precision)
			oss << std::setprecision(precision) << num;
		else
			oss << num;
		return oss.str();
	}
} // namespace EigenLab

#endif // #ifndef EigenLab_H
