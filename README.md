EigenLab
===
[![Build Status](https://travis-ci.org/marcel-goldschen-ohm/EigenLab.svg?branch=master)](https://travis-ci.org/marcel-goldschen-ohm/EigenLab)

C++ header only library for parsing/evaluating matrix math equations which are unknown until run time (e.g. user input) in a format similar to [MATLAB](http://www.mathworks.com/products/matlab). Variables can be defined via the parsed equations as in MATLAB, or alternatively mapped to preallocated data. Matrix math is performed using [Eigen](http://eigen.tuxfamily.org). Supports basic matrix operations, matrix reductions, submatrix indexing as in MATLAB (except indices are 0 based), coefficient-wise operations and coefficient-wise function evaluation. This allows the user to evaluate matrix math equations in an interactive fashion like in MATLAB, or, for example, to fit a data set using an arbitrary function defined by the user at run time. Finally, EigenLab's interface is incredibly simple and easy to use. Check out the Quick Start Examples seciton and see below for a taste of what EigenLab is capable of.

```cpp
EigenLab::ParserXd parser;
EigenLab::ValueXd result;
Eigen::MatrixXd myData(2,2);
myData << 1,2,3,4;
    
parser.var("x").setShared(myData);
result = parser.eval("x * [[1,2];[1,2]]");
std::cout << result.matrix() << std::endl;
    
parser.eval("y = ([0:0.5:10] .* sqrt(4)) + 100");
std::cout << parser.var("y").matrix() << std::endl;
    
result = parser.eval("x(:,1)");
std::cout << result.matrix() << std::endl;
    
parser.eval("x = x .* 2");
std::cout << myData << std::endl;
```

Author: Dr. Marcel Paz Goldschen-Ohm  
Email:  <marcel.goldschen@gmail.com>   
License: MIT  
Copyright (c) 2014 Dr. Marcel Paz Goldschen-Ohm

### Why should I use EigenLab?
When you need to evaluate a previously unknown matrix math equation that is input by the user at run time, and optionally contains variable names referencing preallocated data. Also because EigenLab has a simple and intuitive interface (similar to MATLAB), and is entirely contained within a single header file.

### When NOT to use EigenLab?
When you are evaluating a math equation that is known at compile time.

### How efficient is EigenLab?
EigenLab uses [Eigen](http://eigen.tuxfamily.org) for all matrix operations, and therefore takes advantage of Eigen's incredible matrix math performance. One caveate is that EigenLab evaluates each individual operation within an equation sequentially, and therefore does NOT take advantage of Eigen's expression templates, which results in some wasteful temporary copies of matrix data.

I have not done an extensive test of performance in EigenLab as compared to other parseres such as those mentioned below. However, a simple comparison between EigenLab and [muParser](http://muparser.beltoforion.de) is included in the file `SpeedTestEigenLabVsMuParser.cpp` that will probably give you the basic idea. *To summarize, if you can formulate your expression to utilize matrix/array math, then EigenLab is incredibly fast as it can take advantage of Eigen's efficeint matrix math manipulations. On the other hand, if you are forced to reevaluate an expression many times on successive scalar values (e.g. some sort of series expansion), then another parser such as muParser is likely your best choice in terms of speed.*

### Are there alternatives to EigenLab?
To my knowledge, there are no native C/C++ alternatives to EigenLab which can both parse matrix math equations and provide an incredibly simple and user friendly interface (i.e. like MATLAB).

*Update:* The [ExprTk](http://www.partow.net/programming/exprtk/index.html) project looks interesting and seems like it can do some vector expression operations, as well as a bunch of stuff that EigenLab doesn't even try to do. However, I've not used it in any depth.

If you don't need matrix operations, there are several equation parsers which can evaluate math equations for floating point numbers, such as [GNU libmatheval](http://www.gnu.org/software/libmatheval), [muParser](http://muparser.beltoforion.de) and [MathPresso/DoublePresso](https://github.com/tartakynov/mathpresso). As I have not spent the time to really dig through these libraries, please feel free to correct me if I am wrong about this, or to inform me of another library not listed here.

Finally, it is of course possible to embed in your C++ program an interpreted language such as Python or Lua that could handle any sort of expression evaluation you want.

### INSTALL
---
All of EigenLab is in the single header file `EigenLab.h`. Just include it.

EigenLab requires the header only C++ library [Eigen](http://eigen.tuxfamily.org) (tested with version 3.2.2) for matrix math. Either make sure Eigen is in your searched include paths, or change the include statement in `EigenLab.h` from `#include <Eigen/Dense>` to `#include <YOUR/PATH/TO/Eigen/Dense>`.

### LICENSE
---
The MIT License (MIT)

Copyright (c) 2014 Dr. Marcel Paz Goldschen-Ohm

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

### Quick Start Examples
---
Here's a series of simple examples of how to use EigenLab. For most applications, this is probably all you'll need to get up and running with EigenLab.

```cpp
#include "EigenLab.h"
#include <iostream>

#define PRINT_RESULT \
  std::cout << expr << " = " << resultXd.matrix() << std::endl;

int main(int argc, char** argv)
{
  // Create a parser and a value to hold the
  // result of evaluating an expression.
  EigenLab::Parser<Eigen::MatrixXd> parserXd;
  EigenLab::Value<Eigen::MatrixXd> resultXd;
  std::string expr;
  
  // Evaluate a numeric expression.
  // Note that all values are treated
  // as matrices in EigenLab.
  // e.g. a scalar is just a 1x1 matrix.
  try {
    expr = "2.5 * sqrt((3 - 1)^2)";
    resultXd = parserXd.eval(expr);
    PRINT_RESULT
  } catch(std::exception &e) {
    std::cout << "ERROR : " << e.what() << std::endl << std::endl;
  }
  
  // From here on I'll leave out the
  // try/catch block for readability.
  
  // Evaluate a matrix expression.
  expr = "[[1,2];[3,4]]";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 1 2
  // 3 4
  
  expr = "[1;2]";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 1
  // 2
  
  expr = "[[1,2];[3,4]] * [1;2]";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 5
  // 11
  
  // Evaluate a matrix reduction.
  expr = "trace([[1,2];[3,4]])";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 5
  
  expr = "sum([[1,2];[3,4]])";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 10
         
  // Evaluate a matrix partial reduction.
  // Note that indices are 0-based
  // whereas MATLAB uses 1-based indices.
  expr = "sum([[1,2];[3,4]], 0)";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 4 6
  
  expr = "sum([[1,2];[3,4]], 1)";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 3
  // 7
  
  // Matrix coefficients and submatrices.
  // Note that indices are 0-based
  // whereas MATLAB uses 1-based indices.
  expr = "[[1,2];[3,4]](1,1)";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 4
  
  expr = "[[1,2];[3,4]](:,0)";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 1
  // 3
  
  expr = "[[1,2,3,4];[5,6,7,8]](1,0:2)";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 5 6 7
  
  // Create a variable and then use it.
  parserXd.eval("x = 2");
  std::cout << "x = " << parserXd.var("x").matrix() << std::endl;
  // x = 2
  
  expr = "[1:2:9]";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 1 3 5 7 9
  
  expr = "[1:2:9] * x";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 2 6 10 14 18
  
  // Assign variables to preallocated data
  // and then use them.
  Eigen::MatrixXd A = Eigen::MatrixXd::Ones(2, 2);
  double B[4] = {1, 2, 3, 4};
  double c = 2.5;
  parserXd.var("A").setShared(A);
  parserXd.var("B").setShared(B, 2, 2);
  parserXd.var("c").setShared(& c);
  
  std::cout << "A = " << A << std::endl;
  // A = 1 1
  //     1 1
  
  std::cout << "B = " << parserXd.var("B").matrix() << std::endl;
  // B = 1 3
  //     2 4
  
  std::cout << "c = " << c << std::endl;
  // c = 2.5
  
  expr = "A .* B";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 1 3
  // 2 4
  
  expr = "A * B";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 3 7
  // 3 7
  
  expr = "A * c";
  resultXd = parserXd.eval(expr);
  PRINT_RESULT
  // 2.5 2.5
  // 2.5 2.5
  
  // Assign to a variable (!!! NO resize required !!!).
  // If a variable is a map to shared data,
  // and assigning to it does not require any resizing,
  // then assign directly to the mapped shared data.
  parserXd.eval("A = B");
  std::cout << "A = " << A << std::endl;
  // A = 1 3
  //     2 4
  
  parserXd.eval("c = 8.2");
  std::cout << "c = " << c << std::endl;
  // c = 8.2
  
  // Assign to a variable (!!! resize required !!!).
  // If a variable is a map to shared data,
  // and assigning to it requires resizing it,
  // then do NOT assign directly to the mapped
  // shared data, but instead assign to a new
  // local variable having the same name.
  parserXd.eval("A = [5;10]");
  std::cout << "A = " << A << std::endl;
  // A = 1 3
  //     2 4
  std::cout << parserXd.var("A").matrix() << std::endl;
  // 5
  // 10
  
  return 0;
};
```

### Example EigenLab command shell similar to MATLAB.
---
Here's an interactive command line application that behaves similar to MATLAB's command window.

```cpp
#include "EigenLab.h"
#include <iostream>

int main(int argc, char** argv)
{
  EigenLab::ParserXd parserXd;
  typedef
    std::map<std::string, EigenLab::ValueXd>
    VariableMap;
  std::string input;
  while(1) {
    // Prompt user for input...
    std::cout << ">> ";
    getline(std::cin, input);
    if(input == "quit") {
      // Quit the program.
      break;
    } else if(input == "whos") {
      // Print variable names and their matrix dimensions.
      for(VariableMap::iterator it = parserXd.vars().begin();
          it != parserXd.vars().end();
          it++) {
        std::string varName = it->first;
        EigenLab::ValueXd varValue = it->second;
        std::cout << varName << " ("
          << varValue.matrix().rows() << "x"
          << varValue.matrix().cols() << ")" << std::endl;
      }
    } else if(input.find("clear ") == 0) {
      // Delete list of comma separated variables.
      std::vector<std::string> args =
        EigenLab::ParserXd::split(input.substr(6), ',');
      for(size_t i = 0; i < args.size(); i++)
        parserXd.clearVar(args[i]);
    } else {
      // Evaluate the input expression.
      try {
        if(input.back() == ';') {
          // Evaluate the expression, but don't print the result.
          parserXd.eval(input.substr(0, int(input.size()) - 1));
        } else {
          // Evaluate the expression and print the result.
          std::cout << parserXd.eval(input).matrix() << std::endl;
        }
      } catch(std::exception & e) {
        // Print the error message.
        std::cout << e.what() << std::endl;
      }
    }
  }
  return 0;
}
```

Users Guide
===

### Namespace.
---
Everything is in the namespace `EigenLab`.

### The Value class.
---
EigenLab treats all values as matrices (i.e. a scalar is a 1x1 matrix). The `Value<Derived>` class template is a wrapper that provides a consistent interface to accessing matrix data whether that data is stored locally or points to some shared data that is stored elsewhere (i.e. variables). The template parameter `Derived` denotes the value's underlying matrix data type, which can be any dynamically sized matrix type supported by Eigen (e.g. `Eigen::MatrixXd`).

```cpp
template <typename Derived> Value;
    
typedef Value<Eigen::MatrixXd> ValueXd;
typedef Value<Eigen::MatrixXf> ValueXf;
typedef Value<Eigen::MatrixXi> ValueXi;
```

Under the hood, there are just two data members. One is a local matrix object of type `Derived`, and the other is a map to some shared matrix data with type `Eigen::Map<Derived>`. If we are useing the local matrix data, the map will point to it, so the map is always valid, whereas the local matrix data is only valid if we are currently using it.

The `matrix()` method returns a reference to an `Eigen::Map<Derived>` object which maps to the actual matrix data regardless of where it is stored (locally or shared). **!!! Note that it is up to you to make sure that any shared data being mapped in a Value instance remains valid.**

```cpp
Eigen::Map<Derived> & matrix();
```

The `setLocal(...)` method stores a local copy of the assigned data in the value, whereas `setShared(...)` sets the value's map to point to the shared data being assigned. When assigning another value using the `=` operator, the behavior depends on whether the assigned value holds local data or a map to shared data.

```cpp
// Example Value assignments.
    
EigenLab::Value<Eigen::MatrixXd> valueXd;
EigenLab::Value<Eigen::MatrixXd> anotherValueXd;
double num = 3.5;
double mat3[3] = {0.0, 1.0, 2.0};
Eigen::MatrixXd mat2x2 = Eigen::MatrixXd::Ones(2, 2);
    
// valueXd contains a local copy
// of the 1x1 matrix num.
valueXd.setLocal(num);
    
// valueXd contains a map to
// the shared 1x1 matrix at & num.
valueXd.setShared(& num);
    
// valueXd contains a local copy
// of the 2x2 matrix mat2x2.
valueXd.setLocal(mat2x2);
    
// valueXd contains a map to
// the shared 2x2 matrix at mat2x2.data().
valueXd.setShared(mat2x2);
    
// valueXd contains a map to
// the shared 1x3 matrix at mat3.
valueXd.setShared(mat3, 1, 3);
    
// valueXd contains a map to
// the shared 3x1 matrix at mat3.
valueXd.setShared(mat3, 3, 1);
    
// valueXd contains a local copy
// of the 3x1 matrix at mat3.
valueXd.setLocal(mat3, 3, 1);
    
// anotherValueXd contains a local copy
// of the 2x2 matrix mat2x2.
valueXd.setLocal(mat2x2);
anotherValueXd = valueXd;
    
// anotherValueXd contains a map to
// the shared 2x2 matrix at mat2x2.data().
valueXd.setShared(mat2x2);
anotherValueXd = valueXd;
```

### The Parser class.
---
To evaluate a matrix math expression, just create a `Parser<Derived>` instance and pass the equation string to its `eval()` method, which will return the result of the evaluation as a `Value<Derived>` object. The template parameter `Derived` denotes the matrix data type used for each value in the expression, which can be any dynamically sized matrix type supported by Eigen (e.g. `Eigen::MatrixXd`).

```cpp
template <typename Derived> Parser;
    
typedef Parser<Eigen::MatrixXd> ParserXd;
typedef Parser<Eigen::MatrixXf> ParserXf;
typedef Parser<Eigen::MatrixXi> ParserXi;
    
Value<Derived> eval(const std::string & equation);
```

### Variables.
---
Each `Parser<Derived>` instance contains its own set of variables. Each variable is stored in a `Value<Derived>` wrapper just like any other value, and thus can be either a locally stored variable or a map to some shared data. A reference to a variable's value is retrieved via the `var(name)` method.  *!!! Note that calling `var(name)` will create the variable 'name' if it doesn't already exist. Thus, you can first check if the variable exists using the `hasVar(name)` method.* You can delete a variable using the `clearVar(name)` method.

```cpp
Value<Derived> & var(const std::string name);
bool hasVar(const std::string name);
void clearVar(const std::string name);
```

Because variables are just `Value<Derived>` objects, we can assign to them in the same way as described above for the Value class.

```cpp
// Example variable assignments.
    
EigenLab::Parser<Eigen::MatrixXd> parserXd;
double num = 3.5;
double mat3[3] = {0.0, 1.0, 2.0};
Eigen::MatrixXd mat2x2 = Eigen::MatrixXd::Ones(2, 2);
    
// parserXd.var("x") contains a local copy
// of the 1x1 matrix num.
parserXd.var("x").setLocal(num);
    
// parserXd.var("x") contains a map to
// the shared 1x1 matrix at & num.
parserXd.var("x").setShared(& num);
    
// parserXd.var("x") contains a local copy
// of the 2x2 matrix mat2x2.
parserXd.var("x").setLocal(mat2x2);
    
// parserXd.var("x") contains a map to
// the shared 2x2 matrix at mat2x2.data().
parserXd.var("x").setShared(mat2x2);
    
// parserXd.var("x") contains a map to
// the shared 1x3 matrix at mat3.
parserXd.var("x").setShared(mat3, 1, 3);
    
// parserXd.var("x") contains a local copy
// of the 1x3 matrix at mat3.
parserXd.var("x").setLocal(mat3, 1, 3);
```

### Assigning to variables within parsed equations.
---
During equation evaluation, assignment to variables which contain maps to preallocated shared data behaves in one of two ways depending on whether the assignment requires a resize of the variable's matrix data. If the assigned value has the same dimensions as the variable's shared data, then the shared data will be overwritten with the assigned value. In contrast, if the assigned value does NOT have the same dimensions as the variable's shared data, then the shared data is left untouched and the variable is updated to contain a local copy of the assigned data. *!!! Note that in the latter case, the shared data no longer reflects the updated variable, which only exists locally within the parser object and can be obtained via the `var(name)` method.* Finally, if an undefined variable is assigned to, a new local variable will be created.

```cpp
EigenLab::Parser<Eigen::MatrixXd> parserXd;
EigenLab::Value<Eigen::MatrixXd> resultXd;
    
Eigen::MatrixXd A(2,2);
A << 1,2,3,4; // A = [[1,2];[3,4]]
Eigen::MatrixXd B(2,1);
B << 0,0; // B = [0;0]
parserXd.var("A").setShared(A);
parserXd.var("B").setShared(B);
    
// Assignment to shared data.
parserXd.eval("B = A * [1;2]");
std::cout << B << std::endl;
// [5;11]
	
// Assignment to new local variable
// (shared data in A is unchanged).
parserXd.eval("A = [1;2]");
std::cout << A << std::endl;
// [[1,2];[3,4]]
std::cout << parserXd.var("A").matrix() << std::endl;
// [1;2]
	
// Assignment to previously undefined variable.
parserXd.eval("x = [10,20]");
std::cout << parserXd.var("x").matrix() << std::endl;
// [10,20]
```

### Defining matrices within parsed expressions.
---
Matrices can be written out in a format similar to MATLAB. For example, `[[1,2];[3,4]]` denotes a 2x2 matrix with rows `[1,2]` and `[3,4]` and columns `[1,3]` and `[2,4]`. Brackets around individual rows are optional, so the same matrix could be written as `[1,2;3,4]`.

Ranges are defined as either `start:step:stop` or `start:stop`, where `step` is assumed to be 1 in the latter case. For example, a 3x3 matrix could be defined as `[[1,2,3];[4:6];[7:2:11]]`.

### 0-based indices, coefficeints and submatrices.
---
Matrix coefficients can be accessed as in MATLAB, **except that indices are 0-based, whereas MATLAB uses 1-based indices**. For example, `A(i,j)` denotes the coefficient with row index `i` and column index `j`. If `A` is a vector, then `A(i)` denotes the `i`'th element of `A`. 

Index ranges can be used to define submatrices as in MATLAB. For example, `A(i:r,j:s)` denotes the `p`x`q` submatrix starting at `A(i,j)`, where `p=r-i+1` and `q=s-j+1`. *!!! Note that noncontiguous submatrix blocks are not currenlty allowed.* If `A` is a vector, then `A(i:r)` denotes the length `p` subvector starting at `A(i)`, where `p=r-i+1`. As in MATLAB, an index range of `:` is interpreted as the entire row/column. Thus, `A(:,j)` denotes the `j`th column of `A`.

### Coefficient-wise operations.
---
For the following, `A` and `B` are matrices where either both have the same dimensions, or one of the two has dimension 1x1 (i.e. is a scalar). In the latter case, the operation will be performed between the scalar and each coefficient of the other matrix, with the result having the same dimensions as the matrix.

Coefficient-wise addtion: `A + B` or `A .+ B`  
Coefficient-wise subtraction: `A - B` or `A .- B`  
Coefficient-wise multiplcation: `A .* B`  
Coefficient-wise division: `A ./ B`  
Coefficient-wise power: `A .^ B`  
Coefficient-wise absolute value: `abs(A)`  
Coefficient-wise square root: `sqrt(A)`  
Coefficient-wise exponential: `exp(A)`  
Coefficient-wise natural logarithm: `log(A)`  
Coefficient-wise base 10 logarithm: `log10(A)`  
Coefficient-wise cosine: `cos(A)`  
Coefficient-wise sine: `sin(A)`  
Coefficient-wise tangent: `tan(A)`  
Coefficient-wise arc cosine: `acos(A)`  
Coefficient-wise arc sine: `asin(A)` 

### Matrix operations.
---
For the following, `A` and `B` are matrices.

Matrix multiplcation: `A * B`  
Matrix transpose: `transpose(A)`  
Matrix conjugate: `conjugate(A)`  
Matrix adjoint: `adjoint(A)`  

### Matrix reductions.
---
For the following, `A` is a matrix and the *optional* argument `dim` denotes the dimension along which to evaluate the reduction operation (0=rows, 1=cols). If `dim` is not given, the reducing function is applied to all of the matrice's coefficients.
  
Matrix norm: `norm(A)`  
Matrix trace: `trace(A)`  
Minimum coefficient: `min(A,dim)`  
Maximum coefficient: `max(A,dim)`  
Coefficient with the maximum absolute value: `absmax(A,dim)`  
Mean coefficient: `mean(A,dim)`  
Sum of coefficients: `sum(A,dim)`  
Product of coefficients: `prod(A,dim)`  

### Matrix initializers.
---
Zero matrix: `zeros(N,M)` (`zeros(N)` = `zeros(N,N)`)  
Ones matrix: `ones(N,M)` (`ones(N)` = `ones(N,N)`)  
Identity matrix: `eye(N,M)` (`eye(N)` = `eye(N,N)`)  

### Caching expressions.
---
If you need to re-evaluate the same expression multiple times, you can optionally turn on expression caching to store the parsed expression and reuse it the next time rather than reparsing it from scratch. This will speed up subsequent expression evaluations. !!! WARNING !!! If your parsed expressions contain large matrices that are stored locally within the parser (i.e. not variable names referencing shared data), then your cached expressions can eat up a lot of memory. Thus, you can turn caching on/off for each parser instance using `setCacheExpressions(bool)` and clear your current cach with `clearCachedExpressions()`. `cacheExpressions()` will return a boolean indicating your current caching status.



Version History
===

**August 26, 2014**

Initial beta release.

**October 17, 2014**

Added the `absmax()` function, which was missing in the previous version.

Fixed submatrix index ranges for variables that used either the `end` keyword or the `:` symbol to denote the entire range. These keywords were previously not interpreted correctly for variables.

Added the `test()` function for debugging. This function is meant to test all functionality in EigenLab as a quick check that everything is OK.

Thanks to Ilja Honkonen <ilja.j.honkonen@nasa.gov> for adding expression caching to speed performance when re-evaluating the same expression multiple times, and also for suggesting several type casts to remove type comparison warning messages.

Made expression caching optional in case one wants to avoid caching large matrices.

Fixed subexpressions for numeric ranges. Previously, the expression "[1 + 2 : 3 * 10 - 4]" was interpreted as 1 + [2:3] * 10 - 4. Now it is correclty evaluated as [3:26].

Added simple performance comparison with muParser (see `SpeedTestEigenLabVsMuParser.cpp`). Thanks to Ilja Honkonen <ilja.j.honkonen@nasa.gov> on which this script was based.

**August 21, 2015**

Released version 0.9.0.

**2016**

Thanks to Michael Tesch for adding some checks for handling complex matrices smoothly. Note, however, that although complex matrices can be used via shared variables, the string parser does not yet support inputing complex values.

Released version 1.0.0.

**September 17, 2016**

Added support for zeros(), ones() and eye().
