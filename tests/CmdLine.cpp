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
