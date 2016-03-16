// --*-  Mode: C++; c-basic-offset: 8 -*--
#include <complex>
#include "EigenLab.h"

int main(int argc, const char * argv[])
{
  EigenLab::ParserXd parserXd;
  auto failXd = parserXd.test();
  
  EigenLab::ParserXf parserXf;
  auto failXf = parserXf.test();
  
  EigenLab::ParserXi parserXi;
  auto failXi = parserXi.test();

  EigenLab::Parser<Eigen::Matrix<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic > > parserXcd;
  auto failXcd = parserXcd.test();

  std::cout << "Test summary, number of failures: Xd=" << failXd << " Xf=" << failXf << " Xi=" << failXi << " Xcd=" << failXcd << std::endl;
  return 0;
}
