#include <iostream>
#include <math.h>
#include <vector>

int main(int argc, char*argv[])
{
  int beta = 2;
  int p = 3;
  int L = -1;
  int U = 1;
  double UFL = pow(beta,L);
  double OFL = pow(beta,U)*(beta - pow(beta, 1-p));
  std::cout << "UFL(F) = " << UFL << std::endl;
  std::cout << "OFL(F) = " << OFL << std::endl;

  std::vector<double> normal{1.00, 1.01, 1.10, 1.11};
  int count = 0;
  std::cout << "Enumerate all values of F : " << std::endl;
  for(int i = 1; i >= -1; i--)
  {
    for(int j = 3; j >= 0; j--)
    {
      std::cout << (-1)*normal[j]*pow(beta,i) << ", ";
      count++;
    }
  }
  std::cout << "0" << ", ";
  for(int i = -1; i <= 1; i++)
  {
    for(int j = 0; j <= 3; j++)
    {
      std::cout << normal[j]*pow(beta,i) << ", ";
      count++;
    }
  }
  std::cout << std::endl;
  std::cout << "Cardinality = " << count+1 << std::endl;

  std::cout << "Enumerate all subnormal numbers in F : " << std::endl;
  std::vector<double> subnormal{0.01, 0.10, 0.11};
  for(int i = 1; i >= -1; i--)
  {
    for(int j = 2; j >= 0; j--)
    {
      std::cout << (-1)*subnormal[j]*pow(beta,i) << ", ";
      count++;
    }
  }
  
  for(int i = -1; i <= 1; i++)
  {
    for(int j = 0; j <= 2; j++)
    {
      std::cout << subnormal[j]*pow(beta,i) << ", ";
      count++;
    }
  }
  std::cout << std::endl;
  return 0;
}
