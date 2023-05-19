#include "Approx.h"
#include <iostream>
#include <vector>

int main(int argc, char* argv[])
{
   std::cout << "Approximate for a_0 + a_1 * x + a_2 * x^2" << std::endl;
   std::vector<double> xdata;
  for(double i = 0; i <= 10; i += 0.5)
    xdata.push_back(i);

  std::vector<double> ydata {2.9, 2.7, 4.8, 5.3, 7.1, 7.6, 7.7, 7.6, 9.4, 9.0, 9.6, 10.0, 10.2, 9.7, 8.3, 8.4, 9.0, 8.3, 6.6, 6.7, 4.1};
  
  Approx data(xdata, ydata);
  data.solveNormal();

  return 0;
}
