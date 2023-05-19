#include "HermiteVec.h"
#include <iostream>
#include <cmath>
#include <math.h>

int main(int argc, char* argv[])
{
  std::vector<double> time{0, 3, 5, 8, 13};
  std::vector<double> dist{0, 225, 383, 623, 993};
  std::vector<double> speed{75, 77, 80, 74, 72};
  InterpConditions data(time, dist, speed);

  NewtonInterp result(data); 
  Polynomial positionpoly(result.NewtonFormula());
  Polynomial speedpoly(positionpoly.derivative());
  std::cout << "Position polynomial : " << std::endl;
  positionpoly.print();
  std::cout << std::endl << std::endl;
  std::cout << "Speed polynomial : " << std::endl;
  speedpoly.print();
  std::cout << std::endl << std::endl;
  
  std::cout << "Position when t = 10s : " << positionpoly.eval(10) << std::endl;
  std::cout << "Speed when t = 10s : " << speedpoly.eval(10) << std::endl;
  std::cout << std::endl;
  for(int i = 0; i<=13; i++)
  {
    std::cout << speedpoly.eval(i) << " ";
  }
  std::cout << std::endl;
  return 0;
}
