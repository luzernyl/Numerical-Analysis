#include <iostream>
#include <math.h>

double fa(double x)
{
  return pow(x,8) - 8*pow(x,7) + 28*pow(x,6) - 56*pow(x,5) + 70*pow(x,4) - 56*pow(x,3) + 28*x*x - 8*x + 1;
}

double fb(double x)
{
  return (((((((x-8)*x + 28)*x - 56)*x + 70)*x -56)*x + 28)*x - 8)*x + 1;
}

double fc(double x)
{
  return pow(x-1,8);
}

int main(int argc, char* argv[])
{
  int i = 0;
  double x = 0.99;
  double delta = (1.01-0.99) / 100;
  int count = 0;
  std::cout << "2.49(a) : " << std::endl;
  while (x <= 1.01)
  {
    std::cout << fa(x) << ", ";
    x = x + delta;
  }
  std::cout << std::endl;

  x = 0.99;
  std::cout << "2.49(b) : " << std::endl;
  while (x <= 1.01)
  {
    std::cout << fb(x) << ", ";
    x = x + delta;
  }
  std::cout << std::endl;

  x = 0.99;
  std::cout << "2.49(c) : " << std::endl;
  while (x <= 1.01)
  {
    std::cout << fc(x) << ", ";
    x = x + delta;
  }
  std::cout << std::endl;
  return 0;
}
