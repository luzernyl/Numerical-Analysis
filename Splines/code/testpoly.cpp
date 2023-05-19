#include "Splines.h"
#include "Vector.h"

int main(int argc, char* argv[])
{
  Vec<double,2> v0 = {0,0};
  Vec<double,2> v1 = {1,1};
  Vec<double,2> v2 = {2,2};
  Vec<double,2> v3 = {3,3};
  Vec<double,2> v4 = {4,4};
  std::vector<Vec<double,2> > vect{v0,v1,v2,v3,v4};
  Polynomial<4,Vec<double,2> > poly(vect);
  std::cout << poly << std::endl;

  Polynomial<4,Vec<double,2> > poly1(poly);
  std::cout << poly1 << std::endl;

  std::cout << poly.get_order() << std::endl;
  poly.set_coeff(4,v0);
  std::cout << poly.get_coeff(4) << std::endl;
  std::cout << poly << std::endl;

  poly.set_coeff(4,v4);
  Polynomial<4,Vec<double,2> > poly3(poly+poly1);
  std::cout << poly3 << std::endl;
  std::cout << poly-poly1 << std::endl;

  std::cout << poly*poly1 << std::endl;
  Polynomial<4,Vec<double,2> > poly4;
  poly4 = poly+poly1;
  std::cout << poly4 << std::endl;

  std::cout << poly4.eval(1) << std::endl;
  std::cout << poly4.derivative() << std::endl;
  
  //Works till here

  std::cout << std::endl;

  std::vector<Vec<double,1> > vectt{0,1,2,3,4};
  Polynomial<4,Vec<double,1> > polyy(vectt);
  std::cout << polyy << std::endl;

  Polynomial<4,Vec<double,1> > poly11(polyy);
  std::cout << poly11 << std::endl;

  std::cout << polyy.get_order() << std::endl;
  polyy.set_coeff(4,0);
  std::cout << polyy.get_coeff(4) << std::endl;
  std::cout << polyy << std::endl;

  polyy.set_coeff(4,4);
  Polynomial<4,Vec<double,1> > poly33(polyy+poly11);
  std::cout << poly33 << std::endl;
  std::cout << polyy-poly11 << std::endl;

  std::cout << polyy*poly11 << std::endl;
  Polynomial<4,Vec<double,1> > poly44;
  poly44 = polyy+poly11;
  std::cout << poly44 << std::endl;

  std::cout << poly44.eval(1) << std::endl;
  std::cout << poly44.derivative() << std::endl;
  return 0;
}
