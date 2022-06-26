#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//- global variables

//- test function
//' Return input string
//'
//' @param String
//' @export
// [[Rcpp::export]]
String ReturnInputC(String line){
return line;}
