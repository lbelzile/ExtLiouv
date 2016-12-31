#ifndef EXPONENTMEASURES_HPP
#define EXPONENTMEASURES_HPP

//Functions by Leonid Rousniak and Adrien de Castro

// Wraps R::pnorm in C++
inline double p_norm(double x) {
	return R::pnorm(x,0,1,true,false);
}

// Wraps R::dnorm in C++
inline double d_norm(double x) {
  return R::dnorm(x,0,1,false);
}

// Wraps R::pnorm in C++
inline double p_t(double x, double nu) {
  return R::pt(x,nu,true,false);
}

// Wraps R::dnorm in C++
inline double d_t(double x, double nu) {
  return R::dt(x,nu,false);
}

///////////////////////////////////////////
// Functions for HR model /////////////////
///////////////////////////////////////////
inline double ztr(double z1, double z2, double a) {
  return a/2.0-log(z1/z2)/a;
}

inline double V_BR(double z1, double z2, double a) {
  return p_norm(ztr(z1,z2,a))/z1 + p_norm(ztr(z2,z1,a))/z2;
}

inline double V_BR1(double z1, double z2, double a) {
  return -(p_norm(ztr(z1,z2,a)) + d_norm(ztr(z1,z2,a))/a)/(z1*z1) +
    d_norm(ztr(z2,z1,a))/(a*z1*z2);
}

inline double V_BR2(double z1, double z2, double a) {
  return V_BR1(z2, z1, a);
}

inline double V_BR12(double z1, double z2, double a) {
  return -(d_norm(ztr(z1,z2,a))/z1*(1-ztr(z1,z2,a)/a) +
    d_norm(ztr(z2,z1,a))/z2*(1-ztr(z2,z1,a)/a)) / (a*z1*z2);
}


///////////////////////////////////////////
// Functions for Extremal-t model /////////
///////////////////////////////////////////
inline double exstudf(double z1, double z2, double b, double rho, double nu) {
	return b*(pow(z2/z1,1/nu) - rho);
}

inline double V_EXST(double z1, double z2, double b, double rho, double nu) {
	return p_t(exstudf(z1,z2,b,rho,nu),nu)/z1 + p_t(exstudf(z2,z1,b,rho,nu),nu)/z2;
}

inline double V_EXST1(double z1, double z2, double b, double rho, double nu) {
  return -p_t(exstudf(z1,z2,b,rho,nu),nu)/(z1*z1);
}

inline double V_EXST2(double z1, double z2, double b, double rho, double nu) {
  return V_EXST1(z2,z1,b,rho,nu);
}

inline double V_EXST12(double z1, double z2, double b, double rho, double nu) {
  return -(b/nu)*d_t(exstudf(z1,z2,b,rho,nu),nu)*pow(z2,1.0/nu-1.0)/pow(z1,1.0/nu+2.0);
}

#endif
