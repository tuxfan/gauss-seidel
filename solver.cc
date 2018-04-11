#include <cstddef>
#include <cmath>
#include <iostream>
#include <fstream>

#include <sys/time.h>

struct field_t {

  field_t(size_t xlen, size_t ylen)
    : xlen_(xlen), ylen_(ylen)
  {
    data_ = new double[xlen_*ylen_];
  } // field_t

  ~field_t() {
    delete[] data_;
  }

  double operator ()(size_t i, size_t j) const {
    return data_[j*xlen_ + i];
  } // operator ()

  double & operator ()(size_t i, size_t j) {
    return data_[j*xlen_ + i];
  } // operator ()

private:

  size_t xlen_;
  size_t ylen_;
  double * data_;

}; // struct field_t

const double PI = 3.14159;
const double K = 5.0;
const double L = 5.0;

#define SQR(x) (x)*(x)
const double SQR_K_L_PI = (SQR(K) + SQR(L))*SQR(PI);

int main(int argc, char ** argv) {
  if(argc < 4) {
    std::cerr << "Usage: " << argv[0] << " XLEN YLEN ITA" << std::endl;
    std::exit(1);
  } // if

  size_t XLEN = atoi(argv[1]);
  size_t YLEN = atoi(argv[2]);
  size_t ITA = atoi(argv[3]);

  std::cerr << "Running " << ITA << " iterations on " <<
    XLEN << "x" << YLEN << " mesh" << std::endl;

  field_t u(XLEN, YLEN);
  field_t f(XLEN, YLEN);
  field_t s(XLEN, YLEN);
  field_t l2(XLEN, YLEN);

  double delta = 1.0/double(XLEN-1);

  std::cerr << "delta^2 = " << SQR(delta) << std::endl;

  // Initialize fields
  for(size_t j{0}; j<YLEN; ++j) {

    const double y = j*delta;

    for(size_t i{0}; i<XLEN; ++i) {

      const double x = i*delta;

      f(i,j) = SQR_K_L_PI*sin(K*PI*x)*sin(L*PI*y);
      s(i,j) = sin(K*PI*x)*sin(L*PI*y);

      if((i == 0 || i == XLEN-1) || (j == 0 || j == YLEN-1)) {
        u(i,j) = sin(K*PI*x)*sin(L*PI*y);
      }
      else {
        u(i,j) = 0.0;
      } // if
    } // for
  } // for

  timeval t0;
  timeval tf;
  gettimeofday(&t0, nullptr);

#if defined (GS)

  // Gauss-Seidel
  for(size_t n{0}; n<ITA; ++n) {
    for(size_t j{1}; j<YLEN-1; ++j) {
      for(size_t i{1}; i<XLEN-1; ++i) {
        u(i,j) = 0.25*(delta*delta*f(i,j) +
            u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)); 
      } // for
    } // for
  } // for

#elif defined(RBGS)

  for(size_t n{0}; n<ITA; ++n) {

    // Red points
    for(size_t j{1}; j<YLEN-1; ++j) {
      for(size_t i{(j-1)%2+1}; i<XLEN-1; i+=2) {
        //std::cerr << "(" << i << "," << j << ")" << std::endl;
        u(i,j) = 0.25*(delta*delta*f(i,j) +
            u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)); 
      } // for
    } // for

    // Black points
    for(size_t j{1}; j<YLEN-1; ++j) {
      for(size_t i{j%2+1}; i<XLEN-1; i+=2) {
        //std::cerr << "(" << i << "," << j << ")" << std::endl;
        u(i,j) = 0.25*(delta*delta*f(i,j) +
            u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1)); 
      } // for
    } // for

  } // for

#endif

  gettimeofday(&tf, nullptr);

  double elapsed = (tf.tv_sec - t0.tv_sec) +
    (tf.tv_usec - t0.tv_usec)/1000000.0;

  std::cerr << "elapsed time: " << elapsed << std::endl;

  std::ofstream solution("solution.dat", std::ofstream::out);

  for(size_t j{0}; j<YLEN; ++j) {
    const double y = j*delta;
    for(size_t i{0}; i<XLEN; ++i) {
      const double x = i*delta;
      solution << x << " " << y << " " << u(i,j) << std::endl;
    } // for
  } // for

  std::ofstream residual("residual.dat", std::ofstream::out);

  for(size_t j{0}; j<YLEN; ++j) {
    const double y = j*delta;
    for(size_t i{0}; i<XLEN; ++i) {
      const double x = i*delta;
      l2(i,j) = sqrt(SQR(s(i,j) - u(i,j)));
      residual << x << " " << y << " " << l2(i,j) << std::endl;
    } // for
  } // for

  return 0;
} // main
