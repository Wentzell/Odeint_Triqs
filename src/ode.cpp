
#include <iostream>
#include <complex>
#include <triqs/gfs.hpp>

#include <boost/numeric/odeint.hpp>
#include <arithmetic_tuple.h>

using namespace ReaK; 
using namespace triqs::gfs;
using triqs::clef::placeholder;
using dcomplex = std::complex< double >; 

const double BETA=1;//inverse temperature
const int N=100;//number of Matsubara frequencies

// Small wrapper around used gf types to fix size and target space automatically

class gf_1d_t : public gf<imfreq>
{
   public:
      typedef gf<imfreq> base_t; 

      using base_t::gf; 

      gf_1d_t()
	 : base_t( { { BETA, Fermion, N},{1,1} } )
	 {}
}; 

double norm( const gf_1d_t& gf_obj ) { return 0.0; } //{ dcomplex val = gf_obj[0]; return std::abs(val); } 

class gf_2d_t : public gf<cartesian_product<imfreq,imfreq>>
{
   public:
      typedef gf<cartesian_product<imfreq,imfreq>> base_t; 

      using base_t::gf; 

      gf_2d_t()
	 : base_t( { {{BETA,Fermion,N},{BETA,Boson,N}},{1,1} } )
	 {}
}; 

double norm( const gf_2d_t& gf_obj ) { return 0.0; } //{ dcomplex val = gf_obj[{0,0}]; return std::abs(val); } 

// Define absolute value of gf, necessary for adaptive integration routines

template<typename... T>
gf< T... > abs( gf< T... > gf_obj )
{
   return gf_obj; 
}


// The state type for the Ode solver, tuple of gf's with arithmetic operations

typedef arithmetic_tuple< gf_1d_t, gf_2d_t > state_t; 


// The rhs of x' = f(x) defined as a class 

class rhs_t{
   public:
      void operator()( const state_t &x , state_t &dxdt , const double  t  )
      {
	 std::cout << " Evaluation at scale " << t << std::endl; 

	 // References to current gf objects
	 const auto &X = std::get< 0 >( x );
	 const auto &Y = std::get< 1 >( x );

	 // References to derivative objects
	 auto &dX = std::get<0>( dxdt );
	 auto &dY = std::get<1>( dxdt );

	 // Initialize derivate objects
	 placeholder<0> w_;
	 placeholder<1> nu_;

	 dX(w_) << 1/(w_);
	 dY(w_,nu_) << 1/(w_+nu_);
      }
};

// Norm of state_t, needed for adaptive stepping routines

namespace boost { namespace numeric { namespace odeint {
   template<>
      struct vector_space_norm_inf< state_t >
      {
	 typedef double result_type;
	 double operator()( const state_t &p ) const
	 {
	    using namespace std; 
	    return max( norm( get<0>(p) ), norm( get<1>(p) ) );
	 }
      };
}}}


int main(int /* argc */ , char** /* argv */ )
{
   using namespace boost::numeric::odeint;
   using namespace std; 

   state_t state_vec; 

   auto &X = get< 0 >( state_vec );
   auto &Y = get< 1 >( state_vec );

   // Initialize current state
   placeholder<0> iw_;
   placeholder<1> iW_;
   X(iw_) << 1.0 / iw_;
   Y(iw_,iW_) << 1.0 / ( iw_ + iW_ );

   // Save copies of initial gfs
   gf_1d_t X0( X );  
   gf_2d_t Y0( Y ); 

   // Some tests
   gf<imfreq> A{ { BETA, Fermion, N},{1,1} }; 
   A(iw_) << 1.0 / iw_; 
   abs( A ); 
   abs( X ); 
   abs( state_vec ); 
   //Y + 0.1; 
   vector_space_norm_inf< state_t > norm; 
   norm( state_vec ); 
   state_vec/state_vec; 

   // instantiate rhs object
   rhs_t rhs;

   // Type of adaptive stepper, use vector_space_algebra here!
   typedef runge_kutta_cash_karp54< state_t, double, state_t, double, vector_space_algebra > error_stepper_t; 

   // Constants
   double ERR_ABS = 0.01; 
   double ERR_REL = 0.01; 

   double LAM_START = 0.0; 
   double LAM_FIN = 1.0; 
   double INIT_STEP = 0.1; 

   // Integrate ODE 
   //int steps = integrate_adaptive( make_controlled< error_stepper_t >( ERR_ABS, ERR_REL ), rhs, state_vec, LAM_START, LAM_FIN, INIT_STEP ); 
   //cout << " Final tuple " << state_vec << endl; 
   error_stepper_t stepper; 
   //int steps = integrate_const( stepper, rhs, state_vec, LAM_START, LAM_FIN, INIT_STEP ); 

   // Output results
   cout << " X0 " << X0(1) << endl; 
   cout << " Y0 " << Y0(0,1) << endl; 

   cout << " X0 " << X0[1] << endl; 
   cout << " Y0 " << Y0[{0,1}] << endl; 

   cout << " X " << X(1) << endl; 
   cout << " Y " << Y(0,1) << endl; 

}
