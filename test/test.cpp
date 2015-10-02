
#include <triqs/gfs.hpp>

using namespace triqs::gfs;
using triqs::clef::placeholder;

int main()
{
   double beta=1;
   int nw=100;

   auto g1 = gf<imfreq>{ {beta,Fermion,nw}, {1,1} };
   auto g2 = gf<cartesian_product<imfreq,imfreq>>{ {{beta,Fermion,nw}, {beta,Fermion,nw}}, {1,1} };

   placeholder<0> w_;
   placeholder<1> nu_;

   g1(w_) << 1/(w_);
   g2(w_,nu_) << 1/(w_+nu_-4);

   g1 + 1.0; 
   g2 + 2.0; 

}
