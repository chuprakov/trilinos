#ifndef KOKKOS_TEST_DUALVIEW_HPP
#define KOKKOS_TEST_DUALVIEW_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_Atomic.hpp>
#include <cmath>

namespace Test {

namespace Impl{

// This test runs the random number generators and uses some statistic tests to
// check the 'goodness' of the random numbers:
//    (i)   mean:         the mean is expected to be 0.5*RAND_MAX
//    (ii)  variance:     the variance is 1/3*mean*mean
//    (iii) covariance:   the covariance is 0
//    (iv)  1-tupledistr: the mean, variance and covariance of a 1D Histrogram of random numbers
//    (v)   3-tupledistr: the mean, variance and covariance of a 3D Histrogram of random numbers

#define HIST_DIM3D 24
#define HIST_DIM1D (HIST_DIM3D*HIST_DIM3D*HIST_DIM3D)

struct RandomProperties {
  uint64_t count;
  double mean;
  double variance;
  double covariance;
  double min;
  double max;

  KOKKOS_INLINE_FUNCTION
  RandomProperties() {
    count = 0;
    mean = 0.0;
    variance = 0.0;
    covariance = 0.0;
    min = 1e64;
    max = -1e64;
  }

  KOKKOS_INLINE_FUNCTION
  RandomProperties& operator+=(const RandomProperties& add) {
    count      += add.count;
    mean       += add.mean;
    variance   += add.variance;
    covariance += add.covariance;
    min         = add.min<min?add.min:min;
    max         = add.max>max?add.max:max;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  volatile RandomProperties& operator+=(const volatile RandomProperties& add) volatile {
    count      += add.count;
    mean       += add.mean;
    variance   += add.variance;
    covariance += add.covariance;
    min         = add.min<min?add.min:min;
    max         = add.max>max?add.max:max;
    return *this;
  }
};

template<class GeneratorPool, class Scalar>
struct test_random_functor {
  typedef typename GeneratorPool::generator_type rnd_type;

  typedef RandomProperties value_type;
  typedef typename GeneratorPool::device_type device_type;

  GeneratorPool rand_pool;
  const double mean;
  typedef Kokkos::View<int[HIST_DIM1D],typename GeneratorPool::device_type> type_1d;
  type_1d density_1d;
  typedef Kokkos::View<int[HIST_DIM3D][HIST_DIM3D][HIST_DIM3D],typename GeneratorPool::device_type> type_3d;
  type_3d density_3d;

  test_random_functor(GeneratorPool rand_pool_,type_1d d1d, type_3d d3d):
      rand_pool(rand_pool_),mean(0.5*Kokkos::rand<rnd_type,Scalar>::max()),
      density_1d(d1d),density_3d(d3d) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, RandomProperties& prop) const {
    rnd_type rand_gen = rand_pool.get_state();
    for(int k = 0;k<1024;k++) {
      const Scalar tmp = Kokkos::rand<rnd_type,Scalar>::draw(rand_gen);
      prop.count++;
      prop.mean += tmp;
      prop.variance += (tmp-mean)*(tmp-mean);
      const Scalar tmp2 = Kokkos::rand<rnd_type,Scalar>::draw(rand_gen);
      prop.count++;
      prop.mean += tmp2;
      prop.variance += (tmp2-mean)*(tmp2-mean);
      prop.covariance += (tmp-mean)*(tmp2-mean);
      const Scalar tmp3 = Kokkos::rand<rnd_type,Scalar>::draw(rand_gen);
      prop.count++;
      prop.mean += tmp3;
      prop.variance += (tmp3-mean)*(tmp3-mean);
      prop.covariance += (tmp2-mean)*(tmp3-mean);

      Kokkos::atomic_fetch_add(&density_1d(uint64_t(1.0*HIST_DIM1D*tmp/Kokkos::rand<rnd_type,Scalar>::max())),1);
      Kokkos::atomic_fetch_add(&density_1d(uint64_t(1.0*HIST_DIM1D*tmp2/Kokkos::rand<rnd_type,Scalar>::max())),1);
      Kokkos::atomic_fetch_add(&density_1d(uint64_t(1.0*HIST_DIM1D*tmp3/Kokkos::rand<rnd_type,Scalar>::max())),1);
      Kokkos::atomic_fetch_add(&density_3d(uint64_t(1.0*HIST_DIM3D*tmp/Kokkos::rand<rnd_type,Scalar>::max()),
          uint64_t(1.0*HIST_DIM3D*tmp2/Kokkos::rand<rnd_type,Scalar>::max()),
          uint64_t(1.0*HIST_DIM3D*tmp3/Kokkos::rand<rnd_type,Scalar>::max())),1);
    }
    rand_pool.free_state(rand_gen);
  }
};

template<class DeviceType>
struct test_histogram1d_functor {
  typedef RandomProperties value_type;
  typedef DeviceType device_type;

  typedef Kokkos::View<int[HIST_DIM1D],device_type> type_1d;
  type_1d density_1d;

  double mean;

  test_histogram1d_functor(type_1d d1d, int num_draws):
      density_1d(d1d),  mean(1.0*num_draws/HIST_DIM1D*3) {printf("Mean: %e\n",mean);}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, RandomProperties& prop) const {
      double count = density_1d(i);
      prop.mean += count;
      prop.variance += 1.0*(count-mean) * (count-mean);
      //prop.covariance += 1.0*count*count;
      prop.min = count<prop.min?count:prop.min;
      prop.max = count>prop.max?count:prop.max;
      if(i<HIST_DIM1D-1)
        prop.covariance += (count-mean) * (density_1d(i+1)-mean);
  }
};

template<class DeviceType>
struct test_histogram3d_functor {
  typedef RandomProperties value_type;
  typedef DeviceType device_type;

  typedef Kokkos::View<int[HIST_DIM3D][HIST_DIM3D][HIST_DIM3D],device_type> type_3d;
  type_3d density_3d;

  double mean;

  test_histogram3d_functor(type_3d d3d, int num_draws):
      density_3d(d3d), mean(1.0*num_draws/HIST_DIM1D) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (int i, RandomProperties& prop) const {
      double count = density_3d( i/(HIST_DIM3D*HIST_DIM3D),
                                (i%(HIST_DIM3D*HIST_DIM3D))/HIST_DIM3D,
                                 i%HIST_DIM3D);
      prop.mean += count;
      prop.variance += (count-mean) * (count-mean);
      if(i<HIST_DIM1D-1) {
        double count_next = density_3d( (i+1)/(HIST_DIM3D*HIST_DIM3D),
                                  ((i+1)%(HIST_DIM3D*HIST_DIM3D))/HIST_DIM3D,
                                  (i+1)%HIST_DIM3D);
        prop.covariance += (count-mean) * (count_next-mean);
      }
  }
};



template <class RandomGenerator,class Scalar>
struct test_random_scalar {
  typedef typename RandomGenerator::generator_type rnd_type;

  int pass_mean,pass_var,pass_covar;
  int pass_hist1d_mean,pass_hist1d_var,pass_hist1d_covar;
  int pass_hist3d_mean,pass_hist3d_var,pass_hist3d_covar;
  test_random_scalar(
      typename test_random_functor<RandomGenerator,int>::type_1d& density_1d,
      typename test_random_functor<RandomGenerator,int>::type_3d& density_3d,
      RandomGenerator& pool, unsigned int num_draws) {
    {
      RandomProperties result;

      Kokkos::parallel_reduce(num_draws/1024,test_random_functor<RandomGenerator,Scalar>(pool,density_1d,density_3d),result);

      //printf("Result: %lf %lf %lf\n",result.mean/num_draws/3,result.variance/num_draws/3,result.covariance/num_draws/2);
      double tolerance = 1.6*sqrt(1.0/num_draws);
      double mean_expect = 0.5*Kokkos::rand<rnd_type,Scalar>::max();
      double variance_expect = 1.0/3.0*mean_expect*mean_expect;
      double mean_eps = mean_expect/(result.mean/num_draws/3)-1.0;
      double variance_eps = variance_expect/(result.variance/num_draws/3)-1.0;
      double covariance_eps = result.covariance/num_draws/2/variance_expect;
      pass_mean  = ((-tolerance < mean_eps) &&
                    ( tolerance > mean_eps)) ? 1:0;
      pass_var   = ((-tolerance < variance_eps) &&
                    ( tolerance > variance_eps)) ? 1:0;
      pass_covar = ((-1.4*tolerance < covariance_eps) &&
                    ( 1.4*tolerance > covariance_eps)) ? 1:0;
      printf("Pass: %i %i %e %e %e || %e\n",pass_mean,pass_var,mean_eps,variance_eps,covariance_eps,tolerance);
    }
    {
      RandomProperties result;

      Kokkos::parallel_reduce(HIST_DIM1D,test_histogram1d_functor<typename RandomGenerator::device_type>(density_1d,num_draws),result);

      //printf("Result: %lf %lf %lf\n",result.mean/num_draws/3,result.variance/num_draws/3,result.covariance/num_draws/2);
      double tolerance = 6*sqrt(1.0/HIST_DIM1D);
      double mean_expect = 1.0*num_draws*3/HIST_DIM1D;
      double variance_expect = 1.0*num_draws*3/HIST_DIM1D*(1.0-1.0/HIST_DIM1D);
      double covariance_expect = -1.0*num_draws*3/HIST_DIM1D/HIST_DIM1D;
      double mean_eps = mean_expect/(result.mean/HIST_DIM1D)-1.0;
      double variance_eps = variance_expect/(result.variance/HIST_DIM1D)-1.0;
      double covariance_eps = (result.covariance/HIST_DIM1D - covariance_expect)/mean_expect;
      pass_hist1d_mean  = ((-tolerance < mean_eps) &&
                           ( tolerance > mean_eps)) ? 1:0;
      pass_hist1d_var   = ((-tolerance < variance_eps) &&
                           ( tolerance > variance_eps)) ? 1:0;
      pass_hist1d_covar = ((-tolerance < covariance_eps) &&
                           ( tolerance > covariance_eps)) ? 1:0;
      printf("Density 1D: %e %e %e || %e %e %e || %e %e || %e %e\n",mean_eps,variance_eps,result.covariance/HIST_DIM1D/HIST_DIM1D,tolerance,
          result.min,result.max,
          result.variance/HIST_DIM1D,1.0*num_draws*3/HIST_DIM1D*(1.0-1.0/HIST_DIM1D),
          result.covariance/HIST_DIM1D,-1.0*num_draws*3/HIST_DIM1D/HIST_DIM1D
          );

    }
    {
      RandomProperties result;

      Kokkos::parallel_reduce(HIST_DIM1D,test_histogram3d_functor<typename RandomGenerator::device_type>(density_3d,num_draws),result);

      //printf("Result: %lf %lf %lf\n",result.mean/num_draws/3,result.variance/num_draws/3,result.covariance/num_draws/2);
      double tolerance = 6*sqrt(1.0/HIST_DIM1D);
      double mean_expect = 1.0*num_draws/HIST_DIM1D;
      double variance_expect = 1.0*num_draws/HIST_DIM1D*(1.0-1.0/HIST_DIM1D);
      double covariance_expect = -1.0*num_draws/HIST_DIM1D/HIST_DIM1D;
      double mean_eps = mean_expect/(result.mean/HIST_DIM1D)-1.0;
      double variance_eps = variance_expect/(result.variance/HIST_DIM1D)-1.0;
      double covariance_eps = (result.covariance/HIST_DIM1D - covariance_expect)/mean_expect;
      pass_hist3d_mean  = ((-tolerance < mean_eps) &&
                           ( tolerance > mean_eps)) ? 1:0;
      pass_hist3d_var   = ((-tolerance < variance_eps) &&
                           ( tolerance > variance_eps)) ? 1:0;
      pass_hist3d_covar = ((-tolerance < covariance_eps) &&
                           ( tolerance > covariance_eps)) ? 1:0;
      printf("Density 3D: %e %e %e || %e %e %e\n",mean_eps,variance_eps,result.covariance/HIST_DIM1D/HIST_DIM1D,tolerance,result.min,result.max);

    }

  }
};

template <class RandomGenerator>
void test_random(unsigned int num_draws)
{
  typedef typename RandomGenerator::generator_type rnd_type;
  typename test_random_functor<RandomGenerator,int>::type_1d density_1d("D1d");
  typename test_random_functor<RandomGenerator,int>::type_3d density_3d("D3d");

  RandomGenerator pool(31891);
  test_random_scalar<RandomGenerator,int> test_int(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_int.pass_mean,1);
  ASSERT_EQ( test_int.pass_var,1);
  ASSERT_EQ( test_int.pass_covar,1);
  ASSERT_EQ( test_int.pass_hist1d_mean,1);
  ASSERT_EQ( test_int.pass_hist1d_var,1);
  ASSERT_EQ( test_int.pass_hist1d_covar,1);
  ASSERT_EQ( test_int.pass_hist3d_mean,1);
  ASSERT_EQ( test_int.pass_hist3d_var,1);
  ASSERT_EQ( test_int.pass_hist3d_covar,1);
  deep_copy(density_1d,0);
  deep_copy(density_3d,0);
  test_random_scalar<RandomGenerator,unsigned int> test_uint(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_uint.pass_mean,1);
  ASSERT_EQ( test_uint.pass_var,1);
  ASSERT_EQ( test_uint.pass_covar,1);
  ASSERT_EQ( test_uint.pass_hist1d_mean,1);
  ASSERT_EQ( test_uint.pass_hist1d_var,1);
  ASSERT_EQ( test_uint.pass_hist1d_covar,1);
  ASSERT_EQ( test_uint.pass_hist3d_mean,1);
  ASSERT_EQ( test_uint.pass_hist3d_var,1);
  ASSERT_EQ( test_uint.pass_hist3d_covar,1);
  deep_copy(density_1d,0);
  deep_copy(density_3d,0);
  test_random_scalar<RandomGenerator,int64_t> test_int64(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_int64.pass_mean,1);
  ASSERT_EQ( test_int64.pass_var,1);
  ASSERT_EQ( test_int64.pass_covar,1);
  ASSERT_EQ( test_int64.pass_hist1d_mean,1);
  ASSERT_EQ( test_int64.pass_hist1d_var,1);
  ASSERT_EQ( test_int64.pass_hist1d_covar,1);
  ASSERT_EQ( test_int64.pass_hist3d_mean,1);
  ASSERT_EQ( test_int64.pass_hist3d_var,1);
  ASSERT_EQ( test_int64.pass_hist3d_covar,1);
  deep_copy(density_1d,0);
  deep_copy(density_3d,0);
  test_random_scalar<RandomGenerator,uint64_t> test_uint64(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_uint64.pass_mean,1);
  ASSERT_EQ( test_uint64.pass_var,1);
  ASSERT_EQ( test_uint64.pass_covar,1);
  ASSERT_EQ( test_uint64.pass_hist1d_mean,1);
  ASSERT_EQ( test_uint64.pass_hist1d_var,1);
  ASSERT_EQ( test_uint64.pass_hist1d_covar,1);
  ASSERT_EQ( test_uint64.pass_hist3d_mean,1);
  ASSERT_EQ( test_uint64.pass_hist3d_var,1);
  ASSERT_EQ( test_uint64.pass_hist3d_covar,1);
  deep_copy(density_1d,0);
  deep_copy(density_3d,0);
  test_random_scalar<RandomGenerator,float> test_float(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_float.pass_mean,1);
  ASSERT_EQ( test_float.pass_var,1);
  ASSERT_EQ( test_float.pass_covar,1);
  ASSERT_EQ( test_float.pass_hist1d_mean,1);
  ASSERT_EQ( test_float.pass_hist1d_var,1);
  ASSERT_EQ( test_float.pass_hist1d_covar,1);
  ASSERT_EQ( test_float.pass_hist3d_mean,1);
  ASSERT_EQ( test_float.pass_hist3d_var,1);
  ASSERT_EQ( test_float.pass_hist3d_covar,1);
  deep_copy(density_1d,0);
  deep_copy(density_3d,0);
  test_random_scalar<RandomGenerator,double> test_double(density_1d,density_3d,pool,num_draws);
  ASSERT_EQ( test_double.pass_mean,1);
  ASSERT_EQ( test_double.pass_var,1);
  ASSERT_EQ( test_double.pass_covar,1);
  ASSERT_EQ( test_double.pass_hist1d_mean,1);
  ASSERT_EQ( test_double.pass_hist1d_var,1);
  ASSERT_EQ( test_double.pass_hist1d_covar,1);
  ASSERT_EQ( test_double.pass_hist3d_mean,1);
  ASSERT_EQ( test_double.pass_hist3d_var,1);
  ASSERT_EQ( test_double.pass_hist3d_covar,1);
}
}

} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP
