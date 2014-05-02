#ifndef KOKKOS_TEST_BITSET_HPP
#define KOKKOS_TEST_BITSET_HPP

#include <gtest/gtest.h>
#include <iostream>


namespace Test {

namespace Impl {

template <typename Bitset, bool Set>
struct TestBitset
{
  typedef Bitset bitset_type;
  typedef typename bitset_type::device_type device_type;
  typedef uint32_t value_type;

  bitset_type m_bitset;

  TestBitset( bitset_type const& bitset)
    : m_bitset(bitset)
  {}

  unsigned apply(unsigned collisions)
  {
    device_type::fence();

    unsigned count = 0;
    Kokkos::parallel_reduce( m_bitset.size()*collisions, *this, count);
    return count;
  }


  KOKKOS_INLINE_FUNCTION
  void init( value_type & v ) const { v = 0; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst, const volatile value_type & src ) const
  { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()(uint32_t i, value_type & v) const
  {
    i = i % m_bitset.size();
    if (Set) {
      if (m_bitset.set(i)) {
        if (m_bitset.test(i)) ++v;
      }
    }
    else {
      if (m_bitset.reset(i)) {
        if (!m_bitset.test(i)) ++v;
      }
    }
  }

};

template <typename Bitset>
struct TestBitsetTest
{
  typedef Bitset bitset_type;
  typedef typename bitset_type::device_type device_type;
  typedef uint32_t value_type;

  bitset_type m_bitset;

  TestBitsetTest( bitset_type const& bitset)
    : m_bitset(bitset)
  {}

  unsigned apply()
  {
    device_type::fence();

    unsigned count = 0;
    Kokkos::parallel_reduce( m_bitset.size(), *this, count);
    return count;
  }


  KOKKOS_INLINE_FUNCTION
  void init( value_type & v ) const { v = 0; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst, const volatile value_type & src ) const
  { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()(uint32_t i, value_type & v) const
  {
    if (m_bitset.test( i )) ++v;
  }
};

template <typename Bitset, bool Set>
struct TestBitsetAny
{
  typedef Bitset bitset_type;
  typedef typename bitset_type::device_type device_type;
  typedef uint32_t value_type;

  bitset_type m_bitset;

  TestBitsetAny( bitset_type const& bitset)
    : m_bitset(bitset)
  {}

  unsigned apply()
  {
    device_type::fence();

    unsigned count = 0;
    Kokkos::parallel_reduce( m_bitset.size(), *this, count);
    return count;
  }


  KOKKOS_INLINE_FUNCTION
  void init( value_type & v ) const { v = 0; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst, const volatile value_type & src ) const
  { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()(uint32_t i, value_type & v) const
  {
    bool result = false;
    unsigned attempts = 0;
    uint32_t hint = (i >> 4) << 4;
    while (attempts < m_bitset.max_hint()) {
      if (Set) {
        Kokkos::tie(result, hint) = m_bitset.find_any_unset_near(hint, i);
        if (result && m_bitset.set(hint)) {
          ++v;
          break;
        }
        else if (!result) {
          ++attempts;
        }
      }
      else {
        Kokkos::tie(result, hint) = m_bitset.find_any_set_near(hint, i);
        if (result && m_bitset.reset(hint)) {
          ++v;
          break;
        }
        else if (!result) {
          ++attempts;
        }
      }
    }
  }

};
} // namespace Impl



template <typename Device>
void test_bitset()
{
  typedef Kokkos::Bitset< Device > bitset_type;
  typedef Kokkos::ConstBitset< Device > const_bitset_type;

  //unsigned test_sizes[] = { 0u, 1000u, 1u<<14, 1u<<16, 10000001 };
  unsigned test_sizes[] = { 1000u, 1u<<14, 1u<<16, 10000001 };

  for (int i=0, end = sizeof(test_sizes)/sizeof(unsigned); i<end; ++i) {

    //std::cout << "Bitset " << test_sizes[i] << std::endl;

    bitset_type bitset(test_sizes[i]);

    //std::cout << "  Check inital count " << std::endl;
    // nothing should be set
    {
      Impl::TestBitsetTest< bitset_type > f(bitset);
      uint32_t count = f.apply();
      EXPECT_EQ(0u, count);
      EXPECT_EQ(count, bitset.count());
    }

    //std::cout << "  Check set() " << std::endl;
    bitset.set();
    // everything should be set
    {
      Impl::TestBitsetTest< const_bitset_type > f(bitset);
      uint32_t count = f.apply();
      EXPECT_EQ(bitset.size(), count);
      EXPECT_EQ(count, bitset.count());
    }

    //std::cout << "  Check reset() " << std::endl;
    bitset.reset();
    EXPECT_EQ(0u, bitset.count());

    //std::cout << "  Check set(i) " << std::endl;
    // test setting bits
    {
      Impl::TestBitset< bitset_type, true > f(bitset);
      uint32_t count = f.apply(10u);
      EXPECT_EQ( bitset.size(), bitset.count());
      EXPECT_EQ( bitset.size(), count );
    }

    //std::cout << "  Check reset(i) " << std::endl;
    // test resetting bits
    {
      Impl::TestBitset< bitset_type, false > f(bitset);
      uint32_t count = f.apply(10u);
      EXPECT_EQ( bitset.size(), count);
      EXPECT_EQ( 0u, bitset.count() );
    }


    //std::cout << "  Check find_any_set(i) " << std::endl;
    // test setting any bits
    {
      Impl::TestBitsetAny< bitset_type, true > f(bitset);
      uint32_t count = f.apply();
      EXPECT_EQ( bitset.size(), bitset.count());
      EXPECT_EQ( bitset.size(), count );
    }

    //std::cout << "  Check find_any_unset(i) " << std::endl;
    // test resetting any bits
    {
      Impl::TestBitsetAny< bitset_type, false > f(bitset);
      uint32_t count = f.apply();
      EXPECT_EQ( bitset.size(), count);
      EXPECT_EQ( 0u, bitset.count() );
    }

  }

}

} // namespace Test

#endif //KOKKOS_TEST_BITSET_HPP

