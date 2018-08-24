// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE interpolation tests
#include <boost/test/included/unit_test.hpp>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Interpolation.hpp"
#include "Acts/Utilities/detail/interpolation_impl.hpp"

namespace Acts {

using namespace detail;

namespace Test {

  BOOST_AUTO_TEST_CASE(interpolation_1d)
  {
    typedef std::array<double, 1u> Point;
    typedef std::array<double, 2u> Values;

    Point  low  = {{1.}};
    Point  high = {{2.}};
    Values v    = {{10., 20.}};

    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{0.5}}), low, high, v), 5., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.}}), low, high, v), 10., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.3}}), low, high, v), 13., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.5}}), low, high, v), 15., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.8}}), low, high, v), 18., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2.}}), low, high, v), 20., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2.3}}), low, high, v), 23., 1e-6);
  }

  BOOST_AUTO_TEST_CASE(interpolation_2d)
  {
    typedef std::array<double, 2u> Point;
    typedef std::array<double, 4u> Values;

    Point  low  = {{1., 1.}};
    Point  high = {{2., 3.}};
    Values v    = {{10., 30., 20., 40.}};

    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 1.}}), low, high, v), 10., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 1.}}), low, high, v), 20., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 3.}}), low, high, v), 30., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 3.}}), low, high, v), 40., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.3, 1.}}), low, high, v), 13., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.5, 1.}}), low, high, v), 15., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.8, 1.}}), low, high, v), 18., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.3, 3.}}), low, high, v), 33., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.5, 3.}}), low, high, v), 35., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.8, 3.}}), low, high, v), 38., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 1.7}}), low, high, v), 17., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 2.}}), low, high, v), 20., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 2.5}}), low, high, v), 25., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 1.7}}), low, high, v), 27., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 2.}}), low, high, v), 30., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 2.5}}), low, high, v), 35., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.5, 2.}}), low, high, v), 25., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.3, 1.7}}), low, high, v), 20., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.3, 2.5}}), low, high, v), 28., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.8, 1.7}}), low, high, v), 25., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.8, 2.5}}), low, high, v), 33., 1e-6);
  }

  BOOST_AUTO_TEST_CASE(interpolation_3d)
  {
    typedef std::array<double, 3u> Point;
    typedef std::array<double, 8u> Values;

    Point  low  = {{1., 1., 1.}};
    Point  high = {{2., 3., 4.}};
    Values v    = {{10., 50., 30., 70., 20., 60., 40., 80.}};

    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 1., 1.}}), low, high, v), 10., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 1., 1.}}), low, high, v), 20., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 3., 1.}}), low, high, v), 30., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 3., 1.}}), low, high, v), 40., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 1., 4.}}), low, high, v), 50., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 1., 4.}}), low, high, v), 60., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 3., 4.}}), low, high, v), 70., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 3., 4.}}), low, high, v), 80., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.5, 1., 1.}}), low, high, v), 15., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.5, 3., 1.}}), low, high, v), 35., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 2., 1.}}), low, high, v), 20., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 2., 1.}}), low, high, v), 30., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.5, 1., 4.}}), low, high, v), 55., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.5, 3., 4.}}), low, high, v), 75., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 2., 4.}}), low, high, v), 60., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 2., 4.}}), low, high, v), 70., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 1., 2.5}}), low, high, v), 30., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1., 3., 2.5}}), low, high, v), 50., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 1., 2.5}}), low, high, v), 40., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{2., 3., 2.5}}), low, high, v), 60., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.5, 2., 2.5}}), low, high, v), 360. / 8, 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate(Point({{1.3, 2.1, 1.6}}), low, high, v), 32., 1e-6);
  }

  BOOST_AUTO_TEST_CASE(interpolation_mixed_point_values)
  {
    typedef ActsVectorD<1> Point1;
    typedef std::array<double, 1u> Point2;
    typedef std::vector<double> Point3;
    typedef std::array<double, 2u> Values;

    Point2 low  = {{1.}};
    Point3 high = {2.};
    Values v    = {{10., 20.}};

    Point1 p;
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate((p << 0.5).finished(), low, high, v), 5., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate((p << 1.).finished(), low, high, v), 10., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate((p << 1.3).finished(), low, high, v), 13., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate((p << 1.5).finished(), low, high, v), 15., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate((p << 1.8).finished(), low, high, v), 18., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate((p << 2.).finished(), low, high, v), 20., 1e-6);
    BOOST_CHECK_CLOSE_FRACTION(
        interpolate((p << 2.3).finished(), low, high, v), 23., 1e-6);
  }
}  // namespace Test

}  // namespace Acts