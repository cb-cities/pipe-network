#include "node.h"

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

TEST_CASE("Node is checked", "[Node]") {
  Node test(123, {11.2, 22.3, 33.4});

  REQUIRE(test.id() == 123);
  REQUIRE(test.coordinates()(2) == Approx(33.4).epsilon(1.e-12));
  REQUIRE(test.coord_at_dir(1) == Approx(22.3).epsilon(1.e-12));
  REQUIRE(test.ishead() == false);
  REQUIRE(test.isdischarge() == false);

  SECTION("reading and writting hydraulic head") {
    test.head(110.1);

    REQUIRE(test.head() == Approx(110.1).epsilon(1.e-12));
    REQUIRE(test.ishead() == true);
  }

  SECTION("reading and writting discharge") {
    test.discharge(-23.4);

    REQUIRE(test.discharge() == Approx(-23.4).epsilon(1.e-12));
    REQUIRE(test.isdischarge() == true);
  }

  // test return nconnections
}
