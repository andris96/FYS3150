#include <UnitTest++/UnitTest++.h>
#include <stdlib.h>

#include "Particle.hpp"
#include "PenningTrap.hpp"

TEST(ResetTrap) {
    double B0 = 1;
    double V0 = 1;
    double d = 1; 

    int q = 1;
    double m = 1;

    PenningTrap trap = PenningTrap(B0, V0 ,d);
    trap.fill_with_particles(10, m, q);
    trap.reset();

    CHECK(trap.get_particles().size() == 0);
    CHECK(trap.get_total_system_time() == 0.0);
}

int main() {
    return UnitTest::RunAllTests();
}