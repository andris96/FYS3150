#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "TestPenningTrap.hpp"

TestPenningTrap::TestPenningTrap() = default; // Empty constructor

void TestPenningTrap::runAllTests() {
    test_add_particle();
    test_reset();
    test_is_within_trap();
    test_fill_with_particles();
}

void TestPenningTrap::test_add_particle() {
    double B0 = 1;
    double V0 = 1;
    double d = 1; 

    int q = 1;
    double m = 1;
    arma::vec r = arma::vec(3).fill(0.0);
    arma::vec v = arma::vec(3).fill(0.0);

    Particle p = Particle(q, m, r, v);
    PenningTrap trap = PenningTrap(B0, V0 ,d);
    trap.add_particle(p);

    assert((trap.get_particles().size() == 1) && ("Particle was not added to the trap!"));
}

void TestPenningTrap::test_reset() {
    double B0 = 1;
    double V0 = 1;
    double d = 1; 

    int q = 1;
    double m = 1;
    arma::vec r = arma::vec(3).fill(0.0);
    arma::vec v = arma::vec(3).fill(0.0);

    Particle p = Particle(q, m, r, v);

    PenningTrap trap = PenningTrap(B0, V0 ,d);
    for (int i = 0; i < 10; i++) {
        trap.add_particle(p);
    }
    trap.reset();

    assert((trap.get_particles().size() == 0) && ("Number of particles not 0!\n"));
    assert((trap.get_total_system_time() == 0.0) && ("Total system time not 0!\n"));
}

void TestPenningTrap::test_is_within_trap() {
    double B0 = 1;
    double V0 = 1;
    double d = 1; 

    arma::vec r_within = arma::vec(3).fill(0.577*d);
    arma::vec r_not_within = arma::vec(3).fill(1.3*d);

    PenningTrap trap = PenningTrap(B0, V0 ,d);
    
    assert((trap.is_within_trap(r_within) == true) && ("Distance within trap returned as not within\n"));
    assert((trap.is_within_trap(r_not_within) == false) && ("Distance not within trap returned as within\n"));
}

void TestPenningTrap::test_fill_with_particles() {

    double B0 = 1;
    double V0 = 1;
    double d = 1; 

    double q = 1;
    double m = 1;

    int N = 10;

    PenningTrap trap = PenningTrap(B0, V0 ,d);
    trap.fill_with_particles(q, m, N);

    std::vector<Particle> particles = trap.get_particles();
    assert((particles.size() == N) && ("The trap was not filled with " + std::to_string(N) + " particles!\n"));
    
    for (int i = 0; i < N; i++) {
        assert((arma::norm(particles.at(i).r(), 2) <= d) && ("Particle not initiated within the trap!\n"))
    }
}