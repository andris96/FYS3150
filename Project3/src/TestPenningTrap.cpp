#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "TestPenningTrap.hpp"

TestPenningTrap::TestPenningTrap() = default; // Empty constructor

void TestPenningTrap::runAllTests() {
    test_add_particle();
    test_reset();
    test_is_within_trap();
    test_fill_with_particles();
    test_count_particles();
    test_external_E_field();
    test_external_B_field();
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
    assert((particles.size() == N) && ("The trap was not filled with 10 particles!\n"));
    
    for (int i = 0; i < N; i++) {
        assert((arma::norm(particles.at(i).r(), 2) <= d) && ("Particle not initiated within the trap!\n"));
    }
}

void TestPenningTrap::test_count_particles() {

    double B0 = 1;
    double V0 = 1;
    double d = 1; 

    double q = 1;
    double m = 1;

    arma::vec r_within = arma::vec(3).fill(0.1*d);
    arma::vec v = arma::vec(3).fill(0.0);

    int N = 10;
    Particle p = Particle(q, m, r_within, v);

    PenningTrap trap = PenningTrap(B0, V0 ,d);
    for (int i = 0; i < N; i++) {
        trap.add_particle(p);
    }

    assert((trap.count_particles() ==  N) && ("The particle counter did return 10 as expected!\n"));

}

void TestPenningTrap::test_external_E_field() {

    // Case 1: Particle within trap, no time-dependence

    double B0 = 1; // Not relevant, only for instantiation
    double V0 = 1;
    double d = 1000; 
    double q = 1; // same
    double m = 1; // same
    arma::vec r = arma::vec(3).fill(0.1*d);
    arma::vec v = arma::vec(3).fill(0.0);

    PenningTrap trap = PenningTrap(B0, V0 ,d);

    arma::vec e_field_numerical = trap.external_E_field(r);
    arma::vec e_field_analytical = arma::vec(3).fill(0.0);
    e_field_analytical << 0.0001 << 0.0001 << -0.0002;

    double tolerance = 10e-15;
    assert((arma::approx_equal(e_field_analytical, e_field_numerical, "absdiff", tolerance)) && ("Numerical and analytical result not equal for E-field!\n"));

    // Case 2: Particle within the trap, with time-dependence
    trap.enable_time_dependence(); 

    double t = 1;
    double f = 1;
    double omega = 1;
    trap.set_total_time(t);

    double f_cos_omega_t = f*cos(omega*t);
    
    e_field_numerical = trap.external_E_field(r);
    e_field_analytical *= (1 + f_cos_omega_t);

    assert((arma::approx_equal(e_field_analytical, e_field_numerical, "absdiff", tolerance)) && ("Numerical and analytical result not equal for E-field with time-dependence!\n"));

    // Case 3: Particle outside the trap
    r = arma::vec(3).fill(2*d);

    e_field_numerical = trap.external_E_field(r);
    e_field_analytical = arma::vec(3).fill(0.0);

    assert((arma::approx_equal(e_field_analytical, e_field_numerical, "absdiff", tolerance)) && ("E-field not cancelled outside of the trap!\n"));

}

void TestPenningTrap::test_external_B_field() {

    // Case 1: Particle within the trap
    double B0 = 1; 
    double V0 = 1; // Not relevant, only for instantiation
    double d = 1000; // same
    double q = 1; // same
    double m = 1; // same
    arma::vec r = arma::vec(3).fill(0.1*d);
    arma::vec v = arma::vec(3).fill(0.0);

    PenningTrap trap = PenningTrap(B0, V0 ,d);

    arma::vec B_field_numerical = trap.external_B_field(r);
    arma::vec B_field_analytical = arma::vec(3).fill(0.0);
    B_field_analytical.at(2) = 1.0;
    
    double tolerance = 10e-15;
    assert((arma::approx_equal(B_field_analytical, B_field_numerical, "absdiff", tolerance)) && ("Analytical and numerical B-field not equal!\n"));

    // Case 2: Particle outside the trap
    r = arma::vec(3).fill(2*d);

    B_field_numerical = trap.external_B_field(r);
    B_field_analytical.at(2) = 0.0;

    assert((arma::approx_equal(B_field_analytical, B_field_numerical, "absdiff", tolerance)) && ("B-field not cancelled outside of the trap!\n"));

}