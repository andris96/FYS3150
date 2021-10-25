#ifndef __TestPenningTrap_hpp__
#define __TestPenningTrap_hpp__

#include <armadillo>
#include <assert.h>
#include <string.h>

class TestPenningTrap {

    public:
    TestPenningTrap();

    void runAllTests();

    void test_add_particle();
    void test_reset();
    void test_is_within_trap();
    void test_fill_with_particles();
    void test_count_particles();
    void test_external_E_field();
    
};

#endif