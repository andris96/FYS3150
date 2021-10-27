#ifndef __TestPenningTrap_hpp__
#define __TestPenningTrap_hpp__

#include <armadillo>
#include <assert.h>
#include <string.h>

// This class contains unit tests for the class PenningTrap.
//
// The getters and setters were assumed to be correct, as these
// were straight forward to implement. 
// Not all methods were tested due to time constrains. 
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
    void test_external_B_field();
};

#endif