

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include <iostream>
#include <boost/shared_ptr.hpp>

// To run in parallel
// #include "PetscSetupAndFinalize.hpp"
// When run in serial
#include "FakePetscSetup.hpp"

#include "AbstractCellBasedTestSuite.hpp"


class BaseClass{
public:
    BaseClass(){};
    virtual ~BaseClass(){};
    virtual void testfn()
    {
        std::cout << "Base class" << std::endl;
    }
};

class DerivedClass : public BaseClass{
public:
    virtual void testfn()
    {
        std::cout << "Derived class" << std::endl;
    }
};

class ContainsClass
{
    BaseClass* mpClass;
public:
    ContainsClass(BaseClass* input)
    {
        mpClass = input;
        mpClass->testfn();
    }
};

class TestClasses : public AbstractCellBasedTestSuite
{
public:
    void TestRunClasses() throw(Exception)
    {
        BaseClass baseclass;
        DerivedClass derivedclass;

        ContainsClass containsbase(&baseclass);
        ContainsClass containsderived(&derivedclass);
    }
};