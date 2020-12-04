/*
 * MIT License
 * 
 * Copyright (c) 2020 Michael von Mengershausen
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <cppunit/CompilerOutputter.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TextTestProgressListener.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/XmlOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <stdexcept>
#include <fstream>

// #include "Measure_Test.cpp"

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

// include all test classes
#include "CartesianCoordinate_Test.cpp"
#include "Distribution_Test.cpp"

int main( int argc, char* argv[] )
{
    // Retreive test path from command line first argument. Default to "" which resolve
    // to the top level suite.

    std::string testPath = (argc > 1) ? string(argv[1]) : std::string("");
    
    // declare all tests
    CCartesianCoordinate_Test CartesianCoordinate_Test;
    CDistribution_Test Distribution_Test;

    // Create the event manager and test controller
    CPPUNIT_NS::TestResult controller;

    // Add a listener that colllects test result
    CPPUNIT_NS::TestResultCollector result;
    controller.addListener( &result );        

    // Add a listener that print dots as test run.
    #ifdef WIN32
    CPPUNIT_NS::TextTestProgressListener progress;
    #else
    CPPUNIT_NS::BriefTestProgressListener progress;
    #endif
    controller.addListener( &progress );   

    
    // Add the top suite to the test runner 
    CPPUNIT_NS::TestRunner runner;

    /////////////////////////////////////////////////
    //// add all tests of CCartesianCoordinate_Test
    /////////////////////////////////////////////////
    runner.addTest( new CppUnit::TestCaller<CCartesianCoordinate_Test> ( 
                    "CartesianCoordinate: operators",
                    &CCartesianCoordinate_Test::CheckOperators,
                    &CartesianCoordinate_Test
                    )
                  );
    runner.addTest( new CppUnit::TestCaller<CCartesianCoordinate_Test> ( 
                    "CartesianCoordinate: functions",
                    &CCartesianCoordinate_Test::CheckFunctions,
                    &CartesianCoordinate_Test
                    )
                  );
    /////////////////////////////////////////////////
    //// add all tests of Distribution_Test
    /////////////////////////////////////////////////
    runner.addTest( new CppUnit::TestCaller<CDistribution_Test> ( 
                    "Distribution: DistributionBuildUp",
                    &CDistribution_Test::DistributionBuildUp,
                    &Distribution_Test
                    )
                  );
    runner.addTest( new CppUnit::TestCaller<CDistribution_Test> ( 
                    "Distribution: DistributionOperations",
                    &CDistribution_Test::DistributionOperations,
                    &Distribution_Test
                    )
                  );
    runner.addTest( new CppUnit::TestCaller<CDistribution_Test> ( 
                    "Distribution: DistributionFunctionOnVariable",
                    &CDistribution_Test::DistributionFunctionOnVariable,
                    &Distribution_Test
                    )
                  );
    runner.addTest( new CppUnit::TestCaller<CDistribution_Test> ( 
                    "Distribution: DistributionCalculations",
                    &CDistribution_Test::DistributionCalculations,
                    &Distribution_Test
                    )
                  );
 


    try
    {
        CPPUNIT_NS::stdCOut() << "Running "  <<  testPath;
        runner.run( controller, testPath );

        CPPUNIT_NS::stdCOut() << "\n";

        // Print test in a compiler compatible format.
        CPPUNIT_NS::CompilerOutputter outputter( &result, CPPUNIT_NS::stdCOut() );
        outputter.write(); 

    // Uncomment this for XML output
    //    std::ofstream file( "tests.xml" );
    //    CPPUNIT_NS::XmlOutputter xml( &result, file );
    //    xml.setStyleSheet( "report.xsl" );
    //    xml.write();
    //    file.close();
    }
    catch ( std::invalid_argument &e )  // Test path not resolved
    {
        CPPUNIT_NS::stdCOut()  <<  "\n"  
                                <<  "ERROR: "  <<  e.what()
                                << "\n";
        return 0;
    }

    return result.wasSuccessful() ? 0 : 1;
    }

