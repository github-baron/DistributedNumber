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

#include "../DistributedNumber/ProbabilityDensityDistribution.h"

class CProbabilityDensityDistribution_Test : public CppUnit::TestFixture  {

private:

public:
    void setUp()
    {
    }

    void tearDown() 
    {
    }

    void ProbabilityDensityDistribution_Normalize()
    {
        CProbabilityDensityDistribution pd1;
        
        // normalization after adding the last element
        pd1.Add(-1,3);
        pd1.Add(1,3,true);
        CPPUNIT_ASSERT_MESSAGE( "norm should be 6: " + pd1.OrigIntegral().RawPrint(30), pd1.OrigIntegral() == 6) ; 
        CPPUNIT_ASSERT_MESSAGE( "distri value is normalized: 3 --> 0.5: " + pd1.DistriValue(0).RawPrint(30), pd1.DistriValue(0) == 0.5) ;

        // normalization after calling DistriValue
        pd1.Reset();
        pd1.Add(-1,3);
        pd1.Add(1,3);
        CPPUNIT_ASSERT_MESSAGE( "norm should be 6: " + pd1.OrigIntegral().RawPrint(30), pd1.OrigIntegral() == 6) ; 
        CPPUNIT_ASSERT_MESSAGE( "distri value is normalized: 3 --> 0.5: " + pd1.DistriValue(0).RawPrint(30), pd1.DistriValue(0) == 0.5) ;

        // normalization after calling AbsIntegral
        pd1.Reset();
        pd1.Add(-1,3);
        pd1.Add(1,3);
        CPPUNIT_ASSERT_MESSAGE( "norm should be 6: " + pd1.OrigIntegral().RawPrint(30), pd1.OrigIntegral() == 6) ; 
        CPPUNIT_ASSERT_MESSAGE( "distri normalized: Integral == 1: " + pd1.AbsIntegral().RawPrint(30), pd1.AbsIntegral() == 1) ;
        
        // normalization after adding a new point to the normalized distribution 
        pd1.Add(0,6);
        CPPUNIT_ASSERT_MESSAGE( "norm should be 9: " + pd1.OrigIntegral().RawPrint(30), pd1.OrigIntegral() == 9) ; 
        CPPUNIT_ASSERT_MESSAGE( "distri normalized: Integral == 1: " + pd1.AbsIntegral().RawPrint(30), pd1.AbsIntegral() == 1) ;
        CPPUNIT_ASSERT_MESSAGE( "distri value is normalized: 3 --> 3/norm:\nreceived " + pd1.DistriValue(0).RawPrint(30) + "\nexpected: " +CDigFloat( CDigFloat(6.)/pd1.OrigIntegral()).RawPrint(30) , pd1.DistriValue(0) == CDigFloat(6.)/pd1.OrigIntegral()) ;
    }  

    void ProbabilityDensityDistribution_Coverage()
    {
        CProbabilityDensityDistribution pd1;
        
        // normalization after adding the last element
        pd1.Add(-2,0);
        pd1.Add(-1,3);
        pd1.Add(0,0);
        pd1.Add(1,3);
        pd1.Add(2,0);
        CDigFloat min,max;
        
        pd1.CoverageFromTo(50.,pd1.Distribution().begin()->first,max);
        CPPUNIT_ASSERT_MESSAGE( "max should be 0: " + max.RawPrint(30) , max == 0) ; 

        pd1.CoverageFromTo(50.,prev(pd1.Distribution().end())->first,max,true);
        CPPUNIT_ASSERT_MESSAGE( "max should be 0: " + max.RawPrint(30) , max == 0) ; 
        
        pd1.CoverageIntervall(50.,min,max);
        CPPUNIT_ASSERT_MESSAGE( "min should be -1: " + min.RawPrint(30), min == -1) ;
        CPPUNIT_ASSERT_MESSAGE( "max should be 1: " + max.RawPrint(30) + "\n1: " + CDigFloat(1.).RawPrint(30), max == 1) ;        
        
       }  

    void ProbabilityDensityDistribution_Operators()
    {
        
        CProbabilityDensityDistribution d1,d2, d3;
        
        ////////////////////////////////////
        // squared distributions
        ////////////////////////////////////
        d1.Add(-3,1); // --> (-3, 0.5)
        d1.Add(-1,1); // --> (-1, 0.5)
        d2.Add(1,1);  // --> (+1, 0.5)
        d2.Add(3,1);  // --> (+3, 0.5)
        
        
        // summation results in 3 elements:
        // (-2,0.25) 
        // ( 0,0.25)
        // ( 0,0.25)
        // ( 2,0.25)
        // resulting in:
        // (-2,0.25) --> (-2, 0.25 / 1.5)
        // ( 0,0.50) --> ( 0, 0.50 / 1.5)
        // ( 2,0.25) --> (+2, 0.25 / 1.5)
        d3 = d1+d2;
       // get the neighbours of a variable -1.00001
//         CPPUNIT_ASSERT_MESSAGE("el must be (-2,0.06666666667) but is : " + d3.DistriValue(-2).RawPrint(3) + "\n distri1:\n " + d1.Print(10)+ "\n distri2:\n " + d2.Print(10)+ "\n distri:\n " + d3.Print(10) + "\n integral:\n " + d3.OrigIntegral().RawPrint(10), d3.DistriValue(-2)== 1./6);
//         CPPUNIT_ASSERT_MESSAGE("el must be (-1,0.2666667) but is : " + d3.DistriValue(-1).RawPrint(3) + "\n distri1:\n " + d1.Print(10)+ "\n distri2:\n " + d2.Print(10)+ "\n distri:\n " + d3.Print(10) + "\n integral:\n " + d3.OrigIntegral().RawPrint(10), d3.DistriValue(0)== 1./3);

        
        d3 = d1*d2;
        cout << endl << d3.Print(10) << endl;
        
        d1.Add(-2.5,1);
        d2.Add(1.5,1);
        d3 = d1*d2;
        cout << endl << d3.Print(10) << endl;
        
        d1.Add(-2.6,1);
        d2.Add(1.4,1);
        d3 = d1*d2;
        cout << endl << d3.Print(10) << endl;
      
        // check for same distribution but more points
        for(int i = 1; i<10; i++)
        {
            d1.Add(-3 + i*0.1,1);
            d2.Add(3 - i*0.1,1);
        }
        d3 = d1*d2;
        cout << endl << d3.Print(10) << endl;
      
        
        d1.Add(-2.95,1);
        d2.Add(1.05,1);
        d3 = d1*d2;
        cout << endl << d3.Print(10) << endl;
      
        
//       CPPUNIT_ASSERT_MESSAGE("el must be (-2,0.06666666667) but is : " + d3.DistriValue(-2).RawPrint(3) + "\n distri1:\n " + d1.Print(10)+ "\n distri2:\n " + d2.Print(10)+ "\n distri:\n " + d3.Print(10) + "\n integral:\n " + d3.OrigIntegral().RawPrint(10), d3.DistriValue(-2)== 0.066666666666666);
 
    }
    
 
};
