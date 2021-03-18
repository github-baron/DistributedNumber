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
        
        pd1.CoverageInterval(50.,min,max);
        CPPUNIT_ASSERT_MESSAGE( "min should be -1: " + min.RawPrint(30), min == -1) ;
        CPPUNIT_ASSERT_MESSAGE( "max should be 1: " + max.RawPrint(30) + "\n1: " + CDigFloat(1.).RawPrint(30), max == 1) ;        
        
    }  

    void ProbabilityDensityDistribution_GenerateTriangle()
    {
        CProbabilityDensityDistribution pd1;
        
        pd1.TriangleDistribution(-1,1);
        
        // mid is zero        
        CPPUNIT_ASSERT_MESSAGE( "mean is zero", pd1.Mean() == 0) ;
        
        // integral should be 1
        CPPUNIT_ASSERT_MESSAGE( "integral is 1", pd1.OrigIntegral() == 1) ;        
        
    }

    void ProbabilityDensityDistribution_GenerateConstant()
    {
        CProbabilityDensityDistribution pd1;
        
        pd1.ConstDistribution(-2,1);
        
        // mid is zero        
        CPPUNIT_ASSERT_MESSAGE( "mean is zero", pd1.Mean() == -0.5) ;
        
        // integral should be 1
        CPPUNIT_ASSERT_MESSAGE( "integral is 3", pd1.OrigIntegral() == 3) ;        
        
    }

    void ProbabilityDensityDistribution_GenerateNormal()
    {
        CProbabilityDensityDistribution pd1;
        CDigFloat dfMu = -1, dfSigma = 0.01, dfMinVal = 0.00001;
        
        pd1.NormalDistribution(100,dfMu,dfSigma,dfMinVal);
        
        // check edges
        CPPUNIT_ASSERT_MESSAGE( "left is " + pd1.DistriValue( pd1.firstVariable()).Print(10) + " but should be " + dfMinVal.Print(), pd1.DistriValue( pd1.firstVariable())== dfMinVal) ;
        CPPUNIT_ASSERT_MESSAGE( "right is " +pd1.DistriValue( pd1.lastVariable()).Print(10) + " but should be " + dfMinVal.Print(), pd1.DistriValue( pd1.lastVariable())== dfMinVal) ;
        
        
        // mid is zero        
        CPPUNIT_ASSERT_MESSAGE( "mean is " + pd1.Mean().Print(10) + " but should be " + dfMu.Print(10), pd1.Mean() == dfMu);
    }
 

    void ProbabilityDensityDistribution_StatisticalParameters()
    {
        CProbabilityDensityDistribution pd1;
        
        pd1.TriangleDistribution(-4,-3);

        // check the assumptions: max value at mean is 0.5
        CPPUNIT_ASSERT_MESSAGE( "value at" + pd1.Mean().Print(10) + " is  " + pd1.DistriValue(pd1.Mean()).Print(10) + " but should be 2" , pd1.DistriValue(pd1.Mean()) == 2);
        
        // check the variance
        CPPUNIT_ASSERT_MESSAGE( "variance is " + pd1.Variance().Print(10) + " but should be " + (CDigFloat(4./96)).Print(10), pd1.Variance() == 4./96.);
        
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
        d1.TriangleDistribution(-3,-1);
        d2.TriangleDistribution(1,3);
        d3 = d1*d2;
        
        cout << endl;
        cout << "d3.Variance() " << d3.Variance().Print(10) << endl;
        cout << "d2.Variance() " << d2.Variance().Print(10) << endl;
        cout << "d1.Variance() " << d1.Variance().Print(10) << endl;
        cout << "d1.Mean() " << d1.Mean().Print(10) << endl;
        cout << "d2.Mean() " << d2.Mean().Print(10) << endl;
        cout << "calc " << (d2.Variance() + d1.Variance() + d1.Variance()*d2.Mean()*d2.Mean() + d2.Variance() * d1.Mean() * d1.Mean()).Print(10) << endl;

        CPPUNIT_ASSERT_MESSAGE( "variances should sum: " + d3.Variance().Print(10) + "=" +  d1.Variance().Print(10) + " + " + d2.Variance().Print(10) + " ( " + (d1.Variance() + d2.Variance()).Print(10) + ")"+ " + "  + (d1.Variance()*d2.Mean()*d2.Mean()).Print(10) + " + " + (d2.Variance() * d1.Mean() * d1.Mean()).Print(10), d3.Variance() == d2.Variance() + d1.Variance() + d1.Variance()*d2.Mean()*d2.Mean() + d2.Variance() * d1.Mean() * d1.Mean());

        // invert the multiplication: the result must be the same
        d3 = d2*d1;
        CDigFloat dfVarianceTheoretical = d2.Variance() + d1.Variance() + d1.Variance()*d2.Mean()*d2.Mean() + d2.Variance() * d1.Mean() * d1.Mean();
        cout << endl;
        cout << "d3.Variance() " << d3.Variance().Print(10) << endl;
        cout << "d2.Variance() " << d2.Variance().Print(10) << endl;
        cout << "d1.Variance() " << d1.Variance().Print(10) << endl;
        cout << "d1.Mean() " << d1.Mean().Print(10) << endl;
        cout << "d2.Mean() " << d2.Mean().Print(10) << endl;
        cout << "calc " << dfVarianceTheoretical.Print(10) << endl;

        CPPUNIT_ASSERT_MESSAGE( "variances should sum: " + d3.Variance().Print(10) + "=" +  dfVarianceTheoretical.Print(10), d3.Variance() != dfVarianceTheoretical);

  
    }
    
 
};
