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
        CDigFloat dfMu = -1, dfSigma = 0.01;
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_GenerateNormal","generating normal distri with:");
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_GenerateNormal",string("dfMu = ") + dfMu.RawPrint(10));
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_GenerateNormal",string("dfSigma = ") + dfSigma.RawPrint(10));
        
        pd1.NormalDistribution(dfMu,dfSigma,20000);
        
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_GenerateNormal","comparing calculated values:");
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_GenerateNormal",string("mean           = ") + pd1.Mean().RawPrint(10));
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_GenerateNormal",string("sqrt(variance) = ") + sqrt(pd1.Variance()).RawPrint(10));
        
        // mid is zero        
        CPPUNIT_ASSERT_MESSAGE( "mean is " + pd1.Mean().Print(10) + " but should be " + dfMu.Print(10), pd1.Mean() == dfMu);
        // sigma
        CPPUNIT_ASSERT_MESSAGE( "sqrt(variance) is " + sqrt(pd1.Variance()).Print(10) + " but should be " + dfSigma.Print(10), sqrt(pd1.Variance()) == dfSigma);
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

    void ProbabilityDensityDistribution_Addition()
    {
        
        CProbabilityDensityDistribution d1,d2, d3;
        
        ////////////////////////////////////
        // squared distributions
        ////////////////////////////////////
        d1.Add(-3,1); // --> (-3, 0.5)
        d1.Add(-1,1); // --> (-1, 0.5)
        d2.Add(1,0);  // --> (+1, 0)
        d2.Add(3,1);  // --> (+3, 1)
        d1.NormalDistribution(2,1,1000);
        d2.NormalDistribution(0,10,1000);
        
        // work with high resolution to minimize integration errors
        d1.IntegrationSteps(1000);
                
        ///////////////////////////////////
        // Addition
        ///////////////////////////////////        
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", "starting addition operation");
        d3 = d1 + d2;
        
        //         Mean:
        //         -> d3.Mean == d1.Mean + d2.Mean
//         ostringstream oss;
        CDigFloat dfTheoreticalValue = d1.Mean() + d2.Mean();
//         oss << endl
//             << "****************************************" << endl
//             << "*********** checking addition **********" << endl 
//             << "****************************************" << endl
//             << "**** means sum up: " << endl
//             << "d3.Mean() == d1.Mean()+d2.Mean() : "   << endl 
//             << "d3.Mean = " << d3.Mean().RawPrint(10) << endl  
//             << "d1.Mean = " << d1.Mean().RawPrint(10) << endl 
//             << "d2.Mean = " << d2.Mean().RawPrint(10) << endl 
//             << "d1.Mean + d2.Mean =" << dfTheoreticalValue.RawPrint(10) << endl;
//         LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", oss.str());
//         CPPUNIT_ASSERT_MESSAGE( "means should sum: " + d3.Mean().Print(10) + "=" +  (d1.Mean()+d2.Mean()).RawPrint(10), d3.Mean() == d1.Mean()+d2.Mean());

        //         Variance:
        //         -> d3.Var == d1.Var + d2.Var + 2*d1.Covar(d2)
        CDigFloat dfVariance3 = d3.Variance();
        CDigFloat dfVariance2 = d2.Variance();
        CDigFloat dfVariance1 = d1.Variance();
        dfTheoreticalValue = dfVariance1 + dfVariance2 ;//+ 2* d1.Covariance(d2);
//         oss << endl 
//             << "**** variances sum up (if independent: covar == 0): " << endl
//             << "d3.Var == d1.Var + d2.Var + 2*d1.Covar(d2): " << endl 
//             << "d3.Var = " << d3.Variance().RawPrint(10) << endl 
//             << "d1.Var = " << d1.Variance().RawPrint(10) << endl
//             << "d2.Var = " << d2.Variance().RawPrint(10) << endl
//             << "d1.Covar(d2) = " << d1.Covariance(d2).RawPrint(10) << endl 
//             << "d1.var + d2.var + 2*d1.covar(d2) =" << dfTheoreticalValue.RawPrint(10) << endl;
//         LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", oss.str());
        CPPUNIT_ASSERT_MESSAGE( "means should sum: " + d3.Variance().Print(10) + "=" +  dfTheoreticalValue.RawPrint(10), d3.Variance() == dfTheoreticalValue);
    }
    void ProbabilityDensityDistribution_Subtraction()
    {
        
        CProbabilityDensityDistribution d1,d2, d3;
        
        ////////////////////////////////////
        // squared distributions
        ////////////////////////////////////
        d1.Add(-3,1); // --> (-3, 0.5)
        d1.Add(-1,1); // --> (-1, 0.5)
        d2.Add(1,0);  // --> (+1, 0)
        d2.Add(3,1);  // --> (+3, 1)
        d1.NormalDistribution(2,1,150);
        d2.NormalDistribution(0,10,150);
        
        // work with high resolution to minimize integration errors
        d1.IntegrationSteps(400);
                
        ///////////////////////////////////
        // Subtraction
        ///////////////////////////////////        
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", "starting subtraction operation");
        d3 = d1 - d2;
        
        // Mean:
        // -> d3.Mean == d1.Mean + d2.Mean
        ostringstream oss;
        CDigFloat dfTheoreticalValue = d1.Mean() - d2.Mean();
        oss << endl
            << "****************************************" << endl
            << "*********** checking addition **********" << endl 
            << "****************************************" << endl
            << "**** means sum up: " << endl
            << "d3.Mean() == d1.Mean()-d2.Mean() : "   << endl 
            << "d3.Mean = " << d3.Mean().RawPrint(10) << endl  
            << "d1.Mean = " << d1.Mean().RawPrint(10) << endl 
            << "d2.Mean = " << d2.Mean().RawPrint(10) << endl 
            << "d1.Mean - d2.Mean =" << dfTheoreticalValue.RawPrint(10) << endl;
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", oss.str());
        CPPUNIT_ASSERT_MESSAGE( "means should sum: " + d3.Mean().Print(10) + "=" +  dfTheoreticalValue.RawPrint(10), d3.Mean() == dfTheoreticalValue);
        
        // Variance:
        // -> d3.Var == d1.Var + d2.Var + 2*d1.Covar(d2)
        dfTheoreticalValue = d1.Variance() + d2.Variance();// - 2* d1.Covariance(d2);
        oss << endl 
            << "**** variances sum up (if independent: covar == 0): " << endl
            << "d3.Var == d1.Var + d2.Var - 2*d1.Covar(d2): " << endl 
            << "d3.Var = " << d3.Variance().RawPrint(10) << endl 
            << "d1.Var = " << d1.Variance().RawPrint(10) << endl
            << "d2.Var = " << d2.Variance().RawPrint(10) << endl
            << "d1.var + d2.var - 2*d1.covar(d2) =" << dfTheoreticalValue.RawPrint(10) << endl;
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", oss.str());
        CPPUNIT_ASSERT_MESSAGE( "means should sum: " + d3.Variance().Print(10) + "=" +  dfTheoreticalValue.RawPrint(10), d3.Variance() == dfTheoreticalValue);
    }
    void ProbabilityDensityDistribution_Multiplication()
    {
        
        CProbabilityDensityDistribution d1,d2, d3;
        
        ////////////////////////////////////
        // squared distributions
        ////////////////////////////////////
        d1.Add(-3,1); // --> (-3, 0.5)
        d1.Add(-1,1); // --> (-1, 0.5)
        d2.Add(1,0);  // --> (+1, 0)
        d2.Add(3,1);  // --> (+3, 1)
        d1.NormalDistribution(2,1,500);
        d2.NormalDistribution(0,10,500);
        
        // work with high resolution to minimize integration errors
        d1.IntegrationSteps(1000);      
        ///////////////////////////////////
        // Multiplication
        ///////////////////////////////////
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", "starting multiplication operation");
        d3 = d1 * d2;
        
        // tracing: comment in if you want the data of the distributions
//         LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators",string("factor1 distri :\n")+d1.Print(10,false));        
//         LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators",string("factor2 distri:\n")+d2.Print(10,false));        
//         LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators",string("resulting distri:\n")+d3.Print(10,false));
        // Variance:
        // -> d3.Mean == d2.Var* d1.Var + d1.Var*d2.Mean² + d2.Var*d1.Mean²
        //DEBUG       
        ostringstream oss; 
        CDigFloat dfTheoreticalValue=d1.Mean() * d2.Mean();
        oss << endl 
            << "**********************************************" << endl
            << "*********** checking multiplication **********" << endl 
            << "**********************************************" << endl
            << "**** means multiply: " << endl
            << "d3.Mean() == d1.Mean()*d2.Mean() : "   << endl 
            << "d3.Mean = " << d3.Mean().RawPrint(10) << endl  
            << "d1.Mean = " << d1.Mean().RawPrint(10) << endl 
            << "d2.Mean = " << d2.Mean().RawPrint(10) << endl 
            << "d1.Mean * d2.Mean =" << dfTheoreticalValue.RawPrint(10) << endl;
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", oss.str());
        CPPUNIT_ASSERT_MESSAGE( "means should multiply: " + d3.Mean().Print(10) + "=" +  dfTheoreticalValue.RawPrint(10), d3.Mean() == dfTheoreticalValue);
        
        dfTheoreticalValue=d2.Variance() * d1.Variance() + d1.Variance()*d2.Mean()*d2.Mean() + d2.Variance() * d1.Mean() * d1.Mean();
        oss << endl 
            << "**** variance: " << endl
            << "d3.Var == d2.Var* d1.Var + d1.Var*d2.Mean² + d2.Var*d1.Mean²: " << endl 
            << "d3.Var  = " << d3.Variance().RawPrint(10) << endl 
            << "d1.Mean = " << d1.Mean().RawPrint(10) << endl 
            << "d1.Var  = " << d1.Variance().RawPrint(10) << endl 
            << "d2.Mean = " << d2.Mean().RawPrint(10) << endl 
            << "d2.Var  = " << d2.Variance().RawPrint(10) << endl 
            << "d1.Var * d2.Var =" << (d1.Variance()*d2.Variance()).RawPrint(10) << endl 
            << "d2.Var * d2.Mean² =" << (d2.Variance()*d1.Mean()*d1.Mean()).RawPrint(10) << endl 
            << "d1.Var * d1.Mean² =" << (d1.Variance()*d2.Mean()*d2.Mean()).RawPrint(10) << endl 
            << "d2.Var* d1.Var + d1.Var*d2.Mean² + d2.Var*d1.Mean² =" << dfTheoreticalValue.RawPrint(10) << endl;
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", oss.str());
        CPPUNIT_ASSERT_MESSAGE( "variances should sum: " + d3.Variance().Print(10) + "=" +  dfTheoreticalValue.RawPrint(10), d3.Variance() == dfTheoreticalValue);


    }
    void ProbabilityDensityDistribution_Division()
    {
        
        CProbabilityDensityDistribution d1,d2, d3;
        
        ////////////////////////////////////
        // squared distributions
        ////////////////////////////////////
        d1.Add(-3,1); // --> (-3, 0.5)
        d1.Add(-1,1); // --> (-1, 0.5)
        d2.Add(1,0);  // --> (+1, 0)
        d2.Add(3,1);  // --> (+3, 1)
        d1.NormalDistribution(-10,1,200);
        d2.NormalDistribution(30,5,200); 
        d1.IntegrationSteps(800);
        ///////////////////////////////////
        // Division
        ///////////////////////////////////
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", "starting division operation");
        d3 = d1 / d2;
        CProbabilityDensityDistribution d2Inv; d2Inv = 1/d2;
        
        // tracing: comment in if you want the data of the distributions
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators",string("factor1 distri :\n")+d1.Print(10,false));        
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators",string("factor2 distri:\n")+d2.Print(10,false));        
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators",string("factor2 inv distri:\n")+d2Inv.Print(10,false));        
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators",string("resulting distri:\n")+d3.Print(10,false));
        // Variance:
        // -> d3.Mean == d2.Var* d1.Var + d1.Var*d2.Mean² + d2.Var*d1.Mean²
        //DEBUG     
        ostringstream oss;   
        CDigFloat dfTheoreticalValue=d1.Mean() * d2Inv.Mean();
        oss << endl 
            << "***************************** **********" << endl
            << "*********** checking Division **********" << endl 
            << "****************************************" << endl
            << "**** means : " << endl
            << "d3.Mean()  == d1.Mean()* (1/d2).Mean() : "   << endl 
            << "d3.Mean     = " << d3.Mean().RawPrint(10) << endl  
            << "d1.Mean     = " << d1.Mean().RawPrint(10) << endl 
            << "(1/d2).Mean = " << d2Inv.Mean().RawPrint(10) << endl 
            << "d2.Mean     = " << d2.Mean().RawPrint(10) << endl 
            << "d1.Mean * (1/d2).Mean =" << dfTheoreticalValue.RawPrint(10) << endl;
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", oss.str());
        CPPUNIT_ASSERT_MESSAGE( "means d1 and 1/d2 should multiply: " + d3.Mean().Print(10) + "=" +  dfTheoreticalValue.RawPrint(10), d3.Mean() == dfTheoreticalValue);
        
        dfTheoreticalValue=d2Inv.Variance() * d1.Variance() + d1.Variance()*d2Inv.Mean()*d2Inv.Mean() + d2Inv.Variance() * d1.Mean() * d1.Mean();
        oss << endl 
            << "**** variance: " << endl
            << "d3.Var                 == (1/d2).Var* d1.Var + d1.Var*(1/d2).Mean² + (1/d2).Var*d1.Mean²: " << endl 
            << "d3.Var                  = " << d3.Variance().RawPrint(10) << endl 
            << "d1.Mean                 = " << d1.Mean().RawPrint(10) << endl 
            << "d1.Var                  = " << d1.Variance().RawPrint(10) << endl 
            << "(1/d2).Mean             = " << d2Inv.Mean().RawPrint(10) << endl 
            << "(1/d2).Var              = " << d2Inv.Variance().RawPrint(10) << endl 
            << "d1.Var * d2.Var         =" << (d1.Variance()*d2Inv.Variance()).RawPrint(10) << endl 
            << "(1/d2).Var * d2.Mean²   =" << (d2Inv.Variance()*d1.Mean()*d1.Mean()).RawPrint(10) << endl 
            << "d1.Var * d1.Mean²       =" << (d1.Variance()*d2Inv.Mean()*d2Inv.Mean()).RawPrint(10) << endl 
            << "(1/d2).Var* d1.Var + d1.Var*(1/d2).Mean² + (1/d2).Var*d1.Mean² =" << dfTheoreticalValue.RawPrint(10) << endl;
        LOGTRACE("Distribution_Test::ProbabilityDensityDistribution_Operators", oss.str());
        CPPUNIT_ASSERT_MESSAGE( "variances check: " + d3.Variance().Print(10) + "=" +  dfTheoreticalValue.RawPrint(10), d3.Variance() == dfTheoreticalValue);

    }
    
    void ProbabilityDensityDistribution_GUMSumOfNormal()
    {
        
        CProbabilityDensityDistribution d1,d2,d3,d4, dRes;
        // expectation value = 0
        // sigma = 1;
        int nDistriPoints = 3000;
        int nIntegrationSteps=3000;
        d1.NormalDistribution(0,1,nDistriPoints );
        d1.IntegrationSteps(nIntegrationSteps);
        d2.NormalDistribution(0,1,nDistriPoints );
        d2.IntegrationSteps(nIntegrationSteps);
        d3.NormalDistribution(0,1,nDistriPoints );
        d3.IntegrationSteps(nIntegrationSteps);
        d4.NormalDistribution(0,1,nDistriPoints );
        d4.IntegrationSteps(nIntegrationSteps);
        
//         dTemp1 = d1 + d2;
//         dTemp2 = d3 + d4;
//         dRes = dTemp1 + dTemp2;
        dRes = d1 +d2 + d3 + d4;
        
        CDigFloat dfMean, dfVariance, dfFrom,dfTo;
        dfMean = dRes.Mean();
        dfVariance = dRes.Variance();
        dRes.CoverageInterval(95,dfFrom,dfTo);
        
        dfFrom.Precision(2); dfTo.Precision(2);

        cout << endl;

        cout << "Mean: " << dfMean.Print(10) << endl;
        cout << "Variance: " << dfVariance.Print(10) << endl;
        cout << "dfFrom: " << dfFrom.Print(10) << endl;
        cout << "dfTo: " << dfTo.Print(10) << endl;

        CPPUNIT_ASSERT_MESSAGE( "expected mean is 0 but is \n left: " + dfMean.Print(2), dfMean == 0);
        CPPUNIT_ASSERT_MESSAGE( "expected variance is 4 but is \n left: " + dfVariance.Print(2), dfVariance == 4);
        CPPUNIT_ASSERT_MESSAGE( "expected coverage interval is +-3.92 but is: [" + dfFrom.Print(2) + ", " +  dfTo.Print(2) + "]", dfTo == -3.92 && dfFrom == 3.92);

    }
    void ProbabilityDensityDistribution_GUMSumOfConstant()
    {
        
        CProbabilityDensityDistribution d1,dRes;
        // summing 4 eqaul constant distributions with 
        // <x> =0 and sigma = 1
        // calculation of width:
        // for a constant distri from -a to + a the integral factor is 1 / (2*a)
        // --> variance is 1/3*a² 
        // --> sigma is a/sqrt(3)
        // when sigma is one --> a = sqrt(3)
        d1.ConstDistribution(-sqrt(3),sqrt(3));
        d1.IntegrationSteps(2000);
        
//         dTemp1 = d1 + d2;
//         dTemp2 = d3 + d4;
//         dRes = dTemp1 + dTemp2;
        dRes = d1+d1+d1+d1;
        
        CDigFloat dfMean, dfVariance, dfFrom,dfTo;
        dfMean = dRes.Mean();
        dfVariance = dRes.Variance();
        dRes.CoverageInterval(95,dfFrom,dfTo);
        dfFrom.Precision(2); dfTo.Precision(2);

        cout << endl;

        cout << "Mean: " << dfMean.Print(10) << endl;
        cout << "Variance: " << dfVariance.Print(10) << endl;
        cout << "dfFrom: " << dfFrom.Print(10) << endl;
        cout << "dfTo: " << dfTo.Print(10) << endl;

        CPPUNIT_ASSERT_MESSAGE( "expected mean is 0 but is \n left: " + dfMean.Print(2), dfMean == 0);
        CPPUNIT_ASSERT_MESSAGE( "expected variance is 4 but is \n left: " + dfVariance.Print(2), dfVariance == 4);
        CPPUNIT_ASSERT_MESSAGE( "expected coverage interval is +-3.88 but is: [" + dfFrom.Print(2) + ", " +  dfTo.Print(2) + "]", dfFrom == -3.88 && dfTo == 3.88);

    }
        
    void ProbabilityDensityDistribution_GUMSumOfConstant_2()
    {
        
        CProbabilityDensityDistribution d1,d2,dRes;
        // summing 3 equal constant distributions with 
        // <x> =0 and sigma = 1
        // and one constant distri with:
        // <x> = 0 and sigma = 10
        // calculation of width:
        // expectation value = 0
        // sigma = 1;
        // for a constant distri from -a to + a the integral factor is 1 / (2*a)
        // --> variance is 1/3*a² 
        // --> sigma is a/sqrt(3)
        // when sigma is one --> a = sqrt(3)
        // when simga is ten --> a =10*sqrt(3)
        d1.ConstDistribution(-sqrt(3),sqrt(3));
        d2.ConstDistribution(-10*sqrt(3),10*sqrt(3));
        d1.IntegrationSteps(2000);
        d2.IntegrationSteps(2000);
        dRes = d1 + d1 + d1 + d2;
        
        CDigFloat dfMean, dfVariance, dfFrom,dfTo;
        dfMean = dRes.Mean();
        dfVariance = dRes.Variance();
        dRes.CoverageInterval(95,dfFrom,dfTo);
        dfFrom.Precision(0); dfTo.Precision(0);


        
        CPPUNIT_ASSERT_MESSAGE( "expected mean is 0 but is : " + dfMean.Print(2), dfMean == 0);
// //         CPPUNIT_ASSERT_MESSAGE( "expected variance is 103 but is : " + dfVariance.Print(2), dfVariance == 103);
        CPPUNIT_ASSERT_MESSAGE( "expected coverage interval is +-17. but is: [" + dfFrom.Print(2) + ", " +  dfTo.Print(2) + "]", dfFrom == -17 && dfTo == 17);

    } 
    void ProbabilityDensityDistribution_GUMMassCalibration_Approximation()
    {
        
        CProbabilityDensityDistribution MassCalibDelta,MassRefConv,  MassRefConvDelta; 
        CProbabilityDensityDistribution DensityAir, DensityMassCalib, DensityMassRef;
        CDigFloat dfDensityAirConv = 1.2; // kg/m³
        CDigFloat dfMassNominal = 100000;  // mg
        // model:
        // dm = (m_rc+dm_rc)*(1+(r_a-r_a0)*(1/r_w-1/r_r))-m_nom 
        // EXPLANATION:
        // dm:
        //     the deviation of the calibrated mass (m_w) from the nominal m_nom
        // m_rc:
        //     reference conventional mass (with density r_0) put on the other scale pan 
        //     (Gaussian)
        // dm_rc:
        //     weight for fine tuning reference put on the other scale pan
        //     (Gaussian)
        // r_a:
        //     density of air (at actual conditions)
        // r_a0:
        //     density of air (at conventional conditions)
        // r_w:
        //     density of calibration mass
        // r_r:
        //     density of reference mass
        // m_nom:
        //     nominal reference mass
        // summing 3 equal constant distributions with 
        // <x> =0 and sigma = 1
        // and one constant distri with:
        // <x> = 0 and sigma = 10
        // calculation of width:
        // expectation value = 0
        // sigma = 1;
        // for a constant distri from -a to + a the integral factor is 1 / (2*a)
        // --> variance is 1/3*a² 
        // --> sigma is a/sqrt(3)
        // when sigma is one --> a = sqrt(3)
        // when simga is ten --> a =10*sqrt(3)
        unsigned int uiSteps = 4000;
        MassRefConv.NormalDistribution(100000,0.05,uiSteps); // mg
        MassRefConvDelta.NormalDistribution(1.234,0.02,uiSteps); // mg
        DensityAir.ConstDistribution(1.2-0.1,1.2+0.1); // kg / m³
        DensityMassCalib.ConstDistribution(8000-1000,8000+1000); // kg / m³
        DensityMassRef.ConstDistribution(8000-50,8000+50); // kg / m³

        
        
        MassRefConv.IntegrationSteps(uiSteps);
        DensityMassRef.IntegrationSteps(uiSteps);
        DensityMassCalib.IntegrationSteps(uiSteps);
        DensityAir.IntegrationSteps(uiSteps);
        DensityMassCalib.IntegrationSteps(uiSteps);
        
//         cout << endl << "calculating: MassSum " << endl;
//         
//         CProbabilityDensityDistribution MassSum = MassRefConv + MassRefConvDelta;
//         
//         LOGTRACE("Distribution_Test::GUMMassCalibration","MassRefConv lin pars =  \n" + MassRefConv.PrintLinPars());
//         LOGTRACE("Distribution_Test::GUMMassCalibration","MassRefConvDelta lin pars =  \n" + MassRefConvDelta.PrintLinPars());
//         LOGTRACE("Distribution_Test::GUMMassCalibration","MassSum =  \n" + MassSum.Print(10,false));
//         for(auto iel: MassSum.Distribution())        
//           CPPUNIT_ASSERT_MESSAGE( "element is less zero " + iel.first.RawPrint(10) + ", " + iel.second.RawPrint(10), iel.second >= 0);
// 
//         
//         
        cout << "calculating: DensityAirDiff " << endl;
        CProbabilityDensityDistribution DensityAirDiff = DensityAir-dfDensityAirConv;
        cout << "calculating: inverse DensityMassCalib" << endl;
        CProbabilityDensityDistribution DensityMassCalibInv = 1./DensityMassCalib;
        cout << "calculating: inverse DensityMassRef" << endl;
        CProbabilityDensityDistribution DensityMassRefInv = 1./DensityMassRef;
        cout << "calculating: inverse DensityMassInvDiff" << endl;
        CProbabilityDensityDistribution DensityMassInvDiff = DensityMassCalibInv - DensityMassRefInv ;
//         cout << "calculating: DensityMassDiff " << endl;
//         CProbabilityDensityDistribution DensityMassDiff = DensityMassRef-DensityMassCalib;
//         cout << "calculating: DensityMassProd " << endl;
//         CProbabilityDensityDistribution DensityMassProd = DensityMassCalib*DensityMassRef;
//         
//         MassCalibDelta = DensityMassDiff / DensityMassProd;
//         cout << "multiplying with air density" << endl;
//         MassCalibDelta *=DensityAirDiff;
//         cout << "multiplying mass sum" << endl;
//         MassCalibDelta *=MassSum;
//         cout << "adding mass sum" << endl;
//         MassCalibDelta += MassSum;
//         cout << "subtracting reference mass " << endl;
//         MassCalibDelta -= dfMassNominal;
        
        
        MassCalibDelta = DensityAirDiff * DensityMassInvDiff;
        CProbabilityDensityDistribution MassSum = MassRefConv + MassRefConvDelta;
        MassCalibDelta *= MassSum;
        MassCalibDelta += MassSum;
        MassCalibDelta -= dfMassNominal;
        
        CDigFloat dfMean, dfVariance, dfFrom,dfTo;
        dfMean = MassCalibDelta.Mean(); dfMean.Precision(3);
        dfVariance = MassCalibDelta.Variance();
        MassCalibDelta.CoverageInterval(95,dfFrom,dfTo);
        dfFrom.Precision(4), dfTo.Precision(4);
        CDigFloat dfSigmaCalc =sqrt(dfVariance);
        dfSigmaCalc.Precision(4);
        
//         cout << endl << "-------------------- MassCalibDelta ----------------------" << endl << MassCalibDelta.Print(10,false) << endl;

        
        CPPUNIT_ASSERT_MESSAGE( "expected mean is 1.234 mg but is : " + dfMean.Print(2), dfMean == 1.234);
        CPPUNIT_ASSERT_MESSAGE( "expected variance GUM is 0.0754 but using 0.0763 (mine) but is : " + dfSigmaCalc.Print(3), dfSigmaCalc == 0.0763 /*0.0754 : value from GUM*/);
        CPPUNIT_ASSERT_MESSAGE( "expected coverage interval is [1.0834 / 1.3825] but is: [" + dfFrom.Print(2) + ", " +  dfTo.Print(2) + "]", dfFrom == 1.0834 && dfTo == 1.3825);

    } 
        
  void ProbabilityDensityDistribution_GUMMassCalibration_Exact()
    {
        
        CProbabilityDensityDistribution MassCalibDelta,MassRefConv,  MassRefConvDelta; 
        CProbabilityDensityDistribution DensityAir, DensityMassCalib, DensityMassRef;
        CDigFloat dfDensityAirConv = 1.2; // kg/m³
        CDigFloat dfMassNominal = 100000;  // mg
        // model:
        // dm = (m_rc+dm_rc)*(1- r_a/r_r)/(1-r_a0/r_r)*(1-r_a0/r_w)/(1-r_a/r_w)-m_nom 
        // EXPLANATION:
        // dm:
        //     the deviation of the calibrated mass (m_w) from the nominal m_nom
        // m_rc:
        //     reference conventional mass (with density r_0) put on the other scale pan 
        //     (Gaussian)
        // dm_rc:
        //     weight for fine tuning reference put on the other scale pan
        //     (Gaussian)
        // r_a:
        //     density of air (at actual conditions)
        // r_a0:
        //     density of air (at conventional conditions)
        // r_w:
        //     density of calibration mass
        // r_r:
        //     density of reference mass
        // m_nom:
        //     nominal reference mass
        // summing 3 equal constant distributions with 
        // <x> =0 and sigma = 1
        // and one constant distri with:
        // <x> = 0 and sigma = 10
        // calculation of width:
        // expectation value = 0
        // sigma = 1;
        // for a constant distri from -a to + a the integral factor is 1 / (2*a)
        // --> variance is 1/3*a² 
        // --> sigma is a/sqrt(3)
        // when sigma is one --> a = sqrt(3)
        // when simga is ten --> a =10*sqrt(3)
        unsigned int uiSteps = 500;
        MassRefConv.NormalDistribution(100000,0.05,uiSteps); // mg
        MassRefConvDelta.NormalDistribution(1.234,0.02,uiSteps); // mg
        
//         DensityAir.ConstDistribution(1.2-0.1,1.2+0.1); // kg / m³
//         DensityMassCalib.ConstDistribution(8000-1000,8000+1000); // kg / m³
//         DensityMassRef.ConstDistribution(8000-50,8000+50); // kg / m³        
        
         MassRefConv.IntegrationSteps(uiSteps);
//         MassRefConvDelta .IntegrationSteps(uiSteps);
//         DensityMassRef.IntegrationSteps(uiSteps);
//         DensityMassCalib.IntegrationSteps(uiSteps);
//         DensityAir.IntegrationSteps(uiSteps);
//         DensityMassCalib.IntegrationSteps(uiSteps);
        
        // calculate the mass sum of reference mass and it compensation weight delta reference
        LOGTRACE("Distribution_Test::GUMMassCalibration", string("calculating MassSum with:"));
        LOGTRACE("Distribution_Test::GUMMassCalibration", string("integration steps for MassRefConv \n") + MassRefConv.PrintMetaInfo());
        LOGTRACE("Distribution_Test::GUMMassCalibration", string("integration steps for MassRefConvDelta \n") + MassRefConvDelta.PrintMetaInfo());
        CProbabilityDensityDistribution MassSum = MassRefConv + MassRefConvDelta;
    
       // calculate all the density relations of the form " 1 - r_a0/r_r"
       LOGTRACE("Distribution_Test::GUMMassCalibration", "calculating DensityRelAirMassRef ");
       CProbabilityDensityDistribution DensityRelAirMassRef = 1. - (DensityAir / DensityMassRef);
//         cout << endl << "-------------------- DensityRelAirMassRef ----------------------" << endl << DensityRelAirMassRef.Print(10) << endl;
//        
//        LOGTRACE("Distribution_Test::GUMMassCalibration", "calculating DensityRelAirConvMassRef");
//        CProbabilityDensityDistribution DensityRelAirConvMassRef = 1.- (dfDensityAirConv / DensityMassRef);
////         cout << endl << "-------------------- DensityRelAirConvMassRef ----------------------" << endl << DensityRelAirConvMassRef.Print(10,false) << endl;
//        
//        LOGTRACE("Distribution_Test::GUMMassCalibration", "calculating DensityRelAirMassCalib");
//        CProbabilityDensityDistribution DensityRelAirMassCalib = 1. - (DensityAir / DensityMassCalib);
////         cout << endl << "-------------------- DensityRelAirMassCalib ----------------------" << endl << DensityRelAirMassCalib.Print(10,false) << endl;
//        
//        LOGTRACE("Distribution_Test::GUMMassCalibration", "calculating DensityRelAirConvMassCalib");
//        CProbabilityDensityDistribution DensityRelAirConvMassCalib = 1. - (dfDensityAirConv / DensityMassCalib);
//            
//        // now calculate the formula         
//        LOGTRACE("Distribution_Test::GUMMassCalibration", "calculating MassCalibDelta");
//        MassCalibDelta = MassSum;// * DensityRelAirMassRef / DensityRelAirConvMassRef * DensityRelAirConvMassCalib / DensityRelAirMassCalib - dfMassNominal ;
////         cout << endl << "-------------------- MassCalibDelta ----------------------" << endl << MassCalibDelta.Print(10) << endl;
//        MassCalibDelta *= DensityRelAirMassRef;  
//        cout << endl << "-------------------- MassCalibDelta ----------------------" << endl << MassCalibDelta.Print(10) << endl;
////         MassCalibDelta /= DensityRelAirConvMassRef; 
////         MassCalibDelta *= DensityRelAirConvMassCalib; 
////         MassCalibDelta /= DensityRelAirMassCalib;
////         MassCalibDelta -= dfMassNominal;
//                
//        CDigFloat dfMean, dfVariance, dfFrom,dfTo;
//        dfMean = MassCalibDelta.Mean(); dfMean.Precision(3);
//        dfVariance = MassCalibDelta.Variance();
//        MassCalibDelta.CoverageInterval(95,dfFrom,dfTo);
//        dfFrom.Precision(4), dfTo.Precision(4);
//        CDigFloat dfSigmaCalc =sqrt(dfVariance);
//        dfSigmaCalc.Precision(4);
//        
////         cout << endl << "-------------------- MassCalibDelta ----------------------" << endl << MassCalibDelta.Print(10,false) << endl;
//
//        
//        CPPUNIT_ASSERT_MESSAGE( "expected mean is 1.234 mg but is : " + dfMean.Print(2), dfMean == 1.234);
//        CPPUNIT_ASSERT_MESSAGE( "expected variance GUM is 0.0754 but using 0.0763 (mine) but is : " + dfSigmaCalc.Print(3), dfSigmaCalc == 0.0763 /*0.0754 : value from GUM*/);
//        CPPUNIT_ASSERT_MESSAGE( "expected coverage interval is [1.0834 / 1.3825] but is: [" + dfFrom.Print(2) + ", " +  dfTo.Print(2) + "]", dfFrom == 1.0834 && dfTo == 1.3825);
//
    } 
  
    
};
