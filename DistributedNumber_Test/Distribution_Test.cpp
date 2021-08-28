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

#include "../DistributedNumber/Distribution.h"

class CDistribution_Test : public CppUnit::TestFixture  {

private:

public:
    void setUp()
    {
    }

    void tearDown() 
    {
    }

    void DistributionBuildUp()
    {
        CDistribution d1;
        // build up distribution
        d1.Add(-1,3);
        d1.Add(1,3);
        
        // check size
        CPPUNIT_ASSERT_MESSAGE( "size of distribution should be 2 but is  " + to_string(d1.Distribution().size()) , d1.Distribution().size() == 2) ;
       
        // now check the single elements:
        MapDFDFType::const_iterator iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-1] should be 3 but is : "+ Print(iel,10) , iel->first == -1 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : "+ Print(iel,10) , iel->first == 1 && iel->second == 3) ;
        
        // newly add an element in between
        d1.Add(0.3,6);
       
        MapDFDFType::const_iterator iel2 = d1.Distribution().begin(); 
        CPPUNIT_ASSERT_MESSAGE( "d[-1] should be 3 but is : " + iel2 ->first.RawPrint(30) + "\n" + iel2 ->second.RawPrint(30) , iel2 ->first == -1 && iel2 ->second == 3) ;
        iel2++;
        if(iel2 != d1.Distribution().end())
            CPPUNIT_ASSERT_MESSAGE( "d[0.3] should be 6 but is : " + iel2 ->first.RawPrint(30) + "\n" + iel2 ->second.RawPrint(30) + "\ndistri:\n" + d1.Print(10) , iel2 ->first == 0.3 && iel2 ->second == 6) ;
        
        iel2++;
        if(iel2 != d1.Distribution().end())
            CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : " + iel2 ->first.RawPrint(30) + "\n" + iel2 ->second.RawPrint(30) , iel2 ->first == 1 && iel2 ->second == 3);
        
    }
    void DistributionOperations()
    {
        CDistribution d1,d2;
        d1.Add(-1,3);
        d1.Add( 1,3);
        d2=d1;
        
        // operator *=:
        d1 *= 2;        
        MapDFDFType::const_iterator iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-2] should be 3 but is : " + Print(iel,10) , iel->first == -2 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[2] should be 3 but is : "+ Print(iel,10) , iel->first == 2 && iel->second == 3) ;
        
        // operator *:
        d1 = d2*2;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-2] should be 3 but is : " + Print(iel,10) , iel->first == -2 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[2] should be 3 but is : "+ Print(iel,10) , iel->first == 2 && iel->second == 3) ;
        // needed for next test
        d2 = d1;
        
        // operator /=:
        d1 /= 2;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-1] should be 3 but is : "+ Print(iel,10) , iel->first == -1 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : "+ Print(iel,10) , iel->first == 1 && iel->second == 3) ;
        
        // operator /:
        d1 = d2 / 2;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-1] should be 3 but is : "+ Print(iel,10) , iel->first == -1 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : "+ Print(iel,10) , iel->first == 1 && iel->second == 3) ;
        // needed for next test
        d2 = d1;
        
        // operator +=:
        d1 += 1.5;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[0.5] should be 3 but is : "+ Print(iel,10) , iel->first == 0.5 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[2.5] should be 3 but is : "+ Print(iel,10) , iel->first == 2.5 && iel->second == 3) ;
        
        // operator +:
        d1 = d2 + 1.5;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[0.5] should be 3 but is : "+ Print(iel,10) , iel->first == 0.5 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[2.5] should be 3 but is : "+ Print(iel,10) , iel->first == 2.5 && iel->second == 3) ;
        // needed for next test
        d2 = d1;
        
        // operator -=:
        d1 -= 1.5;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-1] should be 3 but is : "+ Print(iel,10) , iel->first == -1 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : "+ Print(iel,10) , iel->first == 1 && iel->second == 3) ;
        
        // operator -:
        d1 = d2 - 1.5;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-1] should be 3 but is : "+ Print(iel,10) , iel->first == -1 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : "+ Print(iel,10) , iel->first == 1 && iel->second == 3) ;
        
    }
    
    void Distribution_ExternalOperations()
    {
        CDistribution d1, d2;
        d2.Add(-1,1);
        d2.Add( 1,2);
        
        // operator +:
        d1 = 2 + d2;        
        MapDFDFType::const_iterator iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : " + Print(iel,10) , iel->first == 1 && iel->second == 1) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[3] should be 3 but is : "+ Print(iel,10) , iel->first == 3 && iel->second == 2) ;
        
        // operator -:
        d1 = 2 - d2;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : "+ Print(iel,10) , iel->first == 1 && iel->second == 2) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[3] should be 3 but is : "+ Print(iel,10) , iel->first == 3 && iel->second == 1) ;
        
        // operator *:
        d1 = 2 * d2;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-2] should be 3 but is : "+ Print(iel,10) , iel->first == -2 && iel->second == 1) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[2] should be 3 but is : "+ Print(iel,10) , iel->first == 2 && iel->second == 2) ;
        
        // operator * negative number:
        d1 = -2 * d2;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-2] should be 3 but is : "+ Print(iel,10) , iel->first == -2 && iel->second == 2) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[2] should be 3 but is : "+ Print(iel,10) , iel->first == 2 && iel->second == 1) ;
        
        // operator /:
        d1 = 2 / d2;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-2] should be 3 but is : "+ Print(iel,10) , iel->first == -2 && iel->second == 1) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[2] should be 3 but is : "+ Print(iel,10) , iel->first == 2 && iel->second == 2) ;
        
        // operator / with negative number:
        d1 = -2 / d2;        
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-2] should be 3 but is : "+ Print(iel,10) , iel->first == -2 && iel->second == 2) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[2] should be 3 but is : "+ Print(iel,10) , iel->first == 2 && iel->second == 1) ;
    }
    
    void DistributionFunctionOnVariable()
    {
        CDistribution d1;
        d1.Add(-1,3);
        d1.Add(1,3);
        
        // shifting back
        d1.Shift(-3);
        MapDFDFType::const_iterator iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-1] should be 3 but is : "+ Print(iel,10) , iel->first == -4 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : "+ Print(iel,10) , iel->first == -2 && iel->second == 3) ;

        // shifting foreward
        d1.Shift(3);
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-1] should be 3 but is : "+ Print(iel,10) , iel->first == -1 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : "+ Print(iel,10) , iel->first == 1 && iel->second == 3) ;
        
        // scaling: inverting never works !!
        d1.Scale(-1);
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-1] should be 3 but is : "+ Print(iel,10) , iel->first == -1 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : "+ Print(iel,10) , iel->first == 1 && iel->second == 3) ;
        
         // scaling: inverting never works but the scale is applied
        d1.Scale(-2);
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-1] should be 3 but is : "+ Print(iel,10) , iel->first == -2 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : "+ Print(iel,10) , iel->first == 2 && iel->second == 3) ;
 
        // scaling: revert: apply 1/ factor
        d1.Scale(0.5);
        iel = d1.Distribution().begin();
        CPPUNIT_ASSERT_MESSAGE( "d[-1] should be 3 but is : "+ Print(iel,10) , iel->first == -1 && iel->second == 3) ;
        iel++;
        CPPUNIT_ASSERT_MESSAGE( "d[1] should be 3 but is : "+ Print(iel,10) , iel->first == 1 && iel->second == 3) ;

 
    }    

    void Distribution_DistriValue()
    {
        
        CDistribution d1;
        // build up distribution
        d1.Add(-1,3);
        d1.Add(1,3);
        
        // simple check: get value in between : must be 3
        CDigFloat dfValue = d1.DistriValue(-0.1);
        CPPUNIT_ASSERT_MESSAGE("value must be 3 but is : " + dfValue.RawPrint(5), dfValue == 3);
        
        // now add a new value and check again
        d1.Add(0.3, 6);
        dfValue = d1.DistriValue(-0.1);
        CDigFloat x1,x2,y1,y2;
        x1=-1;y1=3;x2=0.3;y2 = 6;
        CDigFloat dfExpected = (y2-y1)/(x2-x1)*(-0.1-x1)+y1;
        CPPUNIT_ASSERT_MESSAGE("value must be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == dfExpected);
        

        
        ///////////////////////////////////
        // adding new interpolated 
        // values --> must not change
        // the interpolated value at -0.1
        //////////////////////////////////
        
        // the value must not change if a new interpolated value is added within this interval
        d1.Add(0.1, d1.DistriValue(0.1));
        CPPUNIT_ASSERT_MESSAGE("after adding interpolated value must still be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == d1.DistriValue(-0.1));
        
        CPPUNIT_ASSERT_MESSAGE("interpolated value at +0.3 should be 6 but is " + d1.DistriValue(0.3).RawPrint(30)  + "\ndistri:\n" + d1.Print(10), d1.DistriValue(0.3) == 6) ;

        
        // the value must not change if a new interpolated value is added within this interval
        d1.Add(-0.1, d1.DistriValue(-0.1));
        CPPUNIT_ASSERT_MESSAGE("after adding interpolated value must still be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == d1.DistriValue(-0.1));

        // the value must not change if a new interpolated value is added within this interval
        d1.Add(-0.2, d1.DistriValue(-0.2));
        CPPUNIT_ASSERT_MESSAGE("after adding interpolated value must still be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == d1.DistriValue(-0.1));

        // the value must not change if a new interpolated value is added within this interval
        d1.Add(-0.3, d1.DistriValue(-0.3));
        CPPUNIT_ASSERT_MESSAGE("after adding interpolated value must still be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == d1.DistriValue(-0.1));

        // the value must not change if a new interpolated value is added within this interval
        d1.Add(-0.4, d1.DistriValue(-0.4));
        CPPUNIT_ASSERT_MESSAGE("after adding interpolated value must still be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == d1.DistriValue(-0.1));

        // the value must not change if a new interpolated value is added within this interval
        d1.Add(-0.5, d1.DistriValue(-0.5));
        CPPUNIT_ASSERT_MESSAGE("after adding interpolated value must still be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == d1.DistriValue(-0.1));

        // the value must not change if a new interpolated value is added within this interval
        d1.Add(-0.6, d1.DistriValue(-0.6));
        CPPUNIT_ASSERT_MESSAGE("after adding interpolated value must still be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == d1.DistriValue(-0.1));

        // the value must not change if a new interpolated value is added within this interval
        d1.Add(0.2, d1.DistriValue(0.2));
        CPPUNIT_ASSERT_MESSAGE("after adding interpolated value must still be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == d1.DistriValue(-0.1));

        // the value must not change if a new interpolated value is added within this interval
        d1.Add(0.3, d1.DistriValue(0.3));
        CPPUNIT_ASSERT_MESSAGE("after adding interpolated value must still be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == d1.DistriValue(-0.1));

        // the value must not change if a new interpolated value is added within this interval
        d1.Add(0.5, d1.DistriValue(0.5));
        CPPUNIT_ASSERT_MESSAGE("after adding interpolated value must still be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == d1.DistriValue(-0.1));

        // the value must not change if a new interpolated value is added within this interval
        d1.Add(0.6, d1.DistriValue(0.6)); 
        CPPUNIT_ASSERT_MESSAGE("after adding interpolated value must still be " + dfExpected.RawPrint(5)  + " but is : " + dfValue.RawPrint(5), dfValue == d1.DistriValue(-0.1));
        
    }
    
    void Distribution_GetInterval()
    {
        CDistribution d1;
        // build up distribution
        d1.Add(-1,3);
        d1.Add(1,3);
        
        // get the neighbours of a variable -1.00001: out of bounds left
        MapDFDFType::const_iterator min, max;
        d1.GetInterval(-1.0001,min,max);
        CPPUNIT_ASSERT_MESSAGE("min must be (-1,3) but is : " + min->first.RawPrint(3) + "," + min->second.RawPrint(3), min->first == -1 && min->second == 3);
        CPPUNIT_ASSERT_MESSAGE("max must be (-1,3) but is : " + max->first.RawPrint(3) + "," + max->second.RawPrint(3), max->first ==  -1 && max->second == 3);
        
        // get the neighbours of a variable -1
        d1.GetInterval(-1,min,max);
        CPPUNIT_ASSERT_MESSAGE("min must be (-1,3) but is : " + min->first.RawPrint(3) + "," + min->second.RawPrint(3), min->first == -1 && min->second == 3);
        CPPUNIT_ASSERT_MESSAGE("max must be (1,3) but is : " + max->first.RawPrint(3) + "," + max->second.RawPrint(3), max->first == 1 && max->second == 3);
        
        // get the neighbours of a variable 0
        d1.GetInterval(0,min,max);
        CPPUNIT_ASSERT_MESSAGE("min must be (-1,3) but is : " + min->first.RawPrint(3) + "," + min->second.RawPrint(3), min->first == -1 && min->second == 3);
        CPPUNIT_ASSERT_MESSAGE("max must be (1,3) but is : " + max->first.RawPrint(3) + "," + max->second.RawPrint(3), max->first ==  1 && max->second == 3);
        
        // get the neighbours of a variable 1
        d1.GetInterval(1,min,max);
        CPPUNIT_ASSERT_MESSAGE("min must be (-1,3) but is : " + min->first.RawPrint(3) + "," + min->second.RawPrint(3), min->first == -1 && min->second == 3);
        CPPUNIT_ASSERT_MESSAGE("max must be (1,3) but is : " + max->first.RawPrint(3) + "," + max->second.RawPrint(3), max->first == 1 && max->second == 3);
        
        // get the neighbours of a variable 1.0000001: out of bounds right
        d1.GetInterval(1.0000001,min,max);
        CPPUNIT_ASSERT_MESSAGE("min and max are the end of the distribution", min == d1.Distribution().end() && max == d1.Distribution().end());
                ///////////////////////////////////////////////
        // change distribution: add (0.3; 6)
        ///////////////////////////////////////////////
        d1.Add(0.3,6);
        
        // get the neighbours of a variable -1.00001
        d1.GetInterval(-1.0001,min,max);
        CPPUNIT_ASSERT_MESSAGE("min must be (-1,3) but is : " + min->first.RawPrint(3) + "," + min->second.RawPrint(3), min->first == -1 && min->second == 3);
        CPPUNIT_ASSERT_MESSAGE("max must be (-1,3) but is : " + max->first.RawPrint(3) + "," + max->second.RawPrint(3), max->first ==  -1 && max->second == 3);
        
        // get the neighbours of a variable -1
        d1.GetInterval(-1,min,max);
        CPPUNIT_ASSERT_MESSAGE("min must be (-1,3) but is : " + min->first.RawPrint(3) + "," + min->second.RawPrint(3), min->first == -1 && min->second == 3);
        CPPUNIT_ASSERT_MESSAGE("max must be (0.3,6) but is : " + max->first.RawPrint(3) + "," + max->second.RawPrint(3), max->first == 0.3 && max->second == 6);

        // get the neighbours of a variable 0
        d1.GetInterval(0,min,max);
        CPPUNIT_ASSERT_MESSAGE("min must be (-1,3) but is : " + min->first.RawPrint(3) + "," + min->second.RawPrint(3), min->first == -1 && min->second == 3);
        CPPUNIT_ASSERT_MESSAGE("max must be (0.3,6) but is : " + max->first.RawPrint(3) + "," + max->second.RawPrint(3), max->first == 0.3 && max->second == 6);
            
        // get the neighbours of a variable 0.3
        d1.GetInterval(0.3,min,max);
        CPPUNIT_ASSERT_MESSAGE("min must be (-1,3) but is : " + min->first.RawPrint(3) + "," + min->second.RawPrint(3), min->first == -1 && min->second == 3);
        CPPUNIT_ASSERT_MESSAGE("max must be (0.3,6) but is : " + max->first.RawPrint(3) + "," + max->second.RawPrint(3), max->first == 0.3 && max->second == 6);

         // get the neighbours of a variable 0.5
         d1.GetInterval(0.5,min,max);
         CPPUNIT_ASSERT_MESSAGE("min must be (0.3,6) but is : " + min->first.RawPrint(3) + "," + min->second.RawPrint(3), min->first == 0.3 && min->second == 6);
         CPPUNIT_ASSERT_MESSAGE("max must be (1,3) but is : " + max->first.RawPrint(3) + "," + max->second.RawPrint(3), max->first ==  1 && max->second == 3);
            
        // get the neighbours of a variable 1
        d1.GetInterval(1,min,max);
        CPPUNIT_ASSERT_MESSAGE("min must be (0.3,6) but is : " + min->first.RawPrint(3) + "," + min->second.RawPrint(3), min->first == 0.3 && min->second == 6);
        CPPUNIT_ASSERT_MESSAGE("max must be (1,3) but is : " + max->first.RawPrint(3) + "," + max->second.RawPrint(3), max->first == 1 && max->second == 3);

        // get the neighbours of a variable 1.0000001
        d1.GetInterval(1.0000001,min,max);
        CPPUNIT_ASSERT_MESSAGE("min and max are the end of the distribution", min == d1.Distribution().end() && max == d1.Distribution().end());

    }
    
    void Distribution_AbsIntegral()
    {      
        
        CDistribution d1;
        // build up distribution
        d1.Add(-1,3);
        d1.Add(1,3);
        
 
        // calculate the absolute integral
        CDigFloat dfIntegral = d1.AbsIntegral();
        CPPUNIT_ASSERT_MESSAGE( "Integral should be 6 but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 6) ;
              
        
        ////////////////////////////////////
        // change distri: adding 0.3 / 6
        ///////////////////////////////////
        d1.Add(0.3,6);
        
        // calculate the absolute integral with new distri
        dfIntegral = d1.AbsIntegral();
        CPPUNIT_ASSERT_MESSAGE( "Integral should be 9 but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 9) ;
        
        // calculate arbitrary integrals of complete distribution elements (i.e. in our case the intervalls
        // [-1, 1]
        // [-1, 0.3]
        // [0.3, 1]
        dfIntegral = d1.AbsIntegral(-1,1);
        CPPUNIT_ASSERT_MESSAGE( "Integral should be 9 but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 9) ;
        dfIntegral = d1.AbsIntegral(-1,0.3);
        CPPUNIT_ASSERT_MESSAGE( "Integral should be 5.85 but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 5.85) ;
        dfIntegral = d1.AbsIntegral(0.3, 1);
        CPPUNIT_ASSERT_MESSAGE( "Integral should be 3.15 but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 3.15) ;

        // halfen the complete distribution intervalls:
        // [-1, 0.3] --> [-0.675, -0.025]
        // [0.3 , 1] --> [0.475, 0.825]
        dfIntegral = d1.AbsIntegral(-0.675,-0.025);
        CPPUNIT_ASSERT_MESSAGE( "Integral should be " + to_string(5.85/2.) + " but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 5.85/2.) ;
        dfIntegral = d1.AbsIntegral(0.475,1-0.175);
        CPPUNIT_ASSERT_MESSAGE( "Integral should be " + to_string(3.15 /2.) + "but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 3.15 / 2.) ;
        
        ////////////////////////////////////////////////
        // adding a new distribution element into the
        // linear interpolation / integration path
        // must not change any of the former integration 
        // and interpolation results
        ////////////////////////////////////////////////
        d1.Add(-0.1, d1.DistriValue(-0.1));
//         d1.Add(-0.2, d1.DistriValue(-0.2));
//         d1.Add(-0.3, d1.DistriValue(-0.3));
//         d1.Add(-0.4, d1.DistriValue(-0.4));
//         d1.Add(-0.5, d1.DistriValue(-0.5));
//         d1.Add(-0.6, d1.DistriValue(-0.6));
//         
//         
//         d1.Add(0.1, d1.DistriValue(0.1));
//         d1.Add(0.2, d1.DistriValue(0.2));
//         d1.Add(0.3, d1.DistriValue(0.3));
//         d1.Add(0.5, d1.DistriValue(0.5));
//         d1.Add(0.6, d1.DistriValue(0.6)); 
        
        // calculate the absolute integral
        dfIntegral = d1.AbsIntegral();
        CPPUNIT_ASSERT_MESSAGE( "Integral should be 9 but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 9) ;
        
        // interpolate value at fixed points: -1 / 0.3 / 1
        CDigFloat dfIntVal = d1.DistriValue(-1);
        CPPUNIT_ASSERT_MESSAGE("interpolated value at -1 should be 3 but is " + dfIntVal.RawPrint(30)  + "\ndistri:\n" + d1.Print(10), dfIntVal == 3) ;
        dfIntVal = d1.DistriValue(0.3);
        dfIntVal = d1.DistriValue(1);
        CPPUNIT_ASSERT_MESSAGE("interpolated value at 1 should be 3 but is " + dfIntVal.RawPrint(30)  + "\ndistri:\n" + d1.Print(10), dfIntVal == 3) ;
        
        // now do the interpolation for outside points:
        dfIntVal = d1.DistriValue(-2);
        CPPUNIT_ASSERT_MESSAGE("interpolated value at -2 should be  0 but is " + dfIntVal.RawPrint(30)  + "\ndistri:\n" + d1.Print(10), dfIntVal == 0);
        dfIntVal = d1.DistriValue(+1.00000001);
        CPPUNIT_ASSERT_MESSAGE("interpolated value at 1.00000001 should be  0 but is " + dfIntVal.RawPrint(30)  + "\ndistri:\n" + d1.Print(10), dfIntVal == 0);
 
        // now interpolate points within the distribution :
        dfIntVal = d1.DistriValue(-0.5);
        CPPUNIT_ASSERT_MESSAGE("interpolated value at 0.5 should be " + to_string(double(1.5/1.3+3.))+ " but is " + dfIntVal.RawPrint(30)  + "\ndistri:\n" + d1.Print(10), dfIntVal == double(3./1.3 * 0.5 +3.));
        dfIntVal = d1.DistriValue(0.2);
        CPPUNIT_ASSERT_MESSAGE("interpolated value at 0.5 should be " + to_string(double(3./1.3 * (0.2+1) +3.))+ " but is " + dfIntVal.RawPrint(30)  + "\ndistri:\n" + d1.Print(10), dfIntVal == double(3./1.3 * (0.2+1) +3.));
        dfIntVal = d1.DistriValue(0.35);
        CPPUNIT_ASSERT_MESSAGE("interpolated value at 0.5 should be " + to_string(double(-3./0.7 * (0.05) +6))+ " but is " + dfIntVal.RawPrint(30)  + "\ndistri:\n" + d1.Print(10), dfIntVal == double(-3./0.7 * (0.05) +6));
        dfIntVal = d1.DistriValue(0.98);
        CPPUNIT_ASSERT_MESSAGE("interpolated value at 0.5 should be " + to_string(double(-3./0.7 * (0.98-0.3) +6))+ " but is " + dfIntVal.RawPrint(30)  + "\ndistri:\n" + d1.Print(10), dfIntVal == double(-3./0.7 * (0.98-0.3) +6));
        
        // calculate arbitrary integrals of complete distribution elements (i.e. in our case the intervalls
        // [-1, 1]
        // [-1, 0.3]
        // [0.3, 1]
        dfIntegral = d1.AbsIntegral(-1,1);
        CPPUNIT_ASSERT_MESSAGE( "Integral should be 9 but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 9) ;
        dfIntegral = d1.AbsIntegral(-1,0.3);
        CPPUNIT_ASSERT_MESSAGE( "Integral should be 5.85 but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 5.85) ;
        dfIntegral = d1.AbsIntegral(0.3, 1);
        CPPUNIT_ASSERT_MESSAGE( "Integral should be 3.15 but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 3.15) ;

        // halfen the complete distribution intervalls:
        // [-1, 0.3] --> [-0.675, -0.025]
        // [0.3 , 1] --> [0.475, 0.825]
        dfIntegral = d1.AbsIntegral(-0.675,-0.025);
        CPPUNIT_ASSERT_MESSAGE( "Integral should be " + to_string(5.85/2.) + " but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 5.85/2.) ;
        dfIntegral = d1.AbsIntegral(0.475,1-0.175);
        CPPUNIT_ASSERT_MESSAGE( "Integral should be " + to_string(3.15 /2.) + "but is: " + dfIntegral.Print() +  "\n distribution: \n" + d1.Print(10) , dfIntegral == 3.15 / 2.) ;
        

    }
       
    void Distribution_Shift()
    {
        CDistribution d1;
        
        d1.Add(-1,3);
        d1.Add(1,3);
        d1.Add(-0.1, d1.DistriValue(-0.1));
        d1.Add(-0.2, d1.DistriValue(-0.2));
        d1.Add(-0.3, d1.DistriValue(-0.3));
        d1.Add(-0.4, d1.DistriValue(-0.4));
        d1.Add(-0.5, d1.DistriValue(-0.5));
        d1.Add(-0.6, d1.DistriValue(-0.6));
        
        
        d1.Add(0.1, d1.DistriValue(0.1));
        d1.Add(0.2, d1.DistriValue(0.2));
        d1.Add(0.3, d1.DistriValue(0.3));
        d1.Add(0.5, d1.DistriValue(0.5));
        d1.Add(0.6, d1.DistriValue(0.6)); 
        
        CDigFloat dfMean = d1.Mean();
        CDigFloat dfShift = 0.012341234;
        
        d1.Shift(-dfShift);
        CDigFloat dfShiftedMean = d1.Mean();
        
        CPPUNIT_ASSERT_MESSAGE( "Mean should be " + (dfMean - dfShift).RawPrint(10) + " but is: " + dfShiftedMean.Print(10),
        dfShiftedMean == dfMean -dfShift);

        
     
    }   
    void Distribution_Coverage()
    {
        CDistribution d1,pd1;
                ////////////////////////////////////
        // coverage calculations 
        ///////////////////////////////////
        d1.Add(-1,3);
        d1.Add(1,3);
               d1.Add(-0.1, d1.DistriValue(-0.1));
        d1.Add(-0.2, d1.DistriValue(-0.2));
        d1.Add(-0.3, d1.DistriValue(-0.3));
        d1.Add(-0.4, d1.DistriValue(-0.4));
        d1.Add(-0.5, d1.DistriValue(-0.5));
        d1.Add(-0.6, d1.DistriValue(-0.6));
        
        
        d1.Add(0.1, d1.DistriValue(0.1));
        d1.Add(0.2, d1.DistriValue(0.2));
        d1.Add(0.3, d1.DistriValue(0.3));
        d1.Add(0.5, d1.DistriValue(0.5));
        d1.Add(0.6, d1.DistriValue(0.6)); 
        
        // shift the distribution to one polarity:
        // avoid zero odd order total integrals

        d1.Shift(-d1.firstVariable());
        
        CDigFloat dfMin, dfMax, dfMaxOld, dfMinOld;
//         cout << endl << endl;
        // check for iterations max. 50
        d1.MaxBinarySearchIterations(50);
        for( int iorder = 0; iorder < 3; iorder++)
        {
            dfMaxOld = 0;
            dfMinOld = 0;
            int nPercStart = 1;
            CDigFloat dfTotalIntegral = d1.AbsIntegral(iorder);
                LOGTRACE("Distribution_Test::Distribution_Coverage",string("total integral =") + dfTotalIntegral.RawPrint(10));
                LOGTRACE("Distribution_Test::Distribution_Coverage",string("int iorder = ") + to_string(iorder));
            CDigFloat dfMean = d1.Median(iorder );
            for(int iperc = nPercStart; iperc < 100; iperc++)
            {
                // check the integral
                CDigFloat dfExpectedIntegral = dfTotalIntegral * iperc / 100.;
                
                LOGTRACE("Distribution_Test::Distribution_Coverage","args:");
                LOGTRACE("Distribution_Test::Distribution_Coverage",string("int iperc = ") + to_string(iperc));
                LOGTRACE("Distribution_Test::Distribution_Coverage",string("expected integral =") +  dfExpectedIntegral.RawPrint(30));
                
                bool bSuccess = d1.CoverageInterval(iperc,dfMin, dfMax, iorder);
                LOGTRACE("Distribution_Test::Distribution_Coverage",string("calculated min =") +  dfMin.RawPrint(30));
                LOGTRACE("Distribution_Test::Distribution_Coverage",string("calculated max =") +  dfMax.RawPrint(30));
                
                // check for success
                CPPUNIT_ASSERT_MESSAGE( "return Succes != true (" + to_string(iorder) + "," + to_string(iperc) +") ", bSuccess==true) ;
      
                // check for new intervall being > than the old one
                if(iperc > nPercStart)
                CPPUNIT_ASSERT_MESSAGE( "old intervall must be less(" + to_string(iorder) + "," + to_string(iperc) +"): " +" \nminOld="+ dfMinOld.Print(30) +  "\nminNew: " + dfMin.RawPrint(30) + "\nmaxOld:" + dfMaxOld.RawPrint(30)+ "\nmaxNew:" + dfMax.RawPrint(30) , dfMinOld > dfMin && dfMaxOld < dfMax) ;
                
                
                CPPUNIT_ASSERT_MESSAGE( "integral must be " + to_string(iperc) + "\% of total integral:\n" +
                "total integral(" + to_string(iorder) + "," + to_string(iperc) +"):" + dfTotalIntegral.RawPrint(10)+ "\n" +
                "calculated from iperc * total: " + dfExpectedIntegral.RawPrint(30) + " \n" +
                "calc. from AbsIntegral:"+ d1.AbsIntegral(dfMin, dfMax,iorder).RawPrint(30) , 
                d1.AbsIntegral(dfMin,dfMax,iorder) == dfExpectedIntegral);
 
                // finally remember old values
                dfMinOld = dfMin; dfMaxOld = dfMax;
                
                // DEBUG
/*                
                cout << "order: " << iorder << endl
                     << "coverage: " << iperc << endl
                     << "mean: " << dfMean.RawPrint(10) << endl
                     << "min: " << dfMin.RawPrint(10) << endl
                     << "max: " << dfMax.RawPrint(10) << endl;*/
                
                
            }   // endfor(int iperc = 0; iperc < 100; iperc++)
        }   // endfor( int iorder = 0; iorder < 2; iorder++)

        
        // normalization after adding the last element
        pd1.Add(-2,0);
        pd1.Add(-1,3);
        pd1.Add(0,0);
        pd1.Add(1,3);
        pd1.Add(2,0);
        CDigFloat min,max;
        
        pd1.CoverageFromTo(50.,pd1.Distribution().begin()->first,max);
        CPPUNIT_ASSERT_MESSAGE( "max should be 0: " + max.RawPrint(30) , max == 0) ; 
        
        pd1.CoverageFromTo(50.,prev(pd1.Distribution().end())->first,max, true);
        CPPUNIT_ASSERT_MESSAGE( "max should be 0: " + max.RawPrint(30) , max == 0) ; 
        
        pd1.CoverageInterval(50.,min,max);
        CPPUNIT_ASSERT_MESSAGE( "min should be -1: " + min.RawPrint(30), min == -1) ;
        CPPUNIT_ASSERT_MESSAGE( "max should be 1: " + max.RawPrint(30) + "\n1: " + CDigFloat(1.).RawPrint(30), max == 1) ; 
        
        
        
        
    }
    
    
 
};
