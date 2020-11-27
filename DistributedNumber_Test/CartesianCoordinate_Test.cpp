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

#include "../DistributedNumber/CartesianCoordinate.h"

class CCartesianCoordinate_Test : public CppUnit::TestFixture  {

private:

public:
    void setUp()
    {
    }

    void tearDown() 
    {
    }

    void CheckOperators()
    {
        
        CCartesianCoordinate c1,c2;
        c1.push_back(1);
        c1.push_back(1);
        c1.push_back(1);
        
        // assignment and multiplication
        c2 = c1*2;
        CPPUNIT_ASSERT_MESSAGE( "c2[0] should be 2" , c2[0] == 2) ;
        CPPUNIT_ASSERT_MESSAGE( "c2[1] should be 2" , c2[1] == 2) ;
        CPPUNIT_ASSERT_MESSAGE( "c2[2] should be 2" , c2[2] == 2) ;
       
        // check inequality operator
        CPPUNIT_ASSERT_MESSAGE( "should be unequal" , c1 != c2) ;
        
        // assignment and division 
        c2/=2;
        CPPUNIT_ASSERT_MESSAGE( "c2[0] should be 1" , c2[0] == 1) ;
        CPPUNIT_ASSERT_MESSAGE( "c2[1] should be 1" , c2[1] == 1) ;
        CPPUNIT_ASSERT_MESSAGE( "c2[2] should be 1" , c2[2] == 1) ;
        
        // equality operator
        CPPUNIT_ASSERT_MESSAGE( "should be equal" , c1 == c2) ;
        
        // check comparison concept: highest coordinate index rules for ">" / "<"
        // compare coordinates:
        // c1 = (0 ; 1 ; 0)
        // c2 = (0 ; 0 ; 1)
        c1[0]=0; c1[1]=1; c1[2]=0;  // coord (0 ; 1 ; 0)
        c2[0]=0; c2[1]=0; c2[2]=1;  // coord (0 ; 0 ; 1)
        CPPUNIT_ASSERT_MESSAGE( "c1 should be < c2" , c1 < c2) ;
        CPPUNIT_ASSERT_MESSAGE( "c2 should be > c1" , c2 > c1) ;
        CPPUNIT_ASSERT_MESSAGE( "c1 should be <= c2" , c1 <= c2) ;
        CPPUNIT_ASSERT_MESSAGE( "c2 should be >= c1" , c2 >= c1) ;
        
        // compare coordinates:
        // c1 = (0 ; 0 ; 1)
        // c2 = (0 ; 1 ; 1)
        c1[0]=0; c1[1]=0; c1[2]=1;  // coord (0 ; 1 ; 0)
        c2[0]=0; c2[1]=1; c2[2]=1;  // coord (0 ; 0 ; 1)
        CPPUNIT_ASSERT_MESSAGE( "c1 should be < c2" , c1 < c2) ;
        CPPUNIT_ASSERT_MESSAGE( "c2 should be > c1" , c2 > c1) ;
        CPPUNIT_ASSERT_MESSAGE( "c1 should be <= c2" , c1 <= c2) ;
        CPPUNIT_ASSERT_MESSAGE( "c2 should be >= c1" , c2 >= c1) ;
        
        // compare coordinates:
        // c1 = (0 ; 0 ; 1)
        // c2 = (1 ; 0 ; 1)
        c1[0]=0; c1[1]=0; c1[2]=1;  // coord (0 ; 1 ; 0)
        c2[0]=1; c2[1]=0; c2[2]=1;  // coord (0 ; 0 ; 1)
        CPPUNIT_ASSERT_MESSAGE( "c1 should be < c2" , c1 < c2) ;
        CPPUNIT_ASSERT_MESSAGE( "c2 should be > c1" , c2 > c1) ;
        CPPUNIT_ASSERT_MESSAGE( "c1 should be <= c2" , c1 <= c2) ;
        CPPUNIT_ASSERT_MESSAGE( "c2 should be >= c1" , c2 >= c1) ;

        
    }
    

    void CheckFunctions()
    {
        
        CCartesianCoordinate c1,c2, c3;
        
        c1.push_back(1);
        c1.push_back(1);
        c1.push_back(1);
        
        // print
        c1.Precision(6);
        CPPUNIT_ASSERT_MESSAGE( "should be equal: (1.000000,1.000000,1.000000)" + c1.Print(), c1.Print() == "(1.000000,1.000000,1.000000)") ;
        c1.Precision(0);
        CPPUNIT_ASSERT_MESSAGE( "should be equal: (1,1,1)" + c1.Print(), c1.Print() == "(1,1,1)") ;
        

        // distance
        CPPUNIT_ASSERT_MESSAGE( "should be equal: sqrt(3.) = " + to_string(sqrt(3.)) + " and " + c1.Distance().Print(), c1.Distance() == sqrt(3.)) ;
        
        // ScalarProduct
        c2 = c1;
        CDigFloat dfScalarProd = ScalarProduct(c2,c1);
        CPPUNIT_ASSERT_MESSAGE( "should be equal: sqrt(3.) = " + to_string(sqrt(3.)) + " and " + dfScalarProd.Print(), dfScalarProd == sqrt(3.)) ;
        
        // vector product 
        CCartesianCoordinate ccVectorProd = VectorProduct(c2,c1);
        c3.push_back(0);
        c3.push_back(0);
        c3.push_back(0);
        CPPUNIT_ASSERT_MESSAGE( "should be equal: (0,0,0) but is:  " + ccVectorProd.Print(), ccVectorProd == c3) ;
        
        c1[0] = 0;
        c2[1] = 0;
        ccVectorProd = VectorProduct(c1,c2);
        // (0,1,1) x (1,0,1) = (1*1-0,1*1-0,0-1) = (1,1,-1)
        c3[0]=1; c3[1]=1; c3[2]=-1;
        CPPUNIT_ASSERT_MESSAGE( "should be equal: (1,1,-1) but is:  " + ccVectorProd.Print(), ccVectorProd == c3) ;
        
        
    }
    
 
};
