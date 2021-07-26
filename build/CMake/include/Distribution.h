
/*
 * MIT License
 * 
 * Copyright (c) 2020 <author>
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

#ifndef CDISTRIBUTION_H
#define CDISTRIBUTION_H

#include "LoggingStrings.h"

// includes 
#include "DigFloat/DigFloat.h"

// macros
#define     MAX_BINARY_SEARCH_ITERATIONS    120
#define     DISTRI_POINTS_WARNING_LIMIT     1000

// typedef
typedef    map<CDigFloat,CDigFloat> MapDFDFType;
typedef    pair<CDigFloat,CDigFloat> PairDFType;




/**
 * @brief This class represents a "quasi continous" distribution defined by discrete  \anchor def-distri-point points which are tuples consisting of the \anchor def-distri-variable distribution variable 
 * and its positive \anchor def-distri-value distribution value. The \anchor def-quasi-continuity "quasi continuity" is achieved by linearly interpolating between the points of the distribution.
 * This class offers:
 * - adding points to the distribution
 * - "quasi conitinously" ...
 *    - ... interpolating distribution values for any given variable
 *    - ... calculating any n'th order weighted integral of this distribution
 *    - ... calculating intervals for a given coverage
 */


/**
 * @todo moving from 1d to nd distribution: replace the representation (CDigFloat) of the \ref def-distri-variable "distribution variable" by CCartesianCoordinates. This will mean
 * a lot of work keeping the \ref def-quasi-continuity "quasi continuity" (e.g. see function CDistribution::DistriValue: generally interpolating values for a nd distribution is much more sophisticated).
 */

class 
#ifdef _WIN32
_WIN_DLL_API
#endif
CDistribution 
{
public:
    
    ///////////////////////////////////
    // construction / destruction
    ///////////////////////////////////
    /**
     * Default constructor
     */
    CDistribution();

    /**
     * Copy constructor
     *
     * @param other CDistribution
     */
    CDistribution(const CDistribution& other);

    /**
     * Destructor
     */
    ~CDistribution();
    
    ///////////////////////////////////
    // assingment operators
    ///////////////////////////////////
    /**
     * Assignment operator
     *
     * @param other CDistribution&
     * @return CDistribution as copy of other
     */
    CDistribution& operator=(const CDistribution& other);

    ///////////////////////////////////
    // comparison operators
    ///////////////////////////////////
    /**
     * @brief equality operator
     *
     * @param other CDistribution&
     * @return bool 
     */
    bool operator==(const CDistribution& other) const;

    /**
     * @brief inequality operator
     *
     * @param[in] other CDistribution&
     * @return bool
     */
    bool operator!=(const CDistribution& other) const;
    
    ///////////////////////////////////
    // arithmetic operators
    ///////////////////////////////////
    /**
     * @brief multiplication assingment operator with factor: is applied on \ref def-distri-value "distribution values"
     *
     * @param[in] Value CDigFloat factor
     * @return CDistribution
     */
    CDistribution& operator*=(const CDigFloat& Value);

    /**
     * @brief division assingment operator with factor: is applied on \ref def-distri-value "distribution values"
     *
     * @param[in] Value CDigFloat divisor
     * @return CDistribution
     */
    CDistribution& operator/=(const CDigFloat& Value);
    
    /**
     * @brief summation assingment operator with factor: is applied on \ref def-distri-value "distribution values"
     *
     * @param[in] Value CDigFloat summand
     * @return CDistribution
     */
    CDistribution& operator+=(const CDigFloat& Value);
    
    /**
     * @brief subtraction assingment operator with factor: is applied on \ref def-distri-value "distribution values"
     *
     * @param[in] Value CDigFloat subtrahend
     * @return CDistribution
     */
    CDistribution& operator-=(const CDigFloat& Value);
    
    ///////////////////////////////////
    // getter 
    ///////////////////////////////////    
    /**
     * @brief return the first valid iterator of the distribution
     * 
     * @return MapDFDFType:const_iterator
     */
    inline MapDFDFType::const_iterator front() { return Distribution().begin();}  
    /**
     * @brief return the last valid iterator of the distribution
     * 
     * @return MapDFDFType:const_iterator
     */
    inline MapDFDFType::const_iterator back() { return prev(Distribution().end());}
    
    /**
     * @brief return the first variable of the distribution
     * 
     * @return CDigFloat
     */
    inline CDigFloat firstVariable() { return front()->first;}  
    /**
     * @brief return the last  varibale of the distribution
     * 
     * @return CDigFloat
     */
    inline CDigFloat lastVariable() { return back()->first;}
    
    
    
    ///////////////////////////////////
    // function on all distribution variables
    ///////////////////////////////////    
   /**
     * @brief shift is applied on \ref def-distri-variable "distribution variables"
     *
     * @param[in] Value CDigFloat shift added on all variables of distribution
     */
   void Shift(const CDigFloat& shift);
   /**
     * @brief scale is applied on \ref def-distri-variable "distribution variables"
     *
     * @param[in] Value CDigFloat scale is multiplies all variables of distribution
     */
    void Scale(const CDigFloat& scale);
    
    ///////////////////////////////////
    // functions
    ///////////////////////////////////    
    /**
     * @brief adds new \ref def-distri-point "point" to the distribution: the distribution value of this point must be positive
     *        otherwise the element will not be added.
     *
     * @param[in] variable CDigFloat as distribution variable of the new point
     * @param[in] value CDigFloat as distribution value of the new point
     */
    bool Add(const CDigFloat& variable, const CDigFloat value); 
    
    /**
     * @brief returns the distribution \ref def-distri-point "points" around the given variable
     *
     * @param[in] variable CDigFloat& variable to search surrounding \ref def-distri-point "points" for
     * @param[out] VariableLeft MapDFDFType::iterator& iterator for left \ref def-distri-point "point" (see ::MapDFDFType)
     * @param[out] VariableRight MapDFDFType::iterator& iterator for right \ref def-distri-point "point" (see ::MapDFDFType)
     */
    void GetInterval(const CDigFloat& variable, MapDFDFType::const_iterator& VariableLeft, MapDFDFType::const_iterator& VariableRight);     
    
    /**
     * @brief returns (interpolated) value at the given variable
     *
     * @param[in] variable CDigFloat variable to interpolate (linearily) the corresponding distribution value
     * @return CDigFloat the corresponding distribution value
     */
    CDigFloat DistriValue(const CDigFloat& variable);        
    
    /**
     * @brief returns the complete variable range
     *
     * @return CDigFloat the variable range
     */
    CDigFloat VariableRange(){ return CDigFloat(prev(Distribution().end())->first) - Distribution().begin()->first;}  
    
    /**
     * @brief returns the offset of the linear function defined by the distribution \ref def-distri-point "points" 
     * around the given variable
     *
     * @param[in] variable CDigFloat& variable to search surrounding \ref def-distri-point "points" for
     * @param[in] Left MapDFDFType::const_iterator* variable to search surrounding \ref def-distri-point "points" for
     * @return the linear offset
     */
    CDigFloat LinearOffset(const CDigFloat& Variable);
    
    /**
     * @brief returns the slope of the linear function defined by the distribution \ref def-distri-point "points" 
     * around the given variable
     * 
     * @param[in] variable CDigFloat& variable to search surrounding \ref def-distri-point "points" for
     * @return the linear slope
     */
    CDigFloat LinearSlope(const CDigFloat& Variable);     
     
    
    /**
     * @brief returns (interpolated) integral of the distribution between the user given limits
     *
     * @param[in] variableLeft CDigFloat for left limit of integral
     * @param[in] variableRight CDigFloat for right limit of integral
     * @param[in] nthOrder int for weighing the integral calculation with variable<SUP>+nthOrder</SUP> 
     * @return CDigFloat as
     */
    CDigFloat AbsIntegral(const CDigFloat& variableLeft, const CDigFloat& variableRight, int nthOrder = 0);     
    
    /**
     * @brief resets distribution 
     *
     */
    void  Reset();  
    
    /**
     * @brief returns distribution as map
     *
     * @return map<CDigFloat, CDigFloat> 
     */
    const MapDFDFType& Distribution(){return mDistribution;}   
    
    /**
     * @brief calculates the integral of the absolute value of this distribution
     *
     * @param[in] nthOrder int for weighing the integral calculation with variable<SUP>+nthOrder</SUP>
     * @return CDigFloat
     */
    CDigFloat AbsIntegral(const int& nthOrder =0);
        
    /**
     * @brief calculates minimal \ref def-distri-value "distribution value"
     *
     * @return CDigFloat minimal distribution value
     */
    CDigFloat Min();
    
    /**
     * @brief calculates maximal \ref def-distri-value "distribution value"
     *
     * @return CDigFloat maximal distribution value
     */
    CDigFloat Max();
        
    /**
     * @brief calculates the nthOrder mean variable weighted by distribution values
     *
     * @param[in] nthOrder int for weighing the integral calculation with variable<SUP>+nthOrder</SUP>
     * @return CDigFloat nthOrder mean variable
     */
    CDigFloat Mean( int nthOrder = 1 ) { return CDigFloat( AbsIntegral(nthOrder) / AbsIntegral(0));}
        
    /**
     * @brief calculates the nthOrder median variable weighted by distribution values
     *
     * @param[in] nthOrder int for weighing the integral calculation with variable<SUP>+nthOrder</SUP>
     * @return CDigFloat nthOrder median variable
     */
    CDigFloat Median( int nthOrder = 0 );
        
    /**
     * @brief calculates the interval for a given coverage (given in %) of the nthOrder integral. The coverage is places symmetrically around the corresponding nthOrder medial (see ::Median)
     *
     * @param[in] dfCoveragePercent CDigFloat between 0 and 100 representing the coverage
     * @param[in] dfStart CDigFloat for fixed starting point of the coverage interval 
     * @param[out] dfTo CDigFloat for end point (to be searched for) of the coverage interval
     * @param[in] bReverse bool flag setting the search direction (false: dfFrom < dfTo)
     * @param[in] nthOrder int setting the order of the integral 
     * @return bool : false if fails (may happen if an impossible setting is supplied)
     */
    bool CoverageFromTo(const CDigFloat& dfCoveragePercent, const CDigFloat& dfFrom, CDigFloat& dfTo, bool bReverse = false, int nthOrder = 0);
        
    /**
     * @brief calculates the interval for a given coverage (given in %) of the nthOrder integral. The coverage is places symmetrically around the corresponding nthOrder medial (see ::Median)
     *
     * @param[in] dfCoveragePercent CDigFloat between 0 and 100 representing the coverage
     * @param[out] dfMin CDigFloat for the minimal value of the coverage interval
     * @param[out] dfMax CDigFloat for the maximal value of the coverage interval
     * @param[in] nthOrder int setting the order of the integral 
     * @return bool : false if fails (only when dfCoveragePercent > 100 or < 0)
     */
    bool CoverageInterval(const CDigFloat& dfCoveragePercent, CDigFloat& dfMin, CDigFloat& dfMax, int nthOrder = 0);
    
    bool OutOfBounds(const CDigFloat& variable) {return variable < Distribution().begin()->first || variable > prev(Distribution().end())->first; }
    
    /**
     * @brief returns string of x-intervall and no of points
     *
     * @param nPrecision int
     * @return string
     */
    string PrintMetaInfo();
    
    /**
     * @brief returns string of all element pairs with given precision
     *
     * @param nPrecision int
     * @return string
     */
    string Print(int nPrecision, bool bWithError = true);
    
    ///////////////////////////////
    // getter / setter
    //////////////////////////////
   /**
     * @brief returns the maximal number used for binary search: is used for ::CoverageFromTo
     * 
     * @return unsigned int
     */
    unsigned int MaxBinarySearchIterations() { return uiMaxBinarySearchIterations;}
   /**
     * @brief sets the maximal number used for binary search: is used for ::CoverageFromTo
     * @param uiIterations unsigned int maximal number of iterations
     */
    void MaxBinarySearchIterations(const unsigned int& uiIterations) { uiMaxBinarySearchIterations = uiIterations;}

protected:    
    
    /**
     * @brief getting integral value of a distribution element (point (variable and value) and the next higher point)
     *
     * @param[in] iel MapDFDFType::const_iterator of given distribtuion point
     */
    CDigFloat _IntegralConsecutiveElements(const MapDFDFType::const_iterator& iel,const int& nthOrder =0);   
    
    /**
     * @brief getting integral value of two distribution points (variable and value) assuming linear connection between these points
     *
     */
    CDigFloat _IntegralOfTwoPoints(const CDigFloat& x1, const CDigFloat& y1, const CDigFloat& x2, const CDigFloat& y2 , const int &nthOrder =0);
    
    /**
     * @brief calculates the value of the n'th order weighed determined integral over a linear function f (x) = slope * x + offset, i.e. Integral ( x^n*f(x)). This function is needed for all integral calculations done within this class. (see ::AbsIntegral, ::CoverageFromTo, ::CoverageInterval)
     *
     * @param[in] dfX1 CDigFloat the lower value to calculate the integral for
     * @param[in] dfX2 CDigFloat the upper value to calculate the integral for
     * @param[in] slope, CDigFloat the slope of the integrated and weighed linear function
     * @param[in] offset CDigFloat the offset of the integrated and weighed linear function
     * @param[in] nthOrder int the order of the weighing of the linear function
     * @return CDigFloat the value of the primitive integral at x
     */
    CDigFloat _nthOrderWeightedIntegral(const CDigFloat& dfX1, const CDigFloat& dfX2,const CDigFloat& slope, const CDigFloat& offset, const int& nthOrder);    
    
    /**
     * @brief returns the offset of the linear function defined by the distribution \ref def-distri-point "points" 
     * for the given interval
     *
     * @param[in] Left MapDFDFType::const_iterator& iterator pointing to the left limit of this distribution
     * @return the linear offset
     */
    CDigFloat _LinearOffset(MapDFDFType::const_iterator& Left);   
    
    /**
     * @brief returns the slope of the linear function defined by the distribution \ref def-distri-point "points" 
     * for the given interval
     *
     * @param[in] Left MapDFDFType::const_iterator& iterator pointing to the left limit of this distribution
     * @return the linear offset
     */
    CDigFloat _LinearSlope(MapDFDFType::const_iterator& Left);
    
    /**
     * @brief initializes member
     *
     */
    void _Init();
    /////////////////////////////////////////
    // protected member
    //////////////////////////////////////
    
    unsigned int uiMaxBinarySearchIterations;    
    MapDFDFType mDistribution;
};


/////////////////////////////////////////
// external functions
//////////////////////////////////////
string 
#ifdef _WIN32
_WIN_DLL_API
#endif
Print(const MapDFDFType::const_iterator& it);

#endif // CDISTRIBUTION_H
