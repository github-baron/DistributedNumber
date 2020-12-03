
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

// includes 
#include "DigFloat/DigFloat.h"

// typedef
typedef    map<CDigFloat,CDigFloat> M_DFDF;
typedef    pair<CDigFloat,CDigFloat> P_DFDF;


/**
 * @brief This class represents a "quasi continous" distribution defined by discrete  \anchor def-distri-point points which are tuples consisting of a variable 
 * and its positive distribution value. The "quasi continuity" is achieved by linearly interpolating between the points of the distribution.
 * This class offers:
 * - adding points to the distribution
 * - "quasi conitinously" ...
 *    - ... interpolating distribution values for any given variable
 *    - ... calculating any n'th order weighted integral of this distribution
 *    - ... calculating intervalls for a given coverage
 */


class CDistribution 
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
    // operators
    ///////////////////////////////////
    /**
     * Assignment operator
     *
     * @param other CDistribution&
     * @return CDistribution as copy of other
     */
    CDistribution& operator=(const CDistribution& other);

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
     * @param[out] VariableLeft M_DFDF::iterator& iterator for left \ref def-distri-point "point" (see ::M_DFDF)
     * @param[out] VariableRight M_DFDF::iterator& iterator for right \ref def-distri-point "point" (see ::M_DFDF)
     */
    void GetIntervall(const CDigFloat& variable, M_DFDF::const_iterator& VariableLeft, M_DFDF::const_iterator& VariableRight);     
    
    /**
     * @brief returns (interpolated) value at the given variable
     *
     * @param[in] variable CDigFloat variable to interpolate (linearily) the corresponding distribution value
     * @return CDigFloat the corresponding distribution value
     */
    CDigFloat DistriValue(const CDigFloat& variable);     
    
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
    const M_DFDF& Distribution(){return mDistribution;}   
    
    /**
     * @brief calculates the integral of the absolute value of this distribution
     *
     * @return CDigFloat
     */
    CDigFloat AbsIntegral(const int& nthOrder =0);
        
    /**
     * @brief calculates the nthOrder mean variable weighted by distribution values
     *
     * @return CDigFloat nthOrder mean variable
     */
    CDigFloat Mean( int nthOrder = 1 ) { return AbsIntegral(nthOrder) / AbsIntegral(0);}
        
    /**
     * @brief calculates the nthOrder median variable weighted by distribution values
     *
     * @return CDigFloat nthOrder median variable
     */
    CDigFloat Median( int nthOrder = 0 );
        
    /**
     * @brief calculates the intervall for a given coverage (given in %) of the nthOrder integral. The coverage is places symmetrically around the corresponding nthOrder medial (see ::Median)
     *
     * @param[in] dfCoveragePercent CDigFloat between 0 and 100 representing the coverage
     * @param[in] dfStart CDigFloat for fixed starting point of the coverage intervall 
     * @param[out] dfTo CDigFloat for end point (to be searched for) of the coverage intervall
     * @param[in] bReverse bool flag setting the search direction (false: dfFrom < dfTo)
     * @param[in] nthOrder int setting the order of the integral 
     * @return bool : false if fails (may happen if an impossible setting is supplied)
     */
    bool CoverageFromTo(const CDigFloat& dfCoveragePercent, const CDigFloat& dfFrom, CDigFloat& dfTo, bool bReverse = false, int nthOrder = 0);
        
    /**
     * @brief calculates the intervall for a given coverage (given in %) of the nthOrder integral. The coverage is places symmetrically around the corresponding nthOrder medial (see ::Median)
     *
     * @param[in] dfCoveragePercent CDigFloat between 0 and 100 representing the coverage
     * @param[out] dfMin CDigFloat for the minimal value of the coverage intervall
     * @param[out] dfMax CDigFloat for the maximal value of the coverage intervall
     * @param[in] nthOrder int setting the order of the integral 
     * @return bool : false if fails (only when dfCoveragePercent > 100 or < 0)
     */
    bool CoverageIntervall(const CDigFloat& dfCoveragePercent, CDigFloat& dfMin, CDigFloat& dfMax, int nthOrder = 0);
    
    /**
     * @brief returns string of all element pairs with given precision
     *
     * @param nPrecision int
     * @return string
     */
    string Print(int nPrecision);

    
protected:    
    
    /**
     * @brief getting integral value of a distribution element (point (variable and value) and the next higher point)
     *
     * @param[in] iel M_DFDF::const_iterator of given distribtuion point
     */
    CDigFloat _IntegralConsecutiveElements(const M_DFDF::const_iterator& iel,const int& nthOrder =0);   
    
    /**
     * @brief getting integral value of two distribution points (variable and value) assuming linear connection between these points
     *
     */
    CDigFloat _IntegralOfTwoPoints(const CDigFloat& x1, const CDigFloat& y1, const CDigFloat& x2, const CDigFloat& y2 , const int &nthOrder =0);
    
    CDigFloat _nthOrderWeightedPrimitiveIntegral(const CDigFloat& x,const CDigFloat& slope, const CDigFloat& offset, const int& nthOrder);
    /**
     * @brief initializes member
     *
     */
    void _Init();
    
    
    M_DFDF mDistribution;
};

#endif // CDISTRIBUTION_H
