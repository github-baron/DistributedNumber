/*
 * 
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
 * 
 * 
 */

#ifndef CPROBABILITYDENSITYDISTRIBUTION_H
#define CPROBABILITYDENSITYDISTRIBUTION_H

// includes
#include "Distribution.h"

/**
 * @brief This class represents a probability density distribution. 
 * It offers :
 * - handling of probability:
 *     - adding non-probability \ref def-distri-value "values" is possible (e.g. countings of the occurrence of a variable)
 *     - normalization will be done when needed (e.g. query for a distribution value will return a probability density value)
 *     - the query for 
 * -
 */
class CProbabilityDensityDistribution :  CDistribution
{
public:
    
    ////////////////////////////////////////////
    // construction / destruction
    ////////////////////////////////////////////
    /**
     * Default constructor
     */
    CProbabilityDensityDistribution();

    /**
     * Copy constructor
     *
     * @param other CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution(const CProbabilityDensityDistribution& other);

    /**
     * Copy constructor
     *
     * @param other CDistribution
     */
    CProbabilityDensityDistribution(const CDistribution& other);

    /**
     * Destructor
     */
    ~CProbabilityDensityDistribution();

    ////////////////////////////////////////////
    // assignment operator
    ////////////////////////////////////////////
    /**
     * @brief Assignment operator
     *
     * @param other CProbabilityDensityDistribution
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution& operator=(const CProbabilityDensityDistribution& other);

    /**
     * @brief Assignment operator
     *
     * @param other CDistribution
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution& operator=(const CDistribution& other);
    
    ////////////////////////////////////////////
    // binary operator
    ////////////////////////////////////////////
   ///////////////////////////////////
    // binary operators
    ///////////////////////////////////
    /**
     * @brief multiplication assingment operator with factor: is applied on \ref def-distri-value "distribution values"
     *
     * @param[in] Value CDigFloat factor
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution& operator*=(const CDigFloat& Value);

    /**
     * @brief division assingment operator with factor: is applied on \ref def-distri-value "distribution values"
     *
     * @param[in] Value CDigFloat divisor
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution& operator/=(const CDigFloat& Value);
    
    /**
     * @brief summation assingment operator with factor: is applied on \ref def-distri-value "distribution values"
     *
     * @param[in] Value CDigFloat summand
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution& operator+=(const CDigFloat& Value);
    
    /**
     * @brief subtraction assingment operator with factor: is applied on \ref def-distri-value "distribution values"
     *
     * @param[in] Value CDigFloat subtrahend
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution& operator-=(const CDigFloat& Value);
    
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

    ////////////////////////////////////////////
    // comparison operator
    ////////////////////////////////////////////
    /**
     * @brief equality operator
     *
     * @param other CProbabilityDensityDistribution
     * @return bool
     */
    bool operator==(const CProbabilityDensityDistribution& other) const;

    /**
     * @todo inequality operator
     *
     * @param other CProbabilityDensityDistribution
     * @return bool
     */
    bool operator!=(const CProbabilityDensityDistribution& other) const;
    
    ////////////////////////////////////////////
    // functions
    ////////////////////////////////////////////   
     /**
     * @brief returns (interpolated) value at the given variable and takes care that the distribution
     * is normalized: if not this will be done
     *
     * @param[in] variable CDigFloat variable to interpolate (linearily) the corresponding distribution value
     * @return CDigFloat the corresponding distribution value
     */
    CDigFloat DistriValue(const CDigFloat& variable);     
    
   /**
     * @brief returns the distribution \ref def-distri-point "points" around the given variable
     *
     * @param[in] variable CDigFloat& variable to search surrounding \ref def-distri-point "points" for
     * @param[out] VariableLeft M_DFDF::iterator& iterator for left \ref def-distri-point "point" (see ::M_DFDF)
     * @param[out] VariableRight M_DFDF::iterator& iterator for right \ref def-distri-point "point" (see ::M_DFDF)
     */
    void GetIntervall(const CDigFloat& variable, M_DFDF::const_iterator& VariableLeft, M_DFDF::const_iterator& VariableRight);   
    
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
    CDigFloat Mean( int nthOrder = 1 );
   
    /**
     * @brief adds new distribution element if distribution is not already normalized (see ::_Normalize).
     *        In case of adding the last element the surrounding elements will be set to zero (... if not done)
     *        and the whole distribution will be normalized (see ::_Normalize).
     *
     * @param variable CDigFloat& distribution variable
     * @param value CDigFloat&
     * @param bLastElement bool if set --> ::_Normalize is called and no further element can be added
     * @return bool
     */
    bool Add(const CDigFloat& variable, const CDigFloat value, bool bLastElement = false); 
  
protected:    
    /**
     * @brief normalize: divides by norm (integral over the distribution) 
     *
     * @param other CDistribution&
     * @return bool
     */
    void _Normalize();     
    /**
     * @brief revert normalizatiion: multiplies by norm (integral over the distribution) 
     *
     * @param other CDistribution&
     * @return bool
     */
    void _DeNormalize();  
    /**
     * @brief initializes member
     *
     */
    void _Init();    
    
    bool bNormalized;
    CDigFloat dfAbsIntegral;

};

#endif // CPROBABILITYDENSITYDISTRIBUTION_H
