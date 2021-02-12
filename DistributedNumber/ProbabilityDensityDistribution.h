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


/**
 * @brief This class represents a probability density distribution. 
 * It offers :
 * - handling of probability:
 *     - adding non-probability \ref def-distri-value "values" is possible (e.g. countings of the occurrence of a variable)
 *     - normalization will be done when needed (e.g. query for a distribution value will return a probability density value)
 *     - the query for 
 * -
 */

#ifndef CPROBABILITYDENSITYDISTRIBUTION_H
#define CPROBABILITYDENSITYDISTRIBUTION_H

// includes
#include "Distribution.h"
#include <iostream>

// typedefs
typedef pair<CDigFloat,CDigFloat> PairDFType;
typedef pair< PairDFType, PairDFType >SubIntLimitsType;
typedef vector< SubIntLimitsType >  VectorSubIntLimitsType;
//  { (targetValue,{ ((minthis,maxthis),(minother,maxother)) }) }
// e.g. [1].second[0].first.first: first integration limit of this distribution for the first sub integration step and the second target value
// e.g. [1].second[0].second.first: first integration limit of other distribution for the first sub integration step and the second target value
// no. of target values : PDD_DEFAULT_INTEGRATION_STEPS
// no. of integration limits for this and other: PDD_DEFAULT_SUBINTEGRATION_STEPS
typedef vector< pair<CDigFloat, VectorSubIntLimitsType > > ConvolutionPlanType;

// macros
#define PDD_DEFAULT_INTEGRATION_STEPS 50
#define PDD_DEFAULT_SUBINTEGRATION_STEPS 10

// enums
enum class ProbDistOp
{
    pdoUnknown = -1,
    pdoPlus,
    pdoFirst = pdoPlus,
    pdoMinus,
    pdoMult,
    pdoDiv,
    pdoLast
};

class 

#ifdef _WIN32
_WIN_DLL_API
#endif
CProbabilityDensityDistribution : public CDistribution
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
    
    ///////////////////////////////////
    // arithmetic operators
    ///////////////////////////////////
    /**
     * @brief summation assignment operator: results in a convolution with the sums of variables as the new \ref def-distri-variable "distribution variable" and the corresponding probability product as the new \ref def-distri-value "distribution values". The resulting distribution has n*m elements when (*this).size() = n and other.size() = m
     *
     * @param[in] other CProbabilityDensityDistribution
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution& operator += (CProbabilityDensityDistribution& other);
    /**
     * @brief summation operator: calls ::operator+= 
     *
     * @param[in] other CProbabilityDensityDistribution
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution operator + (CProbabilityDensityDistribution& other);
    
    /**
     * @brief subtraction assignment operator: results in a convolution with the sums of variables as the new \ref def-distri-variable "distribution variable" and the corresponding probability product as the new \ref def-distri-value "distribution values". The resulting distribution has n*m elements when (*this).size() = n and other.size() = m
     *
     * @param[in] other CProbabilityDensityDistribution
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution& operator -= (CProbabilityDensityDistribution& other);
    /**
     * @brief subtraction operator: calls ::operator-= 
     *
     * @param[in] other CProbabilityDensityDistribution
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution operator - (CProbabilityDensityDistribution& other);

    /**
     * @brief multiplication assignment operator: results in a convolution with the sums of variables as the new \ref def-distri-variable "distribution variable" and the corresponding probability product as the new \ref def-distri-value "distribution values". The resulting distribution has n*m elements when (*this).size() = n and other.size() = m
     *
     * @param[in] other CProbabilityDensityDistribution
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution& operator *= (CProbabilityDensityDistribution& other);
    /**
     * @brief multiplication operator: calls ::operator*= 
     *
     * @param[in] other CProbabilityDensityDistribution
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution operator * (CProbabilityDensityDistribution& other);
    /**
     * @brief multiplication assignment operator: results in a convolution with the sums of variables as the new \ref def-distri-variable "distribution variable" and the corresponding probability product as the new \ref def-distri-value "distribution values". The resulting distribution has n*m elements when (*this).size() = n and other.size() = m
     *
     * @param[in] other CProbabilityDensityDistribution
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution& operator /= (CProbabilityDensityDistribution& other);
    /**
     * @brief division operator: calls ::operator/= 
     *
     * @param[in] other CProbabilityDensityDistribution
     * @return CProbabilityDensityDistribution
     */
    CProbabilityDensityDistribution operator / (CProbabilityDensityDistribution& other);

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
     * @brief returns the actual integral (0th order) of the distribution
     *
     * @return CDigFloat 0th order integral
     */
    CDigFloat OrigIntegral();
   
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
    /**
     * @brief normalize: divides by norm (integral over the distribution) 
     *
     * @param other CDistribution&
     * @return bool
     */
    void Normalize();     
    /**
     * @brief revert normalizatiion: multiplies by norm (integral over the distribution) 
     *
     * @param other CDistribution&
     * @return bool
     */
    void DeNormalize(); 
    
    ////////////////////////////////////////////
    // getter / setter
    ////////////////////////////////////////////   
    void IntegrationSteps(int nSteps) { m_nIntegrationSteps = nSteps;}
    int IntegrationSteps() const { return m_nIntegrationSteps;}
    void SubIntegrationSteps(int nSteps) { m_nSubIntegrationSteps = nSteps;}
    int SubIntegrationSteps() const { return m_nSubIntegrationSteps;}
protected:     
    /**
     * @brief initializes member
     *
     */
    void _Init();   
    
    /**
     * @brief returns the distribution to a corresponding operation between two variables
     * 
     * @return CProbabilityDensityDistribution resulting from corresponding operation
     *
     */     
    CProbabilityDensityDistribution& _GeneralOperatorsFunction(CProbabilityDensityDistribution& other, const ProbDistOp Operation);
    
    /**
     * @brief returns the distribution range for a convolution with this distri and the other 
     * 
     * @param[in] other M_DFDF distribution
     * @param[in] Operation ProbDistOp the convolution operation (e.g. sum , subtration, ..)
     * @param[out] TargetRangeStart CDigFloat start variable of convolution result
     * @param[out] TargetRangeEnd CDigFloat end variable of convolution result
     * 
     */     
    void _GetRangeFromDistributionOperation(CProbabilityDensityDistribution& other, const ProbDistOp Operation, CDigFloat& TargetRangeStart, CDigFloat& TargetRangeEnd);
   
   void _SetConvolutionPlan(CProbabilityDensityDistribution& other, const ProbDistOp Operation);
   void _SetSubIntegrationIntervalls4TargetValue(const ProbDistOp Operation, const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervalls4TargetValue);
   void _SetSubIntegrationIntervalls4TargetValue4Addition(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervalls4TargetValue);
   void _SetSubIntegrationIntervalls4TargetValue4Subtraction(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervalls4TargetValue);
   
    /**
     * @brief calculates the sub integration intervalls needed for multiplication of two random variables (see CProbabilityDensityDistribution::operator*) for a specified resulting random value (called target value)
     * 
     * @param[in] dfActTargetValue CDigFloat target value to find sub integration intervalls for
     * @param[in] Other CProbabilityDensityDistribution the distribution which this distri is multiplied with
     * @param[out] SubIntervalls4TargetValue VectorSubIntLimitsType keeps the sub integration intervalls (see CProbabilityDensityDistribution::SubIntegrationSteps) for this distri and Other needed to calculate the total probability for dfActTargetValue
     * 
     */     
   void _SetSubIntegrationIntervalls4TargetValue4Multiplication(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& Other, VectorSubIntLimitsType& SubIntervalls4TargetValue);
   
    /**
     * @brief calculates the total limits of the sub integration intervalls needed for multiplication of two random variables (see CProbabilityDensityDistribution::operator*) for a specified negative resulting random value (called target value).This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervalls4TargetValue4Multiplication.
     * 
     * @param[in] dfActTargetValue CDigFloat target value to find sub integration intervalls for
     * @param[in] One CProbabilityDensityDistribution the distribution which Other distri is multiplied with
     * @param[in] Other CProbabilityDensityDistribution the distribution which One distri is multiplied with
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalIntervall4NegativeTargetValue4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);
   
    /**
     * @brief calculates the total limits of the sub integration intervalls needed for multiplication of positive ranges of two random variables (see CProbabilityDensityDistribution::operator*) for a specified positive resulting random value (called target value). This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervalls4TargetValue4Multiplication.
     * 
     * @param[in] dfActTargetValue CDigFloat positive target value to find sub integration intervalls for
     * @param[in] One CProbabilityDensityDistribution the distribution which Other distri is multiplied with
     * @param[in] Other CProbabilityDensityDistribution the distribution which One distri is multiplied with
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalIntervall4PositiveTargetValuePositiveVariables4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);   
   
    /**
     * @brief calculates the total limits of the sub integration intervalls needed for multiplication of negative ranges of two random variables (see CProbabilityDensityDistribution::operator*) for a specified positive resulting random value (called target value).This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervalls4TargetValue4Multiplication.
     * 
     * @param[in] dfActTargetValue CDigFloat positive target value to find sub integration intervalls for
     * @param[in] One CProbabilityDensityDistribution the distribution which Other distri is multiplied with
     * @param[in] Other CProbabilityDensityDistribution the distribution which One distri is multiplied with
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalIntervall4PositiveTargetValueNegativeVariables4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);
   
   
   void _SetSubIntegrationIntervalls4TargetValue4Division(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervalls4TargetValue);
   void _GetSubIntegrationTotalIntervall4PositiveTargetValuePositiveVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);
   
   
   void _GetSubIntegrationTotalIntervall4PositiveTargetValueNegativeVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);
   
   void _GetSubIntegrationTotalIntervall4NegativeTargetValuePosNegVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);
   
   
   void _GetSubIntegrationTotalIntervall4NegativeTargetValueNegPosVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);
   
   
   
   void _InitIntegrationLimits(PairDFType& pdtLimits);
   
   
   ////////////////////////////////////////////////////////
   // member variables
   ////////////////////////////////////////////////////////
    
    bool bNormalized; 
    CDigFloat dfAbsIntegral;
    int m_nIntegrationSteps;
    int m_nSubIntegrationSteps;
    ConvolutionPlanType m_ConvolutionPlan;

};

#endif // CPROBABILITYDENSITYDISTRIBUTION_H
