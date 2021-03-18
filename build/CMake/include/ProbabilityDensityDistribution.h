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
typedef pair< PairDFType, PairDFType >SubIntLimitsType;
typedef vector< SubIntLimitsType >  VectorSubIntLimitsType;
//  { (targetValue,{ ((minthis,maxthis),(minother,maxother)) }) }
// e.g. [1].second[0].first.first: first integration limit of this distribution for the first sub integration step and the second target value
// e.g. [1].second[0].second.first: first integration limit of other distribution for the first sub integration step and the second target value
// no. of target values : PDD_DEFAULT_INTEGRATION_STEPS
// no. of integration limits for this and other: PDD_DEFAULT_SUBINTEGRATION_STEPS
typedef vector< pair<CDigFloat, VectorSubIntLimitsType > > ConvolutionPlanType;
typedef vector< PairDFType >  VectorPairDFType;
typedef vector< pair<CDigFloat, VectorPairDFType  > > ConvolutionPlanType2;

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
    pdoCov,
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
    // distribution generation
    ////////////////////////////////////////////
    /**
     * @brief clears existing distri and generates a constant distribution for the given interval
     * 
     * @param[in] dfXMin CDigFloat minimal value of the interval
     * @param[in] dfXMax CDigFloat maximal value of the interval
     */
    void ConstDistribution(const CDigFloat& dfXMin, const CDigFloat& dfXMax);
    /**
     * @brief clears existing distri and generates a symmetrich triangulated distribution given by an interval
     * 
     * @param[in] dfXMin CDigFloat minimal value of the interval
     * @param[in] dfXMax CDigFloat maximal value of the interval
     */
    void TriangleDistribution(const CDigFloat& dfXMin, const CDigFloat& dfXMax);
    
    /**
     * @brief clears existing distri and generates a distribution according to function \f$f(x,\mu, \sigma) = \frac{1}{\sqrt{2\pi\sigma^2}} e^{-\frac{(x-\mu)^2}{2\sigma^2}}\f$ 
     * 
     * @param[in] nPoints int number of points
     * @param[in] dfMu CDigFloat expectation value
     * @param[in] dfSigma CDigFloat width of curve
     * @param[in] dfMinVal CDigFloat minimum value to start calculation with (the lower, the broader the interval for calculating of this distribution)
     */
    void NormalDistribution(const int nPoints, const CDigFloat& dfMu, const CDigFloat& dfSigma, const CDigFloat& dfMinVal = 0.001);
    
     
    /**
     * @brief clears existing distri and calculates covariance distri calculated from this and an Other given distribution
     * 
     * @param[in] Other CProbabilityDensityDistribution other distribution from which the covariance with this distribution is calculated
     */
    void CovarianceDistribution(CProbabilityDensityDistribution& Other);
    
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
     * @brief calculates the nthOrder variance  
     *
     * @return CDigFloat nthOrder mean variable
     */
    CDigFloat Variance( int nthOrder = 1 );
        
    /**
     * @brief calculates the nthOrder mean variable weighted by distribution values
     *
     * @return CDigFloat nthOrder mean variable
     */
    CDigFloat Covariance(CProbabilityDensityDistribution& Other); 
    
   /**
     * @brief returns the distribution \ref def-distri-point "points" around the given variable
     *
     * @param[in] variable CDigFloat& variable to search surrounding \ref def-distri-point "points" for
     * @param[out] VariableLeft M_DFDF::iterator& iterator for left \ref def-distri-point "point" (see ::M_DFDF)
     * @param[out] VariableRight M_DFDF::iterator& iterator for right \ref def-distri-point "point" (see ::M_DFDF)
     */
    void GetInterval(const CDigFloat& variable, MapDFDFType::const_iterator& VariableLeft, MapDFDFType::const_iterator& VariableRight);   
    
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
    /**
     * @brief sets the number of resulting distribution from arithmetic operations with other distributions
     * 
     */
    void IntegrationSteps(int nSteps) { m_nIntegrationSteps = nSteps;}

    /**
     * @brief returns the number of resulting distribution from arithmetic operations with other distributions
     * 
     * @return int 
     * 
     */
    int IntegrationSteps() const { return m_nIntegrationSteps;}
    
    /**
     * @brief sets the number of integration steps for a single target value for a arithmetic operations with other distributions
     * 
     */
   void SubIntegrationSteps(int nSteps) { m_nSubIntegrationSteps = nSteps;}
    
    /**
     * @brief returns the number of integration steps for a single target value for a arithmetic operations with other distributions
     * 
     * @return int
     * 
     */
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
    CProbabilityDensityDistribution& _GeneralOperatorsFunctionNumerical(CProbabilityDensityDistribution& other, const ProbDistOp Operation);
    
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
   
   /**
     * @brief set CProbabilityDensityDistribution::m_ConvolutionPlan which keeps the integration path for an arithmetic operation for all target values of 
     *   the resulting distribution
     * 
     * @param[in] Other CProbabilityDensityDistribution distribution the arithmetic operation is performed with this distri
     * @param[in] Operation ProbDistOp the arithmetic operation (e.g. sum , subtration, ..)
     * 
     */     
   void _SetConvolutionPlan(CProbabilityDensityDistribution& Other, const ProbDistOp Operation);
   
   
   /**
     * @brief calculates the convolution factor for arithmetic operations and covariance calculations  on the distribution 
     * 
     * @param[in] Operation ::ProbDistOp the operation conducted on the distributions
     * @param[in] pdtLimitsOne ::PairDFType sub-integration limits for one distribution (this)
     * @param[in] dfMeanOne CDigFloat the mean of the one (this) distribution
     * 
     * @return CDigFloat as the convolution factor 
     * 
     */     
  CDigFloat _GetIntegrationFactor4Convolution(const ProbDistOp Operation,const PairDFType pdtLimitsOne, const CDigFloat dfMeanOne);
   
   /**
     * @brief set CProbabilityDensityDistribution::m_ConvolutionPlan for a single target value for all kinds of arithmetic operataions
     * 
     * @param[in] Operation ProbDistOp the arithmetic operation (e.g. sum , subtration, ..)
     * @param[in] dfActTargetValue CDigFloat given target value to derive integration intervals fors
     * @param[in] Other CProbabilityDensityDistribution distribution the arithmetic operation is performed with this distri
     * @param[out] SubIntervals4TargetValue VectorSubIntLimitsType keeps pairwise integration intervals for this distri and other distri
     * for a given target value
     * 
     */      
   void _SetSubIntegrationIntervals4TargetValue(const ProbDistOp Operation, const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& Other, VectorSubIntLimitsType& SubIntervals4TargetValue);

   /**
     * @brief set CProbabilityDensityDistribution::m_ConvolutionPlan for a single target value for arithmetic operataion addition
     * 
     * @param[in] Operation ProbDistOp the arithmetic operation (e.g. sum , subtration, ..)
     * @param[in] dfActTargetValue CDigFloat given target value to derive integration intervals fors
     * @param[in] Other CProbabilityDensityDistribution distribution the arithmetic operation is performed with this distri
     * @param[out] SubIntervals4TargetValue VectorSubIntLimitsType keeps pairwise integration intervals for this distri and other distri
     * for a given target value
     * 
     */      
   void _SetSubIntegrationIntervals4TargetValue4Addition(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervals4TargetValue);

   /**
     * @brief set CProbabilityDensityDistribution::m_ConvolutionPlan for a single target value for arithmetic operataion subtraction
     * 
     * @param[in] Operation ProbDistOp the arithmetic operation (e.g. sum , subtration, ..)
     * @param[in] dfActTargetValue CDigFloat given target value to derive integration intervals fors
     * @param[in] Other CProbabilityDensityDistribution distribution the arithmetic operation is performed with this distri
     * @param[out] SubIntervals4TargetValue VectorSubIntLimitsType keeps pairwise integration intervals for this distri and other distri
     * for a given target value
     * 
     */      
   void _SetSubIntegrationIntervals4TargetValue4Subtraction(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervals4TargetValue);
   
    /**
     * @brief calculates the sub integration intervals needed for multiplication of two random variables (see CProbabilityDensityDistribution::operator*) for a specified resulting random value (called target value)
     * 
     * @param[in] dfActTargetValue CDigFloat target value to find sub integration intervals for
     * @param[in] Other CProbabilityDensityDistribution the distribution which this distri is multiplied with
     * @param[out] SubIntervals4TargetValue VectorSubIntLimitsType keeps the sub integration intervals (see CProbabilityDensityDistribution::SubIntegrationSteps) for this distri and Other needed to calculate the total probability for dfActTargetValue
     * 
     */     
   void _SetSubIntegrationIntervals4TargetValue4Multiplication(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& Other, VectorSubIntLimitsType& SubIntervals4TargetValue);
   
    /**
     * @brief calculates the total limits of the sub integration intervals needed for multiplication of two random variables (see CProbabilityDensityDistribution::operator*) for a specified negative resulting random value (called target value).This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Multiplication.
     * 
     * @param[in] dfActTargetValue CDigFloat target value to find sub integration intervals for
     * @param[in] One CProbabilityDensityDistribution the distribution which Other distri is multiplied with
     * @param[in] Other CProbabilityDensityDistribution the distribution which One distri is multiplied with
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalInterval4NegativeTargetValue4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);
   
    /**
     * @brief calculates the total limits of the sub integration intervals needed for multiplication of positive ranges of two random variables (see CProbabilityDensityDistribution::operator*) for a specified positive resulting random value (called target value). This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Multiplication.
     * 
     * @param[in] dfActTargetValue CDigFloat positive target value to find sub integration intervals for
     * @param[in] One CProbabilityDensityDistribution the distribution which Other distri is multiplied with
     * @param[in] Other CProbabilityDensityDistribution the distribution which One distri is multiplied with
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalInterval4PositiveTargetValuePositiveVariables4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);   
   
    /**
     * @brief calculates the total limits of the sub integration intervals needed for multiplication of negative ranges of two random variables (see CProbabilityDensityDistribution::operator*) for a specified positive resulting random value (called target value).This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Multiplication.
     * 
     * @param[in] dfActTargetValue CDigFloat positive target value to find sub integration intervals for
     * @param[in] One CProbabilityDensityDistribution the distribution which Other distri is multiplied with
     * @param[in] Other CProbabilityDensityDistribution the distribution which One distri is multiplied with
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalInterval4PositiveTargetValueNegativeVariables4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);
   
   /**
     * @brief set CProbabilityDensityDistribution::m_ConvolutionPlan for a single target value for arithmetic operataion Division
     * 
     * @param[in] dfActTargetValue CDigFloat given target value to derive integration intervals fors
     * @param[in] Other CProbabilityDensityDistribution distribution which divides this distri
     * @param[out] SubIntervals4TargetValue VectorSubIntLimitsType keeps pairwise integration intervals for this distri and other distri
     * for a given target value
     * 
     */  
   void _SetSubIntegrationIntervals4TargetValue4Division(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervals4TargetValue);

    /**
     * @brief calculates the total limits of the sub integration intervals needed for division of positive ranges of two random variables (see CProbabilityDensityDistribution::operator/) for a specified positive resulting random value (called target value).This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Division.
     * 
     * @param[in] dfActTargetValue CDigFloat positive target value to find sub integration intervals for
     * @param[in] One CProbabilityDensityDistribution the distribution which is dividec by Other distri 
     * @param[in] Other CProbabilityDensityDistribution the distribution which is the divisor of the One distri 
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalInterval4PositiveTargetValuePositiveVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);

    /**
     * @brief calculates the total limits of the sub integration intervals needed for division of negative ranges of two random variables (see CProbabilityDensityDistribution::operator/) for a specified positive resulting random value (called target value).This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Division.
     * 
     * @param[in] dfActTargetValue CDigFloat positive target value to find sub integration intervals for
     * @param[in] One CProbabilityDensityDistribution the distribution which is dividec by Other distri 
     * @param[in] Other CProbabilityDensityDistribution the distribution which is the divisor of the One distri 
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalInterval4PositiveTargetValueNegativeVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);

    /**
     * @brief calculates the total limits of the sub integration intervals needed for division of positive and negative ranges of two random variables (see CProbabilityDensityDistribution::operator/) for a specified negative resulting random value (called target value).This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Division.
     * 
     * @param[in] dfActTargetValue CDigFloat positive target value to find sub integration intervals for
     * @param[in] One CProbabilityDensityDistribution the distribution which is dividec by Other distri 
     * @param[in] Other CProbabilityDensityDistribution the distribution which is the divisor of the One distri 
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalInterval4NegativeTargetValuePosNegVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);

    /**
     * @brief calculates the total limits of the sub integration intervals needed for division of negative and positive ranges of two random variables (see CProbabilityDensityDistribution::operator/) for a specified negative resulting random value (called target value).This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Division.
     * 
     * @param[in] dfActTargetValue CDigFloat positive target value to find sub integration intervals for
     * @param[in] One CProbabilityDensityDistribution the distribution which is dividec by Other distri 
     * @param[in] Other CProbabilityDensityDistribution the distribution which is the divisor of the One distri 
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalInterval4NegativeTargetValueNegPosVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);
   
   
  
   /**
     * @brief set CProbabilityDensityDistribution::m_ConvolutionPlan for a single target value for covariance calculation
     * 
     * @param[in] dfActTargetValue CDigFloat given target value to derive integration intervals fors
     * @param[in] Other CProbabilityDensityDistribution distribution for which the covariance with this distri is calculated
     * @param[out] SubIntervals4TargetValue VectorSubIntLimitsType keeps pairwise integration intervals for this distri and other distri
     * for a given target value
     * 
     */  
    void _SetSubIntegrationIntervals4TargetValue4Covariance(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervals4TargetValue);
     /**
     * @brief calculates the total limits of the sub integration intervals needed for the covariance of two random variables (see CProbabilityDensityDistribution::CovarianceDistribution) for a specified negative resulting random value (called target value).This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Covariance.
     * 
     * @param[in] dfActTargetValue CDigFloat target value to find sub integration intervals for
     * @param[in] One CProbabilityDensityDistribution the distribution with which the covariance with Other distri is calculated
     * @param[in] Other CProbabilityDensityDistribution the distribution with which the covariance with One distri is calculated
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalInterval4NegativeTargetValue4Covariance(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);
   
    /**
     * @brief calculates the total limits of the sub integration intervals needed for the calculation of the covariance of positive ranges of two random variables (see CProbabilityDensityDistribution::CovarianceDistribution) for a specified positive resulting random value (called target value). This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Multiplication.
     * 
     * @param[in] dfActTargetValue CDigFloat positive target value to find sub integration intervals for
     * @param[in] One CProbabilityDensityDistribution the distribution with which the covariance with Other distri is calculated
     * @param[in] Other CProbabilityDensityDistribution the distribution with which the covariance with One distri is calculated
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalInterval4PositiveTargetValuePositiveVariables4Covariance(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);   
   
    /**
     * @brief calculates the total limits of the sub integration intervals needed for the calculation of covariance of negative ranges of two random variables (see CProbabilityDensityDistribution::CovarianceDistribution) for a specified positive resulting random value (called target value).This function is called by CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Multiplication.
     * 
     * @param[in] dfActTargetValue CDigFloat positive target value to find sub integration intervals for
     * @param[in] One CProbabilityDensityDistribution the distribution with which the covariance with Other distri is calculated
     * @param[in] Other CProbabilityDensityDistribution the distribution with which the covariance with One distri is calculated
     * @param[out] pdtOneLimits PairDFType keeps total integration limits (min, max.) of One distri
     * @param[out] pdtOtherLimits PairDFType keeps total integration limits (min, max.) of Other distri
     * 
     */     
   void _GetSubIntegrationTotalInterval4PositiveTargetValueNegativeVariables4Covariance(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits);
  
   /**
    * @brief: initializes value with zero and resets errors
    * 
    * @param[in] pdtLimits ::PairDFType variable to be reset
    */ 
   void _InitIntegrationLimits(PairDFType& pdtLimits);
   
   ////////////////////////////////////////////////////////
   // functions for analytical calculation of parts of integrals
   ////////////////////////////////////////////////////////
   
   CProbabilityDensityDistribution& _GeneralOperatorsFunctionAnalytical(CProbabilityDensityDistribution& other, const ProbDistOp Operation);

   void _SetConvolutionPlan_2(CProbabilityDensityDistribution& Other, const ProbDistOp Operation);

   void _GetSubIntegrationIntervals4TargetValue_2(const ProbDistOp Operation, const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& Other, VectorSubIntLimitsType&  SubIntervals4TargetValue);
   
   void _GetXInterval4Operation(const ProbDistOp Operation, VectorPairDFType& vXIntervals);
   
   void _GetTotalIntegrationInterval4TargetValue(const ProbDistOp& Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const VectorPairDFType& xVIntervals, VectorSubIntLimitsType&  TotalIntervals4TargetValue);
   
   void _GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue(const ProbDistOp Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const VectorSubIntLimitsType& TotalIntervals4TargetValue, VectorSubIntLimitsType& SubIntervals4TargetValue);
   
   CDigFloat _GetComplementaryVariable(const ProbDistOp Operation, const CDigFloat& dfTargetValue, const CDigFloat& dfVariable, bool b4X);

   CDigFloat _GetConvolutionIntegralAnalytical4Operation(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, const CDigFloat& dfTargetValue, const SubIntLimitsType& vXYLimits);
   
   CDigFloat _PrimitiveIntegral4TargetValue4Addition(const CDigFloat dfOffsetX, const CDigFloat dfSlopeX, const CDigFloat dfOffsetY, const CDigFloat dfSlopeY, const CDigFloat dfTargetValue, const CDigFloat dfX);
   
   CDigFloat _PrimitiveIntegral4TargetValue4Subtraction(const CDigFloat dfOffsetX, const CDigFloat dfSlopeX, const CDigFloat dfOffsetY, const CDigFloat dfSlopeY, const CDigFloat dfTargetValue, const CDigFloat dfX);
   
   CDigFloat _PrimitiveIntegral4TargetValue4Multiplication(const CDigFloat dfOffsetX, const CDigFloat dfSlopeX, const CDigFloat dfOffsetY, const CDigFloat dfSlopeY, const CDigFloat dfTargetValue, const CDigFloat dfX);
   
   CDigFloat _PrimitiveIntegral4TargetValue4Division(const CDigFloat dfOffsetX, const CDigFloat dfSlopeX, const CDigFloat dfOffsetY, const CDigFloat dfSlopeY, const CDigFloat dfTargetValue, const CDigFloat dfX);
   
    
   
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
