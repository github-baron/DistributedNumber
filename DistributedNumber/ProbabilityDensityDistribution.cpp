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

#include "ProbabilityDensityDistribution.h"
#include "iostream"
////////////////////////////////////////////
// construction / destruction
////////////////////////////////////////////
CProbabilityDensityDistribution::CProbabilityDensityDistribution()
{
    _Init();
}

CProbabilityDensityDistribution::CProbabilityDensityDistribution(const CProbabilityDensityDistribution& other)
{
    operator=(other);
}

CProbabilityDensityDistribution::CProbabilityDensityDistribution(const CDistribution& other)
{
    operator=(other);
}

CProbabilityDensityDistribution::~CProbabilityDensityDistribution()
{
    _Init();
}
////////////////////////////////////////////
// assignment operators
////////////////////////////////////////////
CProbabilityDensityDistribution& CProbabilityDensityDistribution::operator=(const CProbabilityDensityDistribution& other)
{
    CDistribution::operator=(other);
    bNormalized = other.bNormalized;
    dfAbsIntegral = other.dfAbsIntegral;
    m_nIntegrationSteps = other.IntegrationSteps();
    m_nSubIntegrationSteps = other.SubIntegrationSteps();
    return *this;
}
CProbabilityDensityDistribution& CProbabilityDensityDistribution::operator=(const CDistribution& other)
{
    _Init();
    CDistribution::operator=(other);
    return *this;
}

////////////////////////////////////////////
// arithmetic operators with CProbabilityDensityDistribution
////////////////////////////////////////////
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator+=(CProbabilityDensityDistribution& other)
{
    
//     return _GeneralOperatorsFunction(other, ProbDistOp::pdoPlus );
//     return _Convolution4GeneralOperatorsNew(other, ProbDistOp::pdoPlus );
//     return _GeneralOperatorsFunctionNumerical(other, ProbDistOp::pdoPlus );
    return _GeneralOperatorsFunctionAnalytical(other, ProbDistOp::pdoPlus );
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator+(CProbabilityDensityDistribution& other)
{
    CProbabilityDensityDistribution ProbDistResult(*this);
    ProbDistResult += other;
    return ProbDistResult;
}

CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator-=(CProbabilityDensityDistribution& other)
{
//     return _Convolution4GeneralOperatorsNew(other, ProbDistOp::pdoMinus );
//     return _GeneralOperatorsFunctionNumerical(other, ProbDistOp::pdoMinus );
    return _GeneralOperatorsFunctionAnalytical(other, ProbDistOp::pdoMinus );
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator-(CProbabilityDensityDistribution& other)
{
    CProbabilityDensityDistribution ProbDistResult(*this);
    ProbDistResult -= other;
    return ProbDistResult;
}


CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator*=(CProbabilityDensityDistribution& other)
{
//     return _GeneralOperatorsFunction(other, ProbDistOp::pdoMult );
//     return _Convolution4GeneralOperatorsNew(other, ProbDistOp::pdoMult );
//     return _GeneralOperatorsFunctionNumerical(other, ProbDistOp::pdoMult );
    return _GeneralOperatorsFunctionAnalytical(other, ProbDistOp::pdoMult );
}
CProbabilityDensityDistribution CProbabilityDensityDistribution::operator*(CProbabilityDensityDistribution& other)
{
    CProbabilityDensityDistribution ProbDistResult(*this);
    ProbDistResult *= other;
    return ProbDistResult;
}

CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator/=(CProbabilityDensityDistribution& other)
{
//     return _Convolution4GeneralOperatorsNew(other, ProbDistOp::pdoDiv);
//     return _GeneralOperatorsFunctionNumerical(other, ProbDistOp::pdoDiv);
    return _GeneralOperatorsFunctionAnalytical(other, ProbDistOp::pdoDiv);
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator/(CProbabilityDensityDistribution& other)
{
    CProbabilityDensityDistribution ProbDistResult(*this);
    ProbDistResult /= other;
    return ProbDistResult;
}

////////////////////////////////////////////
// arithmetic operators with CDigFloat
////////////////////////////////////////////
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator+=(const CDigFloat& Value)
{
    // _DeNormalize
    DeNormalize();
    
    CDistribution::operator+=(Value);
    
    // normalize again
    Normalize();
    
    return *this;
}
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator-=(const CDigFloat& Value)
{
    // _DeNormalize
    DeNormalize();
    
    CDistribution::operator-=(Value);
    
    // normalize again
    Normalize();
    
    return *this;
}
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator*=(const CDigFloat& Value)
{
    // simply do nothing: due to normalization the factor will be cancelled    
    return *this;
}
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator/=(const CDigFloat& Value)
{
    // simply do nothing: due to normalization the divisor will be cancelled
    return *this;
}

////////////////////////////////////////////
// comparison operators
////////////////////////////////////////////
bool CProbabilityDensityDistribution::operator==(const CProbabilityDensityDistribution& other) const
{
    // at the moment the distribution is the relevant part:
    // dfAbsIntegral can be calculated any time and is the same if the distributions are the same
    // bNormalized is hidden from user and will be performed when necessary
    return CDistribution::operator==(other);
}

bool CProbabilityDensityDistribution::operator!=(const CProbabilityDensityDistribution& other) const
{
    return !operator==(other);
}
    
////////////////////////////////////////////
// distribution generation
////////////////////////////////////////////

void CProbabilityDensityDistribution::ConstDistribution(const CDigFloat& dfXMin, const CDigFloat& dfXMax)
{
    // do nothing in case xmin >= xmax
    if( dfXMin >= dfXMax)
        return;
    
    // clear the actual distribution
    Reset();
    
    // add the points
    Add(dfXMin,1);
    Add(dfXMax,1);
    
    // normalize
    Normalize();
}

void CProbabilityDensityDistribution::TriangleDistribution(const CDigFloat& dfXMin, const CDigFloat& dfXMax)
{
    // do nothing in case xmin >= xmax
    if( dfXMin >= dfXMax)
        return;
    
    // clear the actual distribution
    Reset();
    
    // add the points 
    Add(dfXMin,0);
    Add((dfXMax+dfXMin)/2.,1);
    Add(dfXMax,0);
    
    // normalize
    Normalize();
    
}

void CProbabilityDensityDistribution::NormalDistribution(const int nPoints, const CDigFloat& dfMu, const CDigFloat& dfSigma, const CDigFloat& dfMinVal)
{
    // check consistent parameters:
    if(dfSigma <= 0 || dfMinVal <= 0)
        return;
    
    // reset
    Reset();
    
    CDigFloat dfFac1 = 2*dfSigma*dfSigma;
    CDigFloat dfXDelta = sqrt(dfFac1*-log(dfMinVal*sqrt(dfFac1*M_PI)));
    CDigFloat dfXInc = 2*dfXDelta/(nPoints-1);
    
    for(int ipt=0;ipt<nPoints; ipt++)
    {
        CDigFloat dfX = dfMu-dfXDelta+ipt*dfXInc;
        Add(dfX,1/sqrt(dfFac1*M_PI)*pow(M_El, -pow(dfX-dfMu,2)/dfFac1));
    }

    Normalize();
}

void CProbabilityDensityDistribution::CovarianceDistribution(CProbabilityDensityDistribution& Other)
{
//     _GeneralOperatorsFunction(Other, ProbDistOp::pdoCov);
    CProbabilityDensityDistribution pddOther(Other);
    pddOther.Shift(-Other.Mean());
    Shift(-Mean());
    operator*=(pddOther);
}

///////////////////////////////////
// function on all distribution variables
///////////////////////////////////    
void CProbabilityDensityDistribution::Shift(const CDigFloat& shift)
{
    // this can simply applied as the normalization is not affected 
    CDistribution::Shift(shift);
    
}
void CProbabilityDensityDistribution::Scale(const CDigFloat& scale)
{
    // scaling needs denormalization and normalization 
    DeNormalize();
    
    CDistribution::Scale(scale);
    
    Normalize();
}


////////////////////////////////////////////
// functions
////////////////////////////////////////////

CDigFloat CProbabilityDensityDistribution::Variance(int nthOrder)
{
    // copy from this
    CProbabilityDensityDistribution pdd4Calc(*this);
    
    // shift by mean: mean is now zero because of Integral (x-<x>) = 0
    pdd4Calc.Shift(-pdd4Calc.Mean());
    
    // return the second moment which is (x-<x>)²
    return pdd4Calc.Mean(2);
}

CDigFloat CProbabilityDensityDistribution::Covariance(CProbabilityDensityDistribution& Other)
{
    // copy from this
    CProbabilityDensityDistribution pdd4Calc(*this);
    
     
    pdd4Calc*=Other;
    
    // return mean: should be now <(x-<x>)*(y-<y>)>
    return pdd4Calc.Mean() - Other.Mean()*Mean();
}

CDigFloat CProbabilityDensityDistribution::DistriValue(const CDigFloat& variable)
{
    Normalize();
    return CDistribution::DistriValue(variable);
}
CDigFloat CProbabilityDensityDistribution::AbsIntegral(const CDigFloat& variableLeft, const CDigFloat& variableRight, int nthOrder)
{
    Normalize();
    return CDistribution::AbsIntegral(variableLeft,variableRight,nthOrder);
}
CDigFloat CProbabilityDensityDistribution::AbsIntegral(const int& nthOrder)
{
    Normalize();
    return CDistribution::AbsIntegral(nthOrder);
}
void CProbabilityDensityDistribution::GetInterval(const CDigFloat& variable, MapDFDFType::const_iterator& VariableLeft, MapDFDFType::const_iterator& VariableRight)
{
    Normalize();
    return CDistribution::GetInterval(variable,VariableLeft,VariableRight);
}
CDigFloat CProbabilityDensityDistribution::Mean(int nthOrder)
{
    Normalize();
    return CDistribution::Mean(nthOrder);
}
CDigFloat CProbabilityDensityDistribution::OrigIntegral()
{
    Normalize();
    return dfAbsIntegral;
}

void CProbabilityDensityDistribution::Reset()
{
    _Init();
}
bool CProbabilityDensityDistribution::Add(const CDigFloat& variable, const CDigFloat value, bool bLastElement /*== false*/)
{
    // first check if the distribution is already normalized
    DeNormalize();
    
    // add the new value
    bool bSuccess = CDistribution::Add(variable,value);
    
    // check if this is the last element: if so --> normalize
    if(bLastElement)
        Normalize();
    
    return bSuccess;
}

void CProbabilityDensityDistribution::DeNormalize()
{
    // do not normalize again ... only once
    if(!bNormalized)
        return;
 
    // revert normalization
    CDistribution::operator*=(dfAbsIntegral);
    bNormalized = false;
    
}

void CProbabilityDensityDistribution::Normalize()
{
    // do not normalize again ... only once
    if(bNormalized)
        return;
    
    dfAbsIntegral = CDistribution::AbsIntegral();
    CDistribution::operator/=(dfAbsIntegral);
    
    // remember we have been normalized
    bNormalized = true;
    
}
void CProbabilityDensityDistribution::_Init()
{
    dfAbsIntegral = 0;
    bNormalized = false;
    IntegrationSteps(PDD_DEFAULT_INTEGRATION_STEPS);
    SubIntegrationSteps(PDD_DEFAULT_SUBINTEGRATION_STEPS);
    CDistribution::_Init();

}

CProbabilityDensityDistribution & CProbabilityDensityDistribution::_GeneralOperatorsFunctionAnalytical(CProbabilityDensityDistribution& Other, const ProbDistOp Operation)
{
    
    // generate plan : setting m_ConvolutionPlan;
    _SetConvolutionPlan_2(Other, Operation);
    
    // declare the result distri
    MapDFDFType TargetDistri;
    
    for(auto itarget: m_ConvolutionPlan)
    {
        TargetDistri[itarget.first] = 0;
        
        // DEBUG
        cout << "itarget.first: " << itarget.first.RawPrint(15) << endl;
        
        // iterate over all sub integration limits
        for(auto intlim: itarget.second)
        {
            
            cout << "x = [ " << intlim.first.first.RawPrint(10) << " ; " << intlim.first.second.RawPrint(10) << " ] " << endl;
            cout << "y = [ " << intlim.second.first.RawPrint(10) << " ; " << intlim.second.second.RawPrint(10) << " ] " << endl;
            
            // add to all other probs of this target value
            TargetDistri[itarget.first] += _GetConvolutionIntegralAnalytical4Operation(Other, Operation, itarget.first, intlim);
            
        }   // endfor(auto interval: itarget.second)        
        
    }   // endfor(auto itarget: m_ConvolutionPlan)
    
    _Init();
    
    
    if(m_ConvolutionPlan.size() > 0)
    {
        for(auto iel: TargetDistri)
            mDistribution[iel.first] = iel.second;
        
        
        // DEBUG
//         cout << "before Normalization: " << endl << Print(10);
        
        
        bNormalized=false;
        Normalize();
   
        // DEBUG     
//         cout << "after Normalization: " << endl << Print(10);
   
        assert(AbsIntegral() == 1);
        
    }   // endif(!(TargetRangeStart == TargetRangeEnd && TargetRangeStart == 0))
    
    return *this;
 
}

CDigFloat CProbabilityDensityDistribution::_GetConvolutionIntegralAnalytical4Operation(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, const CDigFloat& dfTargetValue, const SubIntLimitsType& vXYLimits)
{
    // DEBUG
    cout << "_GetConvolutionIntegralAnalytical4Operation: " << endl;
    
    // init the result
    CDigFloat dfResult = 0;
    
    // derive all necessary arguments needed for the analytical calculation of the 
    // convolution integrals
    CDigFloat dfXMean = (vXYLimits.first.first + vXYLimits.first.second)/2.;
    CDigFloat dfYMean = (vXYLimits.second.first + vXYLimits.second.second)/2.;
    CDigFloat dfXOffset = LinearOffset(dfXMean);
    CDigFloat dfXSlope = LinearSlope(dfXMean);
    CDigFloat dfYOffset = Other.LinearOffset(dfYMean);
    CDigFloat dfYSlope = Other.LinearSlope(dfYMean);
    
    cout << "X = [" << vXYLimits.first.first.RawPrint(10) << " ; " << vXYLimits.first.second.RawPrint(10) << " ] "  << endl;
    cout << "dfXMean = " << dfXMean.RawPrint(10) << endl;
    cout << "dfXOffset = " << dfXOffset.RawPrint(10) << endl;
    cout << "dfXSlope = " << dfXSlope.RawPrint(10) << endl;
    cout << "Y = [" << vXYLimits.second.first.RawPrint(10) << " ; " << vXYLimits.second.second.RawPrint(10) << " ] "  << endl;
    cout << "dfYMean = " << dfYMean.RawPrint(10) << endl;
    cout << "dfYOffset = " << dfYOffset.RawPrint(10) << endl;
    cout << "dfYSlope = " << dfYSlope.RawPrint(10) << endl;
    
    switch( Operation)
    {
        case ProbDistOp::pdoPlus:
            dfResult = _PrimitiveIntegral4TargetValue4Addition(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.second) -
                       _PrimitiveIntegral4TargetValue4Addition(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.first);
            break;
        case ProbDistOp::pdoMinus:
            dfResult = _PrimitiveIntegral4TargetValue4Subtraction(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.second) -
                       _PrimitiveIntegral4TargetValue4Subtraction(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.first);
            break;
        case ProbDistOp::pdoMult:
            dfResult = _PrimitiveIntegral4TargetValue4Multiplication(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.second) -
                       _PrimitiveIntegral4TargetValue4Multiplication(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.first);
            break;
        case ProbDistOp::pdoDiv:
            dfResult = _PrimitiveIntegral4TargetValue4Division(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.second) -
                       _PrimitiveIntegral4TargetValue4Division(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.first);
            break;
        default:
            break;
    }   //endswitch( Operation)
    
    cout << "dfResult = " << dfResult.RawPrint(10) << endl;
//     assert(dfResult > 0);
    return dfResult;
    
}

void CProbabilityDensityDistribution::_SetConvolutionPlan_2(CProbabilityDensityDistribution& Other, const ProbDistOp Operation)
{   
    
    CDigFloat dfTargetStart, dfTargetEnd;
    _GetRangeFromDistributionOperation(Other, Operation, dfTargetStart, dfTargetEnd);
    
    // get the target step
    CDigFloat dfTargetStep = (dfTargetEnd - dfTargetStart)/(CDigFloat(IntegrationSteps()));
    
    // DEBUG
//     cout << "target interval: [" << dfTargetStart.RawPrint(15) << " ; " << dfTargetEnd.RawPrint(15) << " ] " << endl;
    
    // clear the plan first
    m_ConvolutionPlan.clear();
    
    // now start building the plan 
    CDigFloat dfActTargetValue = dfTargetStart;
    while(dfActTargetValue <= dfTargetEnd)
    {
        // this is the vector defining the sub integration limits for a given targe value
        VectorSubIntLimitsType SubIntLimits;
        
        // defining the complete integration range for this and other distribution
        _GetSubIntegrationIntervals4TargetValue_2(Operation, dfActTargetValue, Other, SubIntLimits);
        
        // now add the sub integration limits and its corresponding target value as pair
        m_ConvolutionPlan.push_back(pair<CDigFloat, VectorSubIntLimitsType>(dfActTargetValue, SubIntLimits));
        
        // increment actual target value 
        dfActTargetValue += dfTargetStep;
    }   // endwhile(dfActTargetValue <= dfTargetEnd)
    
}

CDigFloat CProbabilityDensityDistribution::_PrimitiveIntegral4TargetValue4Addition(const CDigFloat dfOffsetX, const CDigFloat dfSlopeX, const CDigFloat dfOffsetY, const CDigFloat dfSlopeY, const CDigFloat dfTargetValue, const CDigFloat dfX)
{  
    CDigFloat dfYConst = dfSlopeY*dfTargetValue+dfOffsetY;
    return -dfSlopeX*dfSlopeY/3.*pow(dfX,3)+
            (dfSlopeX*dfYConst-dfSlopeY*dfOffsetX)/2.*dfX*dfX+
            dfYConst*dfOffsetX*dfX;
}
CDigFloat CProbabilityDensityDistribution::_PrimitiveIntegral4TargetValue4Subtraction(const CDigFloat dfOffsetX, const CDigFloat dfSlopeX, const CDigFloat dfOffsetY, const CDigFloat dfSlopeY, const CDigFloat dfTargetValue, const CDigFloat dfX)
{
    CDigFloat dfYConst = -dfSlopeY*dfTargetValue+dfOffsetY;
    return dfSlopeX*dfSlopeY/3.*pow(dfX,3)+
            (dfSlopeX*dfYConst+dfSlopeY*dfOffsetX)/2.*dfX*dfX+
            dfYConst*dfOffsetX*dfX;
}

CDigFloat CProbabilityDensityDistribution::_PrimitiveIntegral4TargetValue4Multiplication(const CDigFloat dfOffsetX, const CDigFloat dfSlopeX, const CDigFloat dfOffsetY, const CDigFloat dfSlopeY, const CDigFloat dfTargetValue, const CDigFloat dfX)
{
    return (dfSlopeX*dfSlopeY*dfTargetValue + dfOffsetX*dfOffsetY)*log(abs(dfX))*abs(dfX)/dfX -
           dfOffsetX*dfSlopeY*dfTargetValue/abs(dfX)+
           dfOffsetY*dfSlopeX*abs(dfX);
}

CDigFloat CProbabilityDensityDistribution::_PrimitiveIntegral4TargetValue4Division(const CDigFloat dfOffsetX, const CDigFloat dfSlopeX, const CDigFloat dfOffsetY, const CDigFloat dfSlopeY, const CDigFloat dfTargetValue, const CDigFloat dfX)
{
    CDigFloat dfPreFac = (dfX < 0) ? -1 : 1;
    return abs(dfX)*dfX*(
        dfSlopeX * dfSlopeY/dfTargetValue/ 4. * dfX*dfX +
        (dfSlopeX*dfOffsetY + dfSlopeY*dfOffsetX/dfTargetValue)/3. *dfX +
        dfOffsetX*dfOffsetY/2.);
}



void CProbabilityDensityDistribution::_GetSubIntegrationIntervals4TargetValue_2(const ProbDistOp Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, VectorSubIntLimitsType& SubIntervals4TargetValue)
{
    
    /////////////////////////////////////////////////////
    // General purpose ot this function:
    //   Determine a valid total integration interval and
    //   its subintervals
    //
    // Procedure:
    //
    //   I Calculation of the total integration interval:
    //     1. calculate the x-interval:
    //        must be unipolar (positive or negative) for 
    //        multiplication and division
    //     2. derive closed intervals for x and y mapped 
    //        onto each other bijectively by the calculation 
    //        of the complementary variable
    //   
    //    II calcuations of the subintervals with constant
    //       slope and offset for x and y for the alls total
    //       intervals
    //
    // Explanation of  notations:
    //   - x-variable: variable of this distribution
    //   - y-variable: variable of Other distribution
    //   - tv: target value (see argument dfTargetValue)
    //   - complementary variable: e.g. for multiplication
    //     the following relation holds: xt*yt = tv
    //     xt and yt are called complementary. 
    //     "Complementarity" depends on the operation: e.g. 
    //        summation xt+yt = tv
    //        subraction: xt-yt = tv
    //        division: xt/yt = tv
    /////////////////////////////////////////////////////
    
    // keeps the total integration interval(s) for a target value
    VectorSubIntLimitsType TotalIntervals4TargetValue;
    
    // fill the x intervals: must be unipolar for multiplication and division
    VectorPairDFType vpdXIntervals;
    
    // DEBUG
    cout << "_SetIntegrationIntervals4TargetValue" << endl;
    cout << "Operation = " <<  to_string((int)Operation ) << endl;
    cout << "dfTargetValue = " << dfTargetValue.RawPrint(10) << endl;
    
    // I 1. get the x-interval(s) for calculation: must be unipolar
    //      for division and multiplication
    _GetXInterval4Operation(Operation, vpdXIntervals);
        
    // I 2. derive closed intervals for x and y (kept in TotalIntervals4TargetValue) mapped onto each other bijectively
    //      by the calculation of the complementary variable
    _GetTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, vpdXIntervals, TotalIntervals4TargetValue);
    
    
    _GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, TotalIntervals4TargetValue, SubIntervals4TargetValue);
    
    
 
    //DEBUG
    cout << endl << endl;
}

void CProbabilityDensityDistribution::_GetXInterval4Operation(const ProbDistOp Operation, VectorPairDFType& vXIntervals)
{     
    
    /////////////////////////////////////////////////////
    // General purpose ot this function:
    //   return the x - interval for the given operation.
    //   In case of multiplication / division and this 
    //   is a bipolar distribution two unipolar intervals
    //   will be returned. In all other cases the complete
    //   interval of x variable is returned.
    //
    // Procedure:
    //   1. set first variable as first lower limit
    //   2. for mult. / div:
    //      if bipolar:
    //       - finalize first interval setting upper limit
    //         to -epsilon (epsilon is close to zero)
    //       - add first interval to result
    //       - start second interval wiht +epsilon
    //   3. finalize interval setting upper limit to
    //      last value
    //   4. if upper limit > lower limit:
    //      add this to result
    //
    // Explanation:
    //   - bipolar distribution:
    //     a distribution which has non-zero values for
    //     negative AND positive variable values.
    //   - unipolar interval:
    //     an interval of unique sign of its values
    /////////////////////////////////////////////////////
    
    // clear result variable
    vXIntervals.clear();
    
    PairDFType pdXInterval;
    pdXInterval.first = firstVariable();
   // take care for unipolar intervals for multiplication and division
    if(Operation == ProbDistOp::pdoMult || Operation == ProbDistOp::pdoDiv)
    {
        // do we have a negative unipolar interval ..
        if(firstVariable() < 0 && lastVariable() >= 0)
        {
            // finish the negative unipolar interval with something close to zero ...
            pdXInterval.second = -CDigFloat(1).Error()*10;
            
            // ... add it
            vXIntervals.push_back(pdXInterval);
            
            // ... start the second interval with something positive close to zero
            pdXInterval.first = - pdXInterval.second;
            
        }   //endif(firstVariable() < 0 && lastVariable() >= 0)
        
    }   //endif(Operation == ProbDistOp::pdoMult || Operation == ProbDistOp::pdoDiv)
    
    // now finish for all cases
    if( lastVariable() > pdXInterval.first)
    {
        // finalize the interval for different cases:
        // a) summation / subtraction: simple x interval
        // b) division / multiplication:
        //    i) simple x-interval in case it is unipolar (purely negative / positive)
        //    ii) divided x-interval into two unipolar intervals in case it bipolar (positive and negative)
        pdXInterval.second = lastVariable();
        vXIntervals.push_back(pdXInterval);
        
    }   //endif( lastVariable() > xInterval.first)
    
    // DEBUG
    cout << "_GetXInterval4Operation:" << endl;
    for(auto iint: vXIntervals)
        cout << "[ " << iint.first.RawPrint(10) << " ; " << iint.second.RawPrint(10) << " ]" << endl;
    cout << endl;

}

void CProbabilityDensityDistribution::_GetTotalIntegrationInterval4TargetValue(const ProbDistOp& Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const VectorPairDFType& vpdXIntervals, VectorSubIntLimitsType& TotalIntervals4TargetValue)
{
    /////////////////////////////////////////////////////
    // General purpose ot this function:
    //   returns the total integration interval(s) for x
    //   and y needed for the analytical calculation 
    //   of the convolution for a given operation by the functions
    // 
    //     _PrimitiveIntegral4TargetValue4...
    //         ...Addition
    //         ...Subtraction
    //         ...Multiplication
    //         ...Division
    //
    //   Each interval pair for x and y is mapped onto 
    //   each other bijectively by the calculation of 
    //   the complementary variable.
    //
    // Procedure:
    //   iterate over the x-interval(s) given by vpdXIntervals:
    //   determine valid interval for the given 
    //   target value doing the following calculations
    //   starting on each x-interval given by vpdXIntervals:
    //     a) calculate the complementary interval YComp from
    //        the x interval
    //     b) calculate the overlapping interval YOverlap of
    //        YComp with the interval of variables of y
    //     c) calculate the Complementarity interval of
    //        YOverlap (depends on operation: see _GetComplementaryVariable)
    //     d) calculate the overlapping interval of 
    //        YOverlap with the starting x-interval
    //     the resulting interval is the total integration
    //     interval used for a given target value which is
    //     handed over to the functions
    /////////////////////////////////////////////////////
    
    // clear result variable
    TotalIntervals4TargetValue.clear();
    
    // so iterate over the x-intervals and determine the integration limits
    for(auto xint: vpdXIntervals)
    {
        //DEBUG
        cout << "xint.first = " << xint.first.RawPrint(10) << endl;
        cout << "xint.second = " << xint.second.RawPrint(10) << endl;
        
        /////////////////////////////////////////////////
        // calculate the complementary interval YComp from x
        /////////////////////////////////////////////////
        CDigFloat dfYCompMin, dfYCompMax;
        dfYCompMin = _GetComplementaryVariable(Operation, dfTargetValue,xint.first,false);
        dfYCompMax = _GetComplementaryVariable(Operation, dfTargetValue,xint.second,false);
        
        // swapping here avoids a lot of if-cases:
        // the if cases would have to consider the sign of target value, x value, y value, and operation
        if(dfYCompMax < dfYCompMin)
            swap(dfYCompMin, dfYCompMax);        
        
        //DEBUG
        cout << "dfYCompMin = " << dfYCompMin.RawPrint(10) << endl;
        cout << "dfYCompMax = " << dfYCompMax.RawPrint(10) << endl;
        
        
        /////////////////////////////////////////////////
        // calculate the overlap of the complementary 
        // ycomp and the original y interval --> YOverlap
        /////////////////////////////////////////////////
        CDigFloat dfYOverlapMin, dfYOverlapMax;
        dfYOverlapMin = max(dfYCompMin,Other.firstVariable());
        dfYOverlapMax = min(dfYCompMax,Other.lastVariable());  
        
        //DEBUG
        cout << "dfYOverlapMin = " << dfYOverlapMin.RawPrint(10) << endl;
        cout << "dfYOverlapMax = " << dfYOverlapMax.RawPrint(10) << endl;
        
        // the overlap is empty in case dfYOverlapMax <= dfYOverlapMin
        if(dfYOverlapMax <= dfYOverlapMin)
        {
            // DEBUG
            cout << "empty interval: dfYOverlapMax(" << dfYOverlapMax.RawPrint(4) << ")  <= dfYOverlapMin(" << dfYOverlapMin.RawPrint(4) << ")" << endl;
            continue;
        }
        
        /////////////////////////////////////////////////
        // calculate the complementary interval XComp from YOverlap
        /////////////////////////////////////////////////
        CDigFloat dfXCompMin, dfXCompMax;
        dfXCompMin = _GetComplementaryVariable(Operation, dfTargetValue,dfYOverlapMin,true);
        dfXCompMax = _GetComplementaryVariable(Operation, dfTargetValue,dfYOverlapMax,true);
        
        // swapping here avoids a lot of if-cases:
        // the if cases would have to consider the sign of target value, x value, y value, and operation
        if(dfXCompMax < dfXCompMin)
            swap(dfXCompMin, dfXCompMax);        
        
        //DEBUG
        cout << "dfXCompMin = " << dfXCompMin.RawPrint(10) << endl;
        cout << "dfXCompMax = " << dfXCompMax.RawPrint(10) << endl;
        
        /////////////////////////////////////////////////
        // calculate the overlap of the complementary 
        // Xcomp and the original x interval --> XOverlap 
        /////////////////////////////////////////////////
        CDigFloat dfXOverlapMin, dfXOverlapMax;
        dfXOverlapMin = max(dfXCompMin,firstVariable());
        dfXOverlapMax = min(dfXCompMax,lastVariable());
        
        // the overlap is empty in case dfXOverlapMax <= dfXOverlapMin
        if(dfXOverlapMax <= dfXOverlapMin)
        {
            // DEBUG
            cout << "empty interval: dfXOverlapMax (" << dfXOverlapMax.RawPrint(4) << ")  <= dfXOverlapMin(" << dfXOverlapMin.RawPrint(4) << ")" << endl;
            continue;
        }
        //DEBUG
        cout << "dfXOverlapMin = " << dfXOverlapMin.RawPrint(10) << endl;
        cout << "dfXOverlapMax = " << dfXOverlapMax.RawPrint(10) << endl;
        
        cout << "X=[ " <<  dfXOverlapMin.RawPrint(10) << " ; " << dfXOverlapMax.RawPrint(10) << " ] " << endl;
        cout << "Y=[ " <<  dfYOverlapMin.RawPrint(10) << " ; " << dfYOverlapMax.RawPrint(10) << " ] " << endl;
         
        /////////////////////////////////////////////////
        // calculate the overlap of the complementary 
        // Xcomp and the original x interval --> XOverlap
        /////////////////////////////////////////////////
        // assure that the limits are pairwise complementary
        assert( (dfXOverlapMin == _GetComplementaryVariable(Operation,dfTargetValue,dfYOverlapMax,true) && dfXOverlapMax == _GetComplementaryVariable(Operation,dfTargetValue,dfYOverlapMin,true)) || 
        (dfXOverlapMin == _GetComplementaryVariable(Operation,dfTargetValue,dfYOverlapMin,true) && dfXOverlapMax == _GetComplementaryVariable(Operation,dfTargetValue,dfYOverlapMax,true)));
        
        // add for calculation
        TotalIntervals4TargetValue.push_back(SubIntLimitsType(PairDFType(dfXOverlapMin,dfXOverlapMax),PairDFType(dfYOverlapMin,dfYOverlapMax))); 
        
        
    }   //endfor(auto xint: xIntervals)

}

void CProbabilityDensityDistribution::_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue(const ProbDistOp Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const VectorSubIntLimitsType& TotalIntervals4TargetValue, VectorSubIntLimitsType& SubIntervals4TargetValue)
 {   
     // DEBUG
     cout << "_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue:" << endl;
     
    /////////////////////////////////////////////////////
    // General purpose ot this function:
    //   returns the sub integration interval(s) for x
    //   and y basing on the total integration intervals.
    // 
    // Explanation:
    //   The sub integration intervals keep ranges of constant
    //   slope and offset for x AND y distri. This is necessary
    //   because the convolution functions for a given operation 
    //  
    // 
    //     _PrimitiveIntegral4TargetValue4...
    //         ...Addition
    //         ...Subtraction
    //         ...Multiplication
    //         ...Division
    //   expect slope and offset as argument.
    //
    // Procedure:
    //   calcuations of the subintervals with constant
    //   slope and offset for x and y for the given total
    //   interval :
    //   At this point we are not finished yet: the
    //   functions mentioned in 1d need the slope and 
    //   the offset of the corresponding linear part 
    //   of a distribution.
    //   Therefore, the total interval for integration
    //   must be divided into subparts with constant
    //   slope and offset for x and y distri.
    //   To derive the sub intervals the following is
    //   done for each interval: 
    //   1. getting all supporting points of x distribution
    //      lying within the given interval
    //   2. getting all supporting points of y distribution
    //      (as complementary values) lying within the given 
    //      interval
    //   3. sorting all the supporting points: x and complementary
    //      y values
    //   4. deriving the corresponding integration intervals
    //      for x and y from the intervals of the sorted points
    /////////////////////////////////////////////////////
    
    
    // now set the subintervals for each total interval
    for(auto xyint: TotalIntervals4TargetValue)
    {
        //DEBUG
        cout << "x = [" << xyint.first.first.RawPrint(10) << " ; " << xyint.first.second.RawPrint(10) << " ]" << endl;
        cout << "y = [" << xyint.second.first.RawPrint(10) << " ; " << xyint.second.second.RawPrint(10) << " ]" << endl;
        
        //////////////////////////////////////////////
        // get all relevant supporting points of the x
        // distri
        //////////////////////////////////////////////
        vector<CDigFloat> vdfSubPoints;
        
        // get the limits of the total integral
        vdfSubPoints.push_back(xyint.first.first);
        vdfSubPoints.push_back(xyint.first.second);
        
        // iterate this for supporting points within the total interval
        for(auto x: Distribution())
        {
            // DEBUG
            cout << "x=" << x.first.RawPrint(10) << endl;
            
            // break conditions
            if(x.first > xyint.first.second)
                break;
            
            if(x.first > xyint.first.first && x.first <xyint.first.second)
            {
            
                vdfSubPoints.push_back(x.first);
                
                //DEBUG
                cout << "added " << x.first.RawPrint(10) << endl;
                
            }
            
        }   //endfor(auto x: Distribution())        
        
        //////////////////////////////////////////////
        // get all relevant supporting points of the y
        // distri: get the Complementary values
        //////////////////////////////////////////////
        // iterate this for supporting points within the total interval
        for(auto x: Other.Distribution())
        {
            // DEBUG
            cout << "y=" << x.first.RawPrint(10) << endl;
            
            // break conditions
            if(x.first > xyint.second.second)
                break;
            
            if(x.first > xyint.second.first && x.first <xyint.second.second)
            {
                vdfSubPoints.push_back(_GetComplementaryVariable(Operation,dfTargetValue,x.first,true));
                
                //DEBUG
                cout << "added " << _GetComplementaryVariable(Operation,dfTargetValue,x.first,true).RawPrint(10) << endl;
            }
            
        }   //endfor(auto x: Other.Distribution())
        
        //////////////////////////////////////////////
        // sort the supporting points and derive
        // sub intervals for x and y
        //////////////////////////////////////////////
        sort(vdfSubPoints.begin(),vdfSubPoints.end());
        //DEBUG
        for(auto x: vdfSubPoints)
            cout << "SubPoint: " << x.RawPrint(10) << endl;
        
        // iterate over sorted sub points
        for(int idx=0; idx < vdfSubPoints.size()-1; idx++)
        {
            // keeping the limits for x and y
            PairDFType dfXSubLimit, dfYSubLimit;
            
            // set the x limits
            dfXSubLimit.first = vdfSubPoints[idx];
            dfXSubLimit.second = vdfSubPoints[idx+1];
            
            // derive the y limits
            dfYSubLimit.first = _GetComplementaryVariable(Operation, dfTargetValue,dfXSubLimit.first,false);
            dfYSubLimit.second = _GetComplementaryVariable(Operation, dfTargetValue,dfXSubLimit.second,false);
            
            // swapping here avoids a lot of if-cases:
            // the if cases would have to consider the sign of target value, x value, y value, and operation
            if(dfYSubLimit.second < dfYSubLimit.first )
                swap(dfYSubLimit.first, dfYSubLimit.second); 
            
            // now add the pair of limits (xmin,xmax) (ymin, ymax) to result vector
            SubIntervals4TargetValue.push_back(SubIntLimitsType(dfXSubLimit, dfYSubLimit));            
            
            // DEBUG
            cout << "adding:" << endl;
            cout << "x = [" << dfXSubLimit.first.RawPrint(10) << " ; " << dfXSubLimit.second.RawPrint(10) << " ]" << endl;
            cout << "y = [" << dfYSubLimit.first.RawPrint(10) << " ; " << dfYSubLimit.second.RawPrint(10) << " ]" << endl;
            
        }   //endfor(int idx=0; idx < vdfSubPoints.size()-1; idx++)
        
        
    }   //endfor(auto xyint: TotalIntervals4TargetValue)
    
}


CDigFloat CProbabilityDensityDistribution::_GetComplementaryVariable(const ProbDistOp Operation, const CDigFloat& dfTargetValue, const CDigFloat& dfVariable, bool bVarIsYVar)
{
    /////////////////////////////////////////////////////
    // General purpose ot this function:
    //   returns the complementary variable for a given
    //   variable, operation, and target value (tv).
    //
    // Explanation:
    //   Complementarity in this content means that
    //   two variables x and y are connected by an
    //   expression.
    //   This expression is for:
    //   Addition: x+y=tv
    //   Subtraction: x-y=tv
    //   Multiplication: x*y=tv
    //   Division: x/y=tv
    /////////////////////////////////////////////////////
    
    CDigFloat dfComplementaryVariable;
    
    switch(Operation)
    {
        case ProbDistOp::pdoPlus:
            dfComplementaryVariable = dfTargetValue - dfVariable;
            break;
        case ProbDistOp::pdoMinus:
            if( bVarIsYVar )
                dfComplementaryVariable = dfTargetValue + dfVariable;
            else                
                dfComplementaryVariable = dfVariable - dfTargetValue;
            break;
        case ProbDistOp::pdoMult:
            dfComplementaryVariable = dfTargetValue / dfVariable;
            break;
        case ProbDistOp::pdoDiv:
            if(bVarIsYVar)
                dfComplementaryVariable = dfTargetValue * dfVariable;
            else
                dfComplementaryVariable = dfVariable / dfTargetValue ;
            break;
        default:
            break;
            
    }   //endswitch(Operation)
    
    return dfComplementaryVariable;
}


CProbabilityDensityDistribution& CProbabilityDensityDistribution::_GeneralOperatorsFunctionNumerical(CProbabilityDensityDistribution& Other, const ProbDistOp Operation)
{
    // generate plan : setting m_ConvolutionPlan;
    _SetConvolutionPlan(Other, Operation);
    
    // declare the result distri
    MapDFDFType TargetDistri;
    
    // this is only needed for covar operations
    CDigFloat dfMeanThis = Mean();
    
    for(auto itarget: m_ConvolutionPlan)
    {
        TargetDistri[itarget.first] = 0;
        
        // DEBUG
//         cout << "itarget.first: " << itarget.first.RawPrint(15) << endl;
        
        // iterate over all sub integration limits
        for(auto intlim: itarget.second)
        {
            // set the probabilities
            CDigFloat dfProbThis = AbsIntegral(intlim.first.first, intlim.first.second);
            CDigFloat dfProbOther = Other.AbsIntegral(intlim.second.first, intlim.second.second);
            
            // DEBUG
//             cout << "ProbThis(" << intlim.first.first.RawPrint(15) << " - " << intlim.first.second.RawPrint(15) << ") = " << dfProbThis.RawPrint(15) << endl;
//             cout << "ProbOther(" << intlim.second.first.RawPrint(15) << " - " << intlim.second.second.RawPrint(15) << ") = " << dfProbOther.RawPrint(15) << endl;
            
            // add to all other probs of this target value
            TargetDistri[itarget.first] += dfProbThis*dfProbOther;//*_GetIntegrationFactor4Convolution(Operation,intlim.first, dfMeanThis);
            
        }   // endfor(auto interval: itarget.second)        
        
    }   // endfor(auto itarget: m_ConvolutionPlan)
    
    _Init();
    
    
    if(m_ConvolutionPlan.size() > 0)
    {
        for(auto iel: TargetDistri)
            mDistribution[iel.first] = iel.second;
        
        
        // DEBUG
//         cout << "before Normalization: " << endl << Print(10);
        
        
        bNormalized=false;
        Normalize();
   
        // DEBUG     
//         cout << "after Normalization: " << endl << Print(10);
   
        assert(AbsIntegral() == 1);
        
    }   // endif(!(TargetRangeStart == TargetRangeEnd && TargetRangeStart == 0))
    
    return *this;
   
}

CDigFloat CProbabilityDensityDistribution::_GetIntegrationFactor4Convolution(const ProbDistOp Operation, const PairDFType pdtLimitsOne, const CDigFloat dfMeanOne)
{
    CDigFloat dfConvFactor = 1;
    
    switch(Operation)
    {
        case ProbDistOp::pdoPlus:
        case ProbDistOp::pdoMinus:
            // do nothin: 1 is correct
            break;
        case ProbDistOp::pdoMult:
            // is 1/|x|: calculate x as mean :
            // x = (xmin + xmax) / 2
            dfConvFactor = 1./abs(pdtLimitsOne.first + pdtLimitsOne.second)*2.;
            break;
        case ProbDistOp::pdoDiv:
            // is |x|: calculate x as mean :
            // x = (xmin + xmax) / 2
            dfConvFactor = abs(pdtLimitsOne.first + pdtLimitsOne.second)/2.;
            break;
        case ProbDistOp::pdoCov:
            // is 1./|x-<X>|: calculate x as mean :
            // x = (xmin + xmax) / 2
            dfConvFactor = 1./abs((pdtLimitsOne.first + pdtLimitsOne.second)/2. - dfMeanOne);
            break;
        default:
            assert(false);
    }
    
    // DEBUG
    cout << endl ;
    cout << "[" << pdtLimitsOne.first.Print(10) << " , " << pdtLimitsOne.second.Print(10) << " ]" << endl;
    cout << "dfConvFactor : "  << dfConvFactor.Print(10) << endl;
    
    return dfConvFactor;
}


void CProbabilityDensityDistribution::_SetConvolutionPlan(CProbabilityDensityDistribution& other, const ProbDistOp Operation)
{
    CDigFloat dfTargetStart, dfTargetEnd;
    _GetRangeFromDistributionOperation(other, Operation, dfTargetStart, dfTargetEnd);
    
    // get the target step
    CDigFloat dfTargetStep = (dfTargetEnd - dfTargetStart)/(CDigFloat(IntegrationSteps()));
    
    // DEBUG
//     cout << "target interval: [" << dfTargetStart.RawPrint(15) << " ; " << dfTargetEnd.RawPrint(15) << " ] " << endl;
    
    // clear the plan first
    m_ConvolutionPlan.clear();
    
    // now start building the plan 
    CDigFloat dfActTargetValue = dfTargetStart;
    while(dfActTargetValue <= dfTargetEnd)
    {
        // this is the vector defining the sub integration limits for a given targe value
        VectorSubIntLimitsType SubIntLimits;
        
        // defining the complete integration range for this and other distribution
        _SetSubIntegrationIntervals4TargetValue(Operation, dfActTargetValue, other, SubIntLimits);
        
        // now add the sub integration limits and its corresponding target value as pair
        m_ConvolutionPlan.push_back(pair<CDigFloat, VectorSubIntLimitsType>(dfActTargetValue, SubIntLimits));
        
        // increment actual target value 
        dfActTargetValue += dfTargetStep;
    }   // endwhile(dfActTargetValue <= dfTargetEnd)
    
}

void CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue(const ProbDistOp Operation, const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervals4TargetValue)
{
    switch(Operation)
    {
        case ProbDistOp::pdoPlus:
            _SetSubIntegrationIntervals4TargetValue4Addition(dfActTargetValue, other, SubIntervals4TargetValue);
            break;
        case ProbDistOp::pdoMinus:
            _SetSubIntegrationIntervals4TargetValue4Subtraction(dfActTargetValue, other, SubIntervals4TargetValue);
            break;
        case ProbDistOp::pdoMult:
            _SetSubIntegrationIntervals4TargetValue4Multiplication(dfActTargetValue, other, SubIntervals4TargetValue);
            break;
        case ProbDistOp::pdoDiv:
            _SetSubIntegrationIntervals4TargetValue4Division(dfActTargetValue, other, SubIntervals4TargetValue);
            break;
        case ProbDistOp::pdoCov:
            _SetSubIntegrationIntervals4TargetValue4Covariance(dfActTargetValue, other, SubIntervals4TargetValue);
            break;
        default:
            break;
    }   // endswitch(Operation)
}

void CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Addition(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& Other, VectorSubIntLimitsType& SubIntervals4TargetValue)
{
    // init result:
    SubIntervals4TargetValue.clear();
    
    // now set local variables
    VectorSubIntLimitsType SubTotalIntervalRanges4TargetValue;
    PairDFType pdtThisLimits, pdtOtherLimits;
    SubIntLimitsType SingleSubTotalIntervalRange;
    
    
    // DEBUG
//     cout << endl << "_SetSubIntegrationIntervals4TargetValue4Addition:" << endl;
//     cout << "dfActTargetValue = " << dfActTargetValue.RawPrint(15) << endl;
//     cout << "This-range: [" << firstVariable().RawPrint(15) << " ; " << lastVariable().RawPrint(15) << " ]" << endl;
//     cout << "Other-range: [" << Other.firstVariable().RawPrint(15) << " ; " << Other.lastVariable().RawPrint(15) << " ]" << endl;

    // invalidate result as initialization
    _InitIntegrationLimits(pdtThisLimits);
    _InitIntegrationLimits(pdtOtherLimits);
    
    // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue (tv)
    // dfXMin: lowest value of complete X-range 
    // dfXMax: 0 (if X-range covers also positive values) or largest value of complete X-range otherwise 
    // dfYMin: lowset value of complete Y-range
    // dfXMax: 0 (if Y-range covers also positive values) or largest value of complete Y-range otherwise 
    CDigFloat dfXMin = firstVariable(); 
    CDigFloat dfXMax = lastVariable();
    CDigFloat dfYMin = Other.firstVariable();
    CDigFloat dfYMax = Other.lastVariable();

    // calculate the extremal upper and lower limits defined by X/YMax/Min and dfActTargetValue 
    // dfYMinCalc: the lowest possible value for Y
    // dfYMaxCalc: the larges possible value for Y 
    CDigFloat dfYMinCalc = dfActTargetValue - dfXMax;
    CDigFloat dfYMaxCalc = dfActTargetValue - dfXMin;
            

    // DEBUG
//     cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
//     cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
//     cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
//     cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
//     cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
    
    // is calculation not possible ? --> leave 
    if(dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin )
        return;

    // now calculate the limits for Y
    pdtOtherLimits.first = max(dfYMinCalc, dfYMin);
    pdtOtherLimits.second = min(dfYMaxCalc, dfYMax);
    
    // and now for X
    pdtThisLimits.first = dfActTargetValue - pdtOtherLimits.second;
    pdtThisLimits.second = dfActTargetValue -  pdtOtherLimits.first;
    
    // DEBUG
//     cout << "X-Limits: [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << " ]" << endl;
//     cout << "Y-Limits: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl;

        
    // check for x (this distri) being within range
    if(pdtThisLimits.first > dfXMax || pdtThisLimits.second < dfXMin) 
        return;
        
    // add only if everything is ok
//     cout << "this: first < second : " << string((pdtThisLimits.first < pdtThisLimits.second) ? "less" : " ge") << endl;
//     cout << "other: first < second : " << string((pdtOtherLimits.first < pdtOtherLimits.second) ? "less" : "ge") << endl;
    // add only valid results
    if((pdtThisLimits.first >= pdtThisLimits.second) || (pdtOtherLimits.first >= pdtOtherLimits.second))
    {
        //DEBUG
//         cout << "zero interval or inverted order" << endl; 
        return;
    }

    SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
//     cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//     cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
    
        // now derive for all elements (total integration intervals) of SubTotalIntervalRanges4TargetValue (max. 2)
    // the SubIntegrationSteps sub integration intervals and add them to SubIntervals4TargetValue
         // do we have negative values
    bool bFirstNegative = (pdtThisLimits.first<0);
    bool bSecondNegative = (pdtOtherLimits.first<0);
    
    // calculate the factor needed for dividing the integration sub intervals from the total interval
    CDigFloat dfStepInc = (pdtThisLimits.second-pdtThisLimits.first)/SubIntegrationSteps();
    
    // in case dfFacInc might be equal one (due to errors in calculation of that factor ) --> leave out this calculations
    if(dfStepInc == 0)
        return;
    
    // DEBUG
//     cout << "totIv.first:  [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << " ]" << endl;
//     cout << "totIv.second: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl; 
//     cout << "dfStepInc: " << dfStepInc.RawPrint(15) << endl;
    
    // set the actual factor
    CDigFloat dfActAdd = 0;
    
    for(int isub = 0; isub < SubIntegrationSteps(); isub++)
    {
        // the sub integration limits
        PairDFType pdtThis, pdtOther;
        
        // set for first limit for this and second limit for other
        pdtThis.first = pdtThisLimits.first + dfActAdd;
        
        // increment dfActFac
        dfActAdd += dfStepInc;
        
        // set the second limit for this and the first limit for other
        pdtThis.second = pdtThisLimits.first + dfActAdd;
        
        pdtOther.first = dfActTargetValue - pdtThis.second;
        pdtOther.second = dfActTargetValue - pdtThis.first;
        
//         cout << "pdtThis (" << isub << "):  [" << pdtThis.first.RawPrint(15) << " ; " << pdtThis.second.RawPrint(15) << " ]" << endl;
//         cout << "pdtOther (" << isub << "): [" << pdtOther.first.RawPrint(15) << " ; " << pdtOther.second.RawPrint(15) << " ]" << endl; 
            
        // checks: take care for correct order it might happen that within the numerical error
        // equality might be achieved
        assert(pdtThis.first <= pdtThis.second);
        assert(pdtOther.first <= pdtOther.second);
        
        // now add the interval
        SubIntervals4TargetValue.push_back(SubIntLimitsType(pdtThis, pdtOther));
    
    }   // endfor(int isub = 0; isub < SubIntegrationSteps(); isub++)
    
}
void CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Subtraction(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& Other, VectorSubIntLimitsType& SubIntervals4TargetValue)
{
        // init result:
    SubIntervals4TargetValue.clear();
    
    // now set local variables
    VectorSubIntLimitsType SubTotalIntervalRanges4TargetValue;
    PairDFType pdtThisLimits, pdtOtherLimits;
    SubIntLimitsType SingleSubTotalIntervalRange;
    
    
    // DEBUG
//     cout << endl << "_SetSubIntegrationIntervals4TargetValue4Subtraction:" << endl;
//     cout << "dfActTargetValue = " << dfActTargetValue.RawPrint(15) << endl;
//     cout << "This-range: [" << firstVariable().RawPrint(15) << " ; " << lastVariable().RawPrint(15) << " ]" << endl;
//     cout << "Other-range: [" << Other.firstVariable().RawPrint(15) << " ; " << Other.lastVariable().RawPrint(15) << " ]" << endl;

    // invalidate result as initialization
    _InitIntegrationLimits(pdtThisLimits);
    _InitIntegrationLimits(pdtOtherLimits);
    
    // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue (tv)
    // dfXMin: lowest value of complete X-range 
    // dfXMax: 0 (if X-range covers also positive values) or largest value of complete X-range otherwise 
    // dfYMin: lowset value of complete Y-range
    // dfXMax: 0 (if Y-range covers also positive values) or largest value of complete Y-range otherwise 
    CDigFloat dfXMin = firstVariable(); 
    CDigFloat dfXMax = lastVariable();
    CDigFloat dfYMin = Other.firstVariable();
    CDigFloat dfYMax = Other.lastVariable();

    // calculate the extremal upper and lower limits defined by X/YMax/Min and dfActTargetValue 
    // dfYMinCalc: the lowest possible value for Y
    // dfYMaxCalc: the larges possible value for Y 
    CDigFloat dfYMinCalc = dfXMin - dfActTargetValue;
    CDigFloat dfYMaxCalc = dfXMax - dfActTargetValue;
            

    // DEBUG
//     cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
//     cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
//     cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
//     cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
//     cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
    
    // is calculation not possible ? --> leave 
    if(dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin )
        return;

    // now calculate the limits for Y
    pdtOtherLimits.first = max(dfYMinCalc, dfYMin);
    pdtOtherLimits.second = min(dfYMaxCalc, dfYMax);
    
    // and now for X
    pdtThisLimits.first = dfActTargetValue + pdtOtherLimits.first;
    pdtThisLimits.second = dfActTargetValue +  pdtOtherLimits.second;
    
    // DEBUG
//     cout << "X-Limits: [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << " ]" << endl;
//     cout << "Y-Limits: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl;

        
    // check for x (this distri) being within range
    if(pdtThisLimits.first > dfXMax || pdtThisLimits.second < dfXMin) 
        return;
        
    // add only if everything is ok
//     cout << "this: first < second : " << string((pdtThisLimits.first < pdtThisLimits.second) ? "less" : " ge") << endl;
//     cout << "other: first < second : " << string((pdtOtherLimits.first < pdtOtherLimits.second) ? "less" : "ge") << endl;

    // add only valid results
    if((pdtThisLimits.first >= pdtThisLimits.second) || (pdtOtherLimits.first >= pdtOtherLimits.second))
    {
        //DEBUG
//         cout << "zero interval or inverted order" << endl; 
        return;
    }

    SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
    // DEBUG
    //     cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//     cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
    
        // now derive for all elements (total integration intervals) of SubTotalIntervalRanges4TargetValue (max. 2)
    // the SubIntegrationSteps sub integration intervals and add them to SubIntervals4TargetValue
         // do we have negative values
    bool bFirstNegative = (pdtThisLimits.first<0);
    bool bSecondNegative = (pdtOtherLimits.first<0);
    
    // calculate the factor needed for dividing the integration sub intervals from the total interval
    CDigFloat dfStepInc = (pdtThisLimits.second-pdtThisLimits.first)/SubIntegrationSteps();
    
    // in case dfFacInc might be equal one (due to errors in calculation of that factor ) --> leave out this calculations
    if(dfStepInc == 0)
        return;
    
    // DEBUG
//     cout << "totIv.first:  [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << " ]" << endl;
//     cout << "totIv.second: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl; 
//     cout << "dfStepInc: " << dfStepInc.RawPrint(15) << endl;
    
    // set the actual factor
    CDigFloat dfActAdd = 0;
    
    for(int isub = 0; isub < SubIntegrationSteps(); isub++)
    {
        // the sub integration limits
        PairDFType pdtThis, pdtOther;
        
        // set for first limit for this and second limit for other
        pdtThis.first = pdtThisLimits.first + dfActAdd;
        
        // increment dfActFac
        dfActAdd += dfStepInc;
        
        // set the second limit for this and the first limit for other
        pdtThis.second = pdtThisLimits.first + dfActAdd;
        
        pdtOther.first = pdtThis.first - dfActTargetValue;
        pdtOther.second = pdtThis.second - dfActTargetValue;

        // DEBUG
//         cout << "pdtThis (" << isub << "):  [" << pdtThis.first.RawPrint(15) << " ; " << pdtThis.second.RawPrint(15) << " ]" << endl;
//         cout << "pdtOther (" << isub << "): [" << pdtOther.first.RawPrint(15) << " ; " << pdtOther.second.RawPrint(15) << " ]" << endl; 
            
        // checks: take care for correct order it might happen that within the numerical error
        // equality might be achieved
        assert(pdtThis.first <= pdtThis.second);
        assert(pdtOther.first <= pdtOther.second);
        
        // now add the interval
        SubIntervals4TargetValue.push_back(SubIntLimitsType(pdtThis, pdtOther));
    
    }   // endfor(int isub = 0; isub < SubIntegrationSteps(); isub++)

}
void CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Multiplication(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervals4TargetValue)
{
    // init result:
    SubIntervals4TargetValue.clear();
    
    // now set local variables
    VectorSubIntLimitsType SubTotalIntervalRanges4TargetValue;
    PairDFType pdtThisLimits, pdtOtherLimits;
    SubIntLimitsType SingleSubTotalIntervalRange;
    
    
    // DEBUG
//     cout << endl << "_SetSubIntegrationIntervals4TargetValue4Multiplication:" << endl;
//     cout << "dfActTargetValue = " << dfActTargetValue.RawPrint(15) << endl;
//     cout << "This-range: [" << firstVariable().RawPrint(15) << " ; " << lastVariable().RawPrint(15) << " ]" << endl;
//     cout << "Other-range: [" << other.firstVariable().RawPrint(15) << " ; " << other.lastVariable().RawPrint(15) << " ]" << endl;

    
    // for zero we cant do anything:
    // zero is a singularity and will be approximated by surrounding target values near zero
    if(dfActTargetValue == 0)
        return;
    
    if(dfActTargetValue < 0)
    {
        // 1. case: random variable this < 0
        //          random variable other > 0
        // --> total integration ranges for this case:
        // this  [ minThis               ; targetValue / maxOther ]
        // other [ targetValue / minThis ; maxOther ]
        if(firstVariable() < 0 && other.lastVariable() > 0)
        {
            _GetSubIntegrationTotalInterval4NegativeTargetValue4Multiplication(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
             
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;

            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
                
                // DEBUG
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
        }
        
        // add another limits
        if(other.firstVariable() < 0 && lastVariable() > 0)
        {
            _GetSubIntegrationTotalInterval4NegativeTargetValue4Multiplication(dfActTargetValue,other,*this,pdtOtherLimits,pdtThisLimits);
            
            // DEBUG
//                 cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//                 cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;

            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
                
                // DEBUG
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl; 
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl << endl;
            }
            
        }
        
        
    }   // endif(dfActTargetValue < 0)
    else
    {
        // check negative variable ranges case 
        if( firstVariable() < 0 && other.firstVariable() < 0)
        {
            
            _GetSubIntegrationTotalInterval4PositiveTargetValueNegativeVariables4Multiplication(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
            
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;
            
            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

            // DEBUG
//             cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//             cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
            
        }   // endif( firstVariable() < 0 && other.firstVariable() < 0)
        
        // check positive case 
        if( lastVariable() > 0 && other.lastVariable() > 0)
        {
            
            _GetSubIntegrationTotalInterval4PositiveTargetValuePositiveVariables4Multiplication(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
            
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;

            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

                // DEBUG
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
            
        }   // endif( lastVariable() > 0 && other.lastVariable() > 0)
        
    }   // endelse(dfActTargetValue < 0)
    
    // now derive for all elements (total integration intervals) of SubTotalIntervalRanges4TargetValue (max. 2)
    // the SubIntegrationSteps sub integration intervals and add them to SubIntervals4TargetValue
    for(auto totIv: SubTotalIntervalRanges4TargetValue)
    {
        // do we have negative values
        bool bFirstNegative = (totIv.first.first<0);
        bool bSecondNegative = (totIv.second.first<0);
        
        // calculate the factor needed for dividing the integration sub intervals from the total interval
        CDigFloat dfFacInc = pow( totIv.first.second / totIv.first.first , 1./SubIntegrationSteps());
        
        // in case dfFacInc might be equal one (due to errors in calculation of that factor ) --> leave out this calculations
        if(dfFacInc == 1)
            continue;
        
        // DEBUG
//         cout << "totIv.first:  [" << totIv.first.first.RawPrint(15) << " ; " << totIv.first.second.RawPrint(15) << " ]" << endl;
//         cout << "totIv.second: [" << totIv.second.first.RawPrint(15) << " ; " << totIv.second.second.RawPrint(15) << " ]" << endl; 
//         cout << "dfFacInc: " << dfFacInc.RawPrint(15) << endl;
        
        // set the actual factor
        CDigFloat dfActFacThis = 1.;
        
        for(int isub = 0; isub < SubIntegrationSteps(); isub++)
        {
            // the sub integration limits
            PairDFType pdtThis, pdtOther;
            
            // set for first limit for this and second limit for other
            pdtThis.first = dfActFacThis * totIv.first.first;
            pdtOther.first = dfActTargetValue / pdtThis.first;
//             pdtOther.second = totIv.second.second / dfActFacThis;
            
            // increment dfActFac
            dfActFacThis *= dfFacInc;
            
            // set the second limit for this and the first limit for other
            pdtThis.second = dfActFacThis * totIv.first.first;
            pdtOther.second = dfActTargetValue / pdtThis.second;
//             pdtOther.first = totIv.second.second / dfActFacThis;
            
            if( pdtOther.first > pdtOther.second )
                swap(pdtOther.first,pdtOther.second);
            
            // DEBUG
//             cout << "pdtThis (" << isub << "):  [" << pdtThis.first.RawPrint(15) << " ; " << pdtThis.second.RawPrint(15) << " ]" << endl;
//             cout << "pdtOther (" << isub << "): [" << pdtOther.first.RawPrint(15) << " ; " << pdtOther.second.RawPrint(15) << " ]" << endl; 
            
            // checks: take care for correct order 
            assert(pdtThis.first <= pdtThis.second);
            assert(pdtOther.first <= pdtOther.second);
            
            // now add the interval
            SubIntervals4TargetValue.push_back(SubIntLimitsType(pdtThis, pdtOther));
            
        }   // endfor(int isub = 0; isub < SubIntegrationSteps(); isub++)
         
    }   // endfor(auto totIntervals: SubTotalIntervalRanges4TargetValue)
    
}

void CProbabilityDensityDistribution::_GetSubIntegrationTotalInterval4NegativeTargetValue4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType &pdtOneLimits, PairDFType& pdtOtherLimits )
{
    // invalidate result as initialization
    pdtOneLimits.first =0;
    pdtOneLimits.second = 0;
    pdtOtherLimits.first =0;
    pdtOtherLimits.second = 0;
    
    // reset errors in case of re-usage of variables the error might accumulate and lead to strange effects
    pdtOneLimits.first.ResetError();
    pdtOneLimits.second.ResetError();
    pdtOtherLimits.first.ResetError();
    pdtOtherLimits.second.ResetError();
    
    
    // do this only if the One distri has a negative part AND the Other distri a positive part
    if(One.firstVariable() < 0 && Other.lastVariable() > 0)
    {   
        
        // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue
        // dfXMin: lowest value of complete X-range
        // dfXMax: 0 (if X-range covers positive values) or largest value of complete X-range otherwise
        // dfYMin: 0 (if Y-range covers negative values) or lowest value of complete Y-range otherwise
        // dfXMax: largest value of complete range of Y
        CDigFloat dfXMin = One.firstVariable();
        CDigFloat dfXMax = (One.lastVariable() >= 0) ? One.lastVariable().Error()*(-10.) /*close to zero (negative)*/ : One.lastVariable();
        CDigFloat dfYMin = (Other.firstVariable() <= 0) ? Other.firstVariable().Error()*10. /* close to zero (positive) */ : Other.firstVariable();
        CDigFloat dfYMax = Other.lastVariable();
        
        // calculate the extremal upper and lower limits defined by X/YMax/Min and dfActTargetValue 
        // dfXMinCalc: the lowest possible value for X 
        // dfXMaxCalc: the largest possible value for X 
        // dfYMinCalc: the lowest possible value for Y
        // dfYMaxCalc: the larges possible value for Y 
        // REMARK:
        // dfXMax and dfYMin could be zero --> this is considered for the calculation of dfXMinCalc and dfYMaxCalc
        CDigFloat dfXMinCalc = dfActTargetValue / dfYMin;
        CDigFloat dfXMaxCalc = dfActTargetValue / dfYMax;
        CDigFloat dfYMinCalc = dfActTargetValue / dfXMin;
        CDigFloat dfYMaxCalc = dfActTargetValue / dfXMax;
        
        
        // DEBUG
//         cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
//         cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
//         cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
//         cout << "dfXMinCalc: " << dfXMinCalc.RawPrint(15) << endl;
//         cout << "dfXMaxCalc: " << dfXMaxCalc.RawPrint(15) << endl;
//         cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
//         cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(dfXMinCalc > dfXMax || dfXMaxCalc < dfXMin || dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin) 
        {
            //DEBUG
//             cout << "Bereiche schließen sich aus --> keine Berechnung möglich" << endl;
            return;
            
        }   // endif(dfXMinCalc > dfXMax || dfXMaxCalc < dfXMin || dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin) 
        
        // get the min. for the other distri: determined by the division of max. abs. value of the one distri
        pdtOneLimits.first = max(dfXMin, dfXMinCalc);
        pdtOneLimits.second = min(dfXMax,dfXMaxCalc);   
        pdtOtherLimits.first = max(dfYMin, dfYMinCalc);
        pdtOtherLimits.second = min(dfYMax,dfYMaxCalc); 
        
        // DEBUG
//         cout << "X-Limits: [" << pdtOneLimits.first.RawPrint(15) << " ; " << pdtOneLimits.second.RawPrint(15) << " ]" << endl;
//         cout << "Y-Limits: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl;
        
        // now we got both limits of one and other: assure consistency
        // i.e. Xll = tv / YLL
        // AND
        // XUL = tv / YUL
        assert(pdtOneLimits.first == dfActTargetValue / pdtOtherLimits.first);
        assert(pdtOneLimits.second == CDigFloat(dfActTargetValue) / pdtOtherLimits.second);
        
    }   // endif(One.firstVariable() < 0 && Other.lastVariable() > 0)
    
    // check for non-zero / non-inverted ranges and the second (other) starts after the first (one)
    assert(pdtOneLimits.first <=  pdtOneLimits.second);
    assert(pdtOtherLimits.first <= pdtOtherLimits.second);
}

void CProbabilityDensityDistribution::_GetSubIntegrationTotalInterval4PositiveTargetValuePositiveVariables4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType &pdtOneLimits, PairDFType& pdtOtherLimits )
{
    
    // invalidate result as initialization
    pdtOneLimits.first =0;
    pdtOneLimits.second = 0;
    pdtOtherLimits.first =0;
    pdtOtherLimits.second = 0;
    
    // reset errors in case of re-usage of variables the error might accumulate and lead to strange effects
    pdtOneLimits.first.ResetError();
    pdtOneLimits.second.ResetError();
    pdtOtherLimits.first.ResetError();
    pdtOtherLimits.second.ResetError();
    
    // do this only if the One AND Other distri have a positive part 
    if(One.lastVariable() > 0 && Other.lastVariable() > 0)
    {
        
        // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue (tv)
        // dfXMin: 0 (if X-range covers also negative values) or lowest value of complete X-range otherwise 
        // dfXMax: largest value of complete X-range 
        // dfYMin: 0 (if Y-range covers also negative values) or lowest value of complete Y-range otherwise 
        // dfXMax: largest value of complete Y-range
        CDigFloat dfXMin = (One.firstVariable() < 0) ? 0 : One.firstVariable();
        CDigFloat dfXMax = One.lastVariable();
        CDigFloat dfYMin = (Other.firstVariable() < 0) ? 0 : Other.firstVariable();
        CDigFloat dfYMax = Other.lastVariable();
        
        // calculate the extremal upper and lower limits defined by X/YMax/Min and dfActTargetValue 
        // dfXMinCalc: the lowest possible value for X 
        // dfXMaxCalc: the largest possible value for X 
        // dfYMinCalc: the lowest possible value for Y
        // dfYMaxCalc: the larges possible value for Y 
        // REMARK:
        // dfXMin and dfYMin could be zero --> this is considered for the calculation of dfXMaxCalc and dfYMaxCalc
        CDigFloat dfXMinCalc = dfActTargetValue / dfYMax;
        CDigFloat dfXMaxCalc = (dfYMin < (dfActTargetValue/One.lastVariable())) ? One.lastVariable() : dfActTargetValue / dfYMin;
        CDigFloat dfYMinCalc = dfActTargetValue / dfXMax;
        CDigFloat dfYMaxCalc = (dfXMin < (dfActTargetValue/Other.lastVariable())) ? Other.lastVariable() : dfActTargetValue / dfXMin;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(dfXMinCalc > dfXMax || dfXMaxCalc < dfXMin || dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin) 
        {
            return;
            
        }   // endif(dfXMinCalc > dfXMax || dfXMaxCalc < dfXMin || dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin)         
        
        // get the min. for the other distri: determined by the division of max. abs. value of the one distri
        pdtOneLimits.first = max(dfXMin, dfXMinCalc);
        pdtOneLimits.second = min(dfXMax,dfXMaxCalc);        
        pdtOtherLimits.first = max(dfYMin, dfYMinCalc);
        pdtOtherLimits.second = min(dfYMax,dfYMaxCalc);
        
        // DEBUG
//         cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
//         cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
//         cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
//         cout << "dfXMinCalc: " << dfXMinCalc.RawPrint(15) << endl;
//         cout << "dfXMaxCalc: " << dfXMaxCalc.RawPrint(15) << endl;
//         cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
//         cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
//         cout << "X-Limits: [" << pdtOneLimits.first.RawPrint(15) << " ; " << pdtOneLimits.second.RawPrint(15) << " ]" << endl;
//         cout << "Y-Limits: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl;
        
        // now we got both limits of one and other: assure consistency
        // i.e. Xll = tv / YUL
        // AND
        // XUL = tv / YLL
        assert(pdtOneLimits.first == dfActTargetValue / pdtOtherLimits.second);
        assert(pdtOneLimits.second == CDigFloat(dfActTargetValue) / pdtOtherLimits.first);
        
    }
    
    // check for non-zero / non-inverted ranges and the second (other) starts after the first (one)
    assert(pdtOneLimits.first <=  pdtOneLimits.second);
    assert(pdtOtherLimits.first <= pdtOtherLimits.second);
}

void CProbabilityDensityDistribution:: _GetSubIntegrationTotalInterval4PositiveTargetValueNegativeVariables4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
{

    // invalidate result as initialization
    _InitIntegrationLimits(pdtOneLimits);
    _InitIntegrationLimits(pdtOtherLimits);
    
    // do this only if the One AND Other distri have a positive part 
    if(One.firstVariable() < 0 && Other.firstVariable() < 0)
    {
        // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue (tv)
        // dfXMin: lowest value of complete X-range 
        // dfXMax: 0 (if X-range covers also positive values) or largest value of complete X-range otherwise 
        // dfYMin: lowset value of complete Y-range
        // dfXMax: 0 (if Y-range covers also positive values) or largest value of complete Y-range otherwise 
        CDigFloat dfXMin = One.firstVariable(); 
        CDigFloat dfXMax = (One.lastVariable() > 0) ? 0 : One.lastVariable();
        CDigFloat dfYMin = Other.firstVariable();
        CDigFloat dfYMax = (Other.lastVariable() > 0) ? 0 : Other.lastVariable();
        
        // calculate the extremal upper and lower limits defined by X/YMax/Min and dfActTargetValue 
        // dfXMinCalc: the lowest possible value for X 
        // dfXMaxCalc: the largest possible value for X 
        // dfYMinCalc: the lowest possible value for Y
        // dfYMaxCalc: the larges possible value for Y 
        // REMARK:
        // dfXMax and dfYMax could be zero --> this is considered for the calculation of dfXMinCalc and dfYMinCalc
        CDigFloat dfXMinCalc = (dfYMax > (dfActTargetValue / One.firstVariable()) ) ? One.firstVariable() : dfActTargetValue / dfYMax;
        CDigFloat dfXMaxCalc = dfActTargetValue / dfYMin;
        CDigFloat dfYMinCalc = (dfXMax > (dfActTargetValue/ Other.firstVariable())) ? Other.firstVariable() : dfActTargetValue / dfXMax;
        CDigFloat dfYMaxCalc = dfActTargetValue / dfXMin;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(dfXMinCalc > dfXMax || dfXMaxCalc < dfXMin || dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin) 
        {
            return;
            
        }   // endif(dfXMinCalc > dfXMax || dfXMaxCalc < dfXMin || dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin)         
        
        // get the min. for the other distri: determined by the division of max. abs. value of the one distri
        pdtOneLimits.first = max(dfXMin, dfXMinCalc);
        pdtOneLimits.second = min(dfXMax,dfXMaxCalc);        
        pdtOtherLimits.first = max(dfYMin, dfYMinCalc);
        pdtOtherLimits.second = min(dfYMax,dfYMaxCalc);
        
        // DEBUG
//         cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
//         cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
//         cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
//         cout << "dfXMinCalc: " << dfXMinCalc.RawPrint(15) << endl;
//         cout << "dfXMaxCalc: " << dfXMaxCalc.RawPrint(15) << endl;
//         cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
//         cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
//         cout << "X-Limits: [" << pdtOneLimits.first.RawPrint(15) << " ; " << pdtOneLimits.second.RawPrint(15) << " ]" << endl;
//         cout << "Y-Limits: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl;

        // now we got both limits of one and other: assure consistency
        // i.e. Xll = tv / YUL
        // AND
        // XUL = tv / YLL
        assert(pdtOneLimits.first == dfActTargetValue / pdtOtherLimits.second);
        assert(pdtOneLimits.second == CDigFloat(dfActTargetValue) / pdtOtherLimits.first);
        
    }
    
    // check for non-zero / non-inverted ranges and the second (other) starts after the first (one)
    assert(pdtOneLimits.first <=  pdtOneLimits.second);
    assert(pdtOtherLimits.first <= pdtOtherLimits.second);
}


void CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Division(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervals4TargetValue)
{
       // init result:
    SubIntervals4TargetValue.clear();
    
    // now set local variables
    VectorSubIntLimitsType SubTotalIntervalRanges4TargetValue;
    PairDFType pdtThisLimits, pdtOtherLimits;
    SubIntLimitsType SingleSubTotalIntervalRange;
    
    
    // DEBUG
//     cout << endl << "_SetSubIntegrationIntervals4TargetValue4Division:" << endl;
//     cout << "dfActTargetValue = " << dfActTargetValue.RawPrint(15) << endl;
//     cout << "This-range: [" << firstVariable().RawPrint(15) << " ; " << lastVariable().RawPrint(15) << " ]" << endl;
//     cout << "Other-range: [" << other.firstVariable().RawPrint(15) << " ; " << other.lastVariable().RawPrint(15) << " ]" << endl;

    
    // for zero we cant do anything:
    // zero is a singularity and will be approximated by surrounding target values near zero
    if(dfActTargetValue == 0)
        return;
    
    if(dfActTargetValue < 0)
    {
        // 1. case: random variable this < 0
        //          random variable other > 0
        // --> total integration ranges for this case:
        // this  [ minThis               ; targetValue / maxOther ]
        // other [ targetValue / minThis ; maxOther ]
        if(firstVariable() < 0 && other.lastVariable() > 0)
        {
            // DEBUG
//             cout << "calling _GetSubIntegrationTotalInterval4NegativeTargetValueNegPosVariables4Division " << endl;
   
            _GetSubIntegrationTotalInterval4NegativeTargetValueNegPosVariables4Division(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
            
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;
 
            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
            
                // DEBUG
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
        }
        
        // add another limits
        if(lastVariable() > 0 && other.firstVariable() < 0)
        {
            // DEBUG
//             cout << "calling _GetSubIntegrationTotalInterval4NegativeTargetValuePosNegVariables4Division " << endl;
   
            _GetSubIntegrationTotalInterval4NegativeTargetValuePosNegVariables4Division(dfActTargetValue,*this,other,pdtThisLimits,pdtOtherLimits);
            
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;

            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

                // DEBUG                
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl; 
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl << endl;
            }
            
        }
        
        
    }   // endif(dfActTargetValue < 0)
    else
    {
        // check negative variable ranges case 
        if( firstVariable() < 0 && other.firstVariable() < 0)
        {
            // DEBUG
//             cout << "calling _GetSubIntegrationTotalInterval4PositiveTargetValueNegativeVariables4Division " << endl;
            
            _GetSubIntegrationTotalInterval4PositiveTargetValueNegativeVariables4Division(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
   
            // DEBUG             
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;
            
            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

                // DEBUG                
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
            
        }   // endif( firstVariable() < 0 && other.firstVariable() < 0)
        
        // check positive case 
        if( lastVariable() > 0 && other.lastVariable() > 0)
        {
            
            // DEBUG
//             cout << "calling _GetSubIntegrationTotalInterval4PositiveTargetValuePositiveVariables4Division " << endl;
            _GetSubIntegrationTotalInterval4PositiveTargetValuePositiveVariables4Division(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
            
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;
            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

                // DEBUG                
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
            
        }   // endif( lastVariable() > 0 && other.lastVariable() > 0)
        
    }   // endelse(dfActTargetValue < 0)
    
    // now derive for all elements (total integration intervals) of SubTotalIntervalRanges4TargetValue (max. 2)
    // the SubIntegrationSteps sub integration intervals and add them to SubIntervals4TargetValue
    for(auto totIv: SubTotalIntervalRanges4TargetValue)
    {
        // do we have negative values
        bool bFirstNegative = (totIv.first.first<0);
        bool bSecondNegative = (totIv.second.first<0);
        
        // calculate the factor needed for dividing the integration sub intervals from the total interval
        CDigFloat dfFacInc = pow( totIv.first.second / totIv.first.first , 1./SubIntegrationSteps());
        
        // in case dfFacInc might be equal one (due to errors in calculation of that factor ) --> leave out this calculations
        if(dfFacInc == 1)
            continue;
        
        // DEBUG
//         cout << "totIv.first:  [" << totIv.first.first.RawPrint(15) << " ; " << totIv.first.second.RawPrint(15) << " ]" << endl;
//         cout << "totIv.second: [" << totIv.second.first.RawPrint(15) << " ; " << totIv.second.second.RawPrint(15) << " ]" << endl; 
//         cout << "dfFacInc: " << dfFacInc.RawPrint(15) << endl;
        
        // set the actual factor
        CDigFloat dfActFacThis = 1.;
        
        for(int isub = 0; isub < SubIntegrationSteps(); isub++)
        {
            // the sub integration limits
            PairDFType pdtThis, pdtOther;
            
            // set for first limit for this and second limit for other
            pdtThis.first = dfActFacThis * totIv.first.first;
            pdtOther.first = pdtThis.first / dfActTargetValue;
//             pdtOther.second = totIv.second.second / dfActFacThis;
            
            // increment dfActFac
            dfActFacThis *= dfFacInc;
            
            // set the second limit for this and the first limit for other
            pdtThis.second = dfActFacThis * totIv.first.first;
            pdtOther.second = pdtThis.second / dfActTargetValue;
//             pdtOther.first = totIv.second.second / dfActFacThis;
            
            if( pdtOther.first > pdtOther.second )
                swap(pdtOther.first,pdtOther.second);
            
            // DEBUG
//             cout << "pdtThis (" << isub << "):  [" << pdtThis.first.RawPrint(15) << " ; " << pdtThis.second.RawPrint(15) << " ]" << endl;
//             cout << "pdtOther (" << isub << "): [" << pdtOther.first.RawPrint(15) << " ; " << pdtOther.second.RawPrint(15) << " ]" << endl; 
            
            // checks: take care for correct order / sensfull integration limits ... might be close to identical in some cases
            if(pdtThis.first >= pdtThis.second)
                continue;
            
            if(pdtOther.first >= pdtOther.second)
                continue;
            
            // checks: consistency
//             assert(pdtThis.first * pdtOther.second == dfActTargetValue);
//             assert(pdtThis.second * pdtOther.first == dfActTargetValue);
            
            // now add the interval
            SubIntervals4TargetValue.push_back(SubIntLimitsType(pdtThis, pdtOther));
            
        }   // endfor(int isub = 0; isub < SubIntegrationSteps(); isub++)
         
    }   // endfor(auto totIntervals: SubTotalIntervalRanges4TargetValue)
 
}

void CProbabilityDensityDistribution::_GetSubIntegrationTotalInterval4PositiveTargetValueNegativeVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
{
        // invalidate result as initialization
    _InitIntegrationLimits(pdtOneLimits);
    _InitIntegrationLimits(pdtOtherLimits);
    
    // do this only if the One AND Other distri have a positive part 
    if(One.firstVariable() < 0 && Other.firstVariable() < 0)
    {
        // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue (tv)
        // dfXMin: lowest value of complete X-range 
        // dfXMax: 0 (if X-range covers also positive values) or largest value of complete X-range otherwise 
        // dfYMin: lowset value of complete Y-range
        // dfXMax: 0 (if Y-range covers also positive values) or largest value of complete Y-range otherwise 
        CDigFloat dfXMin = One.firstVariable(); 
        CDigFloat dfXMax = (One.lastVariable() >= 0) ? CDigFloat(100).Error()*(-10.) /*close to zero (negative)*/ : One.lastVariable();
        CDigFloat dfYMin = Other.firstVariable();
        CDigFloat dfYMax = (Other.lastVariable() >= 0) ? CDigFloat(100).Error()*(-10.) /*close to zero (negative)*/: Other.lastVariable();
        
        // calculate the extremal upper and lower limits defined by X/YMax/Min and dfActTargetValue 
        // dfYMinCalc: the lowest possible value for Y
        // dfYMaxCalc: the larges possible value for Y 
        CDigFloat dfYMinCalc = dfXMin / dfActTargetValue ;
        CDigFloat dfYMaxCalc = dfXMax / dfActTargetValue ;
               
        // DEBUG
//         cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
//         cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
//         cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
//         cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
//         cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin) 
        {
            _InitIntegrationLimits(pdtOneLimits);
            return;
            
        }   // endif(dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin)
   
        
        // now calculate the limits for Y
        pdtOtherLimits.first = max(dfYMinCalc, dfYMin);
        pdtOtherLimits.second = min(dfYMaxCalc, dfYMax);
        
        // and now for X
        pdtOneLimits.first = dfActTargetValue * pdtOtherLimits.first;
        pdtOneLimits.second = dfActTargetValue * pdtOtherLimits.second;
        
        // DEBUG
//         cout << "X-Limits: [" << pdtOneLimits.first.RawPrint(15) << " ; " << pdtOneLimits.second.RawPrint(15) << " ]" << endl;
//         cout << "Y-Limits: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl;

        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(pdtOneLimits.first > dfXMax || pdtOneLimits.second < dfXMin || pdtOtherLimits.first > dfYMax || pdtOtherLimits.second < dfYMin) 
        {
            _InitIntegrationLimits(pdtOneLimits);
            return;
            
        }   // endif(dfXMinCalc > dfXMax || dfXMaxCalc < dfXMin || dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin)  
         
    }
    
    // check for non-zero / non-inverted ranges and the second (other) starts after the first (one)
    assert(pdtOneLimits.first <=  pdtOneLimits.second);
    assert(pdtOtherLimits.first <= pdtOtherLimits.second);

}
void CProbabilityDensityDistribution::_GetSubIntegrationTotalInterval4PositiveTargetValuePositiveVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
{   
    // invalidate result as initialization
    _InitIntegrationLimits(pdtOneLimits);
    _InitIntegrationLimits(pdtOtherLimits);
    
    // do this only if the One AND Other distri have a positive part 
    if(One.lastVariable() > 0 && Other.lastVariable() > 0)
    {
        // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue (tv)
        // dfXMin: lowest value of complete X-range 
        // dfXMax: 0 (if X-range covers also positive values) or largest value of complete X-range otherwise 
        // dfYMin: lowset value of complete Y-range
        // dfXMax: 0 (if Y-range covers also positive values) or largest value of complete Y-range otherwise 
        CDigFloat dfXMin = (One.firstVariable() <= 0) ? CDigFloat(100).Error()*10. /*close to zero (positive)*/ : One.firstVariable();
        CDigFloat dfXMax = One.lastVariable();
        CDigFloat dfYMin = (Other.firstVariable() <= 0) ? CDigFloat(100).Error()*10. /*close to zero (positive)*/ : Other.firstVariable();
        CDigFloat dfYMax = Other.lastVariable();
        
        // calculate the extremal upper and lower limits defined by X/YMax/Min and dfActTargetValue 
        // dfYMinCalc: the lowest possible value for Y
        // dfYMaxCalc: the larges possible value for Y 
        CDigFloat dfYMinCalc = dfXMin / dfActTargetValue ;
        CDigFloat dfYMaxCalc = dfXMax / dfActTargetValue ;
        
        
        // DEBUG
//         cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
//         cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
//         cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
//         cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
//         cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin) 
        {
            _InitIntegrationLimits(pdtOneLimits);
            return;
            
        }   // endif(dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin)
        
        // now calculate the limits for Y
        pdtOtherLimits.first = max(dfYMinCalc, dfYMin);
        pdtOtherLimits.second = min(dfYMaxCalc, dfYMax);
        
        // and now for X
        pdtOneLimits.first = dfActTargetValue * pdtOtherLimits.first;
        pdtOneLimits.second = dfActTargetValue * pdtOtherLimits.second;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(pdtOneLimits.first > dfXMax || pdtOneLimits.second < dfXMin || pdtOtherLimits.first > dfYMax || pdtOtherLimits.second < dfYMin) 
        {
            _InitIntegrationLimits(pdtOneLimits);
            return;
            
        }   // endif(dfXMinCalc > dfXMax || dfXMaxCalc < dfXMin || dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin)         
       
        // DEBUG
//         cout << "X-Limits: [" << pdtOneLimits.first.RawPrint(15) << " ; " << pdtOneLimits.second.RawPrint(15) << " ]" << endl;
//         cout << "Y-Limits: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl;

         
    }
    
    // check for non-zero / non-inverted ranges and the second (other) starts after the first (one)
    assert(pdtOneLimits.first <=  pdtOneLimits.second);
    assert(pdtOtherLimits.first <= pdtOtherLimits.second);

}

void CProbabilityDensityDistribution::_GetSubIntegrationTotalInterval4NegativeTargetValuePosNegVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
{
        // invalidate result as initialization
    _InitIntegrationLimits(pdtOneLimits);
    _InitIntegrationLimits(pdtOtherLimits);
    
    // do this only if the One AND Other distri have a positive part 
    if(One.lastVariable() > 0 && Other.firstVariable() < 0 )
    {
        // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue (tv)
        // dfXMin: lowest value of complete X-range 
        // dfXMax: 0 (if X-range covers also positive values) or largest value of complete X-range otherwise 
        // dfYMin: lowset value of complete Y-range
        // dfXMax: 0 (if Y-range covers also positive values) or largest value of complete Y-range otherwise 
        CDigFloat dfXMin = (One.firstVariable() <= 0) ? CDigFloat(100).Error()*10. /*close to zero (positive)*/ : One.firstVariable();
        CDigFloat dfXMax = One.lastVariable();
        CDigFloat dfYMin = Other.firstVariable();
        CDigFloat dfYMax = (Other.lastVariable() >= 0) ? CDigFloat(100).Error()*(-10.) /*close to zero (negative)*/ : Other.lastVariable();
        
        // calculate the extremal upper and lower limits defined by X/YMax/Min and dfActTargetValue 
        // dfYMinCalc: the lowest possible value for Y
        // dfYMaxCalc: the larges possible value for Y 
        CDigFloat dfYMinCalc = dfXMax / dfActTargetValue ;
        CDigFloat dfYMaxCalc = dfXMin / dfActTargetValue ;
        
        // DEBUG
//         cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
//         cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
//         cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
//         cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
//         cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin) 
        {
            _InitIntegrationLimits(pdtOneLimits);
            return;
            
        }   // endif(dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin)
        
        // now calculate the limits for Y
        pdtOtherLimits.first = max(dfYMinCalc, dfYMin);
        pdtOtherLimits.second = min(dfYMaxCalc, dfYMax);
        
        // and now for X
        pdtOneLimits.first = dfActTargetValue * pdtOtherLimits.second;
        pdtOneLimits.second = dfActTargetValue * pdtOtherLimits.first;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(pdtOneLimits.first > dfXMax || pdtOneLimits.second < dfXMin || pdtOtherLimits.first > dfYMax || pdtOtherLimits.second < dfYMin) 
        {
            _InitIntegrationLimits(pdtOneLimits);
            return;
            
        }   // endif(dfXMinCalc > dfXMax || dfXMaxCalc < dfXMin || dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin)         
       
        // DEBUG
//         cout << "X-Limits: [" << pdtOneLimits.first.RawPrint(15) << " ; " << pdtOneLimits.second.RawPrint(15) << " ]" << endl;
//         cout << "Y-Limits: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl;

         
    }
    
    // check for non-zero / non-inverted ranges and the second (other) starts after the first (one)
    assert(pdtOneLimits.first <=  pdtOneLimits.second);
    assert(pdtOtherLimits.first <= pdtOtherLimits.second);

}

void CProbabilityDensityDistribution::_GetSubIntegrationTotalInterval4NegativeTargetValueNegPosVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
{
           // invalidate result as initialization
    _InitIntegrationLimits(pdtOneLimits);
    _InitIntegrationLimits(pdtOtherLimits);
    
    // do this only if the One has negative AND Other distri has a positive part 
    if(One.firstVariable() < 0 && Other.lastVariable() > 0 )
    {
        // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue (tv)
        // dfXMin: lowest value of complete X-range 
        // dfXMax: 0 (if X-range covers also positive values) or largest value of complete X-range otherwise 
        // dfYMin: lowset value of complete Y-range
        // dfXMax: 0 (if Y-range covers also positive values) or largest value of complete Y-range otherwise 
        CDigFloat dfXMin = One.firstVariable();
        CDigFloat dfXMax = (One.lastVariable() >= 0) ? CDigFloat(100).Error()*(-10.) /*close to zero (negative)*/ : One.lastVariable();
        CDigFloat dfYMin = (Other.firstVariable() <= 0) ? CDigFloat(100).Error()*10. /*close to zero (positive)*/ : Other.firstVariable();
        CDigFloat dfYMax = Other.lastVariable();
        
        // calculate the extremal upper and lower limits defined by X/YMax/Min and dfActTargetValue 
        // dfYMinCalc: the lowest possible value for Y
        // dfYMaxCalc: the larges possible value for Y 
        CDigFloat dfYMinCalc = dfXMax / dfActTargetValue ;
        CDigFloat dfYMaxCalc = dfXMin / dfActTargetValue ;
        
        // now calculate the limits for Y
        pdtOtherLimits.first = max(dfYMinCalc, dfYMin);
        pdtOtherLimits.second = min(dfYMaxCalc, dfYMax);
        
        // and now for X
        pdtOneLimits.first = dfActTargetValue * pdtOtherLimits.second;
        pdtOneLimits.second = dfActTargetValue * pdtOtherLimits.first;
        
        // DEBUG
//         cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
//         cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
//         cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
//         cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
//         cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(pdtOneLimits.first > dfXMax || pdtOneLimits.second < dfXMin || pdtOtherLimits.first > dfYMax || pdtOtherLimits.second < dfYMin) 
        {
            _InitIntegrationLimits(pdtOneLimits);
            return;
            
        }   // endif(dfXMinCalc > dfXMax || dfXMaxCalc < dfXMin || dfYMinCalc > dfYMax || dfYMaxCalc < dfYMin)         
       
        
        // DEBUG
//         cout << "X-Limits: [" << pdtOneLimits.first.RawPrint(15) << " ; " << pdtOneLimits.second.RawPrint(15) << " ]" << endl;
//         cout << "Y-Limits: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl;

         
    }
    
    // check for non-zero / non-inverted ranges and the second (other) starts after the first (one)
    assert(pdtOneLimits.first <=  pdtOneLimits.second);
    assert(pdtOtherLimits.first <= pdtOtherLimits.second);

}

void CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue4Covariance(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& Other, VectorSubIntLimitsType& SubIntervals4TargetValue)
{
    // init result:
    SubIntervals4TargetValue.clear();
    
    // now set local variables
    VectorSubIntLimitsType SubTotalIntervalRanges4TargetValue;
    PairDFType pdtThisLimits, pdtOtherLimits;
    SubIntLimitsType SingleSubTotalIntervalRange;
    
    // need the means very often
    CDigFloat dfOneMean = Mean();
    CDigFloat dfOtherMean = Other.Mean();
    
    
    // DEBUG
    cout << endl << "_SetSubIntegrationIntervals4TargetValue4Covariance:" << endl;
    cout << "dfActTargetValue = " << dfActTargetValue.RawPrint(15) << endl;
    cout << "This-range: [" << firstVariable().RawPrint(15) << " ; " << lastVariable().RawPrint(15) << " ]" << endl;
    cout << "Other-range: [" << Other.firstVariable().RawPrint(15) << " ; " << Other.lastVariable().RawPrint(15) << " ]" << endl;

    
    // for zero we cant do anything:
    // zero is a singularity and will be approximated by surrounding target values near zero
    if(dfActTargetValue == 0)
        return;
    
    if(dfActTargetValue < 0)
    {
        // 1. case: random variable this < 0
        //          random variable other > 0
        // --> total integration ranges for this case:
        // this  [ minThis               ; targetValue / maxOther ]
        // other [ targetValue / minThis ; maxOther ]
        if( (firstVariable() - dfOneMean) < 0 && (Other.lastVariable()-dfOtherMean) > 0)
        {
            _GetSubIntegrationTotalInterval4NegativeTargetValue4Covariance(dfActTargetValue,*this,Other,pdtThisLimits, pdtOtherLimits);
             
            // DEBUG
            cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
            cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;

            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
                
                // DEBUG
                cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
                cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
        }
        
        // add another limits
        if( (Other.firstVariable()-dfOtherMean) < 0 && (lastVariable()-dfOneMean) > 0)
        {
            _GetSubIntegrationTotalInterval4NegativeTargetValue4Covariance(dfActTargetValue,Other,*this,pdtOtherLimits,pdtThisLimits);
            
            // DEBUG
                cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
                cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;

            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
                
                // DEBUG
                cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl; 
                cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl << endl;
            }
            
        }
        
        
    }   // endif(dfActTargetValue < 0)
    else
    {
        // check negative variable ranges case 
        if( (firstVariable()-dfOneMean) < 0 && (Other.firstVariable()-dfOtherMean) < 0)
        {
            
            _GetSubIntegrationTotalInterval4PositiveTargetValueNegativeVariables4Covariance(dfActTargetValue,*this,Other,pdtThisLimits, pdtOtherLimits);
            
            // DEBUG
            cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
            cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;
            
            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

            // DEBUG
            cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
            cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
            
        }   // endif( firstVariable() < 0 && other.firstVariable() < 0)
        
        // check positive case 
        if( (lastVariable()-dfOneMean) > 0 && (Other.lastVariable()-dfOtherMean) > 0)
        {
            
            _GetSubIntegrationTotalInterval4PositiveTargetValuePositiveVariables4Covariance(dfActTargetValue,*this,Other,pdtThisLimits, pdtOtherLimits);
            
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;

            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervalRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

                // DEBUG
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
            
        }   // endif( lastVariable() > 0 && other.lastVariable() > 0)
        
    }   // endelse(dfActTargetValue < 0)
    
    // now derive for all elements (total integration intervals) of SubTotalIntervalRanges4TargetValue (max. 2)
    // the SubIntegrationSteps sub integration intervals and add them to SubIntervals4TargetValue
    for(auto totIv: SubTotalIntervalRanges4TargetValue)
    {   
        // calculate the factor needed for dividing the integration sub intervals from the total interval
        CDigFloat dfFacInc = pow( totIv.first.second / totIv.first.first , 1./SubIntegrationSteps());
        
        // in case dfFacInc might be equal one (due to errors in calculation of that factor ) --> leave out this calculations
        if(dfFacInc == 1)
            continue;
        
        // DEBUG
        cout << "totIv.first:  [" << totIv.first.first.RawPrint(15) << " ; " << totIv.first.second.RawPrint(15) << " ]" << endl;
        cout << "totIv.second: [" << totIv.second.first.RawPrint(15) << " ; " << totIv.second.second.RawPrint(15) << " ]" << endl; 
        cout << "dfFacInc: " << dfFacInc.RawPrint(15) << endl;
        
        // set the actual factor
        CDigFloat dfActFacThis = 1.;
        
        for(int isub = 0; isub < SubIntegrationSteps(); isub++)
        {
            // the sub integration limits
            PairDFType pdtThis, pdtOther;
            
            // set for first limit for this and second limit for other
            pdtThis.first = dfActFacThis * totIv.first.first;
            pdtOther.first = dfActTargetValue / ( pdtThis.first - dfOneMean) + dfOtherMean;
//             pdtOther.second = totIv.second.second / dfActFacThis;
            
            // increment dfActFac
            dfActFacThis *= dfFacInc;
            
            // set the second limit for this and the first limit for other
            pdtThis.second = dfActFacThis * totIv.first.first;
            pdtOther.second = dfActTargetValue / (pdtThis.second - dfOneMean) + dfOtherMean;;
//             pdtOther.first = totIv.second.second / dfActFacThis;
            
            if( pdtOther.first > pdtOther.second )
                swap(pdtOther.first,pdtOther.second);
            
            // DEBUG
            cout << "pdtThis (" << isub << "):  [" << pdtThis.first.RawPrint(15) << " ; " << pdtThis.second.RawPrint(15) << " ]" << endl;
            cout << "pdtOther (" << isub << "): [" << pdtOther.first.RawPrint(15) << " ; " << pdtOther.second.RawPrint(15) << " ]" << endl; 
            
            // checks: take care for correct order 
            assert(pdtThis.first <= pdtThis.second);
            assert(pdtOther.first <= pdtOther.second);
            
            // now add the interval
            SubIntervals4TargetValue.push_back(SubIntLimitsType(pdtThis, pdtOther));
            
        }   // endfor(int isub = 0; isub < SubIntegrationSteps(); isub++)
         
    }   // endfor(auto totIntervals: SubTotalIntervalRanges4TargetValue)
    
}

void CProbabilityDensityDistribution::_GetSubIntegrationTotalInterval4NegativeTargetValue4Covariance(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
{
    // invalidate result as initialization
    _InitIntegrationLimits(pdtOneLimits);
    _InitIntegrationLimits(pdtOtherLimits);
    
    // need the means very often
    CDigFloat dfOneMean = One.Mean();
    CDigFloat dfOtherMean = Other.Mean();    
    
    // do this only if the One distri has a negative part AND the Other distri a positive part
    if( (One.firstVariable()-dfOneMean) < 0 && (Other.lastVariable()-dfOtherMean) > 0)
    {   
        
        // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue
        // dfXMin: lowest value of complete X-range
        // dfXMax: 0 (if X-range covers positive values) or largest value of complete X-range otherwise
        // dfYMin: 0 (if Y-range covers negative values) or lowest value of complete Y-range otherwise
        // dfXMax: largest value of complete range of Y
        CDigFloat dfXMin = One.firstVariable();
        CDigFloat dfXMax = ( (One.lastVariable() - dfOneMean) >= 0) ? dfOneMean - CDigFloat(1).Error()*10 /*close to mean (less)*/ : One.lastVariable();
        CDigFloat dfYMin = ( (Other.firstVariable() - dfOtherMean) <= 0) ? dfOtherMean + Other.firstVariable().Error()*10. /* close to mean (greater) */ : Other.firstVariable();
        CDigFloat dfYMax = Other.lastVariable();
        
        // calculate the theoretical y limits from x
        CDigFloat dfYMinCalc = dfActTargetValue / (dfXMin-dfOneMean) + dfOtherMean;
        CDigFloat dfYMaxCalc = dfActTargetValue / (dfXMax-dfOneMean) + dfOtherMean;
        
        // now calculate the real limits for Y
        pdtOtherLimits.first = max(dfYMinCalc, dfYMin);
        pdtOtherLimits.second = min(dfYMaxCalc, dfYMax);
        
        // and now for X
        pdtOneLimits.first = dfActTargetValue / (pdtOtherLimits.first - dfOtherMean) + dfOneMean;
        pdtOneLimits.second = dfActTargetValue / (pdtOtherLimits.second- dfOtherMean) + dfOneMean;
        
        // DEBUG
        cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
        cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
        cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
        cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
        cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(pdtOneLimits.first > dfXMax || pdtOneLimits.second < dfXMin || pdtOtherLimits.first > dfYMax || pdtOtherLimits.second < dfYMin) 
        {
            _InitIntegrationLimits(pdtOneLimits);
            return;
            
        }   // endif(pdtOneLimits.first > dfXMax || pdtOneLimits.second < dfXMin || pdtOtherLimits.first > dfYMax || pdtOtherLimits.second < dfYMin)        
       
        
        // DEBUG
        cout << "X-Limits: [" << pdtOneLimits.first.RawPrint(15) << " ; " << pdtOneLimits.second.RawPrint(15) << " ]" << endl;
        cout << "Y-Limits: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl;
        
        
    }   // endif( (One.firstVariable()-One.Mean()) < 0 && (Other.lastVariable()-Other.Mean()) > 0)
    
    // check for non-zero / non-inverted ranges and the second (other) starts after the first (one)
    assert(pdtOneLimits.first <=  pdtOneLimits.second);
    assert(pdtOtherLimits.first <= pdtOtherLimits.second);
}

void CProbabilityDensityDistribution::_GetSubIntegrationTotalInterval4PositiveTargetValueNegativeVariables4Covariance(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
{

    // invalidate result as initialization
    _InitIntegrationLimits(pdtOneLimits);
    _InitIntegrationLimits(pdtOtherLimits);
    
    // need the means very often
    CDigFloat dfOneMean = One.Mean();
    CDigFloat dfOtherMean = Other.Mean();
    
    // do this only if the One AND Other distri have a positive part 
    if( (One.firstVariable()-dfOneMean) < 0 && (Other.firstVariable()-dfOtherMean) < 0)
    {
        // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue (tv)
        // dfXMin: lowest value of complete X-range 
        // dfXMax: 0 (if X-range covers also positive values) or largest value of complete X-range otherwise 
        // dfYMin: lowset value of complete Y-range
        // dfXMax: 0 (if Y-range covers also positive values) or largest value of complete Y-range otherwise 
        CDigFloat dfXMin = One.firstVariable(); 
        CDigFloat dfXMax = ( (One.lastVariable()-dfOneMean) > 0) ? dfOneMean - CDigFloat(1).Error()*10 /*close to mean (less)*/  : One.lastVariable();
        CDigFloat dfYMin = Other.firstVariable();
        CDigFloat dfYMax = ( (Other.lastVariable()-dfOtherMean) > 0) ? dfOtherMean - CDigFloat(1).Error()*10 /*close to mean (less)*/  : Other.lastVariable();
        
        // take care for correct order:
        if(dfXMin > dfXMax)
            swap(dfXMin, dfXMax);
        if(dfYMin > dfYMax)
            swap(dfYMin, dfYMax);
        
        // calculate the theoretical y limits from x
        CDigFloat dfYMinCalc = dfActTargetValue / (dfXMin-dfOneMean) + dfOtherMean;
        CDigFloat dfYMaxCalc = dfActTargetValue / (dfXMax-dfOneMean) + dfOtherMean;
        
        // take care for correct order:
        if(dfYMinCalc > dfYMaxCalc)
            swap(dfYMinCalc, dfYMaxCalc);
        
        // now calculate the real limits for Y
        pdtOtherLimits.first = max(dfYMinCalc, dfYMin);
        pdtOtherLimits.second = min(dfYMaxCalc, dfYMax);
        
        // and now for X
        pdtOneLimits.first = dfActTargetValue / (pdtOtherLimits.second - dfOtherMean) + dfOneMean;
        pdtOneLimits.second = dfActTargetValue / (pdtOtherLimits.first- dfOtherMean) + dfOneMean;
        
        // DEBUG
        cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
        cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
        cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
        cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
        cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(pdtOneLimits.first > dfXMax || pdtOneLimits.second < dfXMin || pdtOtherLimits.first > dfYMax || pdtOtherLimits.second < dfYMin) 
        {
            _InitIntegrationLimits(pdtOneLimits);
            return;
            
        }   // endif(pdtOneLimits.first > dfXMax || pdtOneLimits.second < dfXMin || pdtOtherLimits.first > dfYMax || pdtOtherLimits.second < dfYMin)        
       
    }
    
    // check for non-zero / non-inverted ranges and the second (other) starts after the first (one)
    assert(pdtOneLimits.first <=  pdtOneLimits.second);
    assert(pdtOtherLimits.first <= pdtOtherLimits.second);

}

void CProbabilityDensityDistribution::_GetSubIntegrationTotalInterval4PositiveTargetValuePositiveVariables4Covariance(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
{
    
    // invalidate result as initialization
    _InitIntegrationLimits(pdtOneLimits);
    _InitIntegrationLimits(pdtOtherLimits);
    
    // need the means very often
    CDigFloat dfOneMean = One.Mean();
    CDigFloat dfOtherMean = Other.Mean();    
    
    // do this only if the One AND Other distri have a positive part 
    if( (One.lastVariable()-dfOneMean) > 0 && (Other.lastVariable()-dfOtherMean) > 0)
    {
        
        // set the limits for the relevant  parts of the X-/Y-range used for the consecutive calculation for the given dfActTargetValue (tv)
        // dfXMin: 0 (if X-range covers also negative values) or lowest value of complete X-range otherwise 
        // dfXMax: largest value of complete X-range 
        // dfYMin: 0 (if Y-range covers also negative values) or lowest value of complete Y-range otherwise 
        // dfXMax: largest value of complete Y-range
        CDigFloat dfXMin = ( (One.firstVariable()-dfOneMean) <= 0) ? dfOneMean + CDigFloat(1).Error()*10 /*close to mean (greater)*/ : One.firstVariable();
        CDigFloat dfXMax = One.lastVariable();
        CDigFloat dfYMin = ((Other.firstVariable()-dfOtherMean) < 0) ? dfOtherMean + CDigFloat(1).Error()*10 /*close to mean (greater)*/  : Other.firstVariable();
        CDigFloat dfYMax = Other.lastVariable();
        
        // calculate the theoretical y limits from x
        CDigFloat dfYMinCalc = dfActTargetValue / (dfXMin-dfOneMean) + dfOtherMean;
        CDigFloat dfYMaxCalc = dfActTargetValue / (dfXMax-dfOneMean) + dfOtherMean;
        
        // now calculate the real limits for Y
        pdtOtherLimits.first = max(dfYMinCalc, dfYMin);
        pdtOtherLimits.second = min(dfYMaxCalc, dfYMax);
        
        // and now for X
        pdtOneLimits.first = dfActTargetValue / (pdtOtherLimits.second - dfOtherMean) + dfOneMean;
        pdtOneLimits.second = dfActTargetValue / (pdtOtherLimits.first- dfOtherMean) + dfOneMean;
        
        // DEBUG
        cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
        cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
        cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
        cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
        cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
        
        // no calculation possible if the extremal calculated limits exceed the variable ranges of x , y:
        if(pdtOneLimits.first > dfXMax || pdtOneLimits.second < dfXMin || pdtOtherLimits.first > dfYMax || pdtOtherLimits.second < dfYMin) 
        {
            _InitIntegrationLimits(pdtOneLimits);
            return;
            
        }   // endif(pdtOneLimits.first > dfXMax || pdtOneLimits.second < dfXMin || pdtOtherLimits.first > dfYMax || pdtOtherLimits.second < dfYMin)        
       
        
    }   // endif( (One.lastVariable()-dfOneMean) > 0 && (Other.lastVariable()-dfOtherMean) > 0)
    
    // check for non-zero / non-inverted ranges and the second (other) starts after the first (one)
    assert(pdtOneLimits.first <=  pdtOneLimits.second);
    assert(pdtOtherLimits.first <= pdtOtherLimits.second);
}




void CProbabilityDensityDistribution::_InitIntegrationLimits(PairDFType& pdtLimits)
{
    pdtLimits.first = 0; pdtLimits.first.ResetError();
    pdtLimits.second = 0; pdtLimits.second.ResetError();
}

void CProbabilityDensityDistribution::_GetRangeFromDistributionOperation(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, CDigFloat& TargetRangeStart, CDigFloat& TargetRangeEnd)
{
    CDigFloat ThisStart = firstVariable();
    CDigFloat ThisEnd = lastVariable();
    CDigFloat OtherStart = Other.firstVariable();
    CDigFloat OtherEnd = Other.lastVariable();
    
    // DEBUG
    cout << "This X variable = [" << ThisStart.RawPrint(10) << " ; " << ThisEnd.RawPrint(10) << endl;
    cout << "Other Y variable = [" << OtherStart.RawPrint(10) << " ; " << OtherEnd.RawPrint(10) << endl;
    
    // only needed for covariance calcuations
    CDigFloat dfThisMean = Mean();
    CDigFloat dfOtherMean = Other.Mean();
        
    switch(Operation)
    {
        case ProbDistOp::pdoPlus:
            TargetRangeStart = ThisStart + OtherStart;
            TargetRangeEnd = ThisEnd + OtherEnd;
            break;
            
        case ProbDistOp::pdoMinus:
            TargetRangeStart = ThisStart - OtherEnd;
            TargetRangeEnd = ThisEnd - OtherStart;
            break;
            
        case ProbDistOp::pdoCov :
            // covariance is multiplication shifted by first order mean
            ThisStart -= dfThisMean;
            ThisEnd -= dfThisMean;
            OtherStart -= dfOtherMean;
            OtherEnd -= dfOtherMean;
            
            // the rest is comparable to mult: so continue with mult. procedure
        case ProbDistOp::pdoMult:
            TargetRangeStart = min(min(ThisStart*OtherStart, ThisStart*OtherEnd),
                               min(ThisEnd*OtherStart,ThisEnd*OtherEnd));
            TargetRangeEnd = max(max(ThisStart*OtherStart, ThisStart*OtherEnd),
                               max(ThisEnd*OtherStart,ThisEnd*OtherEnd));
            break;
            
        case ProbDistOp::pdoDiv :
            
            // only set, if divisor is non-zero in every case:
            if(!(OtherStart <= 0 && OtherEnd >=0))
            {
                TargetRangeStart = min(min(ThisStart/OtherStart, ThisStart/OtherEnd),
                                min(ThisEnd/OtherStart,ThisEnd/OtherEnd));
                TargetRangeEnd = max(max(ThisStart/OtherStart, ThisStart/OtherEnd),
                                max(ThisEnd/OtherStart,ThisEnd/OtherEnd));
            }
            else
            {
                // in this case the target will range to infinity: so we take 20 times of the max range
                // symmetrically around zero 
                TargetRangeStart = max(ThisEnd- ThisStart, OtherEnd - OtherStart)*(-10);
                TargetRangeEnd = TargetRangeStart*(-1);
            }
            break;
            
        default:
            // not implemented
            assert(false);
    }   //endswitch(Operation)
    
    // DEBUG
    cout << "TargetRangeStart = " << TargetRangeStart.RawPrint(10) << endl;
    cout << "TargetRangeEnd = " << TargetRangeEnd.RawPrint(10) << endl;
}
