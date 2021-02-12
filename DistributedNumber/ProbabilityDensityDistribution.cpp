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
// arithmetic operator
////////////////////////////////////////////
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator+=(CProbabilityDensityDistribution& other)
{
    
//     return _GeneralOperatorsFunction(other, ProbDistOp::pdoPlus );
//     return _Convolution4GeneralOperatorsNew(other, ProbDistOp::pdoPlus );
    return _GeneralOperatorsFunction(other, ProbDistOp::pdoPlus );
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
    return _GeneralOperatorsFunction(other, ProbDistOp::pdoMinus );
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
    return _GeneralOperatorsFunction(other, ProbDistOp::pdoMult );
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
    return _GeneralOperatorsFunction(other, ProbDistOp::pdoDiv);
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator/(CProbabilityDensityDistribution& other)
{
    CProbabilityDensityDistribution ProbDistResult(*this);
    ProbDistResult /= other;
    return ProbDistResult;
}



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
void CProbabilityDensityDistribution::GetIntervall(const CDigFloat& variable, M_DFDF::const_iterator& VariableLeft, M_DFDF::const_iterator& VariableRight)
{
    Normalize();
    return CDistribution::GetIntervall(variable,VariableLeft,VariableRight);
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
CProbabilityDensityDistribution & CProbabilityDensityDistribution::_GeneralOperatorsFunction(CProbabilityDensityDistribution& other, const ProbDistOp Operation)
{
    // generate plan : setting m_ConvolutionPlan;
    _SetConvolutionPlan(other, Operation);
    
    // declare the result distri
    M_DFDF TargetDistri;
    
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
            CDigFloat dfProbOther = other.AbsIntegral(intlim.second.first, intlim.second.second);
            
            // DEBUG
//             cout << "ProbThis(" << intlim.first.first.RawPrint(15) << " - " << intlim.first.second.RawPrint(15) << ") = " << dfProbThis.RawPrint(15) << endl;
//             cout << "ProbOther(" << intlim.second.first.RawPrint(15) << " - " << intlim.second.second.RawPrint(15) << ") = " << dfProbOther.RawPrint(15) << endl;
            
            // add to all other probs of this target value
            TargetDistri[itarget.first] += dfProbThis*dfProbOther;
            
        }   // endfor(auto intervall: itarget.second)        
        
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

void CProbabilityDensityDistribution::_SetConvolutionPlan(CProbabilityDensityDistribution& other, const ProbDistOp Operation)
{
    CDigFloat dfTargetStart, dfTargetEnd;
    _GetRangeFromDistributionOperation(other, Operation, dfTargetStart, dfTargetEnd);
    
    // get the target step
    CDigFloat dfTargetStep = (dfTargetEnd - dfTargetStart)/(CDigFloat(IntegrationSteps()));
    
    // DEBUG
    cout << "target intervall: [" << dfTargetStart.RawPrint(15) << " ; " << dfTargetEnd.RawPrint(15) << " ] " << endl;
    
    // clear the plan first
    m_ConvolutionPlan.clear();
    
    // now start building the plan 
    CDigFloat dfActTargetValue = dfTargetStart;
    while(dfActTargetValue <= dfTargetEnd)
    {
        // this is the vector defining the sub integration limits for a given targe value
        VectorSubIntLimitsType SubIntLimits;
        
        // defining the complete integration range for this and other distribution
        _SetSubIntegrationIntervalls4TargetValue(Operation, dfActTargetValue, other, SubIntLimits);
        
        // now add the sub integration limits and its corresponding target value as pair
        m_ConvolutionPlan.push_back(pair<CDigFloat, VectorSubIntLimitsType>(dfActTargetValue, SubIntLimits));
        
        // increment actual target value 
        dfActTargetValue += dfTargetStep;
    }   // endwhile(dfActTargetValue <= dfTargetEnd)
    
}

void CProbabilityDensityDistribution::_SetSubIntegrationIntervalls4TargetValue(const ProbDistOp Operation, const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervalls4TargetValue)
{
    switch(Operation)
    {
        case ProbDistOp::pdoPlus:
            _SetSubIntegrationIntervalls4TargetValue4Addition(dfActTargetValue, other, SubIntervalls4TargetValue);
            break;
        case ProbDistOp::pdoMinus:
            _SetSubIntegrationIntervalls4TargetValue4Subtraction(dfActTargetValue, other, SubIntervalls4TargetValue);
            break;
        case ProbDistOp::pdoMult:
            _SetSubIntegrationIntervalls4TargetValue4Multiplication(dfActTargetValue, other, SubIntervalls4TargetValue);
            break;
        case ProbDistOp::pdoDiv:
            _SetSubIntegrationIntervalls4TargetValue4Division(dfActTargetValue, other, SubIntervalls4TargetValue);
            break;
        default:
            break;
    }   // endswitch(Operation)
}

void CProbabilityDensityDistribution::_SetSubIntegrationIntervalls4TargetValue4Addition(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& Other, VectorSubIntLimitsType& SubIntervalls4TargetValue)
{
    // init result:
    SubIntervalls4TargetValue.clear();
    
    // now set local variables
    VectorSubIntLimitsType SubTotalIntervallRanges4TargetValue;
    PairDFType pdtThisLimits, pdtOtherLimits;
    SubIntLimitsType SingleSubTotalIntervallRange;
    
    
    // DEBUG
    cout << endl << "_SetSubIntegrationIntervalls4TargetValue4Addition:" << endl;
    cout << "dfActTargetValue = " << dfActTargetValue.RawPrint(15) << endl;
    cout << "This-range: [" << firstVariable().RawPrint(15) << " ; " << lastVariable().RawPrint(15) << " ]" << endl;
    cout << "Other-range: [" << Other.firstVariable().RawPrint(15) << " ; " << Other.lastVariable().RawPrint(15) << " ]" << endl;

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
    cout << "dfActTargetValue: " << dfActTargetValue.RawPrint(15) << endl;
    cout << "X-range: [" << dfXMin.RawPrint(15) << " ; " << dfXMax.RawPrint(15) << " ]" << endl;
    cout << "Y-range: [" << dfYMin.RawPrint(15) << " ; " << dfYMax.RawPrint(15) << " ]" << endl;
    cout << "dfYMinCalc: " << dfYMinCalc.RawPrint(15) << endl;
    cout << "dfYMaxCalc: " << dfYMaxCalc.RawPrint(15) << endl;
    
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
    cout << "X-Limits: [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << " ]" << endl;
    cout << "Y-Limits: [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << " ]" << endl;

        
    // check for x (this distri) being within range
    if(pdtThisLimits.first > dfXMax || pdtThisLimits.second < dfXMin) 
        return;
        
    // add only if everything is ok
    cout << "this: first < second : " << string((pdtThisLimits.first < pdtThisLimits.second) ? "less" : " ge") << endl;
    cout << "other: first < second : " << string((pdtOtherLimits.first < pdtOtherLimits.second) ? "less" : "ge") << endl;
    // add only valid results
    if((pdtThisLimits.first >= pdtThisLimits.second) || (pdtOtherLimits.first >= pdtOtherLimits.second))
    {
        //DEBUG
        cout << "zero intervall or inverted order" << endl; 
        return;
    }

    SubTotalIntervallRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
    cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
    cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
    
        // now derive for all elements (total integration intervalls) of SubTotalIntervallRanges4TargetValue (max. 2)
    // the SubIntegrationSteps sub integration intervalls and add them to SubIntervalls4TargetValue
         // do we have negative values
    bool bFirstNegative = (pdtThisLimits.first<0);
    bool bSecondNegative = (pdtOtherLimits.first<0);
    
    // calculate the factor needed for dividing the integration sub intervalls from the total intervall
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
        
        // now add the intervall
        SubIntervalls4TargetValue.push_back(SubIntLimitsType(pdtThis, pdtOther));
    
    }   // endfor(int isub = 0; isub < SubIntegrationSteps(); isub++)
    
}
void CProbabilityDensityDistribution::_SetSubIntegrationIntervalls4TargetValue4Subtraction(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& Other, VectorSubIntLimitsType& SubIntervalls4TargetValue)
{
        // init result:
    SubIntervalls4TargetValue.clear();
    
    // now set local variables
    VectorSubIntLimitsType SubTotalIntervallRanges4TargetValue;
    PairDFType pdtThisLimits, pdtOtherLimits;
    SubIntLimitsType SingleSubTotalIntervallRange;
    
    
    // DEBUG
//     cout << endl << "_SetSubIntegrationIntervalls4TargetValue4Subtraction:" << endl;
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
//         cout << "zero intervall or inverted order" << endl; 
        return;
    }

    SubTotalIntervallRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
    // DEBUG
    //     cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//     cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
    
        // now derive for all elements (total integration intervalls) of SubTotalIntervallRanges4TargetValue (max. 2)
    // the SubIntegrationSteps sub integration intervalls and add them to SubIntervalls4TargetValue
         // do we have negative values
    bool bFirstNegative = (pdtThisLimits.first<0);
    bool bSecondNegative = (pdtOtherLimits.first<0);
    
    // calculate the factor needed for dividing the integration sub intervalls from the total intervall
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
        
        // now add the intervall
        SubIntervalls4TargetValue.push_back(SubIntLimitsType(pdtThis, pdtOther));
    
    }   // endfor(int isub = 0; isub < SubIntegrationSteps(); isub++)

}
void CProbabilityDensityDistribution::_SetSubIntegrationIntervalls4TargetValue4Multiplication(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervalls4TargetValue)
{
    // init result:
    SubIntervalls4TargetValue.clear();
    
    // now set local variables
    VectorSubIntLimitsType SubTotalIntervallRanges4TargetValue;
    PairDFType pdtThisLimits, pdtOtherLimits;
    SubIntLimitsType SingleSubTotalIntervallRange;
    
    
    // DEBUG
//     cout << endl << "_SetSubIntegrationIntervalls4TargetValue4Multiplication:" << endl;
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
            _GetSubIntegrationTotalIntervall4NegativeTargetValue4Multiplication(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
             
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;

            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervallRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
                
                // DEBUG
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
        }
        
        // add another limits
        if(other.firstVariable() < 0 && lastVariable() > 0)
        {
            _GetSubIntegrationTotalIntervall4NegativeTargetValue4Multiplication(dfActTargetValue,other,*this,pdtOtherLimits,pdtThisLimits);
            
            // DEBUG
//                 cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//                 cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;

            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervallRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
                
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
            
            _GetSubIntegrationTotalIntervall4PositiveTargetValueNegativeVariables4Multiplication(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
            
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;
            
            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervallRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

            // DEBUG
//             cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//             cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
            
        }   // endif( firstVariable() < 0 && other.firstVariable() < 0)
        
        // check positive case 
        if( lastVariable() > 0 && other.lastVariable() > 0)
        {
            
            _GetSubIntegrationTotalIntervall4PositiveTargetValuePositiveVariables4Multiplication(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
            
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;

            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervallRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

                // DEBUG
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
            
        }   // endif( lastVariable() > 0 && other.lastVariable() > 0)
        
    }   // endelse(dfActTargetValue < 0)
    
    // now derive for all elements (total integration intervalls) of SubTotalIntervallRanges4TargetValue (max. 2)
    // the SubIntegrationSteps sub integration intervalls and add them to SubIntervalls4TargetValue
    for(auto totIv: SubTotalIntervallRanges4TargetValue)
    {
        // do we have negative values
        bool bFirstNegative = (totIv.first.first<0);
        bool bSecondNegative = (totIv.second.first<0);
        
        // calculate the factor needed for dividing the integration sub intervalls from the total intervall
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
            
            // now add the intervall
            SubIntervalls4TargetValue.push_back(SubIntLimitsType(pdtThis, pdtOther));
            
        }   // endfor(int isub = 0; isub < SubIntegrationSteps(); isub++)
         
    }   // endfor(auto totIntervalls: SubTotalIntervallRanges4TargetValue)
    
}

void CProbabilityDensityDistribution::_GetSubIntegrationTotalIntervall4NegativeTargetValue4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType &pdtOneLimits, PairDFType& pdtOtherLimits )
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

void CProbabilityDensityDistribution::_GetSubIntegrationTotalIntervall4PositiveTargetValuePositiveVariables4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType &pdtOneLimits, PairDFType& pdtOtherLimits )
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

void CProbabilityDensityDistribution::_GetSubIntegrationTotalIntervall4PositiveTargetValueNegativeVariables4Multiplication(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
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


void CProbabilityDensityDistribution::_SetSubIntegrationIntervalls4TargetValue4Division(const CDigFloat& dfActTargetValue, CProbabilityDensityDistribution& other, VectorSubIntLimitsType& SubIntervalls4TargetValue)
{
       // init result:
    SubIntervalls4TargetValue.clear();
    
    // now set local variables
    VectorSubIntLimitsType SubTotalIntervallRanges4TargetValue;
    PairDFType pdtThisLimits, pdtOtherLimits;
    SubIntLimitsType SingleSubTotalIntervallRange;
    
    
    // DEBUG
//     cout << endl << "_SetSubIntegrationIntervalls4TargetValue4Division:" << endl;
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
//             cout << "calling _GetSubIntegrationTotalIntervall4NegativeTargetValueNegPosVariables4Division " << endl;
   
            _GetSubIntegrationTotalIntervall4NegativeTargetValueNegPosVariables4Division(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
            
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;
 
            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervallRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));
            
                // DEBUG
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
        }
        
        // add another limits
        if(lastVariable() > 0 && other.firstVariable() < 0)
        {
            // DEBUG
//             cout << "calling _GetSubIntegrationTotalIntervall4NegativeTargetValuePosNegVariables4Division " << endl;
   
            _GetSubIntegrationTotalIntervall4NegativeTargetValuePosNegVariables4Division(dfActTargetValue,*this,other,pdtThisLimits,pdtOtherLimits);
            
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;

            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervallRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

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
//             cout << "calling _GetSubIntegrationTotalIntervall4PositiveTargetValueNegativeVariables4Division " << endl;
            
            _GetSubIntegrationTotalIntervall4PositiveTargetValueNegativeVariables4Division(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
   
            // DEBUG             
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;
            
            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervallRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

                // DEBUG                
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
            
        }   // endif( firstVariable() < 0 && other.firstVariable() < 0)
        
        // check positive case 
        if( lastVariable() > 0 && other.lastVariable() > 0)
        {
            
            // DEBUG
//             cout << "calling _GetSubIntegrationTotalIntervall4PositiveTargetValuePositiveVariables4Division " << endl;
            _GetSubIntegrationTotalIntervall4PositiveTargetValuePositiveVariables4Division(dfActTargetValue,*this,other,pdtThisLimits, pdtOtherLimits);
            
            // DEBUG
//             cout << "this: first != second : " << string((pdtThisLimits.first != pdtThisLimits.second) ? "different" : "equal") << endl;
//             cout << "other: first != second : " << string((pdtOtherLimits.first != pdtOtherLimits.second) ? "different" : "equal") << endl;
            // add only valid results
            if((pdtThisLimits.first != pdtThisLimits.second) && (pdtOtherLimits.first != pdtOtherLimits.second))
            {
                SubTotalIntervallRanges4TargetValue.push_back(SubIntLimitsType(pdtThisLimits, pdtOtherLimits));

                // DEBUG                
//                 cout << "adding : [" << pdtThisLimits.first.RawPrint(15) << " ; " << pdtThisLimits.second.RawPrint(15) << endl;
//                 cout << "adding : [" << pdtOtherLimits.first.RawPrint(15) << " ; " << pdtOtherLimits.second.RawPrint(15) << endl << endl; 
            }
            
        }   // endif( lastVariable() > 0 && other.lastVariable() > 0)
        
    }   // endelse(dfActTargetValue < 0)
    
    // now derive for all elements (total integration intervalls) of SubTotalIntervallRanges4TargetValue (max. 2)
    // the SubIntegrationSteps sub integration intervalls and add them to SubIntervalls4TargetValue
    for(auto totIv: SubTotalIntervallRanges4TargetValue)
    {
        // do we have negative values
        bool bFirstNegative = (totIv.first.first<0);
        bool bSecondNegative = (totIv.second.first<0);
        
        // calculate the factor needed for dividing the integration sub intervalls from the total intervall
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
            
            // now add the intervall
            SubIntervalls4TargetValue.push_back(SubIntLimitsType(pdtThis, pdtOther));
            
        }   // endfor(int isub = 0; isub < SubIntegrationSteps(); isub++)
         
    }   // endfor(auto totIntervalls: SubTotalIntervallRanges4TargetValue)
 
}

void CProbabilityDensityDistribution::_GetSubIntegrationTotalIntervall4PositiveTargetValueNegativeVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
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
void CProbabilityDensityDistribution::_GetSubIntegrationTotalIntervall4PositiveTargetValuePositiveVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
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

void CProbabilityDensityDistribution::_GetSubIntegrationTotalIntervall4NegativeTargetValuePosNegVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
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

void CProbabilityDensityDistribution::_GetSubIntegrationTotalIntervall4NegativeTargetValueNegPosVariables4Division(const CDigFloat dfActTargetValue, CProbabilityDensityDistribution& One, CProbabilityDensityDistribution& Other, PairDFType& pdtOneLimits, PairDFType& pdtOtherLimits)
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


void CProbabilityDensityDistribution::_InitIntegrationLimits(PairDFType& pdtLimits)
{
    pdtLimits.first = 0; pdtLimits.first.ResetError();
    pdtLimits.second = 0; pdtLimits.second.ResetError();
}

void CProbabilityDensityDistribution::_GetRangeFromDistributionOperation(CProbabilityDensityDistribution& other, const ProbDistOp Operation, CDigFloat& TargetRangeStart, CDigFloat& TargetRangeEnd)
{
    CDigFloat ThisStart = firstVariable();
    CDigFloat ThisEnd = lastVariable();
    CDigFloat OtherStart = other.firstVariable();
    CDigFloat OtherEnd = other.lastVariable();
        
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
                TargetRangeStart = max(ThisEnd- ThisStart, OtherEnd - OtherStart)*(-10);
                TargetRangeEnd = TargetRangeStart*(-1);
            }
            break;
            
    }   //endswitch(Operation)
}
