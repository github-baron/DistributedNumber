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
    CDigFloat dfPreFac = 1/sqrt(dfFac1*M_PI);
    
    /////////////////////////////////////////////
    // calculation of the gaussian curve
    // keeping intervals between the x-values
    // constant
    /////////////////////////////////////////////
//     //set the starting point and increment depending on the min. value.
//     CDigFloat dfXDelta = sqrt(dfFac1*-log(dfMinVal*sqrt(dfFac1*M_PI)));
//     CDigFloat dfXInc = 2*dfXDelta/(nPoints-1);
//     
//     for(int ipt=0;ipt<nPoints; ipt++)
//     {
//         CDigFloat dfX = dfMu-dfXDelta+ipt*dfXInc;
//         Add(dfX,dfPreFac*pow(M_El, -pow(dfX-dfMu,2)/dfFac1));
//     }

    
    /////////////////////////////////////////////
    // calculation of the gaussian curve
    // keeping intervals between the y-values constant
    // STEPS:
    // 1. the maximal amplitude is at x= dfMu --> dfPreFac
    // 2. y-steps are dfPreFac / nPoints
    // 3. for each y-step the value of the gaussian curve
    //    decrements about y-step, i.e. for the i-th
    //    step the y - value is 
    //    dfPreFac * (1-i/n)
    //    x_{+-i} is derived by solving the Gaussian
    //    function for x and one obtains:
    //    x_{+-i} = sqrt( 2 * dfSigma * ln(1- i / nPoints) )
    //    
    /////////////////////////////////////////////
    // set npoints to an even number 
    int nPtsEvenHalf = nPoints / 2;
    CDigFloat dfdY = dfPreFac / (nPtsEvenHalf);
    CDigFloat dfY = dfPreFac;
    CDigFloat dfRelPreFac = 1.;
    
    // add the first point: max at x = dfMu
    Add(dfMu,dfY);
    
    for(int ipt=1;ipt<nPtsEvenHalf; ipt++)
    {
        dfY -= dfdY;
        dfRelPreFac = dfY / dfPreFac;
        CDigFloat dfdXPos = sqrt(2*dfSigma*(-log(dfRelPreFac)));
        //DEBUG
//         cout << "ipt = " << to_string(ipt) << endl 
//              << "relPts = " << dfRelPreFac.Print(10) << endl 
//              << "y = " << dfY.Print(10) << endl 
//              << "dx = " << dfdXPos.Print(10) << endl 
//              << "x-i = " << (dfMu - dfdXPos).Print(10) << endl 
//              << "x+i = " << (dfMu + dfdXPos).Print(10) << endl; 
        Add(dfMu + dfdXPos,dfY);
        Add(dfMu - dfdXPos,dfY);
    }

    Normalize();
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
    /////////////////////////////////////////////////////
    // General purpose ot this function:
    //   caclculate <(x-<x>)*(y-<y>)> by solving the 
    //   double integration:
    //   Int(Int( (x-<x>)*px(x)* (y-<y>)*py(y)) dx) dy
    //   where:
    //   px: this distribution (probability density):
    //       px(x) = Xslope * x + XOffset
    //   py: Other distribution (probability density):
    //       px(x) = Yslope * x + YOffset
    //
    // Procedure:
    //   The solution of the upper double integral is mathematically
    //   simple becaus it is a polynome in x and y. The tricky
    //   thing is it consists of a lot of terms: that makes
    //   the solution look pretty complicated.
    //   With the solution:
    //   iterate over x and y intervals with constant
    //   linear slope and offset and sum up the 
    //   covariance
    //
    /////////////////////////////////////////////////////
    
    
    CDigFloat dfCovar = 0;
    
    // we need the mean values
    CDigFloat dfXMean = Mean();
    CDigFloat dfYMean = Other.Mean();
    
    //DEBUG
//     cout << "dfXMean " << dfXMean.RawPrint(10) << endl;
//     cout << "dfYMean " << dfYMean.RawPrint(10) << endl;
    
    // remember old x-value
    CDigFloat dfXOld;
    CDigFloat dfYOld;
    
    // DEBUG
    CDigFloat dfXPOld;
    CDigFloat dfYPOld;
    
    for(auto x: Distribution())
    {
        // set the old x-value
        if( !dfXOld.Valid() )
        {
            // set it ...
            dfXOld = x.first;
            dfXPOld = x.second;
            
            // hop to next x
            continue;
            
        }   //endif( !dfXOld.Valid() )
        
        // set all x-dependend variables needed for integration
        CDigFloat dfXMid = (x.first + dfXOld) / 2.;
        CDigFloat dfXSlope = LinearSlope(dfXMid);
        CDigFloat dfXOffset = LinearOffset(dfXMid);
        CDigFloat dfDx = (x.first - dfXOld);
        CDigFloat dfDx2 = (pow(x.first,2) -pow(dfXOld,2))/2.;
        CDigFloat dfDx3 = (pow(x.first,3) - pow(dfXOld,3))/3.;
        
        // helping variables used for simplification of integration
        CDigFloat dfC = dfXSlope*dfDx3 + (dfXOffset-dfXSlope*dfXMean)*dfDx2 - dfXOffset*dfXMean*dfDx;
        
        // continue 
        if(dfC == 0 )
            continue;
        
        //DEBUG
//         cout << "x = ["<< dfXOld.RawPrint(10) << " ; " << x.first.RawPrint(10) << "]" << endl;
//         cout << "px= ["<< dfXPOld.RawPrint(10) << " ; " << x.second.RawPrint(10) << "]" << endl;
//         cout << "dfXMid = "<< dfXMid.RawPrint(10)<< endl;
//         cout << "dfXOffset = "<< dfXOffset.RawPrint(10)<< endl;
//         cout << "dfXSlope = "<< dfXSlope.RawPrint(10)<< endl;
//         cout << "dfDx = "<< dfDx.RawPrint(10)<< endl;
//         cout << "dfDx2 = "<< dfDx2.RawPrint(10)<< endl;
//         cout << "dfDx3 = "<< dfDx3.RawPrint(10)<< endl;
//         cout << "dfC = "<< dfC.RawPrint(10)<< endl;
//         cout << "dfXSlope*dfDx3  = "<< (dfXSlope*dfDx3 ).RawPrint(10)<< endl;
//         cout << "(dfXOffset-dfXSlope*dfXMean)*dfDx2 = "<< ((dfXOffset-dfXSlope*dfXMean)*dfDx2).RawPrint(10)<< endl;
//         cout << "- dfXOffset*dfXMean*dfDx = "<< (- dfXOffset*dfXMean*dfDx).RawPrint(10)<< endl;
        
        // 
        for(auto y: Other.Distribution())
        {                        
            // set the old y-value
            if( !dfYOld.Valid() )
            {
                // set it ...
                dfYOld = y.first;
                dfYPOld = y.second;
                
                // hop to next y
                continue;
                
            }   //endif( !dfXOld.Valid() )
                
            // set all y-dependend variables needed for integration
            CDigFloat dfYMid = (y.first + dfYOld) / 2.;
            CDigFloat dfDy = (dfYOld - y.first);
            CDigFloat dfDy2 = (pow(y.first,2) - pow(dfYOld,2))/2.;
            CDigFloat dfDy3 = (pow(y.first,3) - pow(dfYOld,3))/3.;
            CDigFloat dfYSlope = Other.LinearSlope(dfYMid);
            CDigFloat dfYOffset =  Other.LinearOffset(dfYMid);
            
            // add to covariance:
            dfCovar += dfC *(dfDy3 + ( dfYOffset - dfYSlope * dfYMean)*dfDy2 - dfYOffset*dfYMean*dfDy);
        
            //DEBUG
//             cout << "y = ["<< dfYOld.RawPrint(10) << " ; " << y.first.RawPrint(10) << "]" << endl;
//             cout << "py= ["<< dfYPOld.RawPrint(10) << " ; " << y.second.RawPrint(10) << "]" << endl;
//             cout << "dfYMid = "<< dfYMid.RawPrint(10)<< endl;
//             cout << "dfYOffset = "<< dfYOffset.RawPrint(10)<< endl;
//             cout << "dfYSlope = "<< dfYSlope.RawPrint(10)<< endl;
//             cout << "dfDy = "<<  dfDy.RawPrint(10)<< endl;
//             cout << "dfDy2 = "<< dfDy2.RawPrint(10)<< endl;
//             cout << "dfDy3 = "<< dfDy3.RawPrint(10)<< endl;
//             cout << "dfCovar " << dfCovar.RawPrint(10) << endl;
            
            
            // increment dfYOld
            dfYOld = y.first;
            
                
        }   //endfor(int iy=0;iy<Other.Distribution().size()-1; iy++)
        
        // invalidate y
        dfYOld = myNAN;
        
        // increment dfXOld
        dfXOld = x.first;
        
    }   //endfor(int ix = 0; ix < Distribution().size()-1; ix++)
    
    // DEBUG
//     cout << endl;
    
    // return mean: should be now <(x-<x>)*(y-<y>)>
    return dfCovar;
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

////////////////////////////////////////////
// internal functions
////////////////////////////////////////////

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

    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("called with args:"));
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("CProbabilityDensityDistribution&:") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical","ProbDistOp: " + GetProbDistOpAsString(Operation));
    
    // generate plan : setting m_ConvolutionPlan;
    _SetConvolutionPlan(Other, Operation);
    
    // declare the result distri
    MapDFDFType TargetDistri;
    
    for(auto itarget: m_ConvolutionPlan)
    {
        TargetDistri[itarget.first] = 0;
        
        LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("init TargetDistri[" + itarget.first.RawPrint(15) + "] = 0"));
        
        // iterate over all sub integration limits
        for(auto intlim: itarget.second)
        {
            // DEBUG
//             cout << "x = [ " << intlim.first.first.RawPrint(10) << " ; " << intlim.first.second.RawPrint(10) << " ] " << endl;
//             cout << "y = [ " << intlim.second.first.RawPrint(10) << " ; " << intlim.second.second.RawPrint(10) << " ] " << endl;
            
            // add to all other probs of this target value
            TargetDistri[itarget.first] += _GetConvolutionIntegralAnalytical4Operation(Other, Operation, itarget.first, intlim);
            
            // logging trace 
            LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("incremented TargetDistri[" + itarget.first.RawPrint(15) + "] = " + TargetDistri[itarget.first].RawPrint(15)));
        
            
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
        // DEBUG
//         if(AbsIntegral() != 1)
//             cout << "AbsIntegral != 1: " << AbsIntegral().RawPrint(30) << endl << endl;
        
        assert(AbsIntegral() == 1);
        
    }   // endif(!(TargetRangeStart == TargetRangeEnd && TargetRangeStart == 0))
    
    return *this;
 
}

CDigFloat CProbabilityDensityDistribution::_GetConvolutionIntegralAnalytical4Operation(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, const CDigFloat& dfTargetValue, const SubIntLimitsType& vXYLimits)
{  
    // init the result
    CDigFloat dfResult = 0;
    
    // derive all necessary arguments needed for the analytical calculation of the 
    // convolution integrals
    // calculate slope and offset for the related interval of this distri
    CDigFloat dfXMean = (vXYLimits.first.first + vXYLimits.first.second)/2.;
    MapDFDFType::const_iterator Left, Right;
    GetInterval(dfXMean, Left, Right);
    CDigFloat dfXOffset = _LinearOffset(Left, Right);
    CDigFloat dfXSlope = _LinearSlope(Left,Right);
    
    // calculate slope and offset for the related interval of other distri
    MapDFDFType::const_iterator LeftOther, RightOther;
    CDigFloat dfYMean = (vXYLimits.second.first + vXYLimits.second.second)/2.;
    Other.GetInterval(dfYMean, LeftOther, RightOther);
    CDigFloat dfYOffset = Other._LinearOffset(LeftOther, RightOther);
    CDigFloat dfYSlope = Other._LinearSlope(LeftOther,RightOther);
    
    // trace logging 
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("preparation of convolution integral calculation"));
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("X = [") + vXYLimits.first.first.RawPrint(10) + " ; " + vXYLimits.first.second.RawPrint(10) + " ]");
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("X-distri = [") + Left->first.RawPrint(10) + " ; " + Right->first.RawPrint(10) + " ]");
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfXMean = ") + dfXMean.RawPrint(10));
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfXOffset = ") + dfXOffset.RawPrint(10));
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfXSlope = ") + dfXSlope.RawPrint(10) );
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("Y = [") + vXYLimits.second.first.RawPrint(10) + " ; " + vXYLimits.second.second.RawPrint(10) + " ]");
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("Y-distri = [") + LeftOther->first.RawPrint(10) + " ; " + RightOther->first.RawPrint(10) + " ]");
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfYMean = ") + dfYMean.RawPrint(10));
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfYOffset = ") + dfYOffset.RawPrint(10));
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfYSlope = ") + dfYSlope.RawPrint(10));
   
     
    // variables for primitive integral values of upper and lower limit
    CDigFloat dfPrimIntValLow, dfPrimIntValUp;
    
    switch( Operation)
    {
        case ProbDistOp::pdoPlus:
            dfPrimIntValUp = _PrimitiveIntegral4TargetValue4Addition(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.second);
            dfPrimIntValLow = _PrimitiveIntegral4TargetValue4Addition(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.first);
            break;
        case ProbDistOp::pdoMinus:
            dfPrimIntValUp = _PrimitiveIntegral4TargetValue4Subtraction(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.second);
            dfPrimIntValLow =  _PrimitiveIntegral4TargetValue4Subtraction(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.first);
            break;
        case ProbDistOp::pdoMult:
            dfPrimIntValUp = _PrimitiveIntegral4TargetValue4Multiplication(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.second);
            dfPrimIntValLow = _PrimitiveIntegral4TargetValue4Multiplication(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.first);
            break;
        case ProbDistOp::pdoDiv:
            dfPrimIntValUp = _PrimitiveIntegral4TargetValue4Division(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.second);
            dfPrimIntValLow = _PrimitiveIntegral4TargetValue4Division(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.first);
            break;
        default:
            break;
    }   //endswitch( Operation)
    
    // calculate the result 
    dfResult = dfPrimIntValUp - dfPrimIntValLow;
    // trace 
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation", "convolution integral is :");
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("result = ") + dfResult.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("upper value = ") + dfPrimIntValUp.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("lower value = ") + dfPrimIntValLow.RawPrint(10));
    
    // error case
    if(dfResult < 0)
    {
        LOGERROR("ErrorLogger","CProbabilityDensityDistribution::_GetConvolutionIntegralAnalytical4Operation: convolution integral is negative :");
        LOGERROR("ErrorLogger",string("result = ") + dfResult.RawPrint(10));
        LOGERROR("ErrorLogger",string("upper value = ") + dfPrimIntValUp.RawPrint(10));
        LOGERROR("ErrorLogger",string("lower value = ") + dfPrimIntValLow.RawPrint(10));
    }
    
    
    assert(dfResult >= 0);
    return dfResult;
    
}

void CProbabilityDensityDistribution::_SetConvolutionPlan(CProbabilityDensityDistribution& Other, const ProbDistOp Operation)
{   
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan", string("called with args:\n") +
    "CProbabilityDensityDistribution&:" + Other.PrintMetaInfo() + "\n" +
    "ProbDistOp: " + GetProbDistOpAsString(Operation) +"\n" +
    "this: " + PrintMetaInfo()
    );
 
    
    CDigFloat dfTargetStart, dfTargetEnd;
    _GetRangeFromDistributionOperation(Other, Operation, dfTargetStart, dfTargetEnd);
    
    // trace 
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("calculated target range = [") + dfTargetStart.RawPrint(15) + ", " + dfTargetEnd.RawPrint(15));
    
    
    // get the target step
    CDigFloat dfTargetStep = (dfTargetEnd - dfTargetStart)/(CDigFloat(IntegrationSteps()));
    
    // trace
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("calculated target step size = ") + dfTargetStep.RawPrint(15) + "( " + to_string(IntegrationSteps()) + " ) ");
    
  
    // clear the plan first
    m_ConvolutionPlan.clear();
    
    // now start building the plan 
    CDigFloat dfActTargetValue = dfTargetStart;
    bool bZeroTargetValueDetected = false;
    bool bCheck4ZeroTargetValue = ( (Operation == ProbDistOp::pdoDiv ) || ( Operation== ProbDistOp::pdoMult));
    while(dfActTargetValue < (dfTargetEnd +dfTargetStep))
    {

        LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("iteration dfActTargetValue = ") + dfActTargetValue.RawPrint(15));
        
        // defining the complete integration range for this and other distribution
        _SetSubIntegrationIntervals4TargetValue(Operation, dfActTargetValue, Other);        
        
        
        // deal with the incrementation of target value
        // in case a zero value was found before ....
        if(bZeroTargetValueDetected)
        {
            
            LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("former dfActTargetValue was zero:"));
            
            // ... add only half a step: now we are in sync again
            dfActTargetValue += dfTargetStep/2.;
            
            
            // no zero anymore 
            bZeroTargetValueDetected = false;
            
        }   // endif(bZeroTargetValueDetected)
        else            
            // increment actual target value 
            dfActTargetValue += dfTargetStep;
        
        LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue set to ") + dfActTargetValue.RawPrint(15));
                
        // split zero dftarget value due to unipolarity trouble for multiplication and division
        if(bCheck4ZeroTargetValue)
            if(dfActTargetValue == 0 && dfActTargetValue < (dfTargetEnd +dfTargetStep))
            {
                LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue zero detected:"));
                
                // decrement half a step 
                dfActTargetValue -= dfTargetStep/2.;
        
                LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue decremented half a step to ") + dfActTargetValue.RawPrint(15));    
                
                // defining the complete integration range for this and other distribution ... again
                _SetSubIntegrationIntervals4TargetValue(Operation, dfActTargetValue, Other); 
                
                // increment about one step ... now we are half a step behind because we split the target value == 0
                dfActTargetValue += dfTargetStep;
                
                LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue incremented one step to ") + dfActTargetValue.RawPrint(15));  
                
                // remember we found a zero and deal with it the next round 
                bZeroTargetValueDetected = true;
                
                
            }   // endif(dfTargetValue.RawValue() == 0)
        
    }   // endwhile(dfActTargetValue < (dfTargetEnd +dfTargetStep))
    
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

/**
 *@details The sub integration intervals are needed to differentiate parts of the total integration interval (see CProbabilityDensityDistribution::_GetTotalIntegrationInterval4TargetValue) where the linear function for x or y distribution changes.
 * 
 */
void CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue(const ProbDistOp Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other)
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
    //       slope and offset for x and y for the all total
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
    
    // keeps sub integration limits
    VectorSubIntLimitsType SubIntervals4TargetValue;
    
    // keeps the total integration interval(s) for a target value
    VectorSubIntLimitsType TotalIntervals4TargetValue;
    
    // fill the x intervals: must be unipolar for multiplication and division
    VectorPairDFType vpdXIntervals;
    
    // DEBUG
//     cout << "_SetIntegrationIntervals4TargetValue" << endl;
//     cout << "Operation = " <<  to_string((int)Operation ) << endl;
//     cout << "dfTargetValue = " << dfTargetValue.RawPrint(30) << endl;
    
    // I 1. get the x-interval(s) for calculation: must be unipolar
    //      for division and multiplication
    _GetXInterval4Operation(Operation, vpdXIntervals);
        
    // I 2. derive closed intervals for x and y (kept in TotalIntervals4TargetValue) mapped onto each other bijectively
    //      by the calculation of the complementary variable
    _GetTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, vpdXIntervals, TotalIntervals4TargetValue);
    
    
    _GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, TotalIntervals4TargetValue, SubIntervals4TargetValue);  
    
    // now add the sub integration limits and its corresponding target value as pair
    m_ConvolutionPlan.push_back(pair<CDigFloat, VectorSubIntLimitsType>(dfTargetValue, SubIntervals4TargetValue)); 
}
/**
 * @details In case of multiplication and division the x interval must be unipolar, i.e. purely negativ or positive. Hence, in case of a bipolar distribution there might be two x-intervals which will integrated separately.
 */
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
//     cout << "_GetXInterval4Operation:" << endl;
//     for(auto iint: vXIntervals)
//         cout << "[ " << iint.first.RawPrint(10) << " ; " << iint.second.RawPrint(10) << " ]" << endl;
//     cout << endl;

}

void CProbabilityDensityDistribution:: _GetTotalIntegrationInterval4TargetValue(const ProbDistOp& Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const VectorPairDFType& vpdXIntervals, VectorSubIntLimitsType& TotalIntervals4TargetValue)
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
//         cout << "xint.first = " << xint.first.RawPrint(10) << endl;
//         cout << "xint.second = " << xint.second.RawPrint(10) << endl;
        
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
//         cout << "dfYCompMin = " << dfYCompMin.RawPrint(10) << endl;
//         cout << "dfYCompMax = " << dfYCompMax.RawPrint(10) << endl;
        
        
        /////////////////////////////////////////////////
        // calculate the overlap of the complementary 
        // ycomp and the original y interval --> YOverlap
        /////////////////////////////////////////////////
        CDigFloat dfYOverlapMin, dfYOverlapMax;
        dfYOverlapMin = max(dfYCompMin,Other.firstVariable());
        dfYOverlapMax = min(dfYCompMax,Other.lastVariable());  
        
        //DEBUG
//         cout << "dfYOverlapMin = " << dfYOverlapMin.RawPrint(10) << endl;
//         cout << "dfYOverlapMax = " << dfYOverlapMax.RawPrint(10) << endl;
        
        // the overlap is empty in case dfYOverlapMax <= dfYOverlapMin
        if(dfYOverlapMax <= dfYOverlapMin)
        {
            // DEBUG
//             cout << "empty interval: dfYOverlapMax(" << dfYOverlapMax.RawPrint(4) << ")  <= dfYOverlapMin(" << dfYOverlapMin.RawPrint(4) << ")" << endl;
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
//         cout << "dfXCompMin = " << dfXCompMin.RawPrint(10) << endl;
//         cout << "dfXCompMax = " << dfXCompMax.RawPrint(10) << endl;
        
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
//             cout << "empty interval: dfXOverlapMax (" << dfXOverlapMax.RawPrint(4) << ")  <= dfXOverlapMin(" << dfXOverlapMin.RawPrint(4) << ")" << endl;
            continue;
        }
        // 
//         cout << "dfXOverlapMin = " << dfXOverlapMin.RawPrint(10) << endl;
//         cout << "dfXOverlapMax = " << dfXOverlapMax.RawPrint(10) << endl;
//         
//         cout << "X=[ " <<  dfXOverlapMin.RawPrint(10) << " ; " << dfXOverlapMax.RawPrint(10) << " ] " << endl;
//         cout << "Y=[ " <<  dfYOverlapMin.RawPrint(10) << " ; " << dfYOverlapMax.RawPrint(10) << " ] " << endl;
         
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
//      cout << "_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue:" << endl;
     
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
//         cout << "x = [" << xyint.first.first.RawPrint(10) << " ; " << xyint.first.second.RawPrint(10) << " ]" << endl;
//         cout << "y = [" << xyint.second.first.RawPrint(10) << " ; " << xyint.second.second.RawPrint(10) << " ]" << endl;
        
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
//             cout << "x=" << x.first.RawPrint(10) << endl;
            
            // break conditions
            if(x.first > xyint.first.second)
                break;
            
            if(x.first > xyint.first.first && x.first <xyint.first.second)
            {
            
                vdfSubPoints.push_back(x.first);
                
                //DEBUG
//                 cout << "added " << x.first.RawPrint(10) << endl;
                
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
//             cout << "y=" << x.first.RawPrint(10) << endl;
            
            // break conditions
            if(x.first > xyint.second.second)
                break;
            
            if(x.first > xyint.second.first && x.first <xyint.second.second)
            {
                vdfSubPoints.push_back(_GetComplementaryVariable(Operation,dfTargetValue,x.first,true));
                
                //DEBUG
//                 cout << "added " << _GetComplementaryVariable(Operation,dfTargetValue,x.first,true).RawPrint(10) << endl;
            }
            
        }   //endfor(auto x: Other.Distribution())
        
        //////////////////////////////////////////////
        // sort the supporting points and derive
        // sub intervals for x and y
        //////////////////////////////////////////////
        sort(vdfSubPoints.begin(),vdfSubPoints.end());
        //DEBUG
//         for(auto x: vdfSubPoints)
//             cout << "SubPoint: " << x.RawPrint(30) << endl;
        
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
//             cout << "adding:" << endl;
//             cout << "x = [" << dfXSubLimit.first.RawPrint(10) << " ; " << dfXSubLimit.second.RawPrint(10) << " ]" << endl;
//             cout << "y = [" << dfYSubLimit.first.RawPrint(10) << " ; " << dfYSubLimit.second.RawPrint(10) << " ]" << endl;
            
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
            // the sign must be kept
            assert(sgn(dfComplementaryVariable)*sgn(dfVariable) == sgn(dfTargetValue));
            break;
        case ProbDistOp::pdoDiv:
            if(bVarIsYVar)
                dfComplementaryVariable = dfTargetValue * dfVariable;
            else
                dfComplementaryVariable = dfVariable / dfTargetValue ;
            // the sign must be kept
            assert(sgn(dfComplementaryVariable)*sgn(dfVariable) == sgn(dfTargetValue));
            break;
        default:
            break;
            
    }   //endswitch(Operation)
    
    return dfComplementaryVariable;
}

void CProbabilityDensityDistribution::_GetRangeFromDistributionOperation(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, CDigFloat& TargetRangeStart, CDigFloat& TargetRangeEnd)
{
    CDigFloat ThisStart = firstVariable();
    CDigFloat ThisEnd = lastVariable();
    CDigFloat OtherStart = Other.firstVariable();
    CDigFloat OtherEnd = Other.lastVariable();
    
    // DEBUG
//     cout << "This X variable = [" << ThisStart.RawPrint(10) << " ; " << ThisEnd.RawPrint(10) << endl;
//     cout << "Other Y variable = [" << OtherStart.RawPrint(10) << " ; " << OtherEnd.RawPrint(10) << endl;
    
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
//     cout << "TargetRangeStart = " << TargetRangeStart.RawPrint(10) << endl;
//     cout << "TargetRangeEnd = " << TargetRangeEnd.RawPrint(10) << endl;
}

////////////////////////////////////////////////////////
// external functions
////////////////////////////////////////////////////////

string GetProbDistOpAsString( ProbDistOp op )
{
    string strOp;
    switch(op)
    {
        case ProbDistOp::pdoUnknown: strOp = "?"; break;
        case ProbDistOp::pdoPlus: strOp = "+"; break;
        case ProbDistOp::pdoMinus: strOp = "-"; break;
        case ProbDistOp::pdoMult: strOp = "*"; break;
        case ProbDistOp::pdoDiv: strOp = "/"; break;
        case ProbDistOp::pdoCov: strOp = "cov"; break;
        default:
            break;            
    }   //endswitch(op)
    
    return strOp;
};
