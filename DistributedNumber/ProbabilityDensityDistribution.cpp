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
    return _GeneralOperatorsFunctionAnalytical3(other, ProbDistOp::pdoPlus );
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator+(CProbabilityDensityDistribution& other)
{
    CProbabilityDensityDistribution ProbDistResult(*this);
    ProbDistResult += other;
    return ProbDistResult;
}

CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator-=(CProbabilityDensityDistribution& other)
{
    return _GeneralOperatorsFunctionAnalytical3(other, ProbDistOp::pdoMinus );
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator-(CProbabilityDensityDistribution& other)
{
    CProbabilityDensityDistribution ProbDistResult(*this);
    ProbDistResult -= other;
    return ProbDistResult;
}


CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator*=(CProbabilityDensityDistribution& other)
{
    return _GeneralOperatorsFunctionAnalytical3(other, ProbDistOp::pdoMult );
}
CProbabilityDensityDistribution CProbabilityDensityDistribution::operator*(CProbabilityDensityDistribution& other)
{
    CProbabilityDensityDistribution ProbDistResult(*this);
    ProbDistResult *= other;
    return ProbDistResult;
}

CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator/=(CProbabilityDensityDistribution& other)
{
    return _GeneralOperatorsFunctionAnalytical3(other, ProbDistOp::pdoDiv);
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

void CProbabilityDensityDistribution::NormalDistribution(const CDigFloat& dfMu, const CDigFloat& dfSigma, const int nPoints /* = 50*/)
{
    LOGTRACE(LS_ProbDist+"NormalDistribution","called with args:");
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("int Points              =") + to_string(nPoints));
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("CDigFloat mean          =") + dfMu.RawPrint(10));
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("CDigFloat sigma         =") + dfSigma.RawPrint(10));
    
    // check consistent parameters:
    if(dfSigma <= 0 )
        return;
    
    // reset
    Reset();
    
    CDigFloat dfFac1 = 2*dfSigma*dfSigma;
    CDigFloat dfPreFac = 1/sqrt(dfFac1*M_PIf64);
    
    LOGTRACE(LS_ProbDist+"NormalDistribution","calculated:");
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("CDigFloat dfFac1 =") + dfFac1.RawPrint(10));
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("CDigFloat dfPreFac (amplitude) =") + dfPreFac.RawPrint(10));
    
    /////////////////////////////////////////////
    // calculation of the gaussian curve
    // keeping intervals between the x-values
    // constant
    /////////////////////////////////////////////
//     //set the starting point and increment depending on the min. value.
//     CDigFloat dfXDelta = sqrt(dfFac1*-log(dfYResolution*sqrt(dfFac1*M_PI)));
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
    nPtsEvenHalf += (nPtsEvenHalf % 2);
    
    // for npts we want to have npts+1 intervals: otherwise the last point is close to zero
    // which leads to a very broad x-intervall (y = 0 --> x = +- infinity)
    CDigFloat dfdY = dfPreFac / (nPtsEvenHalf+1);
    CDigFloat dfY = dfPreFac;
    CDigFloat dfRelPreFac = 1.;
    
    
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("int nPtsEvenHalf =") + to_string(nPtsEvenHalf));
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("int total points =") + to_string(2* nPtsEvenHalf + 1));
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("CDigFloat dfPreFac (amplitude) =") + dfPreFac.RawPrint(10));
    
    if( (2*nPtsEvenHalf + 1) > DISTRI_POINTS_WARNING_LIMIT)
        LOGWARN("WarningLogger", LS_ProbDist+"NormalDistribution: no. of points (" + to_string(2*nPtsEvenHalf+1) + ") exceeds warning limit (" + to_string(DISTRI_POINTS_WARNING_LIMIT) + ")!!");
    
    // add the first point: max at x = dfMu
    Add(dfMu,dfY);
    
    for(int ipt=0;ipt<nPtsEvenHalf; ipt++)
    {
        dfY -= dfdY;
        dfRelPreFac = dfY / dfPreFac;
        CDigFloat dfdXPos = sqrt(dfFac1*(-log(dfRelPreFac)));
        Add(dfMu + dfdXPos,dfY);
        Add(dfMu - dfdXPos,dfY);
        
        LOGTRACE(LS_ProbDist+"NormalDistribution",string("added point = (") + (dfMu-dfdXPos).RawPrint(10) + ", " + dfY.RawPrint(10) + ")" );
        LOGTRACE(LS_ProbDist+"NormalDistribution",string("added point = (") + (dfMu+dfdXPos).RawPrint(10) + ", " + dfY.RawPrint(10) + ")" );
        
    }
    
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("integral (should be 1) = ") + AbsIntegral().RawPrint(10));
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("************************************************************"));
    

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

CDigFloat CProbabilityDensityDistribution::Variance()
{
    LOGTRACE(LS_ProbDist+"Variance", "called");
    
    // copy from this
    CProbabilityDensityDistribution pdd4Calc(*this);
    
    
    LOGTRACE(LS_ProbDist+"Variance", string("original mean = ") + Mean().RawPrint(10));
    LOGTRACE(LS_ProbDist+"Variance", string("actual mean of copy (must be the same) = ") + pdd4Calc.Mean().RawPrint(10));
    
    
    // shift by mean: mean is now zero because of Integral (x-<x>) = 0
    pdd4Calc.Shift(-Mean());
    
    
    LOGTRACE(LS_ProbDist+"Variance", string("control: mean must be zero = ") + pdd4Calc.Mean().RawPrint(10));
    LOGTRACE(LS_ProbDist+"Variance", string("result = ") + pdd4Calc.Mean(2).RawPrint(10));
    
    LOGTRACE(LS_ProbDist+"Variance",string("************************************************************"));
    
    // return the second moment which is (x-<x>)² identical to x²
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

std::__cxx11::string CProbabilityDensityDistribution::PrintLinPars()
{
    ostringstream oss;
//     for(auto iel: mDistribution)
    for(MapDFDFType::iterator iel = mDistribution.begin(); iel != mDistribution.end(); iel++)
    {
        oss << "[" << iel->first.RawPrint(10,false) << " / " << iel->second.RawPrint(10,false) + " ] =";
        oss << "( " << m_LinPars[iel].first.RawPrint(10) << " / ";
        oss << m_LinPars[iel].second.RawPrint(10) << " )" << endl;
        
    }   // endfor(auto iel: Distribution())
    
    return oss.str();
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
    m_LinPars.clear();
    CDistribution::_Init();

}

void CProbabilityDensityDistribution::_SetLinPars()
{
    m_LinPars.clear();
    for(MapDFDFType::const_iterator iel = mDistribution.begin(); iel != mDistribution.end(); iel++)
//     for(auto iel: mDistribution)
        m_LinPars[iel] = PairDFType(_LinearOffset(iel),_LinearSlope(iel));
}


CProbabilityDensityDistribution & CProbabilityDensityDistribution:: _GeneralOperatorsFunctionAnalytical(CProbabilityDensityDistribution& Other, const ProbDistOp Operation)
{

    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("called with args:"));
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("probdist other:") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("operation     : ") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("probdist this :") + PrintMetaInfo());
    
    // generate plan : setting m_ConvolutionPlan;
    _SetConvolutionPlan(Other, Operation);
    
    // set all linear parameters for this and other
    _SetLinPars(); Other._SetLinPars();
    
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("lin par. this (") + Address2String(this)+ ") ...\n" + PrintLinPars() );
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("lin par. other (") + Address2String(&Other)+ ") ...\n"+ Other.PrintLinPars() );
    
    
    // declare the result distri
    MapDFDFType TargetDistri;
    
    for(auto itarget: m_ConvolutionPlan)
    {
        TargetDistri[itarget.first] = 0;
        
        LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("init TargetDistri[" + itarget.first.RawPrint(10) + "] = 0"));
        
        // iterate over all sub integration limits
        for(auto intlim: itarget.second)
        {
            // DEBUG
//             cout << "x = [ " << intlim.first.first.RawPrint(10) << " ; " << intlim.first.second.RawPrint(10) << " ] " << endl;
//             cout << "y = [ " << intlim.second.first.RawPrint(10) << " ; " << intlim.second.second.RawPrint(10) << " ] " << endl;
            
            // add to all other probs of this target value
            TargetDistri[itarget.first] += _GetConvolutionIntegralAnalytical4Operation(Other, Operation, itarget.first, intlim);
            
            // logging trace 
            LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("incremented TargetDistri[" + itarget.first.RawPrint(10) + "] = " + TargetDistri[itarget.first].RawPrint(10)));
        
            
        }   // endfor(auto interval: itarget.second)        
        
    }   // endfor(auto itarget: m_ConvolutionPlan)
    
    // remember some things before init resets all the values 
    int nMaxBinIter = MaxBinarySearchIterations();
    int nIntegrationSteps = IntegrationSteps();
    int nSubIntegrationSteps = SubIntegrationSteps();
    
    // now init: especially the integration step for convolution get lost
    _Init();
    
    // reset to old values 
    MaxBinarySearchIterations(nMaxBinIter);
    IntegrationSteps(nIntegrationSteps);
    SubIntegrationSteps(nSubIntegrationSteps);
    
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
    
    
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("************************************************************"));
    
    return *this;
 
}

CProbabilityDensityDistribution & CProbabilityDensityDistribution::_GeneralOperatorsFunctionAnalytical3(CProbabilityDensityDistribution& Other, const ProbDistOp Operation)
{   
    
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("called with args:"));
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("probdist other:") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("operation     : ") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("probdist this :") + PrintMetaInfo());
    

    // generate plan : setting x-intervals and target values
    VectorPairDFType vpdfXIntervals;
    vector<CDigFloat> vdfTargetValues;
    _PrepareConvolution(Other, Operation, vdfTargetValues, vpdfXIntervals);
    
    // set all linear parameters for this and other
    _SetLinPars(); Other._SetLinPars();
    
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("lin par. this (") + Address2String(this)+ ") ...\n" + PrintLinPars() );
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("lin par. other (") + Address2String(&Other)+ ") ...\n"+ Other.PrintLinPars() );
       
    // prepare threading: need a vector keeping the convolution threads 
    map<unsigned int, std::thread*> vpConvolutionThreads;
    
    // declare the result distri
    MapDFDFType TargetDistri;
    unsigned int uiThreadID = 0;
    
    // need vector of vector of target values 
    vector< vector<CDigFloat>*> vpvdfTargetValuePackages;
    vector< vector<CDigFloat>*> vpvdfResultPackages;
    vector<void*> vpFunctions;
    vector < ProbDistOp*> vpOperationCopies;
    vector < VectorPairDFType*> vpvpdfXIntervalCopies;
    vector< CProbabilityDensityDistribution*> vpdThisCopies;
    vector< CProbabilityDensityDistribution*> vpdOtherCopies;
    unsigned int uiPackageCounter = 0;
    unsigned int uiNoOfThreads = IntegrationSteps()*(Distribution().size()+Other.Distribution().size()) / PDD_SUB_INTERVALL_CALCULATIONS_PER_THREAD;
    if(uiNoOfThreads > vdfTargetValues.size())
        uiNoOfThreads = vdfTargetValues.size();
    if(uiNoOfThreads == 0)
        uiNoOfThreads = 1;
    unsigned int uiTargetValuesPerThread = vdfTargetValues.size() / uiNoOfThreads;
    
    
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("target values per thread = ") + to_string(uiTargetValuesPerThread));
    
    
    for(auto itarget: vdfTargetValues)
    {
        // handle init of new package:
        if(uiPackageCounter == 0)
        {
            vpvdfTargetValuePackages.push_back(new vector<CDigFloat>);
            vpvdfResultPackages.push_back(new vector<CDigFloat>);
            
        }   // endif(uiPackageCounter == 0)
        
        TargetDistri[itarget] = 0;
        
        // hand over target value and result to package vectors 
        vpvdfTargetValuePackages.back()->push_back(itarget);
//         vpvpdfResultPointerPackages.back()->push_back(&(TargetDistri[itarget]));
        vpvdfResultPackages.back()->push_back(0);
        
        // increment package counter 
        uiPackageCounter++;
        
        // handle finished package
        if( uiPackageCounter == uiTargetValuesPerThread || itarget == vdfTargetValues.back() )
        {
            vpdThisCopies.push_back( new CProbabilityDensityDistribution( *this));
            vpdThisCopies.back()->_SetLinPars();
            vpdOtherCopies.push_back(new CProbabilityDensityDistribution( Other));
            vpdOtherCopies.back()->_SetLinPars();
            
            // copy x intervals, too 
            vpvpdfXIntervalCopies.push_back(new VectorPairDFType);
            vpvpdfXIntervalCopies.back()->operator=(vpdfXIntervals);
            
            // copy operators
            vpOperationCopies.push_back(new ProbDistOp(Operation));
            
            // threading call: generate thread objects
            auto f = std::bind(&CProbabilityDensityDistribution::_Convolution4TargetValueThread, vpdThisCopies.back(), vpdOtherCopies.back(), vpOperationCopies.back(), vpvpdfXIntervalCopies.back(), vpvdfTargetValuePackages.back(), vpvdfResultPackages.back());
            
            LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("starting thread ") + to_string(uiThreadID) );
            vpConvolutionThreads[uiThreadID] = new std::thread(f, vpdOtherCopies.back(), vpOperationCopies.back(), vpvpdfXIntervalCopies.back(), vpvdfTargetValuePackages.back(), vpvdfResultPackages.back());
            uiThreadID++;
            
            // resetting package counter 
            uiPackageCounter = 0;
            
        }   // endif( uiPackageCounter == PDD_TARGET_VALUES_PER_THREAD)
        
//         vpConvolutionThreads.push_back(std::thread( CConvolutionThread(),this, &Other, Operation, &vpdfXIntervals, itarget, &(TargetDistri[itarget])));
        // call by member function:
//         vpConvolutionThreads.push_back(std::thread( &CProbabilityDensityDistribution::_Convolution4TargetValue, this, Other, Operation, vpdfXIntervals, itarget, TargetDistri[itarget]));
        
        // non-threading call: simply call internal function
//         _Convolution4TargetValue(Other, Operation, vpdfXIntervals, itarget, TargetDistri[itarget] );
        
    }   // endfor(auto itarget: vdfTargetValues)
    
    // now wait for threads to finish
    for(auto ith: vpConvolutionThreads)
    {
        LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("waiting for thread ") + to_string(ith.first) + " to finish" );
        ith.second->join();
    }
    // now delete threads 
    for(auto ith: vpConvolutionThreads)
    {
        LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("handing over result for thread ") + to_string(ith.first) + "... ");
           
        for(int itv = 0; itv < vpvdfTargetValuePackages[ith.first]->size(); itv++)
        {
            LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("... [") + (*vpvdfTargetValuePackages[ith.first])[itv].RawPrint(10,false) + "] = " + (*vpvdfResultPackages[ith.first])[itv].RawPrint(10)   );
            TargetDistri[(*vpvdfTargetValuePackages[ith.first])[itv]] = (*vpvdfResultPackages[ith.first])[itv];
           
        }   //endfor(int itv = 0; itv < vpvdfTargetValuePackages[ith.first]->size(); itv++)
        LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("deleting thread ") + to_string(ith.first) );
        SecureDeleteObjectPointer( ith.second );
        vpvdfTargetValuePackages[ith.first]->clear();
        vpvdfResultPackages[ith.first]->clear();
        SecureDeleteObjectPointer( vpvdfTargetValuePackages[ith.first]);
        SecureDeleteObjectPointer( vpvdfResultPackages[ith.first]);
        SecureDeleteObjectPointer( vpdThisCopies[ith.first]);
        SecureDeleteObjectPointer( vpdOtherCopies[ith.first]);
    }
    
    
    // remember some things before init resets all the values 
    int nMaxBinIter = MaxBinarySearchIterations();
    int nIntegrationSteps = IntegrationSteps();
    int nSubIntegrationSteps = SubIntegrationSteps();
    
    // now init: especially the integration step for convolution get lost
    _Init();
    
    // reset to old values 
    MaxBinarySearchIterations(nMaxBinIter);
    IntegrationSteps(nIntegrationSteps);
    SubIntegrationSteps(nSubIntegrationSteps);
    
    if(vdfTargetValues.size() > 0)
    {
        for(auto iel: TargetDistri)
            mDistribution[iel.first] = iel.second;
        
        bNormalized=false;
        Normalize();
        
        assert(AbsIntegral() == 1);
        
    }   // endif(!(TargetRangeStart == TargetRangeEnd && TargetRangeStart == 0))
    
    
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("************************************************************"));
    
    return *this;
 
 
}
void CProbabilityDensityDistribution::_Convolution4TargetValueThread(CProbabilityDensityDistribution* pOther, ProbDistOp *pOperation, VectorPairDFType* pvpdfXIntervals, vector<CDigFloat>* pvdfTargetValues, vector<CDigFloat>* pvdfResult)
{
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValueThread", string("called with args:") );
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValueThread", string("this") + Address2String(this) );
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValueThread", string("other") + Address2String(pOther) );
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValueThread", string("operation ") + GetProbDistOpAsString(*pOperation) );
        for(auto x: *pvpdfXIntervals)
            LOGTRACE(LS_ProbDist+"_Convolution4TargetValueThread", string("x in [") + x.first.RawPrint(10) + ", " + x.second.RawPrint(10) + "]" );
        for(auto tv: *pvdfTargetValues)
            LOGTRACE(LS_ProbDist+"_Convolution4TargetValueThread", string("tv = ") + tv.RawPrint(10) );
        
        
        for(int itv=0; itv < pvdfTargetValues->size(); itv++)
        _Convolution4TargetValue(*pOther, *pOperation, (*pvpdfXIntervals), (*pvdfTargetValues)[itv], (*pvdfResult)[itv]);
}
 


void CProbabilityDensityDistribution::_Convolution4TargetValue(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, const VectorPairDFType& vpdfXIntervals, const CDigFloat& dfTargetValue, CDigFloat& dfResult)
{
    
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("called with args:"));
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("operation ") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("dfTargetValue ") + dfTargetValue.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("other distri ") + Other.PrintMetaInfo());
    
    /////////////////////////////////////////////////////
    // General purpose ot this function:
    //   Determine a valid total integration interval and
    //   its subintervals
    //
    // Procedure:
    //
    //   I  derive closed intervals for x and y mapped 
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
    VectorConvPlanElementType SubIntervals4TargetValue;
    
    // keeps the total integration interval(s) for a target value
    VectorSubIntLimitsType TotalIntervals4TargetValue;
    
        
    // I 2. derive closed intervals for x and y (kept in TotalIntervals4TargetValue) mapped onto each other bijectively
    //      by the calculation of the complementary variable
//     _GetTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, vpdXIntervals, TotalIntervals4TargetValue);
    _GetTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, vpdfXIntervals, TotalIntervals4TargetValue);
    
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("calculated complementary intervals with mixed slope and offset: "));
    for(auto ix: TotalIntervals4TargetValue)
    {
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("x in [") + ix.first.first.RawPrint(10) + ", " + ix.first.second.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("y in [") + ix.second.first.RawPrint(10) + ", " + ix.second.second.RawPrint(10));
    }
    
    
    _GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue3(Operation, dfTargetValue, Other, TotalIntervals4TargetValue, SubIntervals4TargetValue);  
    
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("calculated complementary intervals with unique slope and offset ... "));
    for(auto planEl: SubIntervals4TargetValue)
    {
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... x iterator points to [") + planEl.xDistriIter->first.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... x in [") + planEl.xLimits.first.RawPrint(10) + ", " + planEl.xLimits.second.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... y iterator points to [") + planEl.yDistriIter->first.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... y in [") + planEl.yLimits.first.RawPrint(10) + ", " + planEl.yLimits.second.RawPrint(10));
    }
    
    
    // init the result
    dfResult = 0;
    
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("iterating over all sub intervals ...") );
    LOGINFO(LS_ProbDist+"_Convolution4TargetValue",string("no. of sub intervals: ") + to_string(SubIntervals4TargetValue.size()) );

    // now calculate 
    for(auto planEl: SubIntervals4TargetValue)
    {   
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... x iterator points to ")  + planEl.xDistriIter->first.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... x in [") + planEl.xLimits.first.RawPrint(10) + ", " + planEl.xLimits.second.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... y iterator points to ") + planEl.yDistriIter->first.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... y in [") + planEl.yLimits.first.RawPrint(10) + ", " + planEl.yLimits.second.RawPrint(10));

        switch( Operation)
        {
            case ProbDistOp::pdoPlus:
                    dfResult +=  _Integral4TargetValue4Addition2(dfTargetValue,Other, planEl);            
                break;
            case ProbDistOp::pdoMinus:
                    dfResult +=  _Integral4TargetValue4Subtraction2(dfTargetValue, Other, planEl);            
                break;
            case ProbDistOp::pdoMult:
                    dfResult +=  _Integral4TargetValue4Multiplication2(dfTargetValue, Other, planEl);            
                break;
            case ProbDistOp::pdoDiv:
                    dfResult +=  _Integral4TargetValue4Division2(dfTargetValue, Other, planEl);            
                break;
            default:
                break;
        }   //endswitch( Operation)
        
        
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... dfResult incremented ") + dfResult.RawPrint(10));
        
        
    }   // endfor(auto interval: SubIntervals4TargetValue)

}


CProbabilityDensityDistribution & CProbabilityDensityDistribution::_GeneralOperatorsFunctionAnalytical2(CProbabilityDensityDistribution& Other, const ProbDistOp Operation)
{

    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("called with args:"));
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("probdist other:") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("operation     : ") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("probdist this :") + PrintMetaInfo());
    
    // generate plan : setting m_ConvolutionPlan;
    _SetConvolutionPlan2(Other, Operation);
    
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("lin par. this ... \n") + PrintLinPars() );
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("lin par. other ... \n") + Other.PrintLinPars() );
    
    
    // declare the result distri
    MapDFDFType TargetDistri;
    
    for(auto itarget: m_ConvolutionPlan2)
    {
        TargetDistri[itarget.first] = 0;
        
        LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("init TargetDistri[" + itarget.first.RawPrint(10) + "] = 0"));
        
        // iterate over all sub integration limits: these are kept in convolution plan elements
        for(auto planEl: itarget.second)
        {
            
            // add to all other probs of this target value
            TargetDistri[itarget.first] += _GetConvolutionIntegralAnalytical4Operation2(Other, Operation, itarget.first, planEl);
            
            // logging trace 
            LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("incremented TargetDistri[" + itarget.first.RawPrint(10) + "] = " + TargetDistri[itarget.first].RawPrint(10)));
        
            
        }   // endfor(auto interval: itarget.second)        
        
    }   // endfor(auto itarget: m_ConvolutionPlan)
    
    // remember some things before init resets all the values 
    int nMaxBinIter = MaxBinarySearchIterations();
    int nIntegrationSteps = IntegrationSteps();
    int nSubIntegrationSteps = SubIntegrationSteps();
    
    // now init: especially the integration step for convolution get lost
    _Init();
    
    // reset to old values 
    MaxBinarySearchIterations(nMaxBinIter);
    IntegrationSteps(nIntegrationSteps);
    SubIntegrationSteps(nSubIntegrationSteps);
    
    if(m_ConvolutionPlan2.size() > 0)
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
    
    
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("************************************************************"));
    
    return *this;
 
}



CDigFloat CProbabilityDensityDistribution::_GetConvolutionIntegralAnalytical4Operation(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, const CDigFloat& dfTargetValue, const SubIntLimitsType& vXYLimits)
{   
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("called with args:"));
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("other distri:") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("this distri:") + PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("operation   :") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("target value:") + dfTargetValue.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("x-limit: [") + vXYLimits.first.first.RawPrint(10) + ", " +vXYLimits.first.second.RawPrint(10) + " ]" );
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("y-limit: [") + vXYLimits.second.first.RawPrint(10) + ", " +vXYLimits.second.second.RawPrint(10) + " ]" );
    
    // init the result
    CDigFloat dfResult = 0;
    
    // derive all necessary arguments needed for the analytical calculation of the 
    // convolution integrals
    // calculate slope and offset for the related interval of this distri
    CDigFloat dfXMean = (vXYLimits.first.first + vXYLimits.first.second)/2.;
    MapDFDFType::const_iterator Left, Right;
    GetInterval(dfXMean, Left, Right);
    CDigFloat dfXOffset = _LinearOffset(Left);
    CDigFloat dfXSlope = _LinearSlope(Left);
    
    // calculate slope and offset for the related interval of other distri
    MapDFDFType::const_iterator LeftOther, RightOther;
    CDigFloat dfYMean = (vXYLimits.second.first + vXYLimits.second.second)/2.;
    Other.GetInterval(dfYMean, LeftOther, RightOther);
    CDigFloat dfYOffset = Other._LinearOffset(LeftOther);
    CDigFloat dfYSlope = Other._LinearSlope(LeftOther);
    
    // trace logging 
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("preparation of convolution integral calculation"));
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("X integration interval= [") + vXYLimits.first.first.RawPrint(10) + " ; " + vXYLimits.first.second.RawPrint(10) + " ]");
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("enclosing X-distri interval = [") + Left->first.RawPrint(10) + " ; " + Right->first.RawPrint(10) + " ]");
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfXMean = ") + dfXMean.RawPrint(10));
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfXOffset = ") + dfXOffset.RawPrint(10));
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfXSlope = ") + dfXSlope.RawPrint(10) );
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("Y integration interval= [") + vXYLimits.second.first.RawPrint(10) + " ; " + vXYLimits.second.second.RawPrint(10) + " ]");
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("enclosing Y-distri interval = [") + LeftOther->first.RawPrint(10) + " ; " + RightOther->first.RawPrint(10) + " ]");
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfYMean = ") + dfYMean.RawPrint(10));
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfYOffset = ") + dfYOffset.RawPrint(10));
      LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("dfYSlope = ") + dfYSlope.RawPrint(10));
   
     
    // variables for primitive integral values of upper and lower limit
    CDigFloat dfPrimIntValLow, dfPrimIntValUp;
    
    switch( Operation)
    {
        case ProbDistOp::pdoPlus:
//             dfPrimIntValUp = _PrimitiveIntegral4TargetValue4Addition(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.second);
//             dfPrimIntValLow = _PrimitiveIntegral4TargetValue4Addition(dfXOffset, dfXSlope, dfYOffset,dfYSlope,dfTargetValue, vXYLimits.first.first);
                dfPrimIntValUp = _Integral4TargetValue4Addition(dfXOffset,dfXSlope, dfYOffset, dfYSlope, dfTargetValue,vXYLimits.first.first, vXYLimits.first.second);
                dfPrimIntValLow = 0;
            
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
            
        LOGERROR("ErrorLogger",string("preparation of convolution integral calculation"));
        LOGERROR("ErrorLogger",string("X integration interval= [") + vXYLimits.first.first.RawPrint(10) + " ; " + vXYLimits.first.second.RawPrint(10) + " ]");
        LOGERROR("ErrorLogger",string("enclosing X-distri interval = [") + Left->first.RawPrint(10) + " ; " + Right->first.RawPrint(10) + " ]");
        LOGERROR("ErrorLogger",string("dfXMean = ") + dfXMean.RawPrint(10));
        LOGERROR("ErrorLogger",string("dfXOffset = ") + dfXOffset.RawPrint(10));
        LOGERROR("ErrorLogger",string("dfXSlope = ") + dfXSlope.RawPrint(10) );
        LOGERROR("ErrorLogger",string("Y integration interval= [") + vXYLimits.second.first.RawPrint(10) + " ; " + vXYLimits.second.second.RawPrint(10) + " ]");
        LOGERROR("ErrorLogger",string("enclosing Y-distri interval = [") + LeftOther->first.RawPrint(10) + " ; " + RightOther->first.RawPrint(10) + " ]");
        LOGERROR("ErrorLogger",string("dfYMean = ") + dfYMean.RawPrint(10));
        LOGERROR("ErrorLogger",string("dfYOffset = ") + dfYOffset.RawPrint(10));
        LOGERROR("ErrorLogger",string("dfYSlope = ") + dfYSlope.RawPrint(10));
   
    }
    
    
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("************************************************************"));
    assert(dfResult >= 0);
    return dfResult;
    
}
CDigFloat CProbabilityDensityDistribution::_GetConvolutionIntegralAnalytical4Operation2(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, const CDigFloat& dfTargetValue, const ConvPlanElement& planEl)
{
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("called with args:"));
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("other distri:") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("this distri:") + PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("operation   :") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("target value:") + dfTargetValue.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("first value of x-interval: ") + planEl.xDistriIter->first.RawPrint(10) );
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("x-limit: [") + planEl.xLimits.first.RawPrint(10) + ", " + planEl.xLimits.second.RawPrint(10) + " ]" );
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("first value of y-interval: ") + planEl.yDistriIter->first.RawPrint(10) );
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("y-limit: [") + planEl.yLimits.first.RawPrint(10) + ", " + planEl.yLimits.second.RawPrint(10) + " ]" );
    
    // init the result
    CDigFloat dfResult = 0;
    
    switch( Operation)
    {
        case ProbDistOp::pdoPlus:
                dfResult =  _Integral4TargetValue4Addition2(dfTargetValue,Other, planEl);            
            break;
        case ProbDistOp::pdoMinus:
                dfResult =  _Integral4TargetValue4Subtraction2(dfTargetValue,Other, planEl);            
            break;
        case ProbDistOp::pdoMult:
                dfResult =  _Integral4TargetValue4Multiplication2(dfTargetValue,Other, planEl);            
            break;
        case ProbDistOp::pdoDiv:
                dfResult =  _Integral4TargetValue4Division2(dfTargetValue,Other, planEl);            
            break;
        default:
            break;
    }   //endswitch( Operation)
    
    
    // trace 
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation", "convolution integral is :");
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("result = ") + dfResult.RawPrint(10));
    
    // error case
    if(dfResult < 0)
    {
        LOGERROR("ErrorLogger","CProbabilityDensityDistribution::_GetConvolutionIntegralAnalytical4Operation: convolution integral is negative :");
        LOGERROR("ErrorLogger",string("result = ") + dfResult.RawPrint(10));
        LOGERROR("ErrorLogger",string("operator = ") + GetProbDistOpAsString(Operation));
        LOGERROR("ErrorLogger",string("lin. pars. this ... \n") + PrintLinPars());
        LOGERROR("ErrorLogger",string("lin. pars. other ... \n") + Other.PrintLinPars());
            
        LOGERROR("ErrorLogger",string("preparation of convolution integral calculation"));
        LOGERROR("ErrorLogger",string("X integration interval= [") + planEl.xLimits.first.RawPrint(10) + " ; " + planEl.xLimits.second.RawPrint(10) + " ]");
        LOGERROR("ErrorLogger",string("enclosing X-distri interval first value = ") + planEl.xDistriIter->first.RawPrint(10));
        LOGERROR("ErrorLogger",string("dfXOffset = ") + Offset(planEl.xDistriIter).RawPrint(10));
        LOGERROR("ErrorLogger",string("dfXSlope = ") + Slope(planEl.xDistriIter).RawPrint(10) );
        LOGERROR("ErrorLogger",string("Y integration interval= [") + planEl.yLimits.first.RawPrint(10) + " ; " + planEl.yLimits.second.RawPrint(10) + " ]");
        LOGERROR("ErrorLogger",string("enclosing Y-distri interval first value = ") + planEl.yDistriIter->first.RawPrint(10));
        LOGERROR("ErrorLogger",string("dfYOffset = ") + Offset(planEl.yDistriIter).RawPrint(10));
        LOGERROR("ErrorLogger",string("dfYSlope = ") + Slope(planEl.yDistriIter).RawPrint(10) );
   
    }
    
    
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("************************************************************"));
    assert(dfResult >= 0);
    return dfResult;
    
}

void CProbabilityDensityDistribution::_PrepareConvolution(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, vector<CDigFloat>& vdfTargetValues, VectorPairDFType& vpdfXIntervals)
{
    LOGTRACE(LS_ProbDist+"_PrepareConvolution", string("called with args:"));
    LOGTRACE(LS_ProbDist+"_PrepareConvolution", string("probdist other :") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_PrepareConvolution", string("operation      :") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_PrepareConvolution", string("probdist this  :") + PrintMetaInfo());
    
    //////////////////////////////
    // x-intervals
    //   must be unipolar (positive or negative) for 
    //   multiplication and division
    /////////////////////////////
       
    
    // x-interval(s) for calculation: must be unipolar
    // for division and multiplication
    _GetXInterval4Operation(Operation, vpdfXIntervals);
    
    LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("calculated unipolar x-intervals: "));
    for(auto ix: vpdfXIntervals)
        LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("x in [") + ix.first.RawPrint(10) + ", " + ix.second.RawPrint(10));
        
    //////////////////////////////
    // target values:
    // avoid zero for multiplication
    // and division
    /////////////////////////////
       
    // clear and set target values 
    vdfTargetValues.clear();
    
    CDigFloat dfTargetStart, dfTargetEnd;
    _GetRangeFromDistributionOperation(Other, Operation, dfTargetStart, dfTargetEnd);
    
    
    // get the target step
    CDigFloat dfTargetStep = (dfTargetEnd - dfTargetStart)/(CDigFloat(IntegrationSteps()));
    for(int itv = 0; itv < IntegrationSteps(); itv++)
        vdfTargetValues.push_back(dfTargetStart+itv*dfTargetStep);
    
    // take care to avoid zero in in case of multiplication and division
    if( Operation == ProbDistOp::pdoDiv || Operation == ProbDistOp::pdoMult)
    {
        // do we have a zero ?
        vector<CDigFloat>::iterator itTv = find( vdfTargetValues.begin(), vdfTargetValues.end(), 0);
        if( itTv != vdfTargetValues.end())
        {
            // if it is the start of the distri ...
            if(itTv == vdfTargetValues.begin())
                // ... add half a step: avoid bipolarity for this case
                *itTv = (dfTargetStep/2.);
            else
                // ... subtract half a step
                *itTv -= (dfTargetStep/2.);
        
        }   // endif( itTv != vdfTargetValues.end())
    
    }   // endif( Operation == ProbDistOp::pdoDiv || Operation == ProbDistOp::pdoMult)
    
    
    // trace 
    LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("calculated:"));    
    LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("target range = [") + dfTargetStart.RawPrint(10) + ", " + dfTargetEnd.RawPrint(10) + " ]");
    LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("step size = ") + dfTargetStep.RawPrint(10) + "( " + to_string(IntegrationSteps()) + " ) ");
    LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("target values ...") );
    for(auto iel: vdfTargetValues)
            LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("... ") + iel.RawPrint(10) );
    
    LOGTRACE(LS_ProbDist+"_PrepareConvolution", string("********************************************************************"));
}


void CProbabilityDensityDistribution::_SetConvolutionPlan2(CProbabilityDensityDistribution& Other, const ProbDistOp Operation)
{  
    
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan", string("called with args:"));
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan", string("probdist other :") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan", string("operation      :") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan", string("probdist this  :") + PrintMetaInfo());
 
    
    CDigFloat dfTargetStart, dfTargetEnd;
    _GetRangeFromDistributionOperation(Other, Operation, dfTargetStart, dfTargetEnd);
    
    
    // get the target step
    CDigFloat dfTargetStep = (dfTargetEnd - dfTargetStart)/(CDigFloat(IntegrationSteps()));
    
    // trace 
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("calculated:"));    
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("target range = [") + dfTargetStart.RawPrint(10) + ", " + dfTargetEnd.RawPrint(10) + " ]");
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("step size = ") + dfTargetStep.RawPrint(10) + "( " + to_string(IntegrationSteps()) + " ) ");
    
  
    // clear the plan first
    m_ConvolutionPlan2.clear();
    
    // set all linear parameters for this and other
    _SetLinPars(); Other._SetLinPars();
    
    // now start building the plan 
    CDigFloat dfActTargetValue = dfTargetStart;
    bool bZeroTargetValueDetected = false;
    bool bCheck4ZeroTargetValue = ( (Operation == ProbDistOp::pdoDiv ) || ( Operation== ProbDistOp::pdoMult));
    while(dfActTargetValue < (dfTargetEnd +dfTargetStep))
    {

        LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("iteration dfActTargetValue = ") + dfActTargetValue.RawPrint(10));
        
        // defining the complete integration range for this and other distribution
        _SetSubIntegrationIntervals4TargetValue2(Operation, dfActTargetValue, Other);        
        
        
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
        
        LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue set to ") + dfActTargetValue.RawPrint(10));
                
        // split zero dftarget value due to unipolarity trouble for multiplication and division
        if(bCheck4ZeroTargetValue)
            if(dfActTargetValue == 0 && dfActTargetValue < (dfTargetEnd +dfTargetStep))
            {
                LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue zero detected:"));
                
                // decrement half a step 
                dfActTargetValue -= dfTargetStep/2.;
        
                LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue decremented half a step to ") + dfActTargetValue.RawPrint(10));    
                
                // defining the complete integration range for this and other distribution ... again
                _SetSubIntegrationIntervals4TargetValue2(Operation, dfActTargetValue, Other); 
                
                // increment about one step ... now we are half a step behind because we split the target value == 0
                dfActTargetValue += dfTargetStep;
                
                LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue incremented one step to ") + dfActTargetValue.RawPrint(10));  
                
                // remember we found a zero and deal with it the next round 
                bZeroTargetValueDetected = true;
                
                
            }   // endif(dfTargetValue.RawValue() == 0)
        
    }   // endwhile(dfActTargetValue < (dfTargetEnd +dfTargetStep))
    
    
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("************************************************************"));
    
}

 void CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue2(const ProbDistOp Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other)
{
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("called with args:"));
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("operation ") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("dfTargetValue ") + dfTargetValue.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("other distri ") + Other.PrintMetaInfo());
    
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
    VectorConvPlanElementType SubIntervals4TargetValue;
    
    // keeps the total integration interval(s) for a target value
    VectorSubIntLimitsType TotalIntervals4TargetValue;
    
    // TODO: vpdXIntervals does not depend on a target value: 
    // 1. calculate before calling _SetSubIntegrationIntervals4TargetValue2
    // 2. add as argument
    
    // fill the x intervals: must be unipolar for multiplication and division
    VectorPairDFType vpdXIntervals;
    
    // DEBUG
//     cout << "_SetIntegrationIntervals4TargetValue" << endl;
//     cout << "Operation = " <<  to_string((int)Operation ) << endl;
//     cout << "dfTargetValue = " << dfTargetValue.RawPrint(30) << endl;
    
    // I 1. get the x-interval(s) for calculation: must be unipolar
    //      for division and multiplication
    _GetXInterval4Operation(Operation, vpdXIntervals);
    
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("calculated unipolar x-intervals: "));
    for(auto ix: vpdXIntervals)
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("x in [") + ix.first.RawPrint(10) + ", " + ix.second.RawPrint(10));
        
    // I 2. derive closed intervals for x and y (kept in TotalIntervals4TargetValue) mapped onto each other bijectively
    //      by the calculation of the complementary variable
//     _GetTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, vpdXIntervals, TotalIntervals4TargetValue);
    _GetTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, vpdXIntervals, TotalIntervals4TargetValue);
    
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("calculated complementary intervals with mixed slope and offset: "));
    for(auto ix: TotalIntervals4TargetValue)
    {
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("x in [") + ix.first.first.RawPrint(10) + ", " + ix.first.second.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("y in [") + ix.second.first.RawPrint(10) + ", " + ix.second.second.RawPrint(10));
    }
    
    
    _GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue3(Operation, dfTargetValue, Other, TotalIntervals4TargetValue, SubIntervals4TargetValue);  
    
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("calculated complementary intervals with unique slope and offset ... "));
    for(auto ix: SubIntervals4TargetValue)
    {
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("... x iterator points to [") + ix.xDistriIter->first.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("... x in [") + ix.xLimits.first.RawPrint(10) + ", " + ix.xLimits.second.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("... y iterator points to [") + ix.yDistriIter->first.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("... y in [") + ix.yLimits.first.RawPrint(10) + ", " + ix.yLimits.second.RawPrint(10));
    }
    
    // now add the sub integration limits and its corresponding target value as pair
    m_ConvolutionPlan2.push_back(pair<CDigFloat, VectorConvPlanElementType>(dfTargetValue, SubIntervals4TargetValue) ); 
}


void CProbabilityDensityDistribution::_SetConvolutionPlan(CProbabilityDensityDistribution& Other, const ProbDistOp Operation)
{  
    
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan", string("called with args:"));
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan", string("probdist other :") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan", string("operation      :") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan", string("probdist this  :") + PrintMetaInfo());
 
    
    CDigFloat dfTargetStart, dfTargetEnd;
    _GetRangeFromDistributionOperation(Other, Operation, dfTargetStart, dfTargetEnd);
    
    
    // get the target step
    CDigFloat dfTargetStep = (dfTargetEnd - dfTargetStart)/(CDigFloat(IntegrationSteps()));
    
    // trace 
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("calculated:"));    
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("target range = [") + dfTargetStart.RawPrint(10) + ", " + dfTargetEnd.RawPrint(10) + " ]");
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("step size = ") + dfTargetStep.RawPrint(10) + "( " + to_string(IntegrationSteps()) + " ) ");
    
  
    // clear the plan first
    m_ConvolutionPlan.clear();
    
    // now start building the plan 
    CDigFloat dfActTargetValue = dfTargetStart;
    bool bZeroTargetValueDetected = false;
    bool bCheck4ZeroTargetValue = ( (Operation == ProbDistOp::pdoDiv ) || ( Operation== ProbDistOp::pdoMult));
    while(dfActTargetValue < (dfTargetEnd +dfTargetStep))
    {

        LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("iteration dfActTargetValue = ") + dfActTargetValue.RawPrint(10));
        
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
        
        LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue set to ") + dfActTargetValue.RawPrint(10));
                
        // split zero dftarget value due to unipolarity trouble for multiplication and division
        if(bCheck4ZeroTargetValue)
            if(dfActTargetValue == 0 && dfActTargetValue < (dfTargetEnd +dfTargetStep))
            {
                LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue zero detected:"));
                
                // decrement half a step 
                dfActTargetValue -= dfTargetStep/2.;
        
                LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue decremented half a step to ") + dfActTargetValue.RawPrint(10));    
                
                // defining the complete integration range for this and other distribution ... again
                _SetSubIntegrationIntervals4TargetValue(Operation, dfActTargetValue, Other); 
                
                // increment about one step ... now we are half a step behind because we split the target value == 0
                dfActTargetValue += dfTargetStep;
                
                LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("dfActTargetValue incremented one step to ") + dfActTargetValue.RawPrint(10));  
                
                // remember we found a zero and deal with it the next round 
                bZeroTargetValueDetected = true;
                
                
            }   // endif(dfTargetValue.RawValue() == 0)
        
    }   // endwhile(dfActTargetValue < (dfTargetEnd +dfTargetStep))
    
    
    LOGTRACE(LS_ProbDist+"_SetConvolutionPlan",string("************************************************************"));
    
}
CDigFloat CProbabilityDensityDistribution::_Integral4TargetValue4Addition2(const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const ConvPlanElement& planEl)
{    
    CDigFloat dfX1 = planEl.xLimits.first;
    CDigFloat dfX2 = planEl.xLimits.second;
 
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", " called with args:");
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("target value =") + dfTargetValue.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition",string("x iterator points to [") + planEl.xDistriIter->first.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition",string("x in [") + dfX1.RawPrint(10) + ", " + dfX2.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition",string("y iterator points to [") + planEl.yDistriIter->first.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition",string("y in [") + planEl.yLimits.first.RawPrint(10) + ", " + planEl.yLimits.second.RawPrint(10));


    CDigFloat dfOffsetX =  Offset(planEl.xDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("XOffset :") + dfOffsetX.RawPrint(20) );
    CDigFloat dfSlopeX =  Slope(planEl.xDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("XSlope  :") + dfSlopeX.RawPrint(20) );
    CDigFloat dfOffsetY =  Other.Offset(planEl.yDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("YOffset :") + dfOffsetY.RawPrint(20) );
    CDigFloat dfSlopeY =  Other.Slope(planEl.yDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("YSlope  :") + dfSlopeY.RawPrint(20) );
    CDigFloat dfYConst = dfSlopeY*dfTargetValue+dfOffsetY;
    CDigFloat dfDiff1 = dfX2 -dfX1;
    CDigFloat dfDiff2 = (dfX1+dfX2) * dfDiff1;
    CDigFloat dfDiff3 = dfDiff1*(pow(dfX1,2) + pow(dfX2,2) +dfX1*dfX2);

    
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", "**************************************************************************");    
    
    return -dfSlopeX*dfSlopeY/3*dfDiff3+
            (dfSlopeX*dfYConst-dfSlopeY*dfOffsetX)/2.* dfDiff2+
            dfYConst*dfOffsetX*dfDiff1;
}

CDigFloat CProbabilityDensityDistribution::_Integral4TargetValue4Addition(const CDigFloat& dfOffsetX, const CDigFloat& dfSlopeX, const CDigFloat& dfOffsetY, const CDigFloat& dfSlopeY, const CDigFloat& dfTargetValue, const CDigFloat& dfX1, const CDigFloat& dfX2)
{
    CDigFloat dfYConst = dfSlopeY*dfTargetValue+dfOffsetY;
    CDigFloat dfDiff1 = dfX2 -dfX1;
    CDigFloat dfDiff2 = (dfX1+dfX2) * dfDiff1;
    CDigFloat dfDiff3 = dfDiff1*(pow(dfX1,2) + pow(dfX2,2) +dfX1*dfX2);
    CDigFloat dfDiff1Comp = dfX2 -dfX1;
    CDigFloat dfDiff2Comp = pow(dfX2, 2) - pow(dfX1,2);
    CDigFloat dfDiff3Comp = pow(dfX2, 3) - pow(dfX1,3);
    
        
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", " called with args:");
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("target value =") + dfTargetValue.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition",string("x in [") + dfX1.RawPrint(10) + ", " + dfX2.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("XSlope  :") + dfSlopeX.RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("XOffset :") + dfOffsetX.RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("YSlope  :") + dfSlopeY.RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("YOffset :") + dfOffsetY.RawPrint(20) );

    
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", "**************************************************************************");    
    
    
//     LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", " called, comparing two ways of calc.");
//     LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("diff1         :") + dfDiff1.RawPrint(20) );
//     LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("diff1 compare :") + dfDiff1Comp.RawPrint(20) );
//     LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("diff1         :") + dfDiff2.RawPrint(20) );
//     LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("diff1 compare :") + dfDiff2Comp.RawPrint(20) );
//     LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("diff1         :") + dfDiff3.RawPrint(20) );
//     LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("diff1 compare :") + dfDiff3Comp.RawPrint(20) );
//     
//     LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", "**************************************************************************");
    
    
    return -dfSlopeX*dfSlopeY/3*dfDiff3+
            (dfSlopeX*dfYConst-dfSlopeY*dfOffsetX)/2.* dfDiff2+
            dfYConst*dfOffsetX*dfDiff1;
}

CDigFloat CProbabilityDensityDistribution::_Integral4TargetValue4Subtraction2(const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const ConvPlanElement& planEl)
{  
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", " called with args:");
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction",string("x iterator points to [") + planEl.xDistriIter->first.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction",string("x in [") + planEl.xLimits.first.RawPrint(10) + ", " + planEl.xLimits.second.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", string("XSlope  :") + Slope(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", string("XOffset :") + Offset(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction",string("y iterator points to [") + planEl.yDistriIter->first.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction",string("y in [") + planEl.yLimits.first.RawPrint(10) + ", " + planEl.yLimits.second.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", string("YSlope  :") + Other.Slope(planEl.yDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", string("YOffset :") + Other.Offset(planEl.yDistriIter).RawPrint(20) );

    MapDFDFType::const_iterator itX = planEl.xDistriIter;
    MapDFDFType::const_iterator itY = planEl.yDistriIter;
    CDigFloat dfX1 = planEl.xLimits.first;
    CDigFloat dfX2 = planEl.xLimits.second;
    CDigFloat dfOffsetX =  Offset(planEl.xDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", string("XOffset :") + dfOffsetX.RawPrint(20) );
    CDigFloat dfSlopeX =  Slope(planEl.xDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", string("XSlope  :") + dfSlopeX.RawPrint(20) );
    CDigFloat dfOffsetY =  Other.Offset(planEl.yDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", string("YOffset :") + dfOffsetY.RawPrint(20) );
    CDigFloat dfSlopeY =  Other.Slope(planEl.yDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", string("YSlope  :") + dfSlopeY.RawPrint(20) );
    CDigFloat dfYConst = -dfSlopeY*dfTargetValue+dfOffsetY;
    CDigFloat dfDiff1 = dfX2 -dfX1;
    CDigFloat dfDiff2 = (dfX1+dfX2) * dfDiff1;
    CDigFloat dfDiff3 = dfDiff1*(pow(dfX1,2) + pow(dfX2,2) +dfX1*dfX2);
    
    
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", "**************************************************************************");
    
    return dfSlopeX*dfSlopeY/3.* dfDiff3+
            (dfSlopeX*dfYConst+dfSlopeY*dfOffsetX)/2.*dfDiff2+
            dfYConst*dfOffsetX*dfDiff1;

}

CDigFloat CProbabilityDensityDistribution::_Integral4TargetValue4Multiplication2(const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const ConvPlanElement& planEl)
{ 
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", " called with args:");
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication",string("x iterator points to [") + planEl.xDistriIter->first.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication",string("x in [") + planEl.xLimits.first.RawPrint(10) + ", " + planEl.xLimits.second.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("XSlope  :") + Slope(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("XOffset :") + Offset(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication",string("y iterator points to [") + planEl.yDistriIter->first.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication",string("y in [") + planEl.yLimits.first.RawPrint(10) + ", " + planEl.yLimits.second.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("YSlope  :") + Other.Slope(planEl.yDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("YOffset :") + Other.Offset(planEl.yDistriIter).RawPrint(20) );

    MapDFDFType::const_iterator itX = planEl.xDistriIter;
    MapDFDFType::const_iterator itY = planEl.yDistriIter;
    CDigFloat dfX1 = planEl.xLimits.first;
    CDigFloat dfX2 = planEl.xLimits.second;
    CDigFloat dfOffsetX =  Offset(planEl.xDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("XOffset :") + dfOffsetX.RawPrint(20) );
    CDigFloat dfSlopeX =  Slope(planEl.xDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("XSlope  :") + dfSlopeX.RawPrint(20) );
    CDigFloat dfOffsetY =  Other.Offset(planEl.yDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("YOffset :") + dfOffsetY.RawPrint(20) );
    CDigFloat dfSlopeY =  Other.Slope(planEl.yDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("YSlope  :") + dfSlopeY.RawPrint(20) );
    CDigFloat dfDiff1 = abs(dfX2) - abs(dfX1);
    
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", "**************************************************************************");
  
    return (dfSlopeX*dfSlopeY*dfTargetValue + dfOffsetX*dfOffsetY)*log( pow(abs(dfX2),sgn(dfX2)) / pow(abs(dfX1), sgn(dfX1)) ) -
           dfOffsetX*dfSlopeY*dfTargetValue*(-dfDiff1)/abs(dfX1)/abs(dfX2)+
           dfOffsetY*dfSlopeX*dfDiff1;

}

CDigFloat CProbabilityDensityDistribution::_Integral4TargetValue4Division2(const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const ConvPlanElement& planEl)
{    
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", " called with args:");
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division",string("x iterator points to [") + planEl.xDistriIter->first.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division",string("x in [") + planEl.xLimits.first.RawPrint(10) + ", " + planEl.xLimits.second.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", string("XSlope  :") + Slope(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", string("XOffset :") + Offset(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division",string("y iterator points to [") + planEl.yDistriIter->first.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division",string("y in [") + planEl.yLimits.first.RawPrint(10) + ", " + planEl.yLimits.second.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", string("YSlope  :") + Other.Slope(planEl.yDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", string("YOffset :") + Other.Offset(planEl.yDistriIter).RawPrint(20) );

    MapDFDFType::const_iterator itX = planEl.xDistriIter;
    MapDFDFType::const_iterator itY = planEl.yDistriIter;
    CDigFloat dfX1 = planEl.xLimits.first;
    CDigFloat dfX2 = planEl.xLimits.second;
    CDigFloat dfOffsetX =  Offset(planEl.xDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", string("XOffset :") + dfOffsetX.RawPrint(20) );
    CDigFloat dfSlopeX =  Slope(planEl.xDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", string("XSlope  :") + dfSlopeX.RawPrint(20) );
    CDigFloat dfOffsetY =  Other.Offset(planEl.yDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", string("YOffset :") + dfOffsetY.RawPrint(20) );
    CDigFloat dfSlopeY =  Other.Slope(planEl.yDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", string("YSlope  :") + dfSlopeY.RawPrint(20) );
    CDigFloat dfDiff2 =  DifferenceNthOrder(dfX1,dfX2,2);// (dfX1+dfX2) * dfDiff1;
    CDigFloat dfDiff3 =  DifferenceNthOrder(dfX1,dfX2,3);//dfDiff1*(pow(dfX1,2) + pow(dfX2,2) +dfX1*dfX2);
    CDigFloat dfDiff4 =  DifferenceNthOrder(dfX1,dfX2,4); //dfDiff3*(dfX1+dfX2);
    
    // depends on the sign of integration limits
    // this case: signs equal ... 
    if(sgn(dfX1) == sgn(dfX2) )
    {
        // ... and negative
        // --> invert the differences
        if( sgn(dfX1) < 0)
        {
            dfDiff2 *=-1;
            dfDiff3 *=-1;
            dfDiff4 *=-1;
            
        }   // if( sgn(dfX1) < 0)
        
        // ... and positive: do nothing: is already done
        
    }   // endif(sgn(dfX1) == sgn(dfX2) )
    // the signs are unequal ....
    else 
    {
        dfDiff2 = sgn(dfX2)*(pow(dfX2,2) + pow(dfX1,2));
        dfDiff3 = sgn(dfX2)*(pow(dfX2,3) + pow(dfX1,3));
        dfDiff4 = sgn(dfX2)*(pow(dfX2,4) + pow(dfX1,4));
        
    }   // endelseif(sgn(dfX1) == sgn(dfX2) )
    
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", "**************************************************************************");
  
  
    return dfSlopeX * dfSlopeY/dfTargetValue/ 4. * dfDiff4 +
        (dfSlopeX*dfOffsetY + dfSlopeY*dfOffsetX/dfTargetValue)/3. *dfDiff3+
        dfOffsetX*dfOffsetY/2.*dfDiff2;

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
void CProbabilityDensityDistribution::_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue2(const ProbDistOp Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const VectorSubIntLimitsType& TotalIntervals4TargetValue, VectorConvPlanElementType& SubIntervals4TargetValue)
{  
     
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("called with args:"));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("operation:") + GetProbDistOpAsString(Operation));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("target value:") + dfTargetValue.RawPrint(10));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("other distri:") + Other.PrintMetaInfo());
     for(int i = 0; i < TotalIntervals4TargetValue.size(); i++){
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("total interval x[" + to_string(i) + "] = [ ") + TotalIntervals4TargetValue[i].first.first.RawPrint(10) + ", " + TotalIntervals4TargetValue[i].first.second.RawPrint(10) + " ]" );
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("total interval y[" + to_string(i) + "] = [ ") + TotalIntervals4TargetValue[i].second.first.RawPrint(10) + ", " + TotalIntervals4TargetValue[i].second.second.RawPrint(10) + " ]" );
     }
     
     
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
    
    LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("iterating over total intervals ...." ));

    // now set the subintervals for each total interval
    for(auto xyint: TotalIntervals4TargetValue)
    {
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... total interval x = [ ") + xyint.first.first.RawPrint(10) + ", " + xyint.first.second.RawPrint(10) + " ]" );
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... total interval y = [ ") + xyint.second.first.RawPrint(10) + ", " + xyint.second.second.RawPrint(10) + " ]" );

        //////////////////////////////////////////////
        // get all relevant supporting points of the x
        // distri
        //////////////////////////////////////////////
        vector<CDigFloat> vdfSubPoints;
        
        // get the iterator of the first interval for x and y distri
        vdfSubPoints.push_back(xyint.first.first);
        vdfSubPoints.push_back(xyint.first.second);
        
        
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... check for supporting points withing this-distri ..."));
        
        // iterate this for supporting points within the total interval
        for(auto x: Distribution())
        {
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... x = ")+ x.first.RawPrint(10));
        
            // break conditions
            if(x.first > xyint.first.second)
            {
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... leaving loop over x-distri: distri is right apart from limits"));
        
                break;
            }
            
            
            if(x.first > xyint.first.first && x.first <xyint.first.second)
            {
                vdfSubPoints.push_back(x.first);
             
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... adding point = ") + vdfSubPoints.back().RawPrint(10));
   
            }
            
        }   //endfor(auto x: Distribution())        
        
        //////////////////////////////////////////////
        // get all relevant supporting points of the y
        // distri: get the Complementary values
        //////////////////////////////////////////////
        
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... check for supporting points withing other-distri (adding Complementary values) ..."));
        
        
        // iterate this for supporting points within the total interval
        for(auto x: Other.Distribution())
        {
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... x(other) = ")+ x.first.RawPrint(10));
            
            // break conditions
            if(x.first > xyint.second.second)
            {
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... leaving loop over other-distri: distri is right apart from limits"));

                break;
            }

            if(x.first > xyint.second.first && x.first <xyint.second.second)
            {
                vdfSubPoints.push_back(_GetComplementaryVariable(Operation,dfTargetValue,x.first,true));
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... adding point = ") + vdfSubPoints.back().RawPrint(10));
            }
            
        }   //endfor(auto x: Other.Distribution())
        
        //////////////////////////////////////////////
        // sort the supporting points and derive
        // sub intervals for x and y
        //////////////////////////////////////////////
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... sorted supporting points ... ") + vdfSubPoints.back().RawPrint(10));
        sort(vdfSubPoints.begin(),vdfSubPoints.end());
        for(auto x: vdfSubPoints)
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... [") +x.RawPrint(10));



        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... iterating over supporting points ... "));
        // iterate over sorted sub points
        for(int idx=0; idx < vdfSubPoints.size()-1; idx++)
        {
            // keeping the limits for x and y
            ConvPlanElement PlanElement;
            
            
            // set the x limits
//             dfXSubLimit.first = vdfSubPoints[idx];
//             dfXSubLimit.second = vdfSubPoints[idx+1];
            PlanElement.xLimits.first = vdfSubPoints[idx];
            PlanElement.xLimits.second = vdfSubPoints[idx+1];
            
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... this interval [") +PlanElement.xLimits.first.RawPrint(10) + ", " +PlanElement.xLimits.second.RawPrint(10));
            
            // set the iterator x/yDistriIter pointing to the first iterator of the correct interval:
            MapDFDFType::const_iterator itDummy;
            GetInterval(PlanElement.xLimits.first,PlanElement.xDistriIter,itDummy);
            
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... interval iterator x points to ") +PlanElement.xDistriIter->first.RawPrint(10) );
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... with offset ") + Other.Offset(PlanElement.xDistriIter).RawPrint(10));
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... with slope ") + Other.Slope(PlanElement.xDistriIter).RawPrint(10));

            
            // derive the y limits
            PlanElement.yLimits.first = _GetComplementaryVariable(Operation, dfTargetValue,PlanElement.xLimits.first,false);
            PlanElement.yLimits.second = _GetComplementaryVariable(Operation, dfTargetValue,PlanElement.xLimits.second,false);

            
            // swapping here avoids a lot of if-cases:
            // the if cases would have to consider the sign of target value, x value, y value, and operation
            if(PlanElement.yLimits.second < PlanElement.yLimits.first )
                swap(PlanElement.yLimits.first, PlanElement.yLimits.second); 
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... other interval [") +PlanElement.yLimits.first.RawPrint(10) + ", " +PlanElement.yLimits.second.RawPrint(10));
            
            // set the iterator x/yDistriIter pointing to the first iterator of the correct interval:
            Other.GetInterval(PlanElement.yLimits.first,PlanElement.yDistriIter,itDummy);
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... interval iterator y points to ") +PlanElement.yDistriIter->first.RawPrint(10));
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... with offset ") + Other.Offset(PlanElement.yDistriIter).RawPrint(10));
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... with slope ") + Other.Slope(PlanElement.yDistriIter).RawPrint(10));
            
            
            // do some checks:
            // first iterator must be less equal left side of limits, i.e. integration interval
            if( PlanElement.xDistriIter->first > PlanElement.xLimits.first ||
                PlanElement.yDistriIter->first > PlanElement.yLimits.first 
            )
            {
                LOGERROR("ErrorLogger", "*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*");
                LOGERROR("ErrorLogger", LS_ProbDist+"::_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue:"); 
                LOGERROR("ErrorLogger","iterator points outside integration sub interval:");
                LOGERROR("ErrorLogger",string("this interval [") + PlanElement.xLimits.first.RawPrint(10) + ", " + PlanElement.xLimits.second.RawPrint(10));                
                LOGERROR("ErrorLogger",string("interval iterator x points to ") +PlanElement.xDistriIter->first.RawPrint(10) );                
                LOGERROR("ErrorLogger",string("other interval [") +PlanElement.yLimits.first.RawPrint(10) + ", " +PlanElement.yLimits.second.RawPrint(10));
                LOGERROR("ErrorLogger",string("interval iterator y points to ") +PlanElement.yDistriIter->first.RawPrint(10) );
                LOGERROR("ErrorLogger", "*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*");
            }
            
            
            
            // now add the ConvPlanElement to result vector
            SubIntervals4TargetValue.push_back(PlanElement);   
            
        }   //endfor(int idx=0; idx < vdfSubPoints.size()-1; idx++)
        
        
    }   //endfor(auto xyint: TotalIntervals4TargetValue)
    LOGTRACE("DistributedNumber::CProbabilityDensityDistribution::_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("************************************************************************"));

}

void CProbabilityDensityDistribution::_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue3(const ProbDistOp Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const VectorSubIntLimitsType& TotalIntervals4TargetValue, VectorConvPlanElementType& SubIntervals4TargetValue)
{
    
     
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("called with args:"));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("operation:") + GetProbDistOpAsString(Operation));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("target value:") + dfTargetValue.RawPrint(10));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("other distri:") + Other.PrintMetaInfo());
     for(int i = 0; i < TotalIntervals4TargetValue.size(); i++){
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("total interval x[" + to_string(i) + "] = [ ") + TotalIntervals4TargetValue[i].first.first.RawPrint(10) + ", " + TotalIntervals4TargetValue[i].first.second.RawPrint(10) + " ]" );
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("total interval y[" + to_string(i) + "] = [ ") + TotalIntervals4TargetValue[i].second.first.RawPrint(10) + ", " + TotalIntervals4TargetValue[i].second.second.RawPrint(10) + " ]" );
     }
     
     
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
    
    LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("iterating over total intervals ...." ));

    // now set the subintervals for each total interval
    for(auto xyint: TotalIntervals4TargetValue)
    {
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... total interval x = [ ") + xyint.first.first.RawPrint(10) + ", " + xyint.first.second.RawPrint(10) + " ]" );
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... total interval y = [ ") + xyint.second.first.RawPrint(10) + ", " + xyint.second.second.RawPrint(10) + " ]" );

        //////////////////////////////////////////////
        // get all relevant supporting points of the x
        // distri
        //////////////////////////////////////////////
        vector<CDigFloat> vdfSubPoints;
        
        // get the iterator of the first interval for x and y distri
        vdfSubPoints.push_back(xyint.first.first);
        vdfSubPoints.push_back(xyint.first.second);
        
        
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... check for supporting points withing this-distri ..."));
        
        // iterate this for supporting points within the total interval
        for(auto x: Distribution())
        {
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... x = ")+ x.first.RawPrint(10));
        
            // break conditions
            if(x.first > xyint.first.second)
            {
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... leaving loop over x-distri: distri is right apart from limits"));
        
                break;
            }
            
            
            if(x.first > xyint.first.first && x.first <xyint.first.second)
            {
                vdfSubPoints.push_back(x.first);
             
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... adding point = ") + vdfSubPoints.back().RawPrint(10));
   
            }
            
        }   //endfor(auto x: Distribution())        
        
        //////////////////////////////////////////////
        // get all relevant supporting points of the y
        // distri: get the Complementary values
        //////////////////////////////////////////////
        
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... check for supporting points withing other-distri (adding Complementary values) ..."));
        
        
        // iterate this for supporting points within the total interval
        for(auto x: Other.Distribution())
        {
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... x(other) = ")+ x.first.RawPrint(10));
            
            // break conditions
            if(x.first > xyint.second.second)
            {
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... leaving loop over other-distri: distri is right apart from limits"));

                break;
            }

            if(x.first > xyint.second.first && x.first <xyint.second.second)
            {
                vdfSubPoints.push_back(_GetComplementaryVariable(Operation,dfTargetValue,x.first,true));
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... adding point = ") + vdfSubPoints.back().RawPrint(10));
            }
            
        }   //endfor(auto x: Other.Distribution())
        
        //////////////////////////////////////////////
        // sort the supporting points and derive
        // sub intervals for x and y
        //////////////////////////////////////////////
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... sorted supporting points ... ") + vdfSubPoints.back().RawPrint(10));
        sort(vdfSubPoints.begin(),vdfSubPoints.end());
        for(auto x: vdfSubPoints)
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... [") +x.RawPrint(10));
        
        // set the number of steps for simple search being more effective than the binary search:
        unsigned int uiMaxSearchStepsX = (unsigned int)log(CDigFloat(Distribution().size()),2).RawValue();
        unsigned int uiMaxSearchStepsY = (unsigned int)log(CDigFloat(Other.Distribution().size()),2).RawValue();
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... max. simple search steps for x: ") +to_string(uiMaxSearchStepsX));
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... max. simple search steps for y: ") +to_string(uiMaxSearchStepsY));
        
        // iterate over sorted sub points
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... iterating over supporting points ... "));
        for(int idx=0; idx < vdfSubPoints.size()-1; idx++)
        {
            // keeping the limits for x and y
            ConvPlanElement PlanElement;
            
            // set to the values of the last planelement
            if(SubIntervals4TargetValue.size()>0)
                PlanElement = SubIntervals4TargetValue.back();
            else
            {
                PlanElement.xDistriIter = Distribution().begin();
                PlanElement.yDistriIter = Other.Distribution().begin();
            }
            
            // set the x limits
            PlanElement.xLimits.first = vdfSubPoints[idx];
            PlanElement.xLimits.second = vdfSubPoints[idx+1];
            
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... this interval [") +PlanElement.xLimits.first.RawPrint(10) + ", " +PlanElement.xLimits.second.RawPrint(10));
            
            // in case there is a preceeding element: check if the iterator or the consecutive iterator can be used
            // to avoid binary search            
            if(!_GetInterval(PlanElement.xLimits.first,PlanElement.xDistriIter, uiMaxSearchStepsX))
            {
                // set the iterator x/yDistriIter pointing to the first iterator of the correct interval:
                MapDFDFType::const_iterator itDummy;
                GetInterval(PlanElement.xLimits.first,PlanElement.xDistriIter,itDummy);

            }   // endif(bDoBinarySearch)
            
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... interval iterator x points to ") +PlanElement.xDistriIter->first.RawPrint(10) );
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... xOffset(pre calc.)  = ") + Offset(PlanElement.xDistriIter).RawPrint(10) );
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... xSlope(pre calc.)  = ") + Slope(PlanElement.xDistriIter).RawPrint(10) );
            
            // derive the y limits
            PlanElement.yLimits.first = _GetComplementaryVariable(Operation, dfTargetValue,PlanElement.xLimits.first,false);
            PlanElement.yLimits.second = _GetComplementaryVariable(Operation, dfTargetValue,PlanElement.xLimits.second,false);

            
            // swapping here avoids a lot of if-cases:
            // the if cases would have to consider the sign of target value, x value, y value, and operation
            if(PlanElement.yLimits.second < PlanElement.yLimits.first )
                swap(PlanElement.yLimits.first, PlanElement.yLimits.second); 
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... other interval [") +PlanElement.xLimits.first.RawPrint(10) + ", " +PlanElement.xLimits.second.RawPrint(10));
            
            if(!_GetInterval(PlanElement.yLimits.first,PlanElement.yDistriIter, uiMaxSearchStepsY))
            {
                // set the iterator x/yDistriIter pointing to the first iterator of the correct interval:
                MapDFDFType::const_iterator itDummy;
                Other.GetInterval(PlanElement.yLimits.first,PlanElement.yDistriIter,itDummy);

            }   // endif(!_GetInterval(PlanElement.yLimits.first,PlanElement.yDistriIter, uiMaxSearchStepsY))
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... interval iterator y points to ") +PlanElement.yDistriIter->first.RawPrint(10));
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... yOffset(pre calc.)  = ") + Other.Offset(PlanElement.yDistriIter).RawPrint(10) );
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... ySlope(pre calc.)  = ") + Other.Slope(PlanElement.yDistriIter).RawPrint(10) );
            
            
            // do some checks:
            // first iterator must be less equal left side of limits, i.e. integration interval
            if( PlanElement.xDistriIter->first > PlanElement.xLimits.first ||
                PlanElement.yDistriIter->first > PlanElement.yLimits.first 
            )
            {
                LOGERROR("ErrorLogger", "*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*");
                LOGERROR("ErrorLogger", LS_ProbDist+"::_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue:"); 
                LOGERROR("ErrorLogger","iterator points outside integration sub interval:");
                LOGERROR("ErrorLogger",string("this interval [") + PlanElement.xLimits.first.RawPrint(10) + ", " + PlanElement.xLimits.second.RawPrint(10));                
                LOGERROR("ErrorLogger",string("interval iterator x points to ") +PlanElement.xDistriIter->first.RawPrint(10) );                
                LOGERROR("ErrorLogger",string("other interval [") +PlanElement.yLimits.first.RawPrint(10) + ", " +PlanElement.yLimits.second.RawPrint(10));
                LOGERROR("ErrorLogger",string("interval iterator y points to ") +PlanElement.yDistriIter->first.RawPrint(10) );
                LOGERROR("ErrorLogger", "*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*");
            }
            
            
            
            // now add the ConvPlanElement to result vector
            SubIntervals4TargetValue.push_back(PlanElement);   
            
        }   //endfor(int idx=0; idx < vdfSubPoints.size()-1; idx++)
        
        
    }   //endfor(auto xyint: TotalIntervals4TargetValue)
    LOGTRACE("DistributedNumber::CProbabilityDensityDistribution::_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("************************************************************************"));

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
    LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("called with args:"));
    LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("operation:") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("target value :") + dfTargetValue.RawPrint(10));
    LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("other distri :") + Other.PrintMetaInfo());
    for(int i = 0; i<vpdXIntervals.size(); i++)
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("x-interval[") + to_string(i) + "] = [ " + vpdXIntervals[i].first.RawPrint(10) + ", " + vpdXIntervals[i].second.RawPrint(10) + " ]") ;
    
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
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("interval = [ " + xint.first.RawPrint(10) + ", " + xint.second.RawPrint(10) + " ]") );
        
        /////////////////////////////////////////////////
        // calculate the complementary interval YComp from x
        /////////////////////////////////////////////////
        CDigFloat dfYCompMin, dfYCompMax;
        dfYCompMin = _GetComplementaryVariable(Operation, dfTargetValue,xint.first,false);
        dfYCompMax = _GetComplementaryVariable(Operation, dfTargetValue,xint.second,false);
        
        
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("calculated complementary Y:"));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYCompMin:") + dfYCompMin.RawPrint(10));
        LOGTRACE(LS_ProbDist+"::_GetTotalIntegrationInterval4TargetValue",string("dfYCompMax:") + dfYCompMax.RawPrint(10));
        
        // swapping here avoids a lot of if-cases:
        // the if cases would have to consider the sign of target value, x value, y value, and operation
        if(dfYCompMax < dfYCompMin)
            swap(dfYCompMin, dfYCompMax);       
        
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("after swap:"));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYCompMin:") + dfYCompMin.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYCompMax:") + dfYCompMax.RawPrint(10));
                
        /////////////////////////////////////////////////
        // calculate the overlap of the complementary 
        // ycomp and the original y interval --> YOverlap
        /////////////////////////////////////////////////
        CDigFloat dfYOverlapMin, dfYOverlapMax;
        dfYOverlapMin = max(dfYCompMin,Other.firstVariable());
        dfYOverlapMax = min(dfYCompMax,Other.lastVariable());  
        
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("calculated overlapping Y:"));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYOverlapMin:") + dfYOverlapMin.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYOverlapMax:") + dfYOverlapMax.RawPrint(10));
        
        // the overlap is empty in case dfYOverlapMax <= dfYOverlapMin
        if(dfYOverlapMax <= dfYOverlapMin)
        {
            
            LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("skipping this intervall: zero length because Y overlap max <= min"));
            continue;
        }
        
        /////////////////////////////////////////////////
        // calculate the complementary interval XComp from YOverlap
        /////////////////////////////////////////////////
        CDigFloat dfXCompMin, dfXCompMax;
        dfXCompMin = _GetComplementaryVariable(Operation, dfTargetValue,dfYOverlapMin,true);
        dfXCompMax = _GetComplementaryVariable(Operation, dfTargetValue,dfYOverlapMax,true);
        
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("calculated complementary X from complementary Y:"));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfXCompMin:") + dfXCompMin.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfXCompMax:") + dfXCompMax.RawPrint(10));
        
        // swapping here avoids a lot of if-cases:
        // the if cases would have to consider the sign of target value, x value, y value, and operation
        if(dfXCompMax < dfXCompMin)
            swap(dfXCompMin, dfXCompMax);    
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("after swap:"));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfXCompMin:") + dfXCompMin.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfXCompMax:") + dfXCompMax.RawPrint(10));
        
        /////////////////////////////////////////////////
        // calculate the overlap of the complementary 
        // Xcomp and the original x interval --> XOverlap 
        /////////////////////////////////////////////////
        CDigFloat dfXOverlapMin, dfXOverlapMax;
        dfXOverlapMin = max(dfXCompMin,firstVariable());
        dfXOverlapMax = min(dfXCompMax,lastVariable());
        
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("calculated overlapping Y:"));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYOverlapMin:") + dfYOverlapMin.RawPrint(10));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYOverlapMax:") + dfYOverlapMax.RawPrint(10));
   
        
        // the overlap is empty in case dfXOverlapMax <= dfXOverlapMin
        if(dfXOverlapMax <= dfXOverlapMin)
        {
            LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("skipping this intervall: zero length because X overlap max <= min"));
            continue;
        }
        
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
    
    
    LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("************************************************************************"));

}


void CProbabilityDensityDistribution::_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue(const ProbDistOp Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const VectorSubIntLimitsType& TotalIntervals4TargetValue, VectorSubIntLimitsType& SubIntervals4TargetValue)
 {   
     
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("called with args:"));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("operation:") + GetProbDistOpAsString(Operation));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("target value:") + dfTargetValue.RawPrint(10));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("other distri:") + Other.PrintMetaInfo());
     for(int i = 0; i < TotalIntervals4TargetValue.size(); i++){
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("total interval x[" + to_string(i) + "] = [ ") + TotalIntervals4TargetValue[i].first.first.RawPrint(10) + ", " + TotalIntervals4TargetValue[i].first.second.RawPrint(10) + " ]" );
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("total interval y[" + to_string(i) + "] = [ ") + TotalIntervals4TargetValue[i].second.first.RawPrint(10) + ", " + TotalIntervals4TargetValue[i].second.second.RawPrint(10) + " ]" );
     }
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
    
    LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("iterating over total intervals ...." ));

    // now set the subintervals for each total interval
    for(auto xyint: TotalIntervals4TargetValue)
    {
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... total interval x = [ ") + xyint.first.first.RawPrint(10) + ", " + xyint.first.second.RawPrint(10) + " ]" );
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... total interval y = [ ") + xyint.second.first.RawPrint(10) + ", " + xyint.second.second.RawPrint(10) + " ]" );

        //////////////////////////////////////////////
        // get all relevant supporting points of the x
        // distri
        //////////////////////////////////////////////
        vector<CDigFloat> vdfSubPoints;
        
        // get the limits of the total integral
        vdfSubPoints.push_back(xyint.first.first);
        vdfSubPoints.push_back(xyint.first.second);
        
        
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... check for supporting points withing this-distri ..."));
        
        // iterate this for supporting points within the total interval
        for(auto x: Distribution())
        {
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... x = ")+ x.first.RawPrint(10));
        
            // break conditions
            if(x.first > xyint.first.second)
            {
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... leaving loop over x-distri: distri is right apart from limits"));
        
                break;
            }
            
            
            if(x.first > xyint.first.first && x.first <xyint.first.second)
            {
                vdfSubPoints.push_back(x.first);
             
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... adding point = ") + vdfSubPoints.back().RawPrint(10));
   
            }
            
        }   //endfor(auto x: Distribution())        
        
        //////////////////////////////////////////////
        // get all relevant supporting points of the y
        // distri: get the Complementary values
        //////////////////////////////////////////////
        
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... check for supporting points withing other-distri (adding Complementary values) ..."));
        
        
        // iterate this for supporting points within the total interval
        for(auto x: Other.Distribution())
        {
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... x(other) = ")+ x.first.RawPrint(10));
            
            // break conditions
            if(x.first > xyint.second.second)
            {
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... leaving loop over other-distri: distri is right apart from limits"));

                break;
            }

            if(x.first > xyint.second.first && x.first <xyint.second.second)
            {
                vdfSubPoints.push_back(_GetComplementaryVariable(Operation,dfTargetValue,x.first,true));
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... adding point = ") + vdfSubPoints.back().RawPrint(10));
            }
            
        }   //endfor(auto x: Other.Distribution())
        
        //////////////////////////////////////////////
        // sort the supporting points and derive
        // sub intervals for x and y
        //////////////////////////////////////////////
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... sorted supporting points ... ") + vdfSubPoints.back().RawPrint(10));
        sort(vdfSubPoints.begin(),vdfSubPoints.end());
        for(auto x: vdfSubPoints)
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... [") +x.RawPrint(10));



        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... iterating over supporting points ... "));
        // iterate over sorted sub points
        for(int idx=0; idx < vdfSubPoints.size()-1; idx++)
        {
            // keeping the limits for x and y
            PairDFType dfXSubLimit, dfYSubLimit;
            
            // set the x limits
            dfXSubLimit.first = vdfSubPoints[idx];
            dfXSubLimit.second = vdfSubPoints[idx+1];
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... this interval [") +dfXSubLimit.first.RawPrint(10) + ", " +dfXSubLimit.second.RawPrint(10));

            
            // derive the y limits
            dfYSubLimit.first = _GetComplementaryVariable(Operation, dfTargetValue,dfXSubLimit.first,false);
            dfYSubLimit.second = _GetComplementaryVariable(Operation, dfTargetValue,dfXSubLimit.second,false);

            
            // swapping here avoids a lot of if-cases:
            // the if cases would have to consider the sign of target value, x value, y value, and operation
            if(dfYSubLimit.second < dfYSubLimit.first )
                swap(dfYSubLimit.first, dfYSubLimit.second); 
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... other interval [") +dfXSubLimit.first.RawPrint(10) + ", " +dfXSubLimit.second.RawPrint(10));
            
            // now add the pair of limits (xmin,xmax) (ymin, ymax) to result vector
            SubIntervals4TargetValue.push_back(SubIntLimitsType(dfXSubLimit, dfYSubLimit));   
            
        }   //endfor(int idx=0; idx < vdfSubPoints.size()-1; idx++)
        
        
    }   //endfor(auto xyint: TotalIntervals4TargetValue)
    LOGTRACE("DistributedNumber::CProbabilityDensityDistribution::_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("************************************************************************"));

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

CProbabilityDensityDistribution operator/(const CDigFloat& dfValue, CProbabilityDensityDistribution& pdDistri)
{
    CProbabilityDensityDistribution pdResult;
    for(auto iel:pdDistri.Distribution())
        pdResult.Add(dfValue/iel.first,iel.second);
    
    return pdResult;
}


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
