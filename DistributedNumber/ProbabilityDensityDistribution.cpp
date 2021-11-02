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
    LOGTRACE(LS_ProbDist+"ConstructorConstProbDistriRef",string("called with args"));
    LOGTRACE(LS_ProbDist+"ConstructorConstProbDistriRef",string("other distri: \n") + other.PrintMetaInfo());
    *this = other;
}

CProbabilityDensityDistribution::CProbabilityDensityDistribution(const CDistribution& other)
{
    LOGTRACE(LS_ProbDist+"ConstructorConstDistriRef",string("called with args"));
    LOGTRACE(LS_ProbDist+"ConstructorConstDistriRef",string("other distri: \n¸") + other.PrintMetaInfo());
    CDistribution::operator=(other);
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
    LOGTRACE(LS_ProbDist+"OperatorEqual",string("called with args"));
    LOGTRACE(LS_ProbDist+"OperatorEqual",string("other distri: \n") + other.PrintMetaInfo());
    CDistribution::operator=(other);
    bNormalized = other.bNormalized;
    dfAbsIntegral = other.dfAbsIntegral;
    m_nIntegrationSteps = other.IntegrationSteps();
    m_nSubIntegrationSteps = other.SubIntegrationSteps();
    
    LOGTRACE(LS_ProbDist+"OperatorEqual",string("initialize this distri: \n") + PrintMetaInfo()+ "\n" + Print(10));
    // hand over linear pars 
    m_LinPars.clear();
    MapDFDFType::const_iterator itThis=mDistribution.begin();
    LOGTRACE(LS_ProbDist+"OperatorEqual",string("iterating over other distri ...\n") );
    bool bSetNonLinTrafoVar = other.m_NonLinTrafoVariables.size()>0;
   for(MapDFDFType::const_iterator itOther = other.mDistribution.begin(); itOther != other.mDistribution.end(); itOther++)
   {
       
        LOGTRACE(LS_ProbDist+"OperatorEqual",string("... itThis = ") + ::Print(itThis));
        LOGTRACE(LS_ProbDist+"OperatorEqual",string("... itOther = ") + ::Print(itOther));
        
        m_LinPars[itThis] = PairDFType( other.Offset(itOther), other.Slope(itOther));
        LOGTRACE(LS_ProbDist+"OperatorEqual",string("... added = (") + m_LinPars[itThis].first.RawPrint(15) + ", " + m_LinPars[itThis].second.RawPrint(15));
        
        // hand over non-lin trafo variables
        if(bSetNonLinTrafoVar)
            m_NonLinTrafoVariables[itThis] = other.m_NonLinTrafoVariables.at(itOther);
        itThis++;
    }
   
    LOGTRACE(LS_ProbDist+"OperatorEqual",string("lin pars other \n") + other.PrintLinPars());
    LOGTRACE(LS_ProbDist+"OperatorEqual",string("lin pars this \n") + PrintLinPars());
    return *this;
}

CProbabilityDensityDistribution& CProbabilityDensityDistribution::operator=(const CDistribution& other)
{
    LOGTRACE(LS_ProbDist+"OperatorEqualConstDistriRef",string("called with args")); 
    LOGTRACE(LS_ProbDist+"OperatorEqualConstDistriRef",string("other ")+ other.PrintMetaInfo());
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

// CProbabilityDensityDistribution CProbabilityDensityDistribution::operator+(CProbabilityDensityDistribution& other)
// {
//     LOGTRACE(LS_ProbDist+"OperatorPlusProbDist","called with args:");
//     LOGTRACE(LS_ProbDist+"OperatorPlusProbDist",string("other \n") + other.PrintMetaInfo());
//     CProbabilityDensityDistribution ProbDistResult(*this);
//     
//     LOGTRACE(LS_ProbDist+"OperatorPlusProbDist",string("init result from this: \n") + PrintMetaInfo());
//     LOGTRACE(LS_ProbDist+"OperatorPlusProbDist",string("resulting in: \n") + ProbDistResult.PrintMetaInfo());
//     ProbDistResult += other;
//     return ProbDistResult;
// }

CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator-=(CProbabilityDensityDistribution& other)
{
    return _GeneralOperatorsFunctionAnalytical(other, ProbDistOp::pdoMinus );
}

// CProbabilityDensityDistribution CProbabilityDensityDistribution::operator-(CProbabilityDensityDistribution& other)
// {
//     CProbabilityDensityDistribution ProbDistResult(*this);
//     ProbDistResult -= other;
//     return ProbDistResult;
// }


CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator*=(CProbabilityDensityDistribution& other)
{
    return _GeneralOperatorsFunctionAnalytical(other, ProbDistOp::pdoMult );
}
// CProbabilityDensityDistribution CProbabilityDensityDistribution::operator*(CProbabilityDensityDistribution& other)
// {
//     CProbabilityDensityDistribution ProbDistResult(*this);
//     ProbDistResult *= other;
//     return ProbDistResult;
// }

CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator/=(CProbabilityDensityDistribution& other)
{
    return _GeneralOperatorsFunctionAnalytical(other, ProbDistOp::pdoDiv);
}

// CProbabilityDensityDistribution CProbabilityDensityDistribution::operator/(CProbabilityDensityDistribution& other)
// {
//     CProbabilityDensityDistribution ProbDistResult(*this);
//     ProbDistResult /= other;
//     return ProbDistResult;
// }

////////////////////////////////////////////
// arithmetic operators with CDigFloat
////////////////////////////////////////////
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator+=(const CDigFloat& Value)
{

    Shift(Value);
    return *this;
}
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator-=(const CDigFloat& Value)
{
    Shift(-Value);
    return *this;
}
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator*=(const CDigFloat& Value)
{
    CDistribution::operator*=(Value);
    Normalize();
    return *this;
}
CProbabilityDensityDistribution & CProbabilityDensityDistribution::operator/=(const CDigFloat& Value)
{
    CDistribution::operator/=(Value);
    Normalize();    
    return *this;
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator+(const CDigFloat& Value)
{
    CProbabilityDensityDistribution Result(*this);
    Result+=Value;
    return Result;
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator-(const CDigFloat& Value)
{
    CProbabilityDensityDistribution Result(*this);
    Result-=Value;
    return Result;
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator*(const CDigFloat& Value)
{
    CProbabilityDensityDistribution Result(*this);
    Result*=Value;
    return Result;
}

CProbabilityDensityDistribution CProbabilityDensityDistribution::operator/(const CDigFloat& Value)
{
    CProbabilityDensityDistribution Result(*this);
    Result/=Value;
    return Result;
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
    
    // set lin pars manually 
    MapDFDFType::const_iterator it = mDistribution.begin();
    m_LinPars[it]=PairDFType(1./(dfXMax-dfXMin),0.);
    it++;
    m_LinPars[it]=PairDFType(0.,0.);
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
    
    // set lin pars
    _SetLinPars();
    
}

void CProbabilityDensityDistribution::NormalDistribution(const CDigFloat& dfMu, const CDigFloat& dfSigma, const int nPoints /* = 50*/)
{
    LOGTRACE(LS_ProbDist+"NormalDistribution","called with args:");
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("int Points              =") + to_string(nPoints));
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("CDigFloat mean          =") + dfMu.RawPrint(15));
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("CDigFloat sigma         =") + dfSigma.RawPrint(15));
    
    // check consistent parameters:
    if(dfSigma <= 0 )
        return;
    
    // reset
    Reset();
    
    CDigFloat dfFac1 = 2*dfSigma*dfSigma;
    CDigFloat dfPreFac = 1/sqrt(dfFac1*M_PIf64);
    
    LOGTRACE(LS_ProbDist+"NormalDistribution","calculated:");
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("CDigFloat dfFac1 =") + dfFac1.RawPrint(15));
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("CDigFloat dfPreFac (amplitude) =") + dfPreFac.RawPrint(15));
    
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
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("CDigFloat dfPreFac (amplitude) =") + dfPreFac.RawPrint(15));
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("CDigFloat dfDY (amplitude diff.) =") + dfdY.RawPrint(15));
    
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
        
        LOGTRACE(LS_ProbDist+"NormalDistribution",string("added point = (") + (dfMu-dfdXPos).RawPrint(15) + ", " + dfY.RawPrint(15) + ")" );
        LOGTRACE(LS_ProbDist+"NormalDistribution",string("added point = (") + (dfMu+dfdXPos).RawPrint(15) + ", " + dfY.RawPrint(15) + ")" );
        
    }
    
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("integral (should be 1) = ") + AbsIntegral().RawPrint(15));
    

    Normalize();
    
    // DEBUG
    m_LinPars.clear();
    _SetLinPars();

   LOGTRACE(LS_ProbDist+"NormalDistribution",string("lin. pars. set automatically\n") + PrintLinPars());
   
   // now set the linear parameters : init dfY and rel. pre. fac and iterators for up and down
   dfY = dfPreFac;
   dfRelPreFac = 1.;
   MapDFDFType::const_iterator itUp, itDown;
   itUp = Distribution().find(dfMu);
   itDown = Distribution().find(dfMu);
   
   m_LinPars.clear();
   
   LOGTRACE(LS_ProbDist+"NormalDistribution",string("iterating over ") + to_string(nPtsEvenHalf) + " points, setting linear parameters manually ...");
   CDigFloat dfdXPosOld = 0;
   for(int ipt=0;ipt<nPtsEvenHalf; ipt++)
   {
       // calculate the dx from old values
       CDigFloat dfdX = -dfdXPosOld+sqrt(dfdXPosOld*dfdXPosOld+dfFac1*log((dfY)/(dfY-dfdY)));
       
       // now calc. the y - stuff
       dfY -= dfdY;
       dfRelPreFac = dfY / dfPreFac;
       CDigFloat dfdXPos = sqrt(dfFac1*(-log(dfRelPreFac)));
       LOGTRACE(LS_ProbDist+"NormalDistribution",string("... dfdY ") + dfdY.RawPrint(15));
       LOGTRACE(LS_ProbDist+"NormalDistribution",string("... dx ") + (dfdX).RawPrint(15));
       
       
       // lower slope and offset: decrement before
       itDown--;
       m_LinPars[itDown]=PairDFType(itDown->second - dfdY/dfdX*itDown->first,dfdY / dfdX);
       // upper slope and offset: increment after
       // last elements are zero
           m_LinPars[itUp]=PairDFType(itUp->second + dfdY/dfdX*itUp->first, -dfdY / dfdX);
       
       LOGTRACE(LS_ProbDist+"NormalDistribution",string("... added lin par [") + (itUp->first).RawPrint(15) + ", " + (itUp->second).RawPrint(15) + "] = (" + m_LinPars[itUp].first.RawPrint(15) + ", " + m_LinPars[itUp].second.RawPrint(15) );
       
       LOGTRACE(LS_ProbDist+"NormalDistribution",string("... added lin par [") + (itDown->first).RawPrint(15) + ", " + (itDown->second).RawPrint(15) + "] = (" + m_LinPars[itDown].first.RawPrint(15) + ", " + m_LinPars[itDown].second.RawPrint(15) );
       
       
       itUp++;   
       dfdXPosOld = dfdXPos;
       
       
   }
   
   // set the final element on upper side:
   m_LinPars[itUp]=PairDFType(0.,0.);
   
   LOGTRACE(LS_ProbDist+"NormalDistribution",string("lin. pars. set manually\n") + PrintLinPars());
   
    
    
    LOGTRACE(LS_ProbDist+"NormalDistribution",string("************************************************************"));
}

void CProbabilityDensityDistribution::ExponentialDistribution(const CDigFloat& dfA, const CDigFloat& dfF, const unsigned int nPoints)
{
    
}

///////////////////////////////////
// function on all distribution variables
///////////////////////////////////    
void CProbabilityDensityDistribution::Shift(const CDigFloat& shift)
{
    LOGTRACE(LS_ProbDist+"Shift", string("called with args:"));
    LOGTRACE(LS_ProbDist+"Shift", string("shift: ") + shift.RawPrint(15));
    // this can simply applied as the normalization is not affected 
    CDistribution::Shift(shift);
    
    // take care for the linear parameters
    if(m_LinPars.size() == 0)
    {
        // calculate for the first time
        _SetLinPars();
    
    }
    else
    {
        // iterate over the distribution setting the linear offset: 
        // slope is unchanged for shifting operation
        for(MapDFDFType::const_iterator iel = mDistribution.begin(); iel != mDistribution.end(); iel++)
        {
            // the offset must be corrected: the slope stays the same
            // the offset is generally calculated by:
            // offset = p(x) - slope * x
            //
            // the shift is the transformation:
            // x --> x' = x + dx (=shift)
            //
            // the transformed offset is:
            // offset' = p(x') - slope' * x'
            // where
            // p(x') = p(x)
            // slope' = slope
            // x' = x + dx
            //
            // follows:
            // offset' = p(x) - slope *(x + dx)
            // therefore:
            // offset'-offset = p(x) - slope * (x + dx) - [ p(x) - slope * x]
            //                = -slope*dx
            m_LinPars[iel].first -= m_LinPars[iel].second*shift;
            
        }   //endfor(MapDFDFType::const_iterator iel = mDistribution.begin(); iel != mDistribution.end(); iel++)
    }
    
}
void CProbabilityDensityDistribution::Scale(const CDigFloat& scale)
{
    // scaling needs multiplying the x - values ...
    CDistribution::Scale(scale);
    
    // ... and dividing the y - values to keep the integral = 1
    CDigFloat dfAbsScale = abs(scale);
    for(auto iel: mDistribution)
        iel.second /= dfAbsScale;
    
    if(m_LinPars.size() == 0)
    {
        _SetLinPars();
        
    }   //endif(m_LinPars.size() == 0)
    else
    {
        for(MapDFDFType::const_iterator iel = mDistribution.begin(); iel != mDistribution.end(); iel++)
        {
            // slope is divided by scale²:
            // the x interval is multiplied by scale
            // the y interval is multiplied by 1./|scale|
            // therfore : y / x ~ 1/|scale|/scale
            CDigFloat dfSlope = Slope(iel)/dfAbsScale/scale;
            CDigFloat dfOffset = iel->second - dfSlope*iel->first;
            m_LinPars[iel] = PairDFType(dfOffset, dfSlope);
            
        }   //endfor(MapDFDFType::const_iterator iel = mDistribution.begin(); iel != mDistribution.end(); iel++)
        
    }
        
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
    
    
    LOGTRACE(LS_ProbDist+"Variance", string("original mean = ") + Mean().RawPrint(15));
    LOGTRACE(LS_ProbDist+"Variance", string("actual mean of copy (must be the same) = ") + pdd4Calc.Mean().RawPrint(15));
    
    
    // shift by mean: mean is now zero because of Integral (x-<x>) = 0
    pdd4Calc.Shift(-Mean());
    
    
    LOGTRACE(LS_ProbDist+"Variance", string("control: mean must be zero = ") + pdd4Calc.Mean().RawPrint(15));
    LOGTRACE(LS_ProbDist+"Variance", string("result = ") + pdd4Calc.Mean(2).RawPrint(15));
    
    LOGTRACE(LS_ProbDist+"Variance",string("************************************************************"));
    
    // return the second moment which is (x-<x>)² identical to x²
    return pdd4Calc.Mean(2);
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
    for(MapDFDFType::iterator iel = mDistribution.begin(); iel != Distribution().end(); iel++)
        iel->second*=dfAbsIntegral;
    bNormalized = false;
    
}

void CProbabilityDensityDistribution::Normalize(bool bWithIntegralError /* = false*/)
{

    // do not normalize again ... only once
    if(bNormalized)
        return;
    
    dfAbsIntegral = CDistribution::AbsIntegral();
    CDigFloat dfAbsIntegralWithOrWithoutError = dfAbsIntegral;
    if(!bWithIntegralError)
        dfAbsIntegralWithOrWithoutError.ResetError();
    for(MapDFDFType::iterator iel = mDistribution.begin(); iel != Distribution().end(); iel++)
        iel->second/=dfAbsIntegralWithOrWithoutError;
    
    // remember we have been normalized
    bNormalized = true;
    
}

std::__cxx11::string CProbabilityDensityDistribution::PrintLinPars() const
{
    ostringstream oss;
//     for(auto iel: mDistribution)
    for(MapDFDFType::const_iterator iel = mDistribution.begin(); iel != mDistribution.end(); iel++)
    {
        oss << "[" << iel->first.RawPrint(15) << " / " << iel->second.RawPrint(15) + " ] =";
        oss << "( " << m_LinPars.at(iel).first.RawPrint(15) << " / ";
        oss << m_LinPars.at(iel).second.RawPrint(15) << " )" << endl;
        
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
void CProbabilityDensityDistribution::_SetPeripherialMember(const CProbabilityDensityDistribution& Other)
{
    // reset to old values 
    MaxBinarySearchIterations(Other.MaxBinarySearchIterations());
    IntegrationSteps(Other.IntegrationSteps());
    SubIntegrationSteps(Other.SubIntegrationSteps());
}

void CProbabilityDensityDistribution::_SetLinPars()
{
//     if(m_LinPars.size() == 0)
    m_LinPars.clear();
    for(MapDFDFType::const_iterator iel = mDistribution.begin(); iel != mDistribution.end(); iel++)
        m_LinPars[iel] = PairDFType(_LinearOffset(iel),_LinearSlope(iel));
}

CProbabilityDensityDistribution & CProbabilityDensityDistribution::_GeneralOperatorsFunctionAnalytical(CProbabilityDensityDistribution& Other, const ProbDistOp Operation)
{   
    
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("called with args:"));
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("probdist other:\n") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("operation     : ") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("probdist this :\n") + PrintMetaInfo());
    

    // generate plan : setting x-intervals and target values
    // dfTargetValueStep is needed for accurate calculation of the lin. pars.
    VectorPairDFType vpdfXIntervals;
    vector<CDigFloat> vdfTargetValues;
    CDigFloat dfTargetValueStep;
    _PrepareConvolution(Other, Operation, vdfTargetValues, vpdfXIntervals, dfTargetValueStep);
    
    // set all linear parameters for this and other if not already done
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
    unsigned int uiNoOfThreads = PDD_NO_OF_THREADS;
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
            LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical", string("... [") + (*vpvdfTargetValuePackages[ith.first])[itv].RawPrint(10,false) + "] = " + (*vpvdfResultPackages[ith.first])[itv].RawPrint(15)   );
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
        // set distri values
        for(auto iel: TargetDistri)
            mDistribution[iel.first] = iel.second;
        
        bNormalized=false;
        Normalize();
        
        // DEBUG
        m_LinPars.clear();
        _SetLinPars();
       LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("lin. pars. set automatically\n") + PrintLinPars());
       
       // set lin. pars manually with accurate delta X
       m_LinPars.clear();
       for(MapDFDFType::const_iterator iel = Distribution().begin(); iel != Distribution().end(); iel++)
       {
           MapDFDFType::const_iterator ielNext = iel;
           ielNext++;
           if(ielNext == Distribution().end())
               m_LinPars[iel] = PairDFType(0.,0.);
           else
               m_LinPars[iel] = PairDFType(iel->second - (ielNext->second-iel->second)/dfTargetValueStep*iel->first,(ielNext->second-iel->second)/dfTargetValueStep);
       }
       
       LOGTRACE(LS_ProbDist+"_GeneralOperatorsFunctionAnalytical",string("lin. pars. set manually\n") + PrintLinPars());
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
            LOGTRACE(LS_ProbDist+"_Convolution4TargetValueThread", string("x in [") + x.first.RawPrint(15) + ", " + x.second.RawPrint(15) + "]" );
        for(auto tv: *pvdfTargetValues)
            LOGTRACE(LS_ProbDist+"_Convolution4TargetValueThread", string("tv = ") + tv.RawPrint(15) );
        
        
        for(int itv=0; itv < pvdfTargetValues->size(); itv++)
        _Convolution4TargetValue(*pOther, *pOperation, (*pvpdfXIntervals), (*pvdfTargetValues)[itv], (*pvdfResult)[itv]);
}
 


void CProbabilityDensityDistribution::_Convolution4TargetValue(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, const VectorPairDFType& vpdfXIntervals, const CDigFloat& dfTargetValue, CDigFloat& dfResult)
{
    
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("called with args:"));
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("operation ") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("dfTargetValue ") + dfTargetValue.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("other distri \n") + Other.PrintMetaInfo());
    
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
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("x in [") + ix.first.first.RawPrint(15) + ", " + ix.first.second.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("y in [") + ix.second.first.RawPrint(15) + ", " + ix.second.second.RawPrint(15));
    }
    
    
    _GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, TotalIntervals4TargetValue, SubIntervals4TargetValue);  
    
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("calculated complementary intervals with unique slope and offset ... "));
    for(auto planEl: SubIntervals4TargetValue)
    {
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... x iterator points to [") + planEl.xDistriIter->first.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... x in [") + planEl.xLimits.first.RawPrint(15) + ", " + planEl.xLimits.second.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... y iterator points to [") + planEl.yDistriIter->first.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... y in [") + planEl.yLimits.first.RawPrint(15) + ", " + planEl.yLimits.second.RawPrint(15));
    }
    
    
    // init the result
    dfResult = 0;
    
    LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("iterating over all sub intervals ...") );
    LOGINFO(LS_ProbDist+"_Convolution4TargetValue",string("no. of sub intervals: ") + to_string(SubIntervals4TargetValue.size()) );

    // now calculate 
    for(auto planEl: SubIntervals4TargetValue)
    {   
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... x iterator points to ")  + planEl.xDistriIter->first.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... x in [") + planEl.xLimits.first.RawPrint(15) + ", " + planEl.xLimits.second.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... y iterator points to ") + planEl.yDistriIter->first.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... y in [") + planEl.yLimits.first.RawPrint(15) + ", " + planEl.yLimits.second.RawPrint(15));

        CDigFloat dfPartialResult;
        switch( Operation)
        {
            case ProbDistOp::pdoPlus:
                    dfPartialResult =  _Integral4TargetValue4Addition(dfTargetValue,Other, planEl);            
                break;
            case ProbDistOp::pdoMinus:
                    dfPartialResult =  _Integral4TargetValue4Subtraction(dfTargetValue, Other, planEl);            
                break;
            case ProbDistOp::pdoMult:
                    dfPartialResult =  _Integral4TargetValue4Multiplication(dfTargetValue, Other, planEl);            
                break;
            case ProbDistOp::pdoDiv:
                    dfPartialResult =  _Integral4TargetValue4Division(dfTargetValue, Other, planEl);            
                break;
            default:
                break;
        }   //endswitch( Operation)
        
        assert(dfPartialResult >= 0);
        
        // error tracing:
        if( dfPartialResult.RawValue() < 0)
        {
            LOGWARN("WarningLogger", LS_ProbDist+"_Convolution4TargetValue: partial integration raw value is negative:");
            LOGWARN("WarningLogger", string("dfPartialResult = ") + dfPartialResult.RawPrint(15));
            
        }   // endif( dfPartialResult.RawValue() < 0)
        
        // TODO: here we forget the error history using the raw value .... do we need it ?
        dfResult += dfPartialResult.RawValue();
        LOGTRACE(LS_ProbDist+"_Convolution4TargetValue",string("... dfResult incremented ") + dfResult.RawPrint(15));
        
        
    }   // endfor(auto interval: SubIntervals4TargetValue)

}

CDigFloat CProbabilityDensityDistribution::_GetConvolutionIntegralAnalytical4Operation(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, const CDigFloat& dfTargetValue, const ConvPlanElement& planEl)
{
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("called with args:"));
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("other distri:\n") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("this distri:\n") + PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("operation   :") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("target value:") + dfTargetValue.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("first value of x-interval: ") + planEl.xDistriIter->first.RawPrint(15) );
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("x-limit: [") + planEl.xLimits.first.RawPrint(15) + ", " + planEl.xLimits.second.RawPrint(15) + " ]" );
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("first value of y-interval: ") + planEl.yDistriIter->first.RawPrint(15) );
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("y-limit: [") + planEl.yLimits.first.RawPrint(15) + ", " + planEl.yLimits.second.RawPrint(15) + " ]" );
    
    // init the result
    CDigFloat dfResult = 0;
    
    switch( Operation)
    {
        case ProbDistOp::pdoPlus:
                dfResult =  _Integral4TargetValue4Addition(dfTargetValue,Other, planEl);            
            break;
        case ProbDistOp::pdoMinus:
                dfResult =  _Integral4TargetValue4Subtraction(dfTargetValue,Other, planEl);            
            break;
        case ProbDistOp::pdoMult:
                dfResult =  _Integral4TargetValue4Multiplication(dfTargetValue,Other, planEl);            
            break;
        case ProbDistOp::pdoDiv:
                dfResult =  _Integral4TargetValue4Division(dfTargetValue,Other, planEl);            
            break;
        default:
            break;
    }   //endswitch( Operation)
    
    
    // trace 
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation", "convolution integral is :");
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("result = ") + dfResult.RawPrint(15));
    
    // error case
    if(dfResult < 0)
    {
        LOGERROR("ErrorLogger","CProbabilityDensityDistribution::_GetConvolutionIntegralAnalytical4Operation: convolution integral is negative :");
        LOGERROR("ErrorLogger",string("result = ") + dfResult.RawPrint(15));
        LOGERROR("ErrorLogger",string("operator = ") + GetProbDistOpAsString(Operation));
        LOGERROR("ErrorLogger",string("lin. pars. this ... \n") + PrintLinPars());
        LOGERROR("ErrorLogger",string("lin. pars. other ... \n") + Other.PrintLinPars());
            
        LOGERROR("ErrorLogger",string("preparation of convolution integral calculation"));
        LOGERROR("ErrorLogger",string("X integration interval= [") + planEl.xLimits.first.RawPrint(15) + " ; " + planEl.xLimits.second.RawPrint(15) + " ]");
        LOGERROR("ErrorLogger",string("enclosing X-distri interval first value = ") + planEl.xDistriIter->first.RawPrint(15));
        LOGERROR("ErrorLogger",string("dfXOffset = ") + Offset(planEl.xDistriIter).RawPrint(15));
        LOGERROR("ErrorLogger",string("dfXSlope = ") + Slope(planEl.xDistriIter).RawPrint(15) );
        LOGERROR("ErrorLogger",string("Y integration interval= [") + planEl.yLimits.first.RawPrint(15) + " ; " + planEl.yLimits.second.RawPrint(15) + " ]");
        LOGERROR("ErrorLogger",string("enclosing Y-distri interval first value = ") + planEl.yDistriIter->first.RawPrint(15));
        LOGERROR("ErrorLogger",string("dfYOffset = ") + Offset(planEl.yDistriIter).RawPrint(15));
        LOGERROR("ErrorLogger",string("dfYSlope = ") + Slope(planEl.yDistriIter).RawPrint(15) );
   
    }
    
    
    LOGTRACE(LS_ProbDist+"_GetConvolutionIntegralAnalytical4Operation",string("************************************************************"));
    assert(dfResult >= 0);
    return dfResult;
    
}

void CProbabilityDensityDistribution::_PrepareConvolution(CProbabilityDensityDistribution& Other, const ProbDistOp Operation, vector<CDigFloat>& vdfTargetValues, VectorPairDFType& vpdfXIntervals, CDigFloat& dfTargetValueStep)
{
    LOGTRACE(LS_ProbDist+"_PrepareConvolution", string("called with args:"));
    LOGTRACE(LS_ProbDist+"_PrepareConvolution", string("probdist other :\n") + Other.PrintMetaInfo());
    LOGTRACE(LS_ProbDist+"_PrepareConvolution", string("operation      :") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_PrepareConvolution", string("probdist this  :\n") + PrintMetaInfo());
    
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
        LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("x in [") + ix.first.RawPrint(15) + ", " + ix.second.RawPrint(15));
        
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
    LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("setting target values : iteration steps ") + to_string(IntegrationSteps()));
    dfTargetValueStep = (dfTargetEnd - dfTargetStart)/(CDigFloat(IntegrationSteps()));
    for(int itv = 0; itv < IntegrationSteps(); itv++)
        // avoid zero in case of multiplication / division
        vdfTargetValues.push_back(dfTargetStart+itv*dfTargetValueStep);
            
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
                *itTv = (dfTargetValueStep/2.);
            else
                // ... subtract half a step
                *itTv -= (dfTargetValueStep/2.);
        
        }   // endif( itTv != vdfTargetValues.end())
    
    }   // endif( Operation == ProbDistOp::pdoDiv || Operation == ProbDistOp::pdoMult)
    
    
    // trace 
    LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("calculated:"));    
    LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("target range = [") + dfTargetStart.RawPrint(15) + ", " + dfTargetEnd.RawPrint(15) + " ]");
    LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("step size = ") + dfTargetValueStep.RawPrint(15) + "( " + to_string(IntegrationSteps()) + " ) ");
    LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("target values ...") );
    for(auto iel: vdfTargetValues)
            LOGTRACE(LS_ProbDist+"_PrepareConvolution",string("... ") + iel.RawPrint(15) );
    
    LOGTRACE(LS_ProbDist+"_PrepareConvolution", string("********************************************************************"));
}
/**
 *@details The sub integration intervals are needed to differentiate parts of the total integration interval (see CProbabilityDensityDistribution::_GetTotalIntegrationInterval4TargetValue) where the linear function for x or y distribution changes.
 * 
 */
 void CProbabilityDensityDistribution::_SetSubIntegrationIntervals4TargetValue(const ProbDistOp Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other)
{
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("called with args:"));
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("operation ") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("dfTargetValue ") + dfTargetValue.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("other distri \n") + Other.PrintMetaInfo());
    
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
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("x in [") + ix.first.RawPrint(15) + ", " + ix.second.RawPrint(15));
        
    // I 2. derive closed intervals for x and y (kept in TotalIntervals4TargetValue) mapped onto each other bijectively
    //      by the calculation of the complementary variable
//     _GetTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, vpdXIntervals, TotalIntervals4TargetValue);
    _GetTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, vpdXIntervals, TotalIntervals4TargetValue);
    
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("calculated complementary intervals with mixed slope and offset: "));
    for(auto ix: TotalIntervals4TargetValue)
    {
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("x in [") + ix.first.first.RawPrint(15) + ", " + ix.first.second.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("y in [") + ix.second.first.RawPrint(15) + ", " + ix.second.second.RawPrint(15));
    }
    
    
    _GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue(Operation, dfTargetValue, Other, TotalIntervals4TargetValue, SubIntervals4TargetValue);  
    
    LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("calculated complementary intervals with unique slope and offset ... "));
    for(auto ix: SubIntervals4TargetValue)
    {
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("... x iterator points to [") + ix.xDistriIter->first.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("... x in [") + ix.xLimits.first.RawPrint(15) + ", " + ix.xLimits.second.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("... y iterator points to [") + ix.yDistriIter->first.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_SetSubIntegrationIntervals4TargetValue",string("... y in [") + ix.yLimits.first.RawPrint(15) + ", " + ix.yLimits.second.RawPrint(15));
    }
    
    // now add the sub integration limits and its corresponding target value as pair
    m_ConvolutionPlan.push_back(pair<CDigFloat, VectorConvPlanElementType>(dfTargetValue, SubIntervals4TargetValue) ); 
}

CDigFloat CProbabilityDensityDistribution::_Integral4TargetValue4Addition(const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const ConvPlanElement& planEl)
{    
    CDigFloat dfX1 = planEl.xLimits.first;
    CDigFloat dfX2 = planEl.xLimits.second;
 
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", " called with args:");
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("target value =") + dfTargetValue.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition",string("x iterator points to [") + planEl.xDistriIter->first.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition",string("x in [") + dfX1.RawPrint(15) + ", " + dfX2.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition",string("y iterator points to [") + planEl.yDistriIter->first.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition",string("y in [") + planEl.yLimits.first.RawPrint(15) + ", " + planEl.yLimits.second.RawPrint(15));


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
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("YConst :") + dfYConst.RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("Diff1:") + dfDiff1.RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("Diff2 :") + dfDiff2.RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", string("Diff3 :") + dfDiff3.RawPrint(20) );

    
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Addition", "**************************************************************************");    
    
    return -dfSlopeX*dfSlopeY/3*dfDiff3+
            (dfSlopeX*dfYConst-dfSlopeY*dfOffsetX)/2.* dfDiff2+
            dfYConst*dfOffsetX*dfDiff1;
}

CDigFloat CProbabilityDensityDistribution::_Integral4TargetValue4Subtraction(const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const ConvPlanElement& planEl)
{  
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", " called with args:");
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction",string("x iterator points to [") + planEl.xDistriIter->first.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction",string("x in [") + planEl.xLimits.first.RawPrint(15) + ", " + planEl.xLimits.second.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", string("XSlope  :") + Slope(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction", string("XOffset :") + Offset(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction",string("y iterator points to [") + planEl.yDistriIter->first.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Subtraction",string("y in [") + planEl.yLimits.first.RawPrint(15) + ", " + planEl.yLimits.second.RawPrint(15));
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

CDigFloat CProbabilityDensityDistribution::_Integral4TargetValue4Multiplication(const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const ConvPlanElement& planEl)
{ 
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", " called with args:");
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication",string("target value = ") + dfTargetValue.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication",string("x iterator points to ") + planEl.xDistriIter->first.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication",string("x in [") + planEl.xLimits.first.RawPrint(15) + ", " + planEl.xLimits.second.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("XSlope  :") + Slope(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("XOffset :") + Offset(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication",string("y iterator points to ") + planEl.yDistriIter->first.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication",string("y in [") + planEl.yLimits.first.RawPrint(15) + ", " + planEl.yLimits.second.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("YSlope  :") + Other.Slope(planEl.yDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("YOffset :") + Other.Offset(planEl.yDistriIter).RawPrint(20) );

    MapDFDFType::const_iterator itX = planEl.xDistriIter;
    MapDFDFType::const_iterator itY = planEl.yDistriIter;
    CDigFloat dfX1 = planEl.xLimits.first;
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("dfX1 :") + dfX1.RawPrint(20) );
    CDigFloat dfX2 = planEl.xLimits.second;
    
    
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("dfX2 :") + dfX2.RawPrint(20) );
    CDigFloat dfOffsetX =  Offset(planEl.xDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("XOffset :") + dfOffsetX.RawPrint(20) );
    CDigFloat dfSlopeX =  Slope(planEl.xDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("XSlope  :") + dfSlopeX.RawPrint(20) );
    CDigFloat dfOffsetY =  Other.Offset(planEl.yDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("YOffset :") + dfOffsetY.RawPrint(20) );
    CDigFloat dfSlopeY =  Other.Slope(planEl.yDistriIter);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("YSlope  :") + dfSlopeY.RawPrint(20) );
    CDigFloat dfDiff1 = abs(dfX2) - abs(dfX1);
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("abs(dfX2)-abs(dfX1)  :") + dfDiff1.RawPrint(20) );
    CDigFloat dfLogRel = log( pow(abs(dfX2),sgn(dfX2)) / pow(abs(dfX1), sgn(dfX1)) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("ln(|dfX2|^sgn(dfX2)/|dfX1|^sgn(dfX1)  :") + dfLogRel.RawPrint(20) );
    
    
    // DEBUG: check removing difference 
    
    MapDFDFType::const_iterator itXNext = planEl.xDistriIter;
    itXNext++;
    
    CDigFloat dfPX1 = itX->second;
    CDigFloat dfPX2 = itXNext->second;
    CDigFloat dfRelationFactor = dfX2/(itXNext->first-itX->first) - dfX1/((itXNext->first-itX->first));
    
    CDigFloat dfResult2 = -dfSlopeY*dfTargetValue/abs(dfX1)/abs(dfX2)*(-dfDiff1*dfPX1+(dfPX2-dfPX1)*dfRelationFactor*itX->first);
    CDigFloat dfResult3 = dfOffsetY*(dfPX2-dfPX1)*dfRelationFactor;
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("result1  :") + ((dfSlopeX*dfSlopeY*dfTargetValue + dfOffsetX*dfOffsetY)*dfLogRel).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("result2  :") + (-(dfOffsetX*dfSlopeY*dfTargetValue*(-dfDiff1)/abs(dfX1)/abs(dfX2))).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("result2(alter.) :") + dfResult2.RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("result3  :") + (dfOffsetY*dfSlopeX*dfDiff1).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("result3(alter.) :") + dfResult3.RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", string("result  :") + ((dfSlopeX*dfSlopeY*dfTargetValue + dfOffsetX*dfOffsetY)*dfLogRel -
           dfOffsetX*dfSlopeY*dfTargetValue*(-dfDiff1)/abs(dfX1)/abs(dfX2)+
           dfOffsetY*dfSlopeX*dfDiff1).RawPrint(20) );
    
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Multiplication", "**************************************************************************");
  
    return (dfSlopeX*dfSlopeY*dfTargetValue + dfOffsetX*dfOffsetY)*dfLogRel -
           dfOffsetX*dfSlopeY*dfTargetValue*(-dfDiff1)/abs(dfX1)/abs(dfX2)+
           dfOffsetY*dfSlopeX*dfDiff1;

}

CDigFloat CProbabilityDensityDistribution::_Integral4TargetValue4Division(const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const ConvPlanElement& planEl)
{    
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", " called with args:");
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division",string("x iterator points to [") + planEl.xDistriIter->first.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division",string("x in [") + planEl.xLimits.first.RawPrint(15) + ", " + planEl.xLimits.second.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", string("XSlope  :") + Slope(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division", string("XOffset :") + Offset(planEl.xDistriIter).RawPrint(20) );
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division",string("y iterator points to [") + planEl.yDistriIter->first.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_Integral4TargetValue4Division",string("y in [") + planEl.yLimits.first.RawPrint(15) + ", " + planEl.yLimits.second.RawPrint(15));
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


void CProbabilityDensityDistribution::_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue(const ProbDistOp Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const VectorSubIntLimitsType& TotalIntervals4TargetValue, VectorConvPlanElementType& SubIntervals4TargetValue)
{
    
     
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("called with args:"));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("operation:") + GetProbDistOpAsString(Operation));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("target value:") + dfTargetValue.RawPrint(15));
     LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("other distri:\n") + Other.PrintMetaInfo());
     for(int i = 0; i < TotalIntervals4TargetValue.size(); i++){
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("total interval x[" + to_string(i) + "] = [ ") + TotalIntervals4TargetValue[i].first.first.RawPrint(15) + ", " + TotalIntervals4TargetValue[i].first.second.RawPrint(15) + " ]" );
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("total interval y[" + to_string(i) + "] = [ ") + TotalIntervals4TargetValue[i].second.first.RawPrint(15) + ", " + TotalIntervals4TargetValue[i].second.second.RawPrint(15) + " ]" );
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
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... total interval x = [ ") + xyint.first.first.RawPrint(15) + ", " + xyint.first.second.RawPrint(15) + " ]" );
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... total interval y = [ ") + xyint.second.first.RawPrint(15) + ", " + xyint.second.second.RawPrint(15) + " ]" );

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
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... x = ")+ x.first.RawPrint(15));
        
            // break conditions
            if(x.first > xyint.first.second)
            {
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... leaving loop over x-distri: distri is right apart from limits"));
        
                break;
            }
            
            
            if(x.first > xyint.first.first && x.first <xyint.first.second)
            {
                vdfSubPoints.push_back(x.first);
             
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... adding point = ") + vdfSubPoints.back().RawPrint(15));
   
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
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... x(other) = ")+ x.first.RawPrint(15));
            
            // break conditions
            if(x.first > xyint.second.second)
            {
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... leaving loop over other-distri: distri is right apart from limits"));

                break;
            }

            if(x.first > xyint.second.first && x.first <xyint.second.second)
            {
                vdfSubPoints.push_back(_GetComplementaryVariable(Operation,dfTargetValue,x.first,true));
                LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... adding point = ") + vdfSubPoints.back().RawPrint(15));
            }
            
        }   //endfor(auto x: Other.Distribution())
        
        //////////////////////////////////////////////
        // sort the supporting points and derive
        // sub intervals for x and y
        //////////////////////////////////////////////
        LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("... sorted supporting points ... ") + vdfSubPoints.back().RawPrint(15));
        sort(vdfSubPoints.begin(),vdfSubPoints.end());
        for(auto x: vdfSubPoints)
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... [") +x.RawPrint(15));
        
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
            
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... this interval [") +PlanElement.xLimits.first.RawPrint(15) + ", " +PlanElement.xLimits.second.RawPrint(15));
            
            // in case there is a preceeding element: check if the iterator or the consecutive iterator can be used
            // to avoid binary search            
            if(!_GetInterval(PlanElement.xLimits.first,PlanElement.xDistriIter, uiMaxSearchStepsX))
            {
                // set the iterator x/yDistriIter pointing to the first iterator of the correct interval:
                MapDFDFType::const_iterator itDummy;
                GetInterval(PlanElement.xLimits.first,PlanElement.xDistriIter,itDummy);

            }   // endif(bDoBinarySearch)
            
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... interval iterator x points to ") +PlanElement.xDistriIter->first.RawPrint(15) );
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... xOffset(pre calc.)  = ") + Offset(PlanElement.xDistriIter).RawPrint(15) );
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... xSlope(pre calc.)  = ") + Slope(PlanElement.xDistriIter).RawPrint(15) );
            
            // derive the y limits
            PlanElement.yLimits.first = _GetComplementaryVariable(Operation, dfTargetValue,PlanElement.xLimits.first,false);
            PlanElement.yLimits.second = _GetComplementaryVariable(Operation, dfTargetValue,PlanElement.xLimits.second,false);

            
            // swapping here avoids a lot of if-cases:
            // the if cases would have to consider the sign of target value, x value, y value, and operation
            if(PlanElement.yLimits.second < PlanElement.yLimits.first )
                swap(PlanElement.yLimits.first, PlanElement.yLimits.second); 
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... other interval [") +PlanElement.xLimits.first.RawPrint(15) + ", " +PlanElement.xLimits.second.RawPrint(15));
            
            if(!_GetInterval(PlanElement.yLimits.first,PlanElement.yDistriIter, uiMaxSearchStepsY))
            {
                // set the iterator x/yDistriIter pointing to the first iterator of the correct interval:
                MapDFDFType::const_iterator itDummy;
                Other.GetInterval(PlanElement.yLimits.first,PlanElement.yDistriIter,itDummy);

            }   // endif(!_GetInterval(PlanElement.yLimits.first,PlanElement.yDistriIter, uiMaxSearchStepsY))
            
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... interval iterator y points to ") +PlanElement.yDistriIter->first.RawPrint(15));
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... yOffset(pre calc.)  = ") + Other.Offset(PlanElement.yDistriIter).RawPrint(15) );
            LOGTRACE(LS_ProbDist+"_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue",string("...... ySlope(pre calc.)  = ") + Other.Slope(PlanElement.yDistriIter).RawPrint(15) );
            
            
            // do some checks:
            // first iterator must be less equal left side of limits, i.e. integration interval
            if( PlanElement.xDistriIter->first > PlanElement.xLimits.first ||
                PlanElement.yDistriIter->first > PlanElement.yLimits.first 
            )
            {
                LOGERROR("ErrorLogger", "*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*");
                LOGERROR("ErrorLogger", LS_ProbDist+"::_GetSubIntegrationIntervalFromTotalIntegrationInterval4TargetValue:"); 
                LOGERROR("ErrorLogger","iterator points outside integration sub interval:");
                LOGERROR("ErrorLogger",string("this interval [") + PlanElement.xLimits.first.RawPrint(15) + ", " + PlanElement.xLimits.second.RawPrint(15));                
                LOGERROR("ErrorLogger",string("interval iterator x points to ") +PlanElement.xDistriIter->first.RawPrint(15) );                
                LOGERROR("ErrorLogger",string("other interval [") +PlanElement.yLimits.first.RawPrint(15) + ", " +PlanElement.yLimits.second.RawPrint(15));
                LOGERROR("ErrorLogger",string("interval iterator y points to ") +PlanElement.yDistriIter->first.RawPrint(15) );
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
        // do we have a bipolar interval ..
        if(firstVariable() < 0 && lastVariable() >= 0)
        {
            // finish the negative unipolar interval with something close to zero ...
            pdXInterval.second = -CDigFloat(1).Error()*10;
            
            // ... add it
            vXIntervals.push_back(pdXInterval);
            
            // ... start the second interval with something positive close to zero
            pdXInterval.first = - pdXInterval.second;
            
        }   //endif(firstVariable() < 0 && lastVariable() >= 0)
        
        // in case of a unipolar positive interval starting with zero:
        // set it to something above zero
        if(firstVariable() == 0)
            pdXInterval.first = CDigFloat(1).Error()*10;
        
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
//         cout << "[ " << iint.first.RawPrint(15) << " ; " << iint.second.RawPrint(15) << " ]" << endl;
//     cout << endl;

}

void CProbabilityDensityDistribution:: _GetTotalIntegrationInterval4TargetValue(const ProbDistOp& Operation, const CDigFloat& dfTargetValue, CProbabilityDensityDistribution& Other, const VectorPairDFType& vpdXIntervals, VectorSubIntLimitsType& TotalIntervals4TargetValue)
{
    LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("called with args:"));
    LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("operation:") + GetProbDistOpAsString(Operation));
    LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("target value :") + dfTargetValue.RawPrint(15));
    LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("other distri :\n") + Other.PrintMetaInfo());
    for(int i = 0; i<vpdXIntervals.size(); i++)
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("x-interval[") + to_string(i) + "] = [ " + vpdXIntervals[i].first.RawPrint(15) + ", " + vpdXIntervals[i].second.RawPrint(15) + " ]") ;
    
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
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("interval = [ " + xint.first.RawPrint(15) + ", " + xint.second.RawPrint(15) + " ]") );
        
        /////////////////////////////////////////////////
        // calculate the complementary interval YComp from x
        /////////////////////////////////////////////////
        CDigFloat dfYCompMin, dfYCompMax;
        dfYCompMin = _GetComplementaryVariable(Operation, dfTargetValue,xint.first,false);
        dfYCompMax = _GetComplementaryVariable(Operation, dfTargetValue,xint.second,false);
        
        
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("calculated complementary Y:"));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYCompMin:") + dfYCompMin.RawPrint(15));
        LOGTRACE(LS_ProbDist+"::_GetTotalIntegrationInterval4TargetValue",string("dfYCompMax:") + dfYCompMax.RawPrint(15));
        
        // swapping here avoids a lot of if-cases:
        // the if cases would have to consider the sign of target value, x value, y value, and operation
        if(dfYCompMax < dfYCompMin)
            swap(dfYCompMin, dfYCompMax);       
        
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("after swap:"));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYCompMin:") + dfYCompMin.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYCompMax:") + dfYCompMax.RawPrint(15));
                
        /////////////////////////////////////////////////
        // calculate the overlap of the complementary 
        // ycomp and the original y interval --> YOverlap
        /////////////////////////////////////////////////
        CDigFloat dfYOverlapMin, dfYOverlapMax;
        dfYOverlapMin = max(dfYCompMin,Other.firstVariable());
        dfYOverlapMax = min(dfYCompMax,Other.lastVariable());  
        
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("calculated overlapping Y:"));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYOverlapMin:") + dfYOverlapMin.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYOverlapMax:") + dfYOverlapMax.RawPrint(15));
        
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
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfXCompMin:") + dfXCompMin.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfXCompMax:") + dfXCompMax.RawPrint(15));
        
        // swapping here avoids a lot of if-cases:
        // the if cases would have to consider the sign of target value, x value, y value, and operation
        if(dfXCompMax < dfXCompMin)
            swap(dfXCompMin, dfXCompMax);    
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("after swap:"));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfXCompMin:") + dfXCompMin.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfXCompMax:") + dfXCompMax.RawPrint(15));
        
        /////////////////////////////////////////////////
        // calculate the overlap of the complementary 
        // Xcomp and the original x interval --> XOverlap 
        /////////////////////////////////////////////////
        CDigFloat dfXOverlapMin, dfXOverlapMax;
        dfXOverlapMin = max(dfXCompMin,firstVariable());
        dfXOverlapMax = min(dfXCompMax,lastVariable());
        
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("calculated overlapping Y:"));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYOverlapMin:") + dfYOverlapMin.RawPrint(15));
        LOGTRACE(LS_ProbDist+"_GetTotalIntegrationInterval4TargetValue",string("dfYOverlapMax:") + dfYOverlapMax.RawPrint(15));
   
        
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
//     cout << "This X variable = [" << ThisStart.RawPrint(15) << " ; " << ThisEnd.RawPrint(15) << endl;
//     cout << "Other Y variable = [" << OtherStart.RawPrint(15) << " ; " << OtherEnd.RawPrint(15) << endl;
    
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
//     cout << "TargetRangeStart = " << TargetRangeStart.RawPrint(15) << endl;
//     cout << "TargetRangeEnd = " << TargetRangeEnd.RawPrint(15) << endl;
}

std::__cxx11::string CProbabilityDensityDistribution::PrintMetaInfo() const
{
    return CDistribution::PrintMetaInfo() + "\nIntegrationSteps = " + to_string(IntegrationSteps()) +
    "\nMaxBinarySearchIterations = " + to_string(MaxBinarySearchIterations()) +
    "\naddress = " + Address2String(this);
}


////////////////////////////////////////////////////////
// external functions
////////////////////////////////////////////////////////
CProbabilityDensityDistribution operator+(CProbabilityDensityDistribution pdOne, CProbabilityDensityDistribution pdOther)
{
    CProbabilityDensityDistribution pdResult(pdOne);
    pdResult +=pdOther;
    return pdResult;
}

CProbabilityDensityDistribution operator-(CProbabilityDensityDistribution pdOne, CProbabilityDensityDistribution pdOther)
{
    CProbabilityDensityDistribution pdResult(pdOne);
    pdResult -=pdOther;
    return pdResult;
}

CProbabilityDensityDistribution operator*(CProbabilityDensityDistribution pdOne, CProbabilityDensityDistribution pdOther)
{
    CProbabilityDensityDistribution pdResult(pdOne);
    pdResult *=pdOther;
    return pdResult;
}

CProbabilityDensityDistribution operator/(CProbabilityDensityDistribution pdOne, CProbabilityDensityDistribution pdOther)
{
    CProbabilityDensityDistribution pdResult(pdOne);
    pdResult /=pdOther;
    return pdResult;
}

CProbabilityDensityDistribution operator*(const CDigFloat& dfValue, CProbabilityDensityDistribution& pdDistri)
{    
    LOGTRACE(LS_ProbDist+"FloatMultProb",string("is called with args:"));
    LOGTRACE(LS_ProbDist+"FloatMultProb",string("dfValue = ") + dfValue.RawPrint(15));
    LOGTRACE(LS_ProbDist+"FloatMultProb",string("distri = \n") + pdDistri.Print(10));
    
    CProbabilityDensityDistribution pdResult(pdDistri);
    pdResult.Scale(dfValue);
    
    LOGTRACE(LS_ProbDist+"FloatMultProb",string("****************************************************************"));
    return pdResult;

}

CProbabilityDensityDistribution operator/(const CDigFloat& dfValue, CProbabilityDensityDistribution& pdDistri)
{
    
    LOGTRACE(LS_ProbDist+"FloatDivProb",string("is called with args:"));
    LOGTRACE(LS_ProbDist+"FloatDivProb",string("dfValue = ") + dfValue.RawPrint(15));
    LOGTRACE(LS_ProbDist+"FloatDivProb",string("distri = \n") + pdDistri.Print(10));
    
    // be sure of unipolarity:
    assert(pdDistri.firstVariable().RawValue()*pdDistri.lastVariable().RawValue() >= 0);

    // hand over denormalized values
    CDigFloat dfOrigIntegral = pdDistri.OrigIntegral();
    LOGTRACE(LS_ProbDist+"FloatDivProb",string("dfOrigIntegral = \n") + dfOrigIntegral.RawPrint(15));
    CProbabilityDensityDistribution pdResult;
    for(auto iel:pdDistri.Distribution())
        pdResult.Add(dfValue/iel.first,iel.second*dfOrigIntegral);
    
    pdResult._SetPeripherialMember(pdDistri);
    
        
    pdResult.Normalize();
    
    CDigFloat dfProbFactor = dfOrigIntegral/pdResult.OrigIntegral();
    
    LOGTRACE(LS_ProbDist+"FloatDivProb",string("dfProbFactor = \n") + dfProbFactor.RawPrint(15));
    LOGTRACE(LS_ProbDist+"FloatDivProb",string("target distri = \n") + pdResult.Print(10));
    
    // DEBUG
//     pdResult._ResetLinPars();
//     pdResult._SetLinPars();
//     pdResult._SetPeripherialMember(pdDistri);
//    LOGTRACE(LS_ProbDist+"FloatDivProb",string("lin. pars. set automatically\n") + pdResult.PrintLinPars());
   
   /////////////////////////////////////////////////
   // now set the linear parameters: 
   // reuse the parameters from the original distribution pdDistri
   // to minimize error
   /////////////////////////////////////////////////
   // reset linear parameters for result distri
   pdResult._ResetLinPars();
   
   // are the target intervals all inverted :
   // reciproc negative intervals are not when the factor dfValue is positive
   bool bTargetIntervalInverted = (dfValue * pdDistri.firstVariable()) > 0;
   LOGTRACE(LS_ProbDist+"FloatDivProb",string("target intervals are inverted: ") + Bool2String(bTargetIntervalInverted));
   
   // ielSource points to the first point of pdDistri 
   MapDFDFType::const_iterator ielTarget; 
   MapDFDFType::const_iterator ielSourceNext;
   
   LOGTRACE(LS_ProbDist+"FloatDivProb",string("iterating over source distri ..."));
   for(MapDFDFType::const_iterator ielSource = pdDistri.Distribution().begin(); ielSource != pdDistri.Distribution().end(); ielSource++)
   {
    
       // get the iterator to the next element 
       ielSourceNext = ielSource;
       ielSourceNext++;
       
       if(ielSourceNext == pdDistri.Distribution().end())
           ielTarget = pdResult.back();
       else 
           // find the target iterator
           if(bTargetIntervalInverted)
               ielTarget = pdResult.Distribution().find(dfValue/ielSourceNext->first);
           else
               ielTarget = pdResult.Distribution().find(dfValue/ielSource->first);
       
       LOGTRACE(LS_ProbDist+"FloatDivProb",string("... ielSource = (") + ielSource->first.RawPrint(15) + ", " + ielSource->second.RawPrint(15) );
       LOGTRACE(LS_ProbDist+"FloatDivProb",string("... ielSourceNext = (") + ielSourceNext->first.RawPrint(15) + ", " + ielSourceNext->second.RawPrint(15) );
       LOGTRACE(LS_ProbDist+"FloatDivProb",string("... ielTarget = (") + ielTarget->first.RawPrint(15) + ", " + ielTarget->second.RawPrint(15) );
       
       if( ielTarget != pdResult.back() )
       {
           
           LOGTRACE(LS_ProbDist+"FloatDivProb",string("... not at end of target."));
           // set the new slope and offset depending on the new target interval being inverted or not
           if(bTargetIntervalInverted)
           {
               // inverted:
               // new slope = old slope / value * old x1 * old x2
               // new offset = new prob ( x1) - new slope * new x1
               //            = old prob ( old x1) - new slopw * value / old x1
               //            = old prob ( old x1) - old slope * old x2
               pdResult._AddLinPar(ielTarget,
                                   ielSourceNext->second*dfProbFactor-pdDistri.Slope(ielSource)*ielSource->first,
                                   -pdDistri.Slope(ielSource)/dfValue*ielSource->first*ielSourceNext->first*dfProbFactor);
           
               
           }   //endif(bTargetIntervalInverted)
           else
           {
               // not inverted:
               // new slope = -old slope / value * old x1 * old x2
               // new offset = new prob ( x1) - new slope * new x1
               //            = old prob ( old x1) - new slopw * value / old x1
               //            = old prob ( old x1) - old slope * old x2
               pdResult._AddLinPar(ielTarget,
                                   ielSource->second*dfProbFactor+pdDistri.Slope(ielSource)*ielSourceNext->first,
                                   -pdDistri.Slope(ielSource)/dfValue*ielSource->first*ielSourceNext->first*dfProbFactor);
               
           }   //endelseif(bTargetIntervalInverted)
       
       }   //endif( ielTarget != pdResult.back() )
       else
       {
           
           LOGTRACE(LS_ProbDist+"FloatDivProb",string("... at end of target."));
           // at the end set the final point with zero offset and slope        
           pdResult._AddLinPar(ielTarget,0,0);
           
       }   // endelseif( ielTarget != pdResult.back() )
       
   }   // endfor(MapDFDFType::const_iterator ielSource = pdDistri.Distribution().begin(); ielSource != pdDistri.Distribution().end(); ielSource++)
   
   // do we have all elements
   assert(pdResult.m_LinPars.size() == pdResult.Distribution().size());
   
   LOGTRACE(LS_ProbDist+"FloatDivProb",string("lin. pars. set manually\n") + pdResult.PrintLinPars());
    LOGTRACE(LS_ProbDist+"FloatDivProb",string("****************************************************************"));
   
    return pdResult;
}

CProbabilityDensityDistribution operator+(const CDigFloat& dfValue, CProbabilityDensityDistribution pdDistri)
{
    
    LOGTRACE(LS_ProbDist+"FloatPlusProb",string("is called with args:"));
    LOGTRACE(LS_ProbDist+"FloatPlusProb",string("dfValue = ") + dfValue.RawPrint(15));
    LOGTRACE(LS_ProbDist+"FloatPlusProb",string("distri = \n") + pdDistri.Print(10));
    
    CProbabilityDensityDistribution pdResult(pdDistri);
    pdResult.Shift(dfValue);
    
    LOGTRACE(LS_ProbDist+"FloatPlusProb",string("****************************************************************"));
    
    return pdResult;    
}

CProbabilityDensityDistribution operator-(const CDigFloat& dfValue, CProbabilityDensityDistribution pdDistri)
{
    
    LOGTRACE(LS_ProbDist+"FloatMinusProb",string("is called with args:"));
    LOGTRACE(LS_ProbDist+"FloatMinusProb",string("dfValue = ") + dfValue.RawPrint(15));
    LOGTRACE(LS_ProbDist+"FloatMinusProb",string("distri = \n") + pdDistri.Print(10));
    
    CProbabilityDensityDistribution pdResult;
    for(auto iel:pdDistri.Distribution())
        pdResult.Add(dfValue-iel.first,iel.second);
    
    pdResult._SetPeripherialMember(pdDistri);
    pdResult.Normalize();
    
    // DEBUG
//     pdResult._ResetLinPars();
//     pdResult._SetLinPars(); 
//     pdResult._SetPeripherialMember(pdDistri);
//     LOGTRACE(LS_ProbDist+"FloatMinusProb",string("lin. pars. set automatically\n") + pdResult.PrintLinPars());
     
    /////////////////////////////////////////////////
    // now set the linear parameters: 
    // reuse the parameters from the original distribution pdDistri
    // to minimize error
    /////////////////////////////////////////////////
    
    pdResult._ResetLinPars();
    // ielSource points to the last point of pdDistri which has 0 offset and slope
    MapDFDFType::const_iterator ielSource = pdDistri.Distribution().end();
    ielSource--;
    
    LOGTRACE(LS_ProbDist+"FloatMinusProb",string("iterating over result distri ..."));
    for(MapDFDFType::const_iterator iel = pdResult.Distribution().begin(); iel != pdDistri.Distribution().end(); iel++)
    {
        ielSource--;
        LOGTRACE(LS_ProbDist+"FloatMinusProb",string("... iel = (") + iel->first.RawPrint(15) + ", " + iel->second.RawPrint(15) );
        LOGTRACE(LS_ProbDist+"FloatMinusProb",string("... ielSource = (") + ielSource->first.RawPrint(15) + ", " + ielSource->second.RawPrint(15) );
        
        pdResult._AddLinPar(iel, iel->second +pdDistri.Slope(ielSource)*iel->first , -pdDistri.Slope(ielSource));
        
        // at the end set the final point with zero offset and slope
        if( ielSource == pdDistri.Distribution().begin() )
        {   
            // increment this iterator 
            iel++;
            LOGTRACE(LS_ProbDist+"FloatMinusProb",string("... ielSource is at start: incrementing iel "));
            LOGTRACE(LS_ProbDist+"FloatMinusProb",string("... iel = (") + iel->first.RawPrint(15) + ", " + iel->second.RawPrint(15) );
            pdResult._AddLinPar(iel,0,0);
            break;
            
        }   // endif( ielSource == pdDistri.Distribution().begin() )
        
    }   //endfor(MapDFDFType::const_iterator iel = pdDistri.Distribution().begin(); iel != pdDistri.Distribution().end(); iel++)
    
    LOGTRACE(LS_ProbDist+"FloatMinusProb",string("lin. pars. set manually\n") + pdResult.PrintLinPars());
    LOGTRACE(LS_ProbDist+"FloatMinusProb",string("****************************************************************"));
    
    return pdResult;    
}

// CProbabilityDensityDistribution operator+(CProbabilityDensityDistribution pdOne, CProbabilityDensityDistribution pdOther)
// {
//     CProbabilityDensityDistribution pdResult(pdOne);
//     pdResult += pdOther;
//     return pdResult;
// }

CProbabilityDensityDistribution pow(CProbabilityDensityDistribution pdDistri, const int& nExponent, bool bUsePos /*= true*/)
{
    LOGTRACE(LS_ProbDist+"pow",string("is called with args:"));
    LOGTRACE(LS_ProbDist+"pow",string("dfValue = ") + dfValue.RawPrint(15));
    LOGTRACE(LS_ProbDist+"pow",string("distri = \n") + pdDistri.Print(10));
    LOGTRACE(LS_ProbDist+"pow",string("pos flag = \n") + Bool2String(bUsePos));
    
    // hand over original distribution:
    CProbabilityDensityDistribution pdResult(pdDistri);
    
    // handle non-lin trafo:
    // check if there is already a non-linear trafo of the copied distri
    // depending on that the 
    assert(pdResult.m_NonLinTrafoVariables.size() == pdResult.Distribution().size());
    CDigFloat dfVariable = pdResult.m_NonLinTrafoVariables[pdResult.Distribution().begin()];
    for(MapDFDFType::const_iterator iel = pdResult.Distribution().begin(); iel != pdResult.Distribution().end(); iel++)
            pdResult.m_NonLinTrafoVariables[iel]=pow(pdResult._Variable(iel),nExponent);
    
    
//     CDigFloat dfRelShift = 0;
//     for(MapDFDFType::const_iterator iel = pdDistri.Distribution().begin(); iel != pdDistri.Distribution().end(); iel++)
//     {
//         MapDFDFType::const_iterator ielNext = iel;
//         ielNext++;
//         CDigFloat dfVar = pow(iel->first, nExponent);
//         CDigFloat dfProb = iel->second;
//         if(ielNext != pdDistri.Distribution().end())
//         {
//             dfVar += (ielNext->first - iel->first)*dfRelShift;
//             dfProb = pdDistri.DistriValue(dfVar);
//             
//         }   //endif(ielNext != pdDistri.Distribution().end())
//         
//         pdResult.Add(pow(dfVar,nExponent),dfProb);
//         
//     }   //endfor(MapDFDFType::const_iterator iel = pdDistri.Distribution().begin(); iel != pdDistri.Distribution().end(); iel++)

    pdResult._SetPeripherialMember(pdDistri);
    pdResult.Normalize();
    pdResult._SetLinPars();
    
    return pdResult;
    
}

CProbabilityDensityDistribution cut(CProbabilityDensityDistribution pdDistri, const CDigFloat& dfStartUser, const CDigFloat& dfEndUser)
{
    CProbabilityDensityDistribution pdResult;
    
    // use correct order
    CDigFloat dfStart = dfStartUser;
    CDigFloat dfEnd = dfEndUser;
    if( dfStart > dfEnd)
        swap(dfStart, dfEnd);
    
    // default init: no overlap of intervals of this pdDistri and user given interval [dfStart, dfEnd] 
    MapDFDFType::const_iterator itStart = pdDistri.Distribution().end();
    MapDFDFType::const_iterator itEnd = pdDistri.Distribution().end();
    
    // set the start iterater being >= user given end of cut interval
    if(dfStart < pdDistri.lastVariable())
    {
        itStart = pdDistri.Distribution().begin();
        while(itStart->first < dfStart && itStart != itEnd)
            itStart++;
    }
    // set the end iterater being <= user given end of cut interval
    if(dfEnd < pdDistri.lastVariable())
    {
        // prepare: start with last variable of distri
        itEnd--;
        while(itEnd->first > dfEnd && itEnd->first != pdDistri.firstVariable())
            itEnd--;
    }   //endif(dfEnd < pdDistri.lastVariable())
    
    VectorPairDFType LinParCopy;
    for(MapDFDFType::const_iterator itCut = itStart; itCut != itEnd; itCut++)
    {
        pdResult.Add(itCut->first,itCut->second);
        LinParCopy.push_back( PairDFType( pdDistri.Offset(itCut),pdDistri.Slope(itCut)));
        
    }   //endfor(MapDFDFType::const_iterator itCut = itStart; itCut != itEnd; itCut++)
        
    // get linear  parameters by copy and paste
    if(LinParCopy.size()>0)
    {
        // init index for linpars:
        unsigned int uiLinParIndex = 0;
        for(MapDFDFType::const_iterator iel=pdResult.Distribution().begin(); iel != pdResult.Distribution().end(); iel++)
        {
            // add offset and slope
            pdResult._AddLinPar(iel,LinParCopy[uiLinParIndex].first, LinParCopy[uiLinParIndex].second);
            
            // increment lin par index 
            uiLinParIndex++;
            
        }   //endfor(MapDFDFType::const_iterator iel=pdResult.Distribution().begin(); iel != pdResult.Distribution().end(); iel++)
        
    }   //endif(LinParCopy.size()>0)
    
    pdResult.Normalize();
    pdResult._SetPeripherialMember(pdDistri);
    
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
