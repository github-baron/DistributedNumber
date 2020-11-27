
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
#include "CartesianCoordinate.h"

/**
 * @brief This class represents a 1D distribution
 */
class CDistribution
{
public:
    /**
     * Default constructor
     */
    CDistribution();

    /**
     * Copy constructor
     *
     * @param other TODO
     */
    CDistribution(const CDistribution& other);

    /**
     * Destructor
     */
    ~CDistribution();

    /**
     * Assignment operator
     *
     * @param other TODO
     * @return TODO
     */
    CDistribution& operator=(const CDistribution& other);

    /**
     * @todo write docs
     *
     * @param other TODO
     * @return TODO
     */
    bool operator==(const CDistribution& other) const;

    /**
     * @todo write docs
     *
     * @param other TODO
     * @return TODO
     */
    bool operator!=(const CDistribution& other) const;

    bool bNormalized;
    map<double,double> mDistribution;
};

#endif // CDISTRIBUTION_H
