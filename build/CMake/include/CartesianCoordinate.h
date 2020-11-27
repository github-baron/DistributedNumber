/*
 * Copyright (c) 2020 <copyright holder> <email>
 * 
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef CCARTESIANCOORDINATE_H
#define CCARTESIANCOORDINATE_H

// includes
 #include <FloatingMeasure/FloatingMeasure.h>
//#include <Utils/Utils.h>

/**
 * @brief This generic class represents cartesian coordinates
 */

class
#ifdef _WIN32
DistributedNumberDLL_API
#endif
CCartesianCoordinate : public vector<CDigFloat>
{
///////////////////////////////////
// friend functions
///////////////////////////////////
friend CDigFloat ScalarProduct(CCartesianCoordinate& one, CCartesianCoordinate& other);
friend CCartesianCoordinate VectorProduct(CCartesianCoordinate& one, CCartesianCoordinate& other);
    
public:
    ///////////////////////////////////
    // Constructor / Destructor
    ///////////////////////////////////

    /**
     * Default constructor
     */
    CCartesianCoordinate();

    /**
     * Copy constructor
     *
     * @param other CCartesianCoordinate&
     */
    CCartesianCoordinate(const CCartesianCoordinate& other);

    /**
     * Destructor
     */
    ~CCartesianCoordinate();
    
    ///////////////////////////////////
    // operators
    ///////////////////////////////////
  
    /**
     * Assignment operator
     *
     * @param other CCartesianCoordinate&
     * @return CCartesianCoordinate&
     */
    CCartesianCoordinate& operator=(const CCartesianCoordinate& other);

    /**
     * @brief comparison operator : compares ::Coordinates with other
     *
     * @param other CCartesianCoordinate&
     * @return CCartesianCoordinate&
     */
    bool operator==(const CCartesianCoordinate& other) const;

    /**
     * @brief inequality operator: 
     *
     * @param other CCartesianCoordinate&
     * @return bool
     */
    bool operator!=(const CCartesianCoordinate& other) const;
    
    /**
     * @brief operator "<": the highest component of the coordinate rules
     *
     * @param other CCartesianCoordinate&
     * @return bool
     */
    bool operator<(const CCartesianCoordinate& other) const;
    
    /**
     * @brief operator "<=": the highest component of the coordinate rules
     *
     * @param other CCartesianCoordinate&
     * @return bool
     */
    bool operator<=(const CCartesianCoordinate& other) const;
    
    /**
     * @brief operator ">": the highest component of the coordinate rules
     *
     * @param other CCartesianCoordinate&
     * @return bool
     */
    bool operator>(const CCartesianCoordinate& other) const;
    
    /**
     * @brief operator ">=": the highest component of the coordinate rules
     *
     * @param other CCartesianCoordinate&
     * @return bool
     */
    bool operator>=(const CCartesianCoordinate& other) const;

    
    
    ///////////////////////////////////
    // functions
    ///////////////////////////////////
    
    /**
     * @brief get number of dimensions: 
     *
     * @return size_t
     *
    inline size_t size(){return Coordinates.size();}
    */
protected:

    //vector<CDigFloat> Coordinates;
};

/////////////////////////////////////
// external functions
/////////////////////////////////////
   /**
     * @brief scalar product: handles coordinates as origin vectors (starting point is origin) 
     *
     * @param one CCartesianCoordinate&
     * @param other CCartesianCoordinate&
     * @return CDigFloat&
     */
CDigFloat ScalarProduct(CCartesianCoordinate& one, CCartesianCoordinate& other);
 
   /**
     * @brief vector product: handles coordinates as origin vectors (starting point is origin) 
     *
     * @param one CCartesianCoordinate&
     * @param other CCartesianCoordinate&
     * @return CCartesianCoordinate
     */
CCartesianCoordinate VectorProduct(CCartesianCoordinate& one, CCartesianCoordinate& other);

#endif // CCARTESIANCOORDINATE_H
