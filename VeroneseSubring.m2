--! \file VeroneseSubring.m2
--! \brief Package for computing the n-th Veronese subring presentation of a graded ring.
--!
--! This package provides the function \code{veronesePresentation}, which returns a presentation
--! for the n-th Veronese subring of a graded polynomial or quotient ring. (The n-th Veronese subring is
--! the subring formed by taking all elements whose degrees are divisible by n.)
--!
--! An internal helper function, \code{allProductsUpToDegree}, is used to construct the necessary products.
--!
--! \author Jack Westbrook
--! \date 2025-04-04

-*
 * allProductsUpToDegree(L, n)
 *
 * \remark
 * This function is an internal helper that returns a list of all products of elements in the list L
 * that have total degree exactly equal to n. It is not exported.
 *
 * @param L a nonempty list of homogeneous elements in a graded ring.
 * @param n a positive integer specifying the total degree.
 * @return a list of products of elements in L of total degree n.
 *-


allProductsUpToDegree = (L, n) -> (
    if length L == 0 then error "Error in allProductsUpToDegree: The input list L must not be empty.";
    if not (isANumber n and n > 0) then error "Error in allProductsUpToDegree: n must be a positive integer.";
    
    R = ring L#0;        
    d = for f in L list sum degree f; 
    result = {};      

    generate = (i, degSoFar, prodSoFar) -> (
        if i == #L then (
            result = append(result, prodSoFar);
        ) else (
            maxE = floor((n - degSoFar) / d#i);
            if maxE < 0 then (
                generate(i+1, degSoFar, prodSoFar);
            ) else (
                for e from 0 to maxE do (
                    generate(i+1, degSoFar + e*d#i, prodSoFar * (L#i)^e);
                );
            );
        );
    );
    
    generate(0, 0, 1_R);
    return result;
);

-*
 * veronesePresentation(R, n)
 *
 * Returns a presentation for the n-th Veronese subring of the graded ring R.
 *
 * R must be a graded polynomial ring or graded quotient ring, and n must be a positive integer.
 * The presentation is provided as a pair: a ring map from R to a polynomial ring and the quotient
 * of that polynomial ring by the kernel of the map.
 *
 * \note
 * The n-th Veronese subring is the subring of R consisting of all elements whose degrees are divisible by n.
 *
 * @param R a graded polynomial ring or graded quotient ring.
 * @param n a positive integer representing the Veronese degree.
 * @return a list containing the map f and the presentation of the n-th Veronese subring.
 *-

 
veronesePresentation = (R, n) -> (
    if not (isANumber n and (isQuotientRing R or isPolynomialRing R))
        then error "Error in veronesePresentation: R must be a graded polynomial or quotient ring and n must be a positive integer.";
    
    gensList = gens R;
    r = numgens R;
    K = coefficientRing R;
    
    if r == 0 then return R else (
        tempList = {n};
        
        for j to r-1 do (
            if not (length(degree gensList#j) == 1) then
                error "Error in veronesePresentation: All generators of R must have one-dimensional degrees.";
        );
        
        for j to r-1 do (
            tempList = append(tempList, (degree gensList#j)#0);
        );
        m = lcm(tempList);
        s = m // n;
        
        tempList = {};
        for i from 1 to s do (
            myBasis = basis(i*n, R);
            k = numColumns myBasis;
            for j to k-1 do (
                if (length tempList == 0) then (
                    tempList = append(tempList, myBasis_(0,j));
                ) else (
                    tempElt = myBasis_(0, j);
                    linIndList = allProductsUpToDegree(tempList, (degree tempElt)#0);
                    checker = false;
                    for elt in linIndList do (
                        if tempElt == elt then (
                            checker = true;
                        );
                    );
                    -- Note: A more robust check of linear independence might be preferable,
                    -- but this approach works given the behavior of basis.
                    if checker == false then (
                        tempList = append(tempList, tempElt);
                    );
                );
            );
        );
        counter = length tempList;
        basePolyRing = K[x_1 .. x_counter];
        f = map(R, basePolyRing, tempList);
        I = kernel f;
        return {f, basePolyRing/I};
    );
);

-- Export the main function so it is available when the package is loaded.
export veronesePresentation
