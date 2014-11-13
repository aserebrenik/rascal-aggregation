module Aggregation

import IO;
import List;
import Exception;
import util::Math;
import analysis::statistics::Descriptive;

public int count(list[num] nums) = size(nums);

//public real sum(list[real] nums) = ( 0 | it + n | n <- nums);
//public num sum(list[num] nums) = ( 0 | it + n | n <- nums);

public real mean(list[num] nums) = toReal(sum(nums)) / count(nums);

public real coefficientOfVariation(list[num] nums)  = standardDeviation(nums) / mean(nums);

@doc{
The phi function extending the one defined in Mordal, K., Anquetil, N., Laval, J., Serebrenik, A., Vasilescu, B.N. & Ducasse, S. (2013). Software quality metrics aggregation in industry. Journal of Software: Evolution and Process, 25(10), 1117-1135. 
The phi function is used to determine applicability of a number of aggregation indexes as well to control the values of Gini.
}
private bool phi(list[num] nums) {
	if (x <- nums, x < 0) {
		throw ArithmeticException("Negative numbers not supported");
	}
	if (x <- nums, x > 0)
		return true;
	else 
		throw ArithmeticException("At least one number should be positive");	 
}	
private bool atLeastTwo(list[num] nums) {
	if (count(nums) < 2) {
		throw "Undefined for lists shorter then 2";
	}
	return true;
}
@doc{The Gini index, introduced by C. Gini in 1912, ranges between 0 and 1-1/n if the phi-condition holds and over R, otherwise.}
public real gini(list[num] nums) {
	if (atLeastTwo(nums)) {
		n = count(nums);
		s = sum(nums);
		x = (0.0 | it + toReal(abs(xi-xj)) | xi <- nums, xj <- nums);
		return 1.0/(2*n*s) * x;
	}
}

public real theilT(list[num] nums) {
	if (phi(nums) && atLeastTwo(nums)) {
		nums = [x|x<- nums, x > 0]; //otherwise ln will be undefined for x <= 0
		m = mean(nums);
		return 1.0/count(nums) * ( 0.0 | it + (x/m * ln(x/m)) | x <- nums);
	}
}
	
@doc{Theil L, the second Theil's measure is also known as the mean logarithmic deviation (MLD). }
public real MLD(list[num] nums) = theilL(nums);

public real theilL(list[num] nums) {
	if (phi(nums) && atLeastTwo(nums)) {
		nums = [x|x<- nums, x > 0]; //otherwise ln will be undefined for x <= 0
		m = mean(nums);
		return 1.0/count(nums) * ( 0.0 | it + ln(m/x) | x <- nums);
	}
}	

public real generalizedEntropy(list[num] nums, 0.0) = theilL(nums);
public real generalizedEntropy(list[num] nums, 1.0) = theilT(nums);
public default real generalizedEntropy(list[num] nums, real alpha) {
	if (phi(nums) && atLeastTwo(nums)) {
		m = mean(nums);
		return 1.0/(count(nums) * alpha * (alpha - 1)) * 
			( 0.0 | it + (pow(x/m, alpha) - 1) | x <- nums);
	}		
}

@docs{Atkinson's inequality index. We use the same default value for the parameter as R.} 
public real atkinson(list[num] nums, real alpha = 0.5) {
	if (phi(nums) && atLeastTwo(nums)) {
		m = mean(nums);
		n = count(nums);
		if (alpha == 1.0) 
			return 1 - ((1.0 | it * pow(x, toRat(1/n)) | x <- nums)/m); 
		else return 1 - 1.0/m *
				pow((1.0/n * ( 0.0 | it + pow(x, 1.0 - alpha) | x <- nums)),
				1.0/(1.0-alpha));
	}
}
@doc{Hoover (also known as the Ricci–Schutz coefficient, or the Robin Hood index). RS in R.}
public real RS(list[num] nums) = hoover(nums);
public real hoover(list[num] nums) {
	if (atLeastTwo(nums)) {
		m = mean(nums);
		if (m == 0) {
			throw ArithmeticException("Hoover inequality index cannot be computed when the mean is zero.");	 
		}
		s = sum(nums);
		return 1.0/(2*s) * (0.0 | it + abs(xi-m) | xi <- nums);
	}
}
@doc{The default value for the beta parameter in Kolm is the same as in R.
The index has been introduced in Serge-Christophe Kolm: Unequal inequalities I. 
Journal of Economic Theory 12(1976):416–442}
public real kolm(list[num] nums, real beta = 1.0) {
	if (beta <= 0.0) {
		throw ArithmeticException("Kolm index is defined only for beta greater than zero.");
	}
	if (atLeastTwo(nums)) {
		m = mean(nums);
		return 1.0/beta * 
			ln(1.0/count(nums) *
				(0.0 | it + exp(beta*(m-x)) | x<- nums));
	}
}

@doc{The Simpson and Blau diversity indices assume that the values being aggregated are percentages, i.e., 
their sum is 1 and each one of these values is between 0 and 1. We normalize the 
values by dividing them by sum(nums).} 
public real simpson(list[num] nums) {
	if (size(nums) < 1) {
		throw "Diversity indices cannot be computed for empty collections of values"; 
	}		
	if (phi(nums))	{
		s = sum(nums);
		t =  (0.0 | it + x*x | x <- nums);
		return toReal(1.0/(s*s) * t);
	}	
}
public real blau(list[num] nums) = 1.0 - simpson(nums);

@doc{Default values in Squale}
public real squaleSoft(list[num] nums) = squale(nums, 3.0);
public real squaleMedium(list[num] nums) = squale(nums, 9.0);
public real squaleHard(list[num] nums) = squale(nums, 30.0);
public real squale(list[num] nums, real lambda) {
	if (lambda <= 0){
			throw ArithmeticException("Squale inequality index is not defined for lambda less than 0");	 
	}
	if (lambda == 1){
			throw ArithmeticException("Squale inequality index is not defined for lambda == 1.");	 
	}	
	if (atLeastTwo(nums)) {
		t = (0.0 | it + pow(lambda, toReal(-x)) | x <- nums);
		return -log(1.0/count(nums) * t, lambda);
	}
}

// tests
test bool correctTheilT() {
	return theilT([1,2,3]) == 0.087208023941279197605; 
}

test bool theilTNeverTooSmall(list[num] nums) {
	if (size(nums) < 2) return true;
	nums = makeSmallerThan(nums, 500);
	nums = [abs(n) | n <- nums];
	return theilT(nums) >= 0;
}
test bool theilTNeverTooLarge(list[num] nums) {
	if (size(nums) < 2) return true;
	nums = makeSmallerThan(nums, 500);
	nums = [abs(n) | n <- nums];

	return theilT(nums) < ln(count(nums));
}
test bool theilNotSensitiveToZeroes(list[num] nums) {
	if (size(nums) < 2) return true;
	nums = makeSmallerThan(nums, 500);
	nums = [abs(n) | n <- nums];

	return theilT(nums) == theilT(nums+[0,0,0]);
}
test bool giniNeverTooSmall(list[num] nums) {
	if (size(nums) < 2) return true;
	nums = makeSmallerThan(nums, 500);
	nums = [abs(n) | n <- nums];
	return gini(nums) >= 0;
}
test bool giniNeverTooLarge(list[num] nums) {
	if (size(nums) < 2) return true;
	nums = makeSmallerThan(nums, 500);
	nums = [abs(n) | n <- nums];
	return gini(nums) < (1 - 1.0/(count(nums)));
}
test bool hooverIsInRange(list[num] nums) {
	if (size(nums) < 2) return true;
	nums = makeSmallerThan(nums, 500);
	nums = [abs(n) | n <- nums];
	if (phi(nums))
		return (hoover(nums) >= 0 && hoover(nums) <= 1);
}
test bool kolmIsNonneg(list[num] nums) {
	if (size(nums) < 2) return true;
	nums = makeSmallerThan(nums, 500);
	return kolm(nums) >= 0;
}
test bool squaleRangeLo(list[num] nums, real lambda) {
	if (size(nums) < 2) return true;
	if (lambda <= 1) return true;
	nums = makeSmallerThan(nums, 500);
	lambda = makeSmallerThan(lambda, 5);
	return min(nums) <= squale(nums, lambda); 
}
test bool squaleRangeHi(list[num] nums, real lambda) {
	if (size(nums) < 2) return true;
	if (lambda <= 1) return true;
	nums = makeSmallerThan(nums, 500);
	lambda = makeSmallerThan(lambda, 5);
	return mean(nums) >= squale(nums, lambda); 
}
test bool squaleRangeInv(list[num] nums, real lambda, real c) {
	if (size(nums) < 2) return true;
	if (lambda == 1) return true;
	if (lambda < 0) return true;
	nums = makeSmallerThan(nums, 500);
	lambda = makeSmallerThan(lambda, 10);
	c = makeSmallerThan(c, 5);
	cnums = [n+c| n<- nums];
	return eq(squale(cnums, lambda), squale(nums, lambda) + c); 
}
test bool squaleKolm(list[num] nums, real lambda) {
	if (size(nums) < 2) return true;
	if (max(nums) > 500) return true; //otherwise the computation of Kolm can take too long
//	if (lambda == 1) return true;
//	if (lambda < 0) return true;
	lambda = abs(lambda);
	if (lambda <= 1) { // for 0 < lambda < 1 ln(lambda) is negative and Kolm is undefined
		lambda = 1 + lambda;
	}
	nums = makeSmallerThan(nums, 500);
	lambda = makeSmallerThan(lambda, 10);
	return eq(squale(nums, lambda) + kolm(nums, beta = ln(lambda)), mean(nums)); 
}
test bool simpsonLo(list[num] nums) {
	if (size(nums) < 1) return true;
	nums = [abs(n) | n <- nums];
	nums = makeSmallerThan(nums, 500);
	nums = prepareDiversityIndex(nums);
	return simpson(nums) >= 0.0;
}
test bool simpsonHi(list[num] nums) {
	if (size(nums) < 1) return true;
	nums = [abs(n) | n <- nums];
	nums = makeSmallerThan(nums, 500);
	nums = prepareDiversityIndex(nums);
	return simpson(nums) <= 1.0;
}

test bool generalizedEntropyCoV(list[num] nums) {
	if (size(nums) < 1) return true;
	nums = [abs(n) | n <- nums];
	nums = makeSmallerThan(nums, 500);
	cov = coefficientOfVariation(nums);
	return eq(generalizedEntropy(nums, alpha = 2.0), 0.5 * cov * cov);
}


@doc{Diversity indices assume that the values are percentages}
private list[num] prepareDiversityIndex(list[num] nums) {
	if (size(nums) >= 1) {
		s = sum(nums);
		nums = [toReal(n/s) | n <- nums];
		return nums;
	} 
}


private (&T<:int) makeSmallerThan(&T <: int n, int limit) = n % limit;
private (&T<:real) makeSmallerThan(&T <: real n, int limit) {
	if (abs(n) < limit) {
		return n;
	}
	f = trunc(n);
	r = n - f;
	return (f % limit) + r;
}
private (&T<:rat) makeSmallerThan(&T <: rat n, int limit) {
	if (abs(n) < limit) {
		return n;
	}
	return toRat(1, denominator(n));
}


list[num] makeSmallerThan(list[num] nums, int limit) 
	= [ makeSmallerThan(n, limit) | n <- nums];
	
bool eq(num a, num b) {
	error = 1 / pow(10, min(scale(a), scale(b)) - 1);
	return abs(a-b) <= error;
}
bool leq(num a, num b) = a < b ? true : eq(a,b);
	
