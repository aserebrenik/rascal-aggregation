module Aggregation

import IO;
import List;
import Exception;
import util::Math;

public int count(list[num] nums) = size(nums);

//public real sum(list[real] nums) = ( 0 | it + n | n <- nums);
public num sum(list[num] nums) = ( 0 | it + n | n <- nums);

public real mean(list[num] nums) = toReal(sum(nums)) / count(nums);

private bool theilCanBeComputed(list[num] nums) {
	if (count(nums) < 2) {
		throw "Undefined for lists shorter then 2";
	}
	if (x <- nums, x < 0) {
		throw ArithmeticException("Negative numbers not supported");
	}
	if (x <- nums, x > 0)
		return true;
	else 
		throw ArithmeticException("At least one number should be positive");	 
}	
public real theilT(list[num] nums) {
	if (theilCanBeComputed(nums)) {
		m = mean(nums);
		//println([ x/m | x <- nums]);
		return 1.0/count(nums) * ( 0.0 | it + (x/m * ln(x/m)) | x <- nums);
	}
}	

public real theilL(list[num] nums) {
	if (theilCanBeComputed(nums)) {
		m = mean(nums);
		return 1.0/count(nums) * ( 0.0 | it + ln(m/x) | x <- nums);
	}
}	

public real generalizedEntropy(list[num] nums, 0.0) = theilL(nums);
public real generalizedEntropy(list[num] nums, 1.0) = theilT(nums);
// alpha usually ranges between 0 and 1
public default real generalizedEntropy(list[num] nums, real alpha) {
	m = mean(nums);
	return 1.0/(count(nums) * alpha * (alpha - 1)) * 
		( 0.0 | it + pow(m/x, alpha) | x <- nums);
}
test bool correctTheilT() {
	return theilT([1,2,3]) == 0.087208023941279197605; 
}

test bool theilTNeverTooSmall(list[num] nums) {
	if (size(nums) < 2) return true;
	nums = [abs(n) | n <- nums];
	return theilT(nums) >= 0;
}
test bool theilTNeverTooLarge(list[num] nums) {
	if (size(nums) < 2) return true;
	nums = [abs(n) | n <- nums];

	return theilT(nums) < ln(count(nums));
}
