module Aggregation

import List;
import util::Math;

public int count(list[num] nums) = size(nums);

//public real sum(list[real] nums) = ( 0 | it + n | n <- nums);
public num sum(list[num] nums) = ( 0 | it + n | n <- nums);

public real mean(list[num] nums) = toReal(sum(nums)) / count(nums);

public real theilT(list[num] nums) {
	m = mean(nums);
	return 1.0/count(nums) * ( 0.0 | it + (x/m * ln(x/m)) | x <- nums);
}

public real theilL(list[num] nums) {
	m = mean(nums);
	return 1.0/count(nums) * ( 0.0 | it + ln(m/x) | x <- nums);
}

public real generalizedEntropy(list[num] nums, 0.0) = theilL(nums);
public real generalizedEntropy(list[num] nums, 1.0) = theilT(nums);
// alpha usually ranges between 0 and 1
public default real generalizedEntropy(list[num] nums, real alpha) {
	m = mean(nums);
	return 1.0/(count(nums) * alpha * (alpha - 1)) * 
		( 0.0 | it + pow(m/x, alpha) | x <- nums);
}

