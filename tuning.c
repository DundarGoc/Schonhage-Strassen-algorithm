#include "include/SSA.h"

tuningModulo_t tuningModulo[] =
{
	{ // bits = 0
	},
	{ // bits = 1
	},
	{ // bits = 2
		132, // KS1 -> KS2 multiplication threshold
		1127, // KS2 -> KS4 multiplication threshold
		10000000, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 3
		132, // KS1 -> KS2 multiplication threshold
		2755, // KS2 -> KS4 multiplication threshold
		10000000, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 4
		80, // KS1 -> KS2 multiplication threshold
		1053, // KS2 -> KS4 multiplication threshold
		10000000, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 5
		86, // KS1 -> KS2 multiplication threshold
		1077, // KS2 -> KS4 multiplication threshold
		32051, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 6
		78, // KS1 -> KS2 multiplication threshold
		985, // KS2 -> KS4 multiplication threshold
		10000000, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 7
		70, // KS1 -> KS2 multiplication threshold
		412, // KS2 -> KS4 multiplication threshold
		20868, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 8
		56, // KS1 -> KS2 multiplication threshold
		385, // KS2 -> KS4 multiplication threshold
		20868, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 9
		57, // KS1 -> KS2 multiplication threshold
		264, // KS2 -> KS4 multiplication threshold
		17119, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 10
		62, // KS1 -> KS2 multiplication threshold
		173, // KS2 -> KS4 multiplication threshold
		15505, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 11
		47, // KS1 -> KS2 multiplication threshold
		158, // KS2 -> KS4 multiplication threshold
		15505, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 12
		43, // KS1 -> KS2 multiplication threshold
		144, // KS2 -> KS4 multiplication threshold
		12307, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 13
		43, // KS1 -> KS2 multiplication threshold
		144, // KS2 -> KS4 multiplication threshold
		14044, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 14
		48, // KS1 -> KS2 multiplication threshold
		101, // KS2 -> KS4 multiplication threshold
		10785, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 15
		38, // KS1 -> KS2 multiplication threshold
		173, // KS2 -> KS4 multiplication threshold
		12720, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 16
		35, // KS1 -> KS2 multiplication threshold
		124, // KS2 -> KS4 multiplication threshold
		8560, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 17
		33, // KS1 -> KS2 multiplication threshold
		102, // KS2 -> KS4 multiplication threshold
		9451, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 18
		31, // KS1 -> KS2 multiplication threshold
		78, // KS2 -> KS4 multiplication threshold
		7022, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 19
		33, // KS1 -> KS2 multiplication threshold
		78, // KS2 -> KS4 multiplication threshold
		8560, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 20
		29, // KS1 -> KS2 multiplication threshold
		70, // KS2 -> KS4 multiplication threshold
		5393, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 21
		27, // KS1 -> KS2 multiplication threshold
		72, // KS2 -> KS4 multiplication threshold
		7022, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 22
		31, // KS1 -> KS2 multiplication threshold
		66, // KS2 -> KS4 multiplication threshold
		5217, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 23
		29, // KS1 -> KS2 multiplication threshold
		80, // KS2 -> KS4 multiplication threshold
		7022, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 24
		29, // KS1 -> KS2 multiplication threshold
		66, // KS2 -> KS4 multiplication threshold
		5760, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 25
		25, // KS1 -> KS2 multiplication threshold
		66, // KS2 -> KS4 multiplication threshold
		6154, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 26
		29, // KS1 -> KS2 multiplication threshold
		57, // KS2 -> KS4 multiplication threshold
		4280, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 27
		25, // KS1 -> KS2 multiplication threshold
		54, // KS2 -> KS4 multiplication threshold
		3511, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 28
		27, // KS1 -> KS2 multiplication threshold
		43, // KS2 -> KS4 multiplication threshold
		3511, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 29
		25, // KS1 -> KS2 multiplication threshold
		51, // KS2 -> KS4 multiplication threshold
		4280, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 30
		21, // KS1 -> KS2 multiplication threshold
		47, // KS2 -> KS4 multiplication threshold
		3877, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 31
		24, // KS1 -> KS2 multiplication threshold
		39, // KS2 -> KS4 multiplication threshold
		4280, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 32
		24, // KS1 -> KS2 multiplication threshold
		39, // KS2 -> KS4 multiplication threshold
		3877, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 33
		21, // KS1 -> KS2 multiplication threshold
		33, // KS2 -> KS4 multiplication threshold
		3511, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 34
		21, // KS1 -> KS2 multiplication threshold
		31, // KS2 -> KS4 multiplication threshold
		3511, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 35
		21, // KS1 -> KS2 multiplication threshold
		31, // KS2 -> KS4 multiplication threshold
		3511, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 36
		21, // KS1 -> KS2 multiplication threshold
		27, // KS2 -> KS4 multiplication threshold
		3511, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 37
		17, // KS1 -> KS2 multiplication threshold
		30, // KS2 -> KS4 multiplication threshold
		3511, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 38
		19, // KS1 -> KS2 multiplication threshold
		25, // KS2 -> KS4 multiplication threshold
		3511, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 39
		19, // KS1 -> KS2 multiplication threshold
		27, // KS2 -> KS4 multiplication threshold
		3511, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 40
		19, // KS1 -> KS2 multiplication threshold
		27, // KS2 -> KS4 multiplication threshold
		3511, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 41
		17, // KS1 -> KS2 multiplication threshold
		25, // KS2 -> KS4 multiplication threshold
		3077, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 42
		17, // KS1 -> KS2 multiplication threshold
		24, // KS2 -> KS4 multiplication threshold
		3180, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 43
		16, // KS1 -> KS2 multiplication threshold
		19, // KS2 -> KS4 multiplication threshold
		2697, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 44
		17, // KS1 -> KS2 multiplication threshold
		24, // KS2 -> KS4 multiplication threshold
		3077, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 45
		16, // KS1 -> KS2 multiplication threshold
		21, // KS2 -> KS4 multiplication threshold
		2697, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 46
		17, // KS1 -> KS2 multiplication threshold
		16, // KS2 -> KS4 multiplication threshold
		3077, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 47
		16, // KS1 -> KS2 multiplication threshold
		19, // KS2 -> KS4 multiplication threshold
		2697, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 48
		17, // KS1 -> KS2 multiplication threshold
		14, // KS2 -> KS4 multiplication threshold
		2697, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 49
		14, // KS1 -> KS2 multiplication threshold
		15, // KS2 -> KS4 multiplication threshold
		2697, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 50
		14, // KS1 -> KS2 multiplication threshold
		16, // KS2 -> KS4 multiplication threshold
		2697, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 51
		14, // KS1 -> KS2 multiplication threshold
		17, // KS2 -> KS4 multiplication threshold
		2363, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 52
		14, // KS1 -> KS2 multiplication threshold
		17, // KS2 -> KS4 multiplication threshold
		2363, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 53
		13, // KS1 -> KS2 multiplication threshold
		14, // KS2 -> KS4 multiplication threshold
		1876, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 54
		14, // KS1 -> KS2 multiplication threshold
		19, // KS2 -> KS4 multiplication threshold
		2363, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 55
		13, // KS1 -> KS2 multiplication threshold
		12, // KS2 -> KS4 multiplication threshold
		2071, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 56
		13, // KS1 -> KS2 multiplication threshold
		10, // KS2 -> KS4 multiplication threshold
		2363, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 57
		13, // KS1 -> KS2 multiplication threshold
		15, // KS2 -> KS4 multiplication threshold
		1876, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 58
		13, // KS1 -> KS2 multiplication threshold
		10, // KS2 -> KS4 multiplication threshold
		1876, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 59
		13, // KS1 -> KS2 multiplication threshold
		13, // KS2 -> KS4 multiplication threshold
		1539, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 60
		13, // KS1 -> KS2 multiplication threshold
		13, // KS2 -> KS4 multiplication threshold
		1539, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 61
		13, // KS1 -> KS2 multiplication threshold
		12, // KS2 -> KS4 multiplication threshold
		1756, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 62
		13, // KS1 -> KS2 multiplication threshold
		19, // KS2 -> KS4 multiplication threshold
		2363, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 63
		13, // KS1 -> KS2 multiplication threshold
		19, // KS2 -> KS4 multiplication threshold
		2363, // KS4 -> FFT multiplication threshold
	},
	{ // bits = 64
		13, // KS1 -> KS2 multiplication threshold
		23, // KS2 -> KS4 multiplication threshold
		4280, // KS4 -> FFT multiplication threshold
	},
};
