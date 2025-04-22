#ifndef UTIL
#define UTIL

// #include "Vector.hpp"

typedef enum {
	BI_D, 	// 2D
	TRI_D	// 3D
} FileFormat;

typedef enum {
	NUT,
	EGF,
	CHEM
} OutputMode;

typedef enum {
	NECROTIC,
	QUIESCENT,      // G0
	PROLIFERATIVE2, // S-G2-M
	HYPOXIC,
	APOPTOTIC,
	PROLIFERATIVE1, // G1
	NORMOXIC,
	MIGRATE,
	KILLED,
	LAST_CELL_STATE = KILLED
} CellState;

#endif /* end of include guard: UTIL */
