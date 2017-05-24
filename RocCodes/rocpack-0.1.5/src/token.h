#ifndef _TOKEN_H
#define _TOKEN_H

#include <cstdio>

typedef struct {
	int type;
	int line;
	union {
		float fval;
		char  *str;
	} value;
} token_t;

enum token_type
{
	NUMBER    = 258, IDENTIFIER = 259, STRING = 260,
	CREATE    = 261, BOUNDARY   = 262, SHAPE  = 263,
	TRANSLATE = 264, ROTATE     = 265, SCALE  = 266,
	SET       = 267, UNSET      = 268, IMPORT = 269,
	TAG       = 270, MATERIAL   = 271, COLOR  = 272,
	SIZE      = 273, UNIFORM    = 274, GAMMA  = 275,
	WEIBULL   = 276, GAUSSIAN   = 277, CAUCHY = 278,
	LOGNORMAL = 279
};

token_t next_token(FILE* input, int& line);
char* token_description(int type);
void debug_token(token_t token);

#endif
