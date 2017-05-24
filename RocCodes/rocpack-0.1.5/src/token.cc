#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <cstring>

#include <err.h>
#include <sysexits.h>

#include "token.h"

static const int keyword_type[] =
{
	CREATE   , BOUNDARY, SHAPE,
	TRANSLATE, ROTATE  , SCALE,
	SET      , UNSET   , IMPORT,
	TAG      , MATERIAL, COLOR,
	SIZE     , UNIFORM , GAMMA,
	WEIBULL  , GAUSSIAN, CAUCHY,
	LOGNORMAL
};

static const char *const keyword_name[] =
{
	"create", "boundary", "shape",
	"translate", "rotate", "scale",
	"set", "unset", "import",
	"tag", "material", "color",
	"size", "uniform", "gamma",
	"weibull", "gaussian", "cauchy",
	"lognormal", NULL
};

token_t next_token(FILE* input, int& line)
{
	int c;
	token_t token;

    start:

	/* skip white space */

    while (isspace(c = getc(input)))
		if (c == '\n') line++;

    /* skip comments */

    if (c == '#') {
        while ((c = getc(input)) != '\n')
            if (c == EOF) {
				token.type = EOF;
				return token;
			}
		line++;
		goto start;
	}

    /* number */

    if (isdigit(c) || c == '-' || c == '+' || c == '.') {
        ungetc (c, input);

		if (fscanf(input, "%f", &token.value.fval) != 1)
			errx(EX_DATAERR, "syntax error: line %d: invalid number format\n", line);

		token.type = NUMBER;
		return token;
    }

	/* string */

	if (c == '"') {
		int i = 0;
        static int length = 0;
		static char *string = 0;

		if (length == 0)
			length = 32, string = (char *) malloc (length + 1);

		c = getc(input);

		do {
			if (i == length) {
				length *= 2;
				string = (char *) realloc (string, length + 1);
			}
			string[i++] = c, c = getc(input);
		} while(c != '"');

		string[i] = '\0';
		token.type = STRING;
		token.value.str = string;
		return token;
	}

    /* keyword or identifier */

    if (isalpha(c)) {
        int i = 0;
        static char *name = 0;
        static int length = 0;

        if (length == 0)
            length = 32, name = (char *) malloc (length + 1);

        do {
            if (i == length) {
                length *= 2;
                name = (char *) realloc (name, length + 1);
            }

            name[i++] = c, c = getc(input);
        } while (isalnum (c) || c == '_');

        ungetc(c, input);

        name[i] = '\0';

        /* keyword */

		for (i = 0; keyword_name[i] != NULL; i++)
			if(strcmp(name, keyword_name[i]) == 0) {
				token.type = keyword_type[i];
				return token;
			}

        /* identifier */

		token.type = IDENTIFIER;
		token.value.str = name;
        return token;
    }

    /* any other character is a token by itself */

	token.type = c;
    return token;
}

char* token_description(int type)
{
	static char desc[64];

	if (type == NUMBER)
		snprintf(desc, sizeof(desc), "a number");
	else if (type == STRING)
		snprintf(desc, sizeof(desc), "a string");
	else if (type == IDENTIFIER)
		snprintf(desc, sizeof(desc), "an identifier");
	else if (type < 257)
		snprintf(desc, sizeof(desc), "a '%c' character token", (char) type);
	else if (type > 257) {
		int i;
		for (i = 0; keyword_name[i] != NULL; i++)
			if(type == keyword_type[i]) {
				snprintf(desc, sizeof(desc), "a '%s' keyword", keyword_name[i]);
				break;
			}

		if (keyword_name[i] == NULL)
			snprintf(desc, sizeof(desc), "an unknown token (%d)", type);
	}

	return desc;
}

void debug_token(token_t token)
{
	if (token.type == NUMBER)
		fprintf(stderr, "a number: %f\n", token.value.fval);
	else if (token.type == STRING)
		fprintf(stderr, "a string: '%s'\n", token.value.str);
	else if (token.type == IDENTIFIER)
		fprintf(stderr, "an identifier: '%s'\n", token.value.str);
	else if (token.type < 257)
		fprintf(stderr, "a character token: '%c' (%d)\n",
			(char) token.type, token.type);
	else if (token.type > 257) {
		int i;
		for (i = 0; keyword_name[i] != NULL; i++)
			if(token.type == keyword_type[i]) {
				fprintf(stderr, "a keyword: %s\n", keyword_name[i]);
				return;
			}

		if (keyword_name[i] == NULL)
			fprintf(stderr, "an unknown token (%d)\n", token.type);
	}
}
