#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "kvec.h"
#include "kstring.h"

#include "ksort.h"
KSORT_INIT_GENERIC(double)

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 65536)

#define SEP_SPACE 256
#define SEP_CSV   257

typedef kvec_t(uint32_t) vec32_t;
typedef kvec_t(uint64_t) vec64_t;
typedef kvec_t(double) vecdbl_t;

/************************
 *** Generic routines ***
 ************************/

int ttk_parse_cols(vec64_t *cols, const char *str, int reorder)
{
	int32_t beg, end, x;
	char *p = (char*)str;
	cols->n = 0;
	while ((x = strtol(p, &p, 10)) > 0) { // parse the field string
		beg = end = x;
		if (*p == '-')
			end = (x = strtol(p + 1, &p, 10)) >= beg? x : INT32_MAX;
		kv_push(uint64_t, *cols, (uint64_t)(beg-1)<<32 | end);
		if (*p) ++p; // skip ','
		else break;
	}
	if (*p) {
		cols->n = 0;
		return -1;
	}
	if (!reorder) {
		uint64_t *i;
		int j, k;
		for (i = cols->a + 1; i < cols->a + cols->n; ++i) // insertion sort
			if (*i < *(i - 1)) {
				uint64_t *j, tmp = *i;
				for (j = i; j > cols->a && tmp < *(j-1); --j) *j = *(j - 1);
				*j = tmp;
			}
		for (j = k = 1; j < cols->n; ++j) { // merge overlapping col regions
			if (cols->a[j]>>32 <= (uint32_t)cols->a[k-1])
				cols->a[k-1] = cols->a[k-1]>>32<<32 | (uint32_t)cols->a[j];
			else cols->a[k++] = cols->a[j];
		}
		cols->n = k;
	}
	return 0;
}

static inline void ttk_split(vec64_t *buf, int sep, int len, const char *str)
{
	int i, b;
	buf->n = 0;
	if (sep >= 0 && sep < 256) {
		for (i = b = 0; i <= len; ++i) // mark columns
			if (i == len || str[i] == sep) {
				kv_push(uint64_t, *buf, (uint64_t)b<<32 | i);
				b = i + 1;
			}
	} else if (sep == SEP_SPACE) {
		for (i = b = 0; i <= len; ++i) // mark columns
			if (i == len || isspace(str[i])) {
				kv_push(uint64_t, *buf, (uint64_t)b<<32 | i);
				b = i + 1;
			}
	} else if (sep == SEP_CSV) {
		int state = 0;
		for (i = b = 0; i <= len; ++i) {
			if (i == len || (state == 0 && str[i] == ',')) {
				kv_push(uint64_t, *buf, (uint64_t)b<<32 | i);
				b = i + 1;
			} else if (state == 1 && str[i] == '\\') {
				if (i < len - 1) ++i;
			} else if (str[i] == '"') {
				state ^= 1;
			}
		}
	}
}

int ttk_parse_sep(const char *s)
{
	int sep;
	if (strcmp(s, "isspace") == 0 || strcmp(s, "space") == 0) sep = SEP_SPACE;
	else if (strcmp(s, "csv") == 0) sep = SEP_CSV;
	else if (strlen(s) == 1) sep = s[0];
	else {
		fprintf(stderr, "[E::%s] multiple delimitors are not supported\n", __func__);
		return -1;
	}
	return sep;
}

/***********
 *** cut ***
 ***********/

int main_cut(int argc, char *argv[])
{
	vec64_t cols = {0,0,0}, buf = {0,0,0};
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0}, out = {0,0,0};
	int dret, sep = '\t', c, reorder = 0, skip_char = -1;
	const char *fields = 0;

	while ((c = getopt(argc, argv, "rd:f:S:")) >= 0) {
		if (c == 'r') reorder = 1;
		else if (c == 'S') skip_char = optarg[0];
		else if (c == 'f') fields = optarg;
		else if (c == 'd') {
			if ((sep = ttk_parse_sep(optarg)) < 0) return 1;
		}
	}
	
	if (argc == optind && isatty(fileno(stdin))) {
		fprintf(stderr, "\nUsage: tabtk cut [options] [file.txt]\n\n");
		fprintf(stderr, "Options: -d CHAR     delimitor, a single CHAR or 'space' for both SPACE and TAB or 'csv' [TAB]\n");
		fprintf(stderr, "         -S CHAR     keep full lines starting with CHAR [null]\n");
		fprintf(stderr, "         -f STR      fields to cut; format identical to Unix cut [null]\n");
		fprintf(stderr, "         -r          reorder fields\n\n");
		return 1;
	}

	if (fields == 0) {
		fprintf(stderr, "[E::%s] no list of fields is specified.\n", __func__);
		return 2;
	}

	if (ttk_parse_cols(&cols, fields, reorder) < 0) {
		fprintf(stderr, "[E::%s] wrong fields format\n", __func__);
		free(cols.a);
		return 1;
	}

	fp = optind < argc && strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0) >= 0) {
		int i;
		if (skip_char >= 0 && str.s[0] == skip_char) {
			puts(str.s);
			continue;
		}
		ttk_split(&buf, sep, str.l, str.s);
		for (i = 0, out.l = 0; i < cols.n; ++i) { // print columns
			int32_t j, beg = cols.a[i]>>32, end = (int32_t)cols.a[i];
			for (j = beg; j < end && j < buf.n; ++j) {
				uint64_t x = buf.a[j];
				if (out.l) kputc('\t', &out);
				kputsn(&str.s[x>>32], (uint32_t)x - (x>>32), &out);
			}
		}
		puts(out.s? out.s : "");
	}
	ks_destroy(ks);
	gzclose(fp);

	free(str.s); free(out.s); free(cols.a); free(buf.a);
	return 0;
}

/*****************
 *** intersect ***
 *****************/

#include "khash.h"
KHASH_SET_INIT_STR(strset)

int main_isct(int argc, char *argv[])
{
	vec64_t cols1 = {0,0,0}, cols2 = {0,0,0}, buf = {0,0,0};
	char *fields1 = 0, *fields2 = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0}, key = {0,0,0};
	int c, dret, sep = '\t', is_comp = 0, skip_char = -1;
	khash_t(strset) *h;

	while ((c = getopt(argc, argv, "c1:2:d:S")) >= 0) {
		if (c == '1') fields1 = optarg;
		else if (c == '2') fields2 = optarg;
		else if (c == 'c') is_comp = 1;
		else if (c == 'S') skip_char = optarg[0];
		else if (c == 'd') {
			if ((sep = ttk_parse_sep(optarg)) < 0) return 1;
		}
	}

	if (optind + 1 > argc || (optind + 2 > argc && isatty(fileno(stdin)))) {
		fprintf(stderr, "\nUsage:   tabtk isct [options] <loaded.txt> <streamed.txt>\n\n");
		fprintf(stderr, "Options: -1 STR    field(s) of the loaded file [1]\n");
		fprintf(stderr, "         -2 STR    field(s) of the streamed file [same as -1]\n");
		fprintf(stderr, "         -d CHAR   delimitor [TAB]\n");
		fprintf(stderr, "         -S CHAR   skip lines starting with CHAR [null]\n");
		fprintf(stderr, "         -c        print lines not present in loaded.txt\n");
		fputc('\n', stderr);
		return 1;
	}

	if (fields1 == 0) kv_push(uint64_t, cols1, 0ULL<<32|1);
	else ttk_parse_cols(&cols1, fields1, 0);
	if (fields2 == 0) fields2 = fields1;
	if (fields2 == 0) kv_push(uint64_t, cols2, 0ULL<<32|1);
	else ttk_parse_cols(&cols2, fields2, 0);

	h = kh_init(strset);
	fp = gzopen(argv[optind], "r");
	ks = ks_init(fp);
	while (ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0) >= 0) {
		int i, absent;
		khint_t k;
		if (str.l == 0) continue;
		if (skip_char >= 0 && str.s[0] == skip_char) continue;
		ttk_split(&buf, sep, str.l, str.s);
		for (i = 0, key.l = 0; i < cols1.n; ++i) {
			int32_t j, beg = cols1.a[i]>>32, end = (int32_t)cols1.a[i];
			for (j = beg; j < end && j < buf.n; ++j) {
				uint64_t x = buf.a[j];
				if (key.l) kputc('\t', &key);
				kputsn(&str.s[x>>32], (uint32_t)x - (x>>32), &key);
			}
		}
		k = kh_put(strset, h, key.s, &absent);
		if (absent) kh_key(h, k) = strdup(key.s);
	}
	ks_destroy(ks);
	gzclose(fp);

	fp = optind + 1 < argc && strcmp(argv[optind+1], "-")? gzopen(argv[optind+1], "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0) >= 0) {
		int i, present;
		if (str.l == 0) continue;
		if (skip_char >= 0 && str.s[0] == skip_char) continue;
		ttk_split(&buf, sep, str.l, str.s);
		for (i = 0, key.l = 0; i < cols2.n; ++i) {
			int32_t j, beg = cols2.a[i]>>32, end = (int32_t)cols2.a[i];
			for (j = beg; j < end && j < buf.n; ++j) {
				uint64_t x = buf.a[j];
				if (key.l) kputc('\t', &key);
				kputsn(&str.s[x>>32], (uint32_t)x - (x>>32), &key);
			}
		}
		present = (kh_get(strset, h, key.s) != kh_end(h));
		if (present != is_comp) puts(str.s);
	}
	ks_destroy(ks);
	gzclose(fp);

	free(str.s); free(key.s); free(buf.a); free(cols1.a); free(cols2.a);
	kh_destroy(strset, h);
	return 0;
}

/***********
 *** num ***
 ***********/

int main_num(int argc, char *argv[])
{
	int c, col = 0, in_ram = 0, dret, show_more = 0, skip_char = -1;
	uint64_t n = 0;
	double qt = -1, min = DBL_MAX, max = DBL_MIN, avg;
	long double sum = 0.;
	vecdbl_t a = {0,0,0};
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};

	while ((c = getopt(argc, argv, "Qc:q:S:")) >= 0) {
		if (c == 'c') col = atol(optarg) - 1;
		else if (c == 'Q') show_more = in_ram = 1;
		else if (c == 'q') qt = atof(optarg), in_ram = 1;
		else if (c == 'S') skip_char = optarg[0];
	}
	if (argc == optind && isatty(fileno(stdin))) {
		fprintf(stderr, "\nUsage:   tabtk num [options] [file.txt]\n\n");
		fprintf(stderr, "Options: -c INT     column number [1]\n");
		fprintf(stderr, "         -q FLOAT   only compute quantile, negative to disable [-1]\n");
		fprintf(stderr, "         -S CHAR    skip lines starting with CHAR [null]\n");
		fprintf(stderr, "         -Q         output quartiles, stdandard deviation and skewness\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Notes: number, mean, min, max[, std.dev, skewness, 25%%-percentile, median, 75%%, 2.5%%, 97.5%%, 1%%, 99%%, 0.5%%, 99.5%%]\n\n");
		return 1;
	}
	fp = optind < argc && strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0) >= 0) {
		int i, beg;
		double x;
		char *p;
		if (skip_char >= 0 && str.s[0] == skip_char) continue;
		for (i = beg = c = 0; i <= str.l; ++i) // mark columns
			if (isspace(str.s[i]) || i == str.l) {
				if (c++ == col) break;
				beg = i + 1;
			}
		if (i > str.l) continue; // not enough fields
		x = strtod(&str.s[beg], &p);
		if (p == &str.s[beg]) continue; // conversion failed
		++n; sum += x;
		min = min < x? min : x;
		max = max > x? max : x;
		if (in_ram) kv_push(double, a, x);
	}
	if (n == 0) {
		fprintf(stderr, "[E::%s] no data are read\n", __func__);
		return 1;
	}
	avg = sum / n;
	if (qt < 0. || qt > 1.) {
		printf("%llu\t%g\t%g\t%g", (unsigned long long)n, avg, min, max);
		if (show_more) {
			long double sum2 = 0., sum3 = 0.;
			double q[3], CI5[2], CI2[2], CI1[2];
			uint64_t i;
			if (n > 1) {
				double g1, tmp;
				for (i = 0; i < n; ++i) {
					double t = (a.a[i] - avg) * (a.a[i] - avg);
					sum2 += t;
					sum3 += t * (a.a[i] - avg);
				}
				tmp = sqrt(sum2 / n);
				printf("\t%g", sqrt(sum2 / (n - 1)));
				g1 = (sum3 / n) / (tmp * tmp * tmp);
				if (n > 2) printf("\t%g", sqrt((double)n * (n - 1)) / (n - 2) * g1);
				else printf("\tNaN");
			} else printf("\tNaN");
			q[0] = ks_ksmall(double, a.n, a.a, (int)(ceil(n * .25) + .499));
			q[1] = ks_ksmall(double, a.n, a.a, (int)(ceil(n * .50) + .499));
			q[2] = ks_ksmall(double, a.n, a.a, (int)(ceil(n * .75) + .499));
			CI5[0] = ks_ksmall(double, a.n, a.a, (int)(ceil(n * .025) + .499));
			CI5[1] = ks_ksmall(double, a.n, a.a, (int)(ceil(n * .975) + .499));
			CI2[0] = ks_ksmall(double, a.n, a.a, (int)(ceil(n * .01) + .499));
			CI2[1] = ks_ksmall(double, a.n, a.a, (int)(ceil(n * .99) + .499));
			CI1[0] = ks_ksmall(double, a.n, a.a, (int)(ceil(n * .005) + .499));
			CI1[1] = ks_ksmall(double, a.n, a.a, (int)(ceil(n * .995) + .499));
			printf("\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g", q[0], q[1], q[2], CI5[0], CI5[1], CI2[0], CI2[1], CI1[0], CI1[1]);
		}
	} else {
		double q;
		q = ks_ksmall(double, a.n, a.a, (int)(ceil(n * qt) + .499));
		printf("%g", q);
	}
	putchar('\n');
	ks_destroy(ks);
	gzclose(fp);
	free(a.a); free(str.s);
	return 0;
}

/************
 *** grep ***
 ************/

#include "regexp9.h"

static inline int grep1(Reprog *p, char *s, int start, int end)
{
	int ret, c = s[end];
	s[end] = 0;
	ret = regexec9(p, &s[start], 0, 0);
	s[end] = c;
	return ret;
}

int main_grep(int argc, char *argv[])
{
	vec64_t cols = {0,0,0}, buf = {0,0,0};
	Reprog *p;
	gzFile fp;
	kstream_t *ks;
	const char *fields = 0;
	int show_lineno = 0, match = 0, i, j, dret, c, sep = '\t';
	long lineno = 0;
	kstring_t str = {0,0,0};

	while ((c = getopt(argc, argv, "f:d:n")) >= 0) {
		if (c == 'f') fields = optarg;
		else if (c == 'n') show_lineno = 1;
		else if (c == 'd') {
			if ((sep = ttk_parse_sep(optarg)) < 0) return 1;
		}
	}
	if (argc == optind || (argc == optind + 1 && isatty(fileno(stdin)))) {
		fprintf(stderr, "\nUsage: tabtk grep [options] <pattern> [file.txt]\n\n");
		fprintf(stderr, "Options: -d CHAR    delimitor, a single CHAR or 'space' for both SPACE and TAB or 'csv' [TAB]\n");
		fprintf(stderr, "         -f STR     fields [null]\n");
		fprintf(stderr, "         -n         output matching line number and column number\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (fields && ttk_parse_cols(&cols, fields, 0) < 0) {
		fprintf(stderr, "[E::%s] wrong fields format\n", __func__);
		free(cols.a);
		return 1;
	}

	p = regcomp9(argv[optind]);
	fp = optind+1 < argc && strcmp(argv[optind], "-")? gzopen(argv[optind+1], "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0) >= 0) {
		++lineno;
		ttk_split(&buf, sep, str.l, str.s);
		if (fields == 0) {
			for (j = 0; j < buf.n; ++j)
				if (grep1(p, str.s, buf.a[j]>>32, (uint32_t)buf.a[j])) {
					if (show_lineno) printf("%ld\t%d\t", lineno, j+1);
					puts(str.s);
					match = 1;
					break;
				}
		} else {
			for (i = 0; i < cols.n; ++i) {
				int32_t start = cols.a[i]>>32, end = (int32_t)cols.a[i];
				for (j = start; j < end && j < buf.n; ++j)
					if (grep1(p, str.s, buf.a[j]>>32, (uint32_t)buf.a[j])) {
						if (show_lineno) printf("%ld\t%d\t", lineno, j+1);
						puts(str.s);
						match = 1;
						break;
					}
				if (j < end && j < buf.n) break;
			}
		}
	}
	ks_destroy(ks);
	gzclose(fp);
	free(p);
	return !match; // for grep, if no lines match, the return code is 1
}

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   tabtk-r%d <command> [arguments]\n\n", 9);
	fprintf(stderr, "Command: cut         Unix cut with optional column reordering\n");
	fprintf(stderr, "         num         summary statistics on a single numerical column\n");
	fprintf(stderr, "         isct        intersect two files\n");
	fprintf(stderr, "         grep        field-aware grep (slower than Unix grep)\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int ret = 1;
	if (argc == 1) return usage();
	if (strcmp(argv[1], "cut") == 0) ret = main_cut(argc-1, argv+1);
	else if (strcmp(argv[1], "num") == 0) ret = main_num(argc-1, argv+1);
	else if (strcmp(argv[1], "isct") == 0 || strcmp(argv[1], "intersect") == 0) ret = main_isct(argc-1, argv+1);
	else if (strcmp(argv[1], "grep") == 0) ret = main_grep(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized commad '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return ret;
}
