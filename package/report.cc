#include <stdio.h>

#include "report.h"
#include "rect.h"

report::report() : sum(0) {}
void report::print() { printf("Number of pairs = %d\n", sum); }
void report::found(const struct RectInfo *a, const struct RectInfo *b) { sum++; }

