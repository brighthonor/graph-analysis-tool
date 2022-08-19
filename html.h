#ifndef HTML_H
#define HTML_H

#include "plot.h"
#include "stat.h"

namespace snu {
// open file, add title & set styles
FILE *openHtml(const char *name, bool directed);

// add stat to html
void addStatToHtml(FILE *fp, Stat *stat, bool directed);

// add plot, close file
void closeHtml(FILE *fp, Plot &plot);
}  // namespace snu

#endif  // HTML_H
