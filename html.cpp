#include "html.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace snu {
FILE *openHtml(const char *name, bool directed) {
    std::string real_name = name + std::string(".html");
    FILE *fp = fopen(real_name.c_str(), "w");

    fprintf(fp,
            "\
        <!DOCTYPE html>\
        <html>\
            <head>\
                <meta charset = \"utf-8\">\
                <title>SNU Graph Analysis Tool</title>\
                <style type=\"text/css\">\
                    h1{\
                        text-indent: 1em;\
                    }\
                    h2{\
                        text-indent: 2em;\
                    }\
                    h3{\
                        text-indent: 4em;\
                    }\
                </style>\
            </head>\
            <body>\
                <h1><p>SNU Graph Analysis Tool</p></h1>");

    fprintf(fp,
            "\
                <h2>\
                    Graph Information\
                </h2>\
                <h3>\
                    <p> graph name: %s </p>\
                    <p> graph type: %s simple graph </p>\
                </h3>\
                    ",
            name, directed ? "directed" : "undirected");

    return fp;
}

void addStatToHtml(FILE *fp, Stat *stat, bool directed) {
    if (stat && stat->getSuccess()) {
        stat->writeToHTML(fp, directed);
    }
}

void closeHtml(FILE *fp, Plot &plot) {
    if (plot.makeplot)
        fprintf(fp,
                "\
                <h2>\
                    Statistics Image\
                </h2>\
                <h3>\
                    <p> label-vertex image </p>\
                    <img src=\"pyplot/%s_label-vertex.png\" width=\"400\" alt=\"label-vertex image\">\
                    <img src=\"pyplot/%s_label-vertex_log.png\" width=\"400\" alt=\"log scale label-vertex image\">\
                </h3>\
                <h3>\
                    <p> degree image </p>\
                    <img src=\"pyplot/%s_degree.png\" width=\"400\" alt=\"degree image\">\
                    <img src=\"pyplot/%s_degree_log.png\" width=\"400\" alt=\"log scale degree image\">\
                </h3>\
            ",
                plot.name.c_str(), plot.name.c_str(), plot.name.c_str(), plot.name.c_str());

    fprintf(fp, "</body></html>");
    fclose(fp);
}
}  // namespace snu
