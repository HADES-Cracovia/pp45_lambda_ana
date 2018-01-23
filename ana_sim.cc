#include "KTifiniAnalysis.h"
#include "pp45_Lambda.h"

int main(int argc, char **argv)
{
    pp45_Lambda	ana_Lambda_pp45("Lambda analysis for pp45", "T");

    KTifiniAnalysis tifini(argc, argv, KT::Sim, &ana_Lambda_pp45);

    tifini.exec();

    return 0;
}
