#include "KTifiniAnalysis.h"
#include "pp35_Lambda.h"

int main(int argc, char **argv)
{
    pp35_Lambda	ana_Lambda_pp35("Lambda analysis for pp35", "T");

    KTifiniAnalysis tifini(argc, argv, KT::Sim, &ana_Lambda_pp35);

    tifini.exec();

    return 0;
}
