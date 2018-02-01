#include <getopt.h>

#include "KTifiniAnalysis.h"
#include "Tifini3Config.h"

#ifdef HYDRA1COMP
#include "ef_lambda_pp35.h"
#else
#include "ef_lambda_pp45.h"
#endif

int main(int argc, char **argv)
{
    int c;

    int flag_verbose = 0;
    int flag_selection = 0;

    while (1)
    {
        static struct option long_options[] =
            {
            /* These options set a flag. */
            {"verbose", no_argument,    &flag_verbose, 1},
            {"brief",   no_argument,    &flag_verbose, 0},
            /* These options don't set a flag.
                We distinguish them by their indices. */
            {"help",    no_argument,    0, 'h'},
            {"lambda",  no_argument,    &flag_selection, 0},
            {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;
    
        c = getopt_long (argc, argv, "h", long_options, &option_index);
    
        /* Detect the end of the options. */
        if (c == -1)
            break;
    
        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                break;
    
            case 'h':
            case '?':
                /* getopt_long already printed an error message. */
//                 Usage();
                break;
    
            default:
                abort();
        }
    }

    KAbstractAnalysis * ana;
    switch (flag_selection)
    {
        case 0:     // lambda
#ifdef HYDRA1COMP
            ana = new ef_lambda_pp35("Lambda analysis for pp35", "T");
#else
            ana = new ef_lambda_pp45("Lambda analysis for pp45", "T");
#endif
    }

    KTifiniAnalysis tifini(argc, argv, KT::Exp, ana);

    tifini.exec();

    delete ana;

    return 0;
}
