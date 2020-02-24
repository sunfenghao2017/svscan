#ifndef SVCREG_H
#define SVCREG_H

#include <string>
#include <util.h>
#include "cgranges.h"

/** class to fetch creg for svscan */
struct SVCreg{
    std::string glist; ///< gene list, first field must be gene name
    std::string ibed;  ///< input bed to extract creg, 5th field must be gene name
    std::string obed;  ///< output bed of merged extrated regs

    /** get creg */
    void getCreg();
};

#endif
