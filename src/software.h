#ifndef SOFTWARE_H
#define SOFTWARE_H

#include <string>

/** class to store software environment */
struct Software{
    std::string version = "0.0.0"; ///< software version
    std::string cmd;               ///< software execution command
    std::string cwd;               ///< software execution directory
    std::string cmp;               ///< software update time

    /** Software constructor */
    Software(){}

    /** Software destructor */
    ~Software(){}
};

#endif
