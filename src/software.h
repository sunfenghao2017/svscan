#ifndef SOFTWARE_H
#define SOFTWARE_H

#include <string>
#include <chrono>
#include <unistd.h>
#include "timer.h"

/** class to store software environment */
struct Software{
    std::string ver; ///< software version
    std::string cmd; ///< software execution command
    std::string cwd; ///< software execution directory
    std::string dur; ///< software execution time
    std::string cmp; ///< software update time
    Timer* timer;    ///< timer to recording execution time

    /** Software constructor */
    Software(){
        ver = "0.1.1";
        cmp = std::string(__TIME__) + " " + std::string(__DATE__);
        char cpath[FILENAME_MAX];
        getcwd(cpath, FILENAME_MAX);
        cwd = cpath;
        timer = new Timer();
    }

    /** Software destructor */
    ~Software(){
        if(timer){
            delete timer;
        }
    }

    /** update cmd of Software
     * @param argc number of arguments feed to main
     * @param argv array of arguments feed to main
     */
    inline void update(int argc, char** argv){
        for(int i = 0; i < argc; ++i){
            cmd.append(argv[i]);
            cmd.append(" ");
        }
    }

    /** update execution time of Software */
    inline void end(){
        dur = timer->toStr();
    }

    /** get string representation of execution time
     * @return string representation of execution time
     */
    inline std::string getExecutionTime(){
        return timer->toStr();
    }

    /** get information of Software 
     * @return update and version info of software
     */
    inline std::string getSoftInfo(){
        return "updated: " + cmp + "\nversion: " + ver + "\n";
    }
};

#endif
