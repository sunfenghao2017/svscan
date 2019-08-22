#ifndef UTIL_H
#define UTIL_H

#ifdef WINDOWS
#include <direct.h>
#else
#include <unistd.h>
#endif

#include <string>
#include <cerrno>
#include <cstdio>
#include <cctype>
#include <vector>
#include <mutex>
#include <numeric>
#include <fstream>
#include <utility>
#include <iostream>
#include <algorithm>
#include <functional>
#include <dirent.h>
#include <sys/stat.h>

/** utility to operate on strings and directories */
namespace util{
    /** whether a string starts with some substring
     * @param str whole string
     * @param pre substring
     * @return true if str starts with pre
     */
    inline bool startsWith(const std::string& str, const std::string& pre){
        if(str.length() < pre.length()){
            return false;
        }else{
            return std::equal(pre.begin(), pre.end(), str.begin());
        }
    }

    /** whether a string ends with some substring
     * @param str whole string
     * @param suf substring
     * @return true if str ends with suf
     */
    inline bool endsWith(const std::string& str, const std::string& suf){
        if(str.length() < suf.length()){
            return false;
        }else{
            return std::equal(suf.rbegin(), suf.rend(), str.rbegin());
        }
    }

    /** get rid of the leading and ending white space characters of a string
     * @param str string to be stripped in both ends
     * @return a string with white spaces stripped in both ends
     */
    inline std::string strip(const std::string& str){
        std::string::size_type ipos = str.find_first_not_of(" \t\n\v\f\r");
        if(ipos == std::string::npos){
            return "";
        }
        std::string::size_type epos = str.find_last_not_of(" \t\n\v\f\r");
        if(epos == ipos){
            return str.substr(ipos);
        }else{
            return str.substr(ipos, epos - ipos + 1);
        }
    }

    /** get rid of the left leading white space characters of a string
     * @param str string to be stripped from front
     * @return a string with left leading white spaces stripped
     */
    inline std::string lstrip(const std::string& str){
        std::string::size_type pos = str.find_first_not_of(" \t\n\v\f\r");
        if(pos == std::string::npos){
            return "";
        }
        return str.substr(pos);
    }

    /** get rid of the trailling white space characters of a string
     * @param str string to be stripped from back
     * @return a string with right ending white spaces stripped
     */
    inline std::string rstrip(const std::string& str){
        std::string::size_type pos = str.find_last_not_of(" \t\n\v\f\r");
        if(pos == std::string::npos){
           return "";
        }
        return str.substr(0, pos + 1);
    } 

    /** split a string by predefined seperator into a vector
     * @param str string
     * @param vec vector to store the split results
     * @param sep seperators, can contain a series of seperators
     */
    inline void split(const std::string& str, std::vector<std::string>& vec, std::string sep = " "){
        vec.clear();
        std::string::size_type las, cur;
        las = cur = str.find_first_not_of(sep);
        while((las = str.find_first_not_of(sep, las)) != std::string::npos){
            cur = str.find_first_of(sep, las);
            if(cur != std::string::npos){
                vec.push_back(str.substr(las, cur - las));
            }else{
                vec.push_back(str.substr(las));
                break;
            }
            las = cur;
        }
    }
    /** join a list of strings by an seperator
     * @param vec vector to store the split results
     * @param ret joined string
     * @param sep seperator used to join strings
     */
    inline void join(const std::vector<std::string>& vec, std::string& ret, const std::string& sep = " "){
        if(vec.empty()) return;
        if(vec.size() == 1){
            ret = vec[0];
            return;
        }
        for(uint32_t i = 0; i < vec.size() - 1; ++i){
            ret.append(vec[i] + sep);
        }
        ret.append(vec[vec.size() - 1]);
    }
    
    /** join a list of strings by an seperator
     * @param vec vector to store the split results
     * @param sep seperator used to join strings
     * @return joined string
     */
    inline std::string join(const std::vector<std::string>& vec, const std::string& sep = " "){
        std::string ret;
        if(vec.empty()) return ret;
        if(vec.size() == 1){
            ret = vec[0];
            return ret;
        }
        for(uint32_t i = 0; i < vec.size() - 1; ++i){
            ret.append(vec[i] + sep);
        }
        ret.append(vec[vec.size() - 1]);
        return ret;
    }

    /** replace a substr apearing in a string with another string
     * @param str string 
     * @param pat substr of string to be replaced
     * @param des string to be used to replaced with pat
     * @return a string with each pat replaced by des
     */
    inline std::string replace(const std::string& str, const std::string& pat, const std::string& des){
        std::string ret;
        std::string::size_type las = 0, cur = 0;
        while((cur = str.find(pat, cur)) != std::string::npos){
            ret.append(str.substr(las, cur - las));
            ret.append(des);
            cur += pat.length();
            las = cur;
        }
        if(las != std::string::npos){
            ret.append(str.substr(las));
        }
        return ret;
    }

    /** get a reverse sequence of str
     * @param str input string sequence
     * @return reverse sequence of str
     */
    inline std::string reverse(const std::string& str){
        std::string ret = str;
        std::reverse(ret.begin(), ret.end());
        return ret;
    }

    /** get the basename of a path string
     * @param path name
     * @return basename of path
     */
    inline std::string basename(const std::string& path){
        std::string fpath = util::replace(path, "~", std::string(std::getenv("HOME")) + "/");
        if(fpath.find_first_of(" \t\n\v\f\r") != std::string::npos){
            return "";
        }
        std::string::size_type pos1 = fpath.find_last_of("/\\");
        if(pos1 == std::string::npos){
            return fpath;
        }
        std::string::size_type pos2 = fpath.find_last_not_of("/\\");
        if(pos2 == fpath.size() - 1){
            return fpath.substr(pos1 + 1, pos2 - pos1);
        }
        std::string::size_type pos3 = fpath.find_last_of("/\\", pos2);
        if(pos3 == std::string::npos){
            return fpath.substr(0, pos2 + 1);
        }
        return fpath.substr(pos3 + 1, pos2 - pos3);
    }

    /** get the dirname of a path string
     * @param path name
     * @return dirname of path
     */
    inline std::string dirname(const std::string& path){
        std::string fpath = util::replace(path, "~", std::string(std::getenv("HOME")) + "/");
        std::string::size_type pos = fpath.find_last_of("/\\");
        if(pos == std::string::npos){
#ifdef _WIN32
            return ".\\";
#else
            return "./";
#endif
        }
        if(pos == fpath.size() - 1){
            std::string::size_type pos1 = fpath.find_last_not_of("/\\");
            if(pos1 == std::string::npos){
#ifdef _WIN32
                return "\\";
#else
                return "/";
#endif
            }
            std::string::size_type pos2 = fpath.find_last_of("/\\", pos1);
            if(pos2 == std::string::npos){
#ifdef _WIN32
                return ".\\";
#else
                return "./";
#endif
            }else{
                return fpath.substr(0, pos2 + 1);
            }
        }else{
            return fpath.substr(0, pos + 1);
        }
    }

    /** get absolute path of a path
     * @param path string of path
     * @return absolute path string
     */
    inline std::string abspath(const std::string& path){
        std::string newPath = util::strip(path);
        if(newPath.length() == 0){
            return "";
        }
        char cpath[FILENAME_MAX];
#ifdef _WIN32
        _realpath(newPath.c_str(), cpath);
#else
        realpath(newPath.c_str(), cpath);
#endif
        return std::string(cpath);
    }

    /** get current working directory
     * @return current working directory
     */
    inline std::string cwd(void){
        char cpath[FILENAME_MAX];
#ifdef _WIN32
        _getcwd(cpath, FILENAME_MAX);
#else
        getcwd(cpath, FILENAME_MAX);
#endif
        return std::string(cpath);
    }

    /** join dirname and basename into a path
     * @param dirname dirname of the path
     * @param basename basename of the path
     * @return full path
     */
    inline std::string joinpath(const std::string& dirname, const std::string& basename){
#ifdef _WIN32
        return dirname + "\\" + basename;
#else
        return dirname + "/" + basename;
#endif
    }

    /** check a string is a regular file or not
     * @param path string of a file/directory
     * @return true if path is an existing regular file
     */
    inline bool isfile(const std::string& path){
#ifdef _WIN32
        struct _stat info;
        if(_stat(path.c_str(), &info) != 0){
            return false;
        }
        return (info.st_mode & _S_IFREG);
#else
        struct stat info;
        if(stat(path.c_str(), &info) !=0){
            return false;
        }
        return (info.st_mode & S_IFREG);
#endif
    }

    /** check a string is a directory or not
     * @param path string of a file/directory
     * @return true if path is an existing path
     */
    inline bool isdir(const std::string& path){
#ifdef _WIN32
        struct _stat info;
        if(_stat(path.c_str(), &info) != 0){
            return false;
        }
        return (info.st_mode & _S_IFDIR);
#else
        struct stat info;
        if(stat(path.c_str(), &info) != 0){
            return false;
        }
        return (info.st_mode & S_IFDIR);
#endif
    }
    
    /** check a file/directory exists or not
     * @param path string of a file/directory
     * @return true if path exists
     */
    inline bool exists(const std::string& path){
        return util::isdir(path) || util::isfile(path);
    }

    /** exit and print string to std::cerr
     * @param msg string to print to std::cerr
     */
    inline void errorExit(const std::string& msg){
        std::cerr << "ERROR: " << msg << std::endl;
        exit(-1);
    }
    
    /** check a file existence status, if not exists, exit 
     * @param path string of file/directory
     */
    inline void validFile(const std::string& path){
        if(util::isdir(path)){
            util::errorExit("this is not a file path!");
        }
        if(!util::isfile(path)){
            util::errorExit("file does not exist");
        }
    }

    /** make directories recursively
     * @param path path of directory to be created
     * @return true if make directories successfully
     */
    inline bool makedir(const std::string& path){
        std::string fpath = util::replace(path, "~", std::string(std::getenv("HOME")) + "/");
#ifdef _WIN32
        int ret = _mkdir(fpath.c_str());
#else
        mode_t mode = 0755;
        int ret = mkdir(fpath.c_str(), mode);
#endif
        if(ret == 0){
            return true;
        }
        switch(errno){
            case ENOENT:
                {
                    if(!util::makedir(util::dirname(fpath))){
                        return false;
                    }
                }
#ifdef _WIN32
                return 0 == _mkdir(fpath.c_str());
#else
                return 0 == mkdir(fpath.c_str(), mode);
#endif
            case EEXIST:
                return util::isdir(fpath);
            default:
                return false;
        }
    }

    /** remove non-alpha characters from a string
     * @param str string to be filtered
     * @return a string without non-alpha characters
     */
    inline std::string getAlpha(const std::string& str){
        std::string ret;
        std::copy_if(str.cbegin(), str.cend(), std::back_inserter(ret), (int(*)(int))std::isalpha);
        return ret;
    }

    /** remove invalid sequence characters from a string
     * @param str string to be filtered
     * @param upper convert result to uppercase if true
     */
    inline void getValid(std::string& str, bool upper = false){
        uint16_t total = 0;
        for(uint16_t i = 0; i < str.size(); ++i){
            if(std::isalpha(str[i]) || str[i] == '-' || str[i] == '*'){
                str[total++] = (upper ? std::toupper(str[i]) : str[i]);
            }
        }
        str.resize(total);
    }

    /** make a string each character uppercased
     * @param str string to be uppercased
     */
    inline void str2upper(std::string& str){
        std::transform(str.begin(), str.end(), str.begin(), (int (*)(int))std::toupper);
    }

    /** make a string each character lowercased
     * @param str string to be lowercased
     */
    inline void str2lower(std::string& str){
        std::transform(str.begin(), str.end(), str.begin(), (int (*)(int))std::tolower);
    }

    /** get hamming distance of two strings
     * @param str1 string 1
     * @param str2 string 2
     * @return hamming distance of string 1 and string 2
     */
    inline int hamming(const std::string& str1, const std::string& str2){
        int diff = std::abs((int)(str1.size() - str2.size()));
        for(uint16_t i = 0; i < std::min(str1.size(), str2.size()); ++i){
            diff += (str1[i] == str2[i] ? 0 : 1);
        }
        return diff;
    }

    /** convert number to 33 based score character
     * @param num number of score
     * @return 33 based score character
     */
    inline char num2qual(int num){
        if(num > 127 - 33){
            num = 127 - 33;
        }
        if(num < 0){
            num = 0;
        }
        char c = num + 33;
        return c;
    }

    /** convert a quality uint8_t to phred33 score string
     * @param a pointer to uint8_t likely array
     * @param l length of array
     * @return string representation of phred33 based score
     */
    template<typename T>
    inline std::string qual2str(T *a, uint16_t l){
        std::string qstr(l, '\0');
        for(uint16_t i = 0; i < l; ++i){
            qstr[i] = (char)(33 + a[i]);
        }
        return qstr;
    }

    /** get complement base of a nucleotide base
     * @param base nucleotide base character
     * @return the uppercased complementary nucleotide base
     */
    inline char complement(char base){
        switch(base){
            case 'A': case 'a':
                return 'T';
            case 'T': case 't':
                return 'A';
            case 'C': case 'c':
                return 'G';
            case 'G': case 'g':
                return 'C';
            default:
                return 'N';
        }
    }

    /** get reverse completement sequence of a nucleotide sequence
     * @param seq a nucleotide sequence
     * @return the reverse completement sequence of seq
     */
    inline std::string reverseComplement(const std::string& seq){
        std::string retSeq(seq.length(), '\0');
        for(int32_t i = retSeq.length() - 1; i >= 0; --i){
            retSeq[i] = complement(seq[seq.length() - 1 - i]);
        }
        return retSeq;
    }
    
    /** reverse completement an sequence of a nucleotide sequence
     * @param seq a nucleotide sequence to be reverse complemented
     */
    inline void reverseComplement(std::string& seq){
        int32_t i = seq.size() - 1;
        int32_t j = 0;
        char tmpChr = '\0';
        while(j <= i){
            tmpChr = complement(seq[i]);
            seq[i] = complement(seq[j]);
            seq[j] = tmpChr;
            ++j;
            --i;
        }
    }
    
    /** get forward completment sequene of a nucleotide sequence
     * @param seq a nucleotide sequence
     * @return the forward completement sequence of seq
     */
    inline std::string forwardComplement(const std::string& seq){
        std::string retSeq(seq.length(), '\0');
        for(uint32_t i = 0; i < retSeq.length(); ++i){
            retSeq[i] = complement(seq[i]);
        }
        return retSeq;
    }

    /** write a log message to std::cerr in a thread-safe way
     * @param s log message 
     * @param logmtx reference to a std::mutex object
     */
    inline void loginfo(const std::string& s, std::mutex& logmtx){
        std::lock_guard<std::mutex> l(logmtx);
        time_t tt = time(NULL);
        tm* t = std::localtime(&tt);
        char date[60] = {0};
        std::sprintf(date, "[%d-%02d-%02d %02d:%02d:%02d] ",
                t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
                t->tm_hour, t->tm_min, t->tm_sec);
        std::cerr << date << s << std::endl;
    }
    /** write a log message to std::cerr directly
     * @param s log message 
     */
    inline void loginfo(const std::string& s){
        time_t tt = time(NULL);
        tm* t = std::localtime(&tt);
        char date[60] = {0};
        std::sprintf(date, "[%d-%02d-%02d %02d:%02d:%02d] ",
                t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
                t->tm_hour, t->tm_min, t->tm_sec);
        std::cerr << date << s << std::endl;
    }

    /** get current time
     * @return year-mm-dd hh-mm-ss of current time
     */
    inline std::string currentTime(void){
        time_t tt = time(NULL);
        tm* t = std::localtime(&tt);
        char date[60] = {0};
        std::sprintf(date, "%d-%02d-%02d %02d:%02d:%02d",
                t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
                t->tm_hour, t->tm_min, t->tm_sec);
        return date;
    }


    /** make a list from file by line 
     * @param filename input file
     * @param ret vector to store stripped line
     */ 
    inline void makeListFromFileByLine(const std::string& filename, std::vector<std::string>& ret){
        std::ifstream fr(filename);
        std::string line;
        while(std::getline(fr, line)){
            util::strip(line);
            ret.push_back(line);
        }
    }
    /** test whether an element is in a vector
     * @param v vector
     * @param e element
     * @return true if e in v
     */
    template<typename T>
    inline bool inVector(const std::vector<T>& v, const T& e){
        return std::find(v.cbegin(), v.cend(), e) != v.cend();
    }

    /** convert a vector of strings to integers
     * @param vs vector of strings
     * @param vi vector of integers
     */
    template<typename T>
    inline void strvec2intvec(const std::vector<std::string>& vs, std::vector<T>& vi){
        vi.clear();
        std::transform(vs.begin(), vs.end(), std::back_inserter(vi), std::bind(std::atoi, std::bind([&](std::string str){return str.c_str();}, std::placeholders::_1)));
    }

    /** convert two vector to a vector of pairs
     * @param v1 vector 1
     * @param v2 vector 2
     * @param vp pair vector
     */
    template<typename T>
    inline void vec2pairvec(const std::vector<T>& v1, const std::vector<T>& v2, std::vector<std::pair<T, T>>& vp){
        vp.clear();
        int maxlen = std::min(v1.size(), v2.size());
        for(int i = 0; i < maxlen; ++i){
            vp.push_back(std::make_pair(v1[i], v2[i]));
        }
    }

    /** list contents in an directory
     * @param path string of path
     * @param vname vector to store contents under path
     */
    inline void listDir(const std::string& path, std::vector<std::string>& vname){
        std::string dirpath = util::replace(path, "~", std::string(std::getenv("HOME")) + "/");
        DIR *dir;
        struct dirent *ent;
        if((dir = opendir(dirpath.c_str())) != NULL){
            while((ent = readdir(dir)) != NULL){
                vname.push_back(std::string(ent->d_name));
            }
            closedir(dir);
        }
    }

    /** print an array to std::cout 
     * @param a pointer to an array
     * @param l length of the array
     * @param n array name to show
     */
    template<typename T>
    inline void showArray(const T* a, uint16_t l, const std::string& n){
        std::cout << n << ":";
        for(uint16_t i = 0; i < l; ++i){
            std::cout << a[i] << " ";
        }
        std::cout << std::endl;
    }

    /** count different neighbor pairs in a string
     * @param str a string
     * @return different neighbor pairs in this string
     */
    inline int neighborDiffCount(const std::string& str){
        int diff = 0;
        for(uint32_t i = 0; i < str.length() - 1; ++i){
            if(str[i] != str[i+1]){
                ++diff;
            }
        }
        return diff;
    }
    
    /** get mismatch ratio of two uncleotide sequence
     * @param s1 nucleotide sequence
     * @param s2 nucleotide sequence
     * @return mismatch ratio to the shorter string
     */
    inline float mismatchRatio(const std::string& s1, const std::string& s2){
        uint32_t minLen = std::min(s1.length(), s2.length());
        if(minLen == 0 || s1.find_first_of("ATCG") == std::string::npos || s2.find_first_of("ATCG") == std::string::npos){
            return 1.0;
        }
        if(minLen == 1){
            return s1[0] == s2[0];
        }
        size_t beg = 0;
        for(beg = 0; beg < minLen; ++beg){
            if(s1[beg] == 'N' || s2[beg] == 'N'){
                ++beg;
            }
        }
        size_t end = minLen - 1;
        for(end = minLen - 1; end > beg; --end){
            if(s1[end] == 'N' || s2[end] == 'N'){
                --end;
            }
        }
        int32_t mismatch = 0;
        for(uint32_t i = beg; i <= end; ++i){
            if(s1[i] == 'N' || s2[i] == 'N'){
                continue;
            }
            if(s1[i] != s2[i]){
                ++mismatch;
            }
        }
        return float(mismatch)/(end - beg + 1);
    }

    /** get median of a list of values
     * @param v a vector of values
     * @return median of values in v
     */
    template<typename T>
    inline double median(std::vector<T>& v){
        std::nth_element(v.begin(), v.begin() + v.size()/2, v.end()); 
        return v[v.size()/2];
    }

    /** test whether a string contains only AaCcGgTt (ie. is DNA)
     * @param seq string of characters
     * @return true if seq contains only AaCcGgTt
     */
    inline bool isDNA(const std::string& seq){
        return seq.find_first_not_of("AaTtCcGg") == std::string::npos;
    }

    /** count number of characters in a string
     * @param str full string
     * @param chrs characters to be found in str
     */
    inline int32_t countChrs(const std::string& str, const std::string& chrs){
        int32_t count = 0;
        std::string::size_type pos = 0;
        while((pos = str.find_first_of(chrs, pos)) != std::string::npos){
            ++pos;
            ++count;
        }
        return count;
    }
}

#endif
