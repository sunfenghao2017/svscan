#ifndef FILEREADER_HPP
#define FILEREADER_HPP

#include <zlib.h>
#include <cstdio>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <iostream>

/** Class to hold a file reader in the format .gz or plain text */
class FileReader{
    std::string mFileName;  ///< name of file
    gzFile mGzipFile;       ///< gzFile to store opened mZipped file handler
    FILE* mFile;            ///< FILE pointer to store opened plain file handler
    bool mZipped;           ///< the file is gzipped if true
    char* mBuf;             ///< mBuffer to store a chunk of characters read from mGzipFile or mFile
    int mBufDataLen;        ///< the length of characters read into the mBuffer after the last read
    int mBufUsedLen;        ///< the length of characters already consumed in the mBuffer 
    bool mStdinMode;        ///< read from stdin if true
    bool mNoLineBreakAtEnd; ///< the file has no '\n' as a line break at the last line if true
    int32_t mReadBufSize;   ///< the mBuffer size used to read
    
    public:
    /** Construct a file reader with filename
     * @param filename Name of the file
     */
    FileReader(const std::string& filename){
        mFileName = filename;
        mGzipFile = NULL;
        mFile = NULL;
        mStdinMode = false;
        mReadBufSize = (1 << 20);
        mBuf = new char[mReadBufSize];
        mBufDataLen = 0;
        mBufUsedLen = 0;
        mNoLineBreakAtEnd = false;
        init();
    }
    
    /** FileReader Destructor */
    ~FileReader(){
        close();
        delete mBuf;
        mBuf = nullptr;
    }
   
    /** Tell whether the file is zipped or not
     * @return true if the fasq file is zipped 
     */
    inline bool isZipped(){
        return mZipped;
    }

    /** Get the number of Bytes read from(write to) the file in the meantime, store it in bytesRead
     *  Get the total number of Bytes in the file, store it in bytesTotal
     */
    inline void getBytes(size_t& bytesRead, size_t& bytesTotal){
        if(mZipped){
            bytesRead = gzoffset(mGzipFile);
        }else{
            bytesRead = std::ftell(mFile);
        }
        // use another ifstream without affecting the current reader
        std::ifstream is(mFileName);
        is.seekg(0, is.end);
        bytesTotal = is.tellg();
    }

    /** tell whether the file is a zipped file
     * @param filename file name
     * @return true if the file is a zipped file(with suffix ".gz")
     */
    inline static bool isZippedFile(const std::string& filename){
        if(filename.length() < 4){
            return false;
        }
        if(filename.substr(filename.length() - 3) == ".gz"){
            return true;
        }
        return false;
    }
    
    /** read one line into line from buffer
     * @param line strin to store result
     * @return true if read successful
     */
    inline bool getline(std::string& line){
        if(mBufUsedLen >= mBufDataLen && eof()){
            return false;
        }
        if(mZipped && mGzipFile == NULL){
            return false;
        }
        line = getlineFromBuffer();
        return true;
    }

private:
    /** get just one line from the mBuf, update the mBuf if needed
     */
    inline std::string getlineFromBuffer(){
        int start = mBufUsedLen;
        int end = start;
        // look for '\r' or '\n' until the end of mBuf
        while(end < mBufDataLen){
            if(mBuf[end] != '\r' && mBuf[end] != '\n'){
                ++end;
            }else{
                break;
            }
        }
        // if '\r' or '\n' found(this line well contained in this mBuf)
        // or this is the last mBuf of file
        if(end < mBufDataLen || mBufDataLen < mReadBufSize){
            int len = end - start;
            std::string line(mBuf + start, len);
            ++end;
            if(end < mBufDataLen - 1 && mBuf[end - 1] == '\r' && mBuf[end] == '\n'){
                ++end;
            }
            mBufUsedLen = end;
            return line;
        }
        // if '\r' or '\n' not found && this is not the last mBuf of file
        // then this line is not contained in this mBuf, we should read new mBuf
        std::string str(mBuf + start, mBufDataLen - start);
        while(true){
            readToBuf();
            start = 0;
            end = 0;
            // look for '\r' or '\n' until the end of mBuf
            while(end < mBufDataLen){
                if(mBuf[end] != '\r' && mBuf[end] != '\n'){
                    ++end;
                }else{
                    break;
                }
            }
            // if '\r' or '\n' found(this line well contained in this mBuf)
            // or this is the last mBuf of file
            if(end < mBufDataLen || mBufDataLen < mReadBufSize){
                int len = end - start;
                str.append(mBuf + start, len);
                ++end;
                if(end < mBufDataLen - 1 && mBuf[end - 1] == '\r' && mBuf[end] == '\n'){
                    ++end;
                }
                mBufUsedLen = end;
                return str;
            }
            // if '\r' or '\n' not found && this is not the last mBuf of file
            // then this line is not contained in this mBuf, we should read new mBuf
            str.append(mBuf + start, mBufDataLen);
        }
        return std::string();
    }
    
    /** Tell whether the file has no '\n' as a line break at the last line
     * @return true if the file has no '\n' as a line break at the last line
     */ 
    inline bool hasNoLineBreakAtEnd(){
        return mNoLineBreakAtEnd;
    }

    /** Tell whether the FileReader has reach the endof file
     * @return true if eof reached
     */
    inline bool eof(){
        if(mZipped){
            return gzeof(mGzipFile);
        }else{
            return std::feof(mFile);
        }
    }

    /** initialize the FileReader:
     * 1, open file and store file handler into mGzipFile or mFile or read from stdin
     * 2, set the starting position for the next read on compressed file stream file to the beginning of file 
     * 3, update the file format mZipped
     * 4, call readToBuf() to try to fill the mBuf from first reading
     */ 
    inline void init(){
        if(isZippedFile(mFileName)){
            mGzipFile = gzopen(mFileName.c_str(), "r");
            mZipped = true;
            gzrewind(mGzipFile);
        }else{
            if(mFileName == "/dev/stdin"){
                mFile = stdin;
            }else{
                mFile = std::fopen(mFileName.c_str(), "rb");
            }
            if(mFile == NULL){
                std::cerr << "Failed to open file: " <<  mFileName << std::endl;
                std::exit(1);
            }
            mZipped = false;
        }
        readToBuf();
    }
    
    /** close the FileReader:
     * 1, close file handler
     * 2, set file handler to NULL
     */
    inline void close(){
        if(mZipped && mGzipFile){
            gzclose(mGzipFile);
            mGzipFile = NULL;
        }else if(mFile){
            std::fclose(mFile);
            mFile = NULL;
        }else{
            return;
        }
    }

    /** trim \n, \r or \r\n in the tail of the line */
    inline void clearLineBreaks(char* line){
        int len = std::strlen(line);
        if(len >= 2){
            if(line[len - 1] == '\n' || line[len - 1] == '\r'){
                line[len - 1] = '\0';
                if(line[len - 2] == '\r'){
                    line[len - 2] = '\0';
                }
            }
        }
    }
    
    /** try a reading action to fill the mBuf,
     * update mBufDataLen to Bytes read during this reading action
     * reset mBufUsedLen to zero
     * if read the last line(mBuf is not filled), update mNoLineBreakAtEnd
     */
   inline void readToBuf(){
        if(isZipped()){
            mBufDataLen = gzread(mGzipFile, mBuf, mReadBufSize);
            if(mBufDataLen == -1){
                std::cerr << "Error to read gzip file" << std::endl;
                std::exit(1);
            }
        }else{
            mBufDataLen = std::fread(mBuf, 1, mReadBufSize, mFile);
        }
        mBufUsedLen = 0;
        if(mBufDataLen < mReadBufSize){
            if(mBuf[mBufDataLen - 1] != '\n'){
                mNoLineBreakAtEnd = true;
            }
        }
    }
};

#endif
