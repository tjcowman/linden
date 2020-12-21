/**
 * @author Tyler Cowman
 * 
 * Convenience class for parsing rectangular, delimited text files. Returns instances of the matrix class
 * which is implemented as a vector of vectors.
 */

#ifndef TSVPARSER_H
#define TSVPARSER_H

#include "Matrix.h"

#include <iostream>
#include <vector>
#include <string>


struct range
{
    
    bool inRange(unsigned int i)
    {
        if(V.size() == 0)
            return true;
        else if(index>=V.size())
            return false;
        
        bool rVal = (i >= V[index] && i<= V[index+1]);
        
        if(i>=V[index+1])
            index+=2;
        
        return rVal;
    }
    
    void reset()
    {
        index=0;
    }
    
    std::vector<int> V;
    unsigned int index;
};

class TsvParser
{
    public:
        TsvParser();
        TsvParser(char delimiter);
        
        void setColumns(std::vector<int> columns);
        void setRows(std::vector<int> rows);
        void setDelimiter(char delimiter);
        
        /*TODO doesnt end when ranges provided and out of bounds, doesnt crash but loops unecessarily*/
        Matrix<std::string> parse(std::string filename);
        Matrix<std::string> parse(std::istream & stream);
        
        Matrix<char> parseByteWise(std::string filename);
        Matrix<char> parseByteWiseDelimited(std::string filename, int columnBegin);

        void formattedOutput(Matrix<std::string> M, std::ostream & out, char delimiter, char quoteL, char quoteR );
        void formattedOutput(std::vector<std::string> V, std::ostream & out, char delimiter, char quoteL, char quoteR );
        
    
    private:
        
       
        
        char delimiter_;
        //std::vectors of ranges (begin end) of columns and rows, if empty parse all
        std::vector<int> columns_;
        std::vector<int> rows_;
        Matrix<std::string> parsedData_;
        
};


#endif //TSVPARSER_H
