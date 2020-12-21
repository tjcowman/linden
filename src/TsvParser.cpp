#include "TsvParser.h"

TsvParser::TsvParser()
{
    delimiter_ = '\t';
    columns_ = std::vector<int>();
    rows_ = std::vector<int>();
}

TsvParser::TsvParser(char delimiter)
{
    delimiter_ = delimiter;
    columns_ = std::vector<int>();
    rows_ = std::vector<int>();
}

void TsvParser::setColumns(std::vector<int> columns)
{
    columns_ = columns;
}

void TsvParser::setRows(std::vector<int> rows)
{
    rows_ = rows;
}

void TsvParser::setDelimiter(char delimiter)
{
    delimiter_ = delimiter;
}

Matrix<std::string> TsvParser::parse(std::string filename)
{
    std::ifstream FILE(filename.c_str());
    if(!FILE.is_open())
    {
        std::cerr<<"file "<<filename<<" not opened"<<"\n";
        exit(1);
    }
    
    Matrix<std::string> M;
    
    std::string L;
    unsigned int lineCounter = 0;
    
    struct range rowRange{rows_, 0};
    struct range columnRange{columns_, 0};

    while( true)
    {
        getline(FILE, L);
        if(FILE.eof()) break;
        
        std::vector<std::string> parsedRow;
 
        if(rowRange.inRange(lineCounter))
        {
            unsigned int wordCounter = 0;
            unsigned int pos1 = 0;
            unsigned int pos2 = 0;
           
            while((pos2)<=L.size() /*pos2 != std::string::npos*/)
            {
                pos2 = L.find(delimiter_, pos1);
                
                if(columnRange.inRange(wordCounter))
                    parsedRow.push_back(L.substr(pos1, pos2-pos1));
                
                pos1 = pos2+1;
                wordCounter++;
                 
            }

            M.addRow(parsedRow);
            columnRange.reset();
            
        }
        ++lineCounter;
    }
    FILE.close();
    return M;
}

Matrix<std::string> TsvParser::parse(std::istream & stream)
{
    Matrix<std::string> M;
    
    std::string L;
    unsigned int lineCounter = 0;
    
    struct range rowRange{rows_, 0};
    struct range columnRange{columns_, 0};

    while( true)
    {
        getline(stream, L);
        if(stream.eof()) break;
        
        std::vector<std::string> parsedRow;
 
        if(rowRange.inRange(lineCounter))
        {
            unsigned int wordCounter = 0;
            unsigned int pos1 = 0;
            unsigned int pos2 = 0;
           
            while((pos2)<=L.size() /*pos2 != std::string::npos*/)
            {
                pos2 = L.find(delimiter_, pos1);
                
                if(columnRange.inRange(wordCounter))
                    parsedRow.push_back(L.substr(pos1, pos2-pos1));
                
                pos1 = pos2+1;
                wordCounter++;
                 
            }

            M.addRow(parsedRow);
            columnRange.reset();
            
        }
        ++lineCounter;
    }

    return M;
}

Matrix<char> TsvParser::parseByteWise(std::string filename)
{
    std::fstream FILE(filename.c_str());
    if(!FILE.is_open())
    {
        std::cerr<<"file "<<filename<<" not opened"<<"\n";
        exit(1);
    }
    
    Matrix<char> M;
    
    std::string L;
    unsigned int lineCounter = 0;
    
    struct range rowRange{rows_, 0};
    struct range columnRange{columns_, 0};
    
    
    while( true)
    {
        getline(FILE, L);
        if(FILE.eof()) break;
        
        std::vector<char> parsedRow;
        
        if(rowRange.inRange(lineCounter))
        {
            for(int i=0; i<L.size(); ++i)
            {
                if(columnRange.inRange(i))
                    parsedRow.push_back(L[i]);
            }
            
            M.addRow(parsedRow);
            columnRange.reset();
        }
        ++lineCounter;
        
    }
    FILE.close();
    return M;
}

Matrix<char> TsvParser::parseByteWiseDelimited(std::string filename, int columnBegin)
{
    std::fstream FILE(filename.c_str());
    if(!FILE.is_open())
    {
        std::cerr<<"file "<<filename<<" not opened"<<"\n";
        exit(1);
    }
    
    Matrix<char> M;
    
    std::string L;
    unsigned int lineCounter = 0;
    
    struct range rowRange{rows_, 0};
    struct range columnRange{columns_, 0};
    
    
    while( true)
    {
        getline(FILE, L);
            if(FILE.eof()) break;
        
        int dataBegins=0;
        for(int i=0; i<columnBegin; ++i)
        {
            dataBegins = L.find(delimiter_, dataBegins+1);
        }
        
        std::vector<char> parsedRow;
        
        if(rowRange.inRange(lineCounter))
        {
            for(int i=dataBegins+1; i<L.size(); i+=2)
            {
                if(columnRange.inRange((i-dataBegins)/2))
                    parsedRow.push_back(L[i]);
            }

            M.addRow(parsedRow);
            columnRange.reset();
        }
        ++lineCounter;
        
    }
    FILE.close();
    return M;
}

void TsvParser::formattedOutput(Matrix<std::string> M, std::ostream & out, char delimiter, char quoteL, char quoteR )
{
    for(int i=0; i<M.dim(0); ++i)
    {
        for(int j=0; j<M.dim(1)-1; ++j)
        {
            out<<quoteL<<M.a(i,j)<<quoteR<<delimiter;
        }
        out<<quoteL<<M.a(i, M.dim(1)-1)<<quoteR;
        out<<std::endl;
    }        
}

void TsvParser::formattedOutput(std::vector<std::string> V, std::ostream & out, char delimiter, char quoteL, char quoteR )
{
    for(int i=0; i<V.size()-1; ++i)
    {
        out<<quoteL<<V[i]<<quoteR<<delimiter;
    }
    out<<quoteL<<V[V.size()-1]<<quoteR;
    out<<std::endl;
           
}
