//
// Created by liufei on 2015/12/8.
//

#ifndef ERRORDETECTION_ERRORDETECTION_H
#define ERRORDETECTION_ERRORDETECTION_H

#include "Image.h"
#include "GCP.h"
#include <string>
#include <set>

class ErrorDetection{
public:
    /* _____ interface _____ */
    int process(std::string ioe_fileName, std::string gcp_fileName,
                std::string pixel_fileName, std::string output_fileName);

private:
    /* _____ read data _____ */
    int ReadTestFile(std::string ioe_fileName,
                     std::string gcp_fileName, std::string pixel_fileName);

    /* _____ get common pixel _____ */
    int GetCommonPixel();

    /* _____ relative orientation and error detection, main processing step _____ */
    int DetectError(std::string output_fileName);

private:
    const double m0 = 2.8e-3;   //unit: mm
    const double photographyScale = 12000;
    
    Image img1,img2;
    std::vector<int> cPixels;   //common pixels
    std::set<int> ePoints;      //error points
    std::vector<GCP> GCPs;
};


#endif //ERRORDETECTION_ERRORDETECTION_H
