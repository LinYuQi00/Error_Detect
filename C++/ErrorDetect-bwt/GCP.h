//
// Created by liufei on 2015/12/8.
//

#ifndef ERRORDETECTION_GCP_H
#define ERRORDETECTION_GCP_H

#include <vector>

class GCP{
public:
    GCP();
    GCP(int id, double x, double y, double z, int d = 3);

    int getID() const {return ID;}
    double getX() const {return X;}
    double getY() const {return Y;}
    double getZ() const {return Z;}
    int getDemension() const {return demension;}
private:
    int ID;
    double X, Y, Z;
    int demension;

    struct GCPpixel{
        int imgID;
        double x, y;
        GCPpixel(): imgID{-1}, x{0}, y{0}{}
    };
    std::vector<GCPpixel> pixels;
};


#endif //ERRORDETECTION_GCP_H
