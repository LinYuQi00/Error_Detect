//
// Created by liufei on 15-12-10.
//

#ifndef ERRORDETECTION_IMAGE_H
#define ERRORDETECTION_IMAGE_H

#include <unordered_map>
#include "ImagePixel.h"

struct Inner{
    double f;
    double x0, y0;
    Inner(): f{0}, x0{0}, y0{0}{}
};

struct Exterior{
    double Xs, Ys, Zs;
    double phi, omega, kappa;
    Exterior(): Xs{0}, Ys{0}, Zs{0}, phi{0}, omega{0}, kappa{0}{}
    Exterior(double _Xs, double _Ys, double _Zs, double _phi, double _omega, double _kappa):
            Xs{_Xs}, Ys{_Ys}, Zs{_Zs}, phi{_phi}, omega{_omega}, kappa{_kappa}{}
};

class Image {
public:
    Image();
    Image(int _ID, Inner &_in, Exterior &_ext, std::unordered_map<int,ImagePixel> &_pixels);

    Image(const Image &img) = delete;
    Image& operator=(const Image &img) = delete;

    /* _____ change values _____ */
    void setID(const int _ID) {ID = _ID;}
    void setIOE(const Inner &_in) { in = _in;}
    void setEOE(const Exterior &_ext) {ext = _ext;}

    /* _____ get values _____ */
    int getPixelNum() const {return pixelNum;}
    int getID() const {return ID;}
    Inner& getInner() {return in;}
    Exterior& getExt() {return ext;}

    ImagePixel& operator[](const int id) { return pixels[id];}
    ImagePixel& getPixel(const int id) {return pixels[id];}
    std::unordered_map<int,ImagePixel>::iterator begin() {return pixels.begin();}
    std::unordered_map<int,ImagePixel>::iterator end() {return pixels.end();};

    /* _____ insert pixels _____ */
    int insert(const ImagePixel &pixel);

    /* _____ found point in pixels _____ */
    int find(const int id);

private:
    int ID;
    Inner in;
    Exterior ext;
    int pixelNum;
    std::unordered_map<int,ImagePixel> pixels;
};


#endif //ERRORDETECTION_IMAGE_H
