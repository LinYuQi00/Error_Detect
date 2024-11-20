//
// Created by liufei on 2015/12/8.
//

#ifndef ERRORDETECTION_IMAGEPIXEL_H
#define ERRORDETECTION_IMAGEPIXEL_H


class ImagePixel{
public:
    ImagePixel();
    ImagePixel(int id, double _x, double _y);

public:
    int ID;
    double x, y;
};


#endif //ERRORDETECTION_IMAGEPIXEL_H
