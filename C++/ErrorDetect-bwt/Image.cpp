//
// Created by liufei on 15-12-10.
//

#include "Image.h"

Image::Image(): ID{-1}, pixelNum{0} {

}

Image::Image(int _ID, Inner &_in, Exterior &_ext, std::unordered_map<int,ImagePixel> &_pixels): ID{_ID} {
    in = _in;
    ext = _ext;
    pixels = move(_pixels);
    pixelNum = static_cast<int>(pixels.size());
}

int Image::insert(const ImagePixel &pixel){
    pixels.insert({pixel.ID, pixel});
    ++pixelNum;
    return 0;
}

int Image::find(const int id){
    //std::unordered_map<int, ImagePixel>::iterator
    auto got = pixels.find(id);

    if(got == pixels.end())
        return 0;
    else
        return got->first;

    return 0;
}
