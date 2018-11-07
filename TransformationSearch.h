#ifndef P_TRANSFORMATIONSEARCH_H
#define P_TRANSFORMATIONSEARCH_H

#include <mutex>
#include <atomic>
#include <thread>
#include <complex>
#include <DensityMap.h>
#include <fftw3.h>
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>

#include "BestK.h"
#include "DensityMap.h"
#include "FFTCommon.h"


/**
 * A class that performs a transformation search using given known angles
 * and FFT for translation search. Each transformation is consisted of 6
 * degrees of freedom: 3 rotations and 3 translations.
 */
class TransformationSearch {
public:
    TransformationSearch(DensityMap& map_lo, ChemMolecule& mol_hi, double resolution,
                        int nx, int ny, int nz);
    ~TransformationSearch();

    void getAnglesFromFile(const std::string &filename);

    void findBestTransformation(int startRow, std::vector<BestK> & bestKVec, int threadsNum);

    int fftwTranslationalSearch(const DensityMap &mapHi, BestK &bestK, Rotation3 &angles);

    BestK startThreaded();

private:
    DensityMap mapLo;  // a density map in low resolution
    ChemMolecule molHi;  // a molecule object in high resolution
    double resolution;  // the resolution of the density map
    std::atomic<int> nx, ny, nz; // the dimensions of the density maps.
    std::mutex mutexBestKVec;
    std::mutex mutexFFT;
    std::vector<Rotation3> anglesVec; // a vector with rotations.
};


#endif //P_TRANSFORMATIONSEARCH_H
