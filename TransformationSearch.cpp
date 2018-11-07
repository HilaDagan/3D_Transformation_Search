
#include "TransformationSearch.h"

#define K 5 // the number of the results we want to be stored.
#define THREADS_NUM 8

using namespace std;

/**
 * Constructor.
 * @param mapLo - a density map in low resolution (from EM)
 * @param molHi - a molecule object in high resolution (obtained from pdb file)
 */
TransformationSearch::TransformationSearch(DensityMap& mapLo, ChemMolecule& molHi,
                                           double resolution, int nx, int ny, int nz) :
        mapLo(mapLo), molHi(molHi), resolution(resolution), nx(nx), ny(ny), nz(nz) {};

/**
 * Destructor.
 */
TransformationSearch::~TransformationSearch() = default;

/**
 * Parse the angles from the given file.
 * @return a vector of angles ordered in groups of three.
 */
void TransformationSearch::getAnglesFromFile(const std::string &filename) {
    std::ifstream afile(filename.c_str());
    if (!afile.is_open()) {
        std::cerr << "Error in open angles file" << std::endl;
        exit(1);
    }
    std::string line;
    while (!afile.eof()) {
        getline(afile, line);
        if (!line.empty()) {
            std::vector<std::string> angleList;
            trim_if(line, boost::is_any_of("\t "));
            boost::split(angleList, line, boost::is_any_of(" \t"), boost::token_compress_on);
            if (angleList.size() != 4) {
                std::cerr << "Error in file format" << std::endl;
                exit(1);
            } else {
                //Instantiate a rotation using 4 quaternion parameters.
                Rotation3 newAngle(boost::lexical_cast<Rotation3::real>(angleList[0]),
                                 boost::lexical_cast<Rotation3::real>(angleList[1]),
                                 boost::lexical_cast<Rotation3::real>(angleList[2]),
                                 boost::lexical_cast<Rotation3::real>(angleList[3]));
                anglesVec.push_back(newAngle);
            }
        }
    }
}

/**
 * Iterate over all the possible translation using FFT, for the given rotated map high.
 * @param mapHi - a density map in a high resolution, after rotated by 'angles'.
 * @param bestK - contains the best K transformations (translation + rotation)
 * @param angles - a vector contains the angles of the rotation (in RAD). Used for the bestK information.
 */
int TransformationSearch::fftwTranslationalSearch(const DensityMap &mapHi, BestK &bestK, Rotation3 &angles)
{
    // STEP 1
    // convert each map to complex numbers
    unsigned long size =  (long) mapLo.getNumberOfVoxels();
    std::vector<std::complex<double>> mapLoComp(size), mapHiComp(size); // complex
    FFTCommon::cast_values(mapLo.getData(), mapLoComp);
    FFTCommon::cast_values(mapHi.getData(), mapHiComp);

    // apply FFT on each map
    std::vector<std::complex<double> > mapLoOut(size), mapHiOut(size);
    mutexFFT.lock();
    FFTCommon::fftw_forward(mapLoComp, mapLoOut, nx, ny, nz);
    FFTCommon::fftw_forward(mapHiComp, mapHiOut, nx, ny, nz);
    mutexFFT.unlock();

    //STEP 2: computation of the cross-correlation
    std::vector<std::complex<double>> cross_corr_arr(size);
    unsigned long ntrans = (unsigned long) nx * ny * nz; // the number of translations
    for(unsigned long i=0; i< ntrans; i++) {
        if(std::norm(mapLoOut[i]) == 0.0 || std::norm(mapHiOut[i] == 0.0))
            cross_corr_arr[i] = 0;
        else
            cross_corr_arr[i] = (mapLoOut[i] * std::conj(mapHiOut[i])); // (abs(F1[i][j]) * abs(F2[i][j]));
    }

    //STEP 3: transform backward
    std::vector<std::complex<double>> inv_r(size);
    mutexFFT.lock();
    FFTCommon::fftw_backward(cross_corr_arr, inv_r, nx, ny, nz);
    mutexFFT.unlock();

    // convert to doubles (taking norm)
    std::vector<double> ir(size);
    FFTCommon::complex_to_double(inv_r, ir);

    //STEP 4: find K maximum results
    double voxSize = mapLo.getVoxelSize();
    for (long i=0; i<ir.size(); i++) {
        long x = i % nx;
        long zny_plus_y = i / nx;
        long y = zny_plus_y % ny;
        long z = zny_plus_y / ny;

        //make sure that x,y and z are in the right range:
        if (x > nx/2.0) x = -(nx - x);
        if (y > ny/2.0) y = -(ny - y);
        if (z > nz/2.0) z = -(nz - z);

        auto * trans = new transformation(ir[i], x * voxSize,
                                          y * voxSize, z * voxSize, angles);
        bestK.push(trans);
    }
    return 0;
}

/**
 * Find the K best transformations, i.e rotation and translation.
 * This function is executed by each thread.
 * @param index - the thread index. Used for determining the angles to work on for each thread.
 * @param bestKVec - a vector that contains the BestK transformations of each thread.
 * @param threadsNum - the total number of threads.
 */
void TransformationSearch::findBestTransformation(int index, vector<BestK> & bestKVec, int threadsNum){
    BestK bestK = BestK(K, true);
    Vector3 translation = Vector3(0.0, 0.0, 0.0);

    // the number of angles (rows in the file) each thread should operate on.
    int anglesNum = anglesVec.size() / threadsNum;
    // For each rotation, get the best transformation using FFT
    for (int j = 0; j < anglesNum; ++j) {
        // apply the current rotation on the molecule
        ChemMolecule rotatedMolHi = molHi; // do not change the original molecule

        // each thread works on different rotations according to its index.
        rotatedMolHi.rigidTrans(RigidTrans3(anglesVec[j + index * anglesNum], translation));

        // generate a map from the rotated molecule.
        DensityMap rotatedMapHi(mapLo);
        rotatedMapHi.simulateMap(rotatedMolHi, resolution); // create a map density according to the given mol

        // find the best transformation
        fftwTranslationalSearch(rotatedMapHi, bestK, anglesVec[j + index * anglesNum]);
    }
    mutexBestKVec.lock();
    bestKVec.push_back(bestK);
    mutexBestKVec.unlock();
}

/**
 * Manage the program using multi-threading.
 * @return the best K transformations
 */
BestK TransformationSearch::startThreaded(){
    std::vector<std::thread> threads;
    std::vector<BestK> bestKVec;  // a vector that contains the BestK transformations of each thread.
    // create the threads
    for(int i = 0; i < THREADS_NUM; i++) {
        std::thread t(&TransformationSearch::findBestTransformation, this, i, std::ref(bestKVec), THREADS_NUM);
        threads.push_back(std::move(t));
    }

    for(auto&& thread: threads)
        if(thread.joinable())
            thread.join();

//    Code for debugging
//    std::cerr << "===================" << std::endl;
//    for (int j = 0; j < bestKVec.size(); ++j) {
//        std::cerr << "bestKVec size: " <<  bestKVec[j].size() << std::endl;
//        for (int k = 0; k < bestKVec[j].size(); ++k) {
//            std::cerr << "score: " << bestKVec[j][k]->score << std::endl;
//            std::cerr << "x " << bestKVec[j][k]->x << " y " << bestKVec[j][k]->y << " z "
//                      << bestKVec[j][k]->z << " ANGLES: " << bestKVec[j][k]->angles << std::endl;
//        }
//    }
//    std::cerr << "===================" << std::endl;

    // Unify the bestK objects of all the threads into one object.
    BestK bk(5, false);
    for (int j = 0; j < bestKVec.size(); ++j) {
        for (int k = 0; k < bestKVec[j].size(); ++k) {
            bk.push(bestKVec[j][k]);
        }
    }
    // Needed so the destructor won't try to delete BestK elements,
    // which are also in the bestKVec
    for (int j = 0; j < bestKVec.size(); ++j) {
        bestKVec[j].clear();
    }
    return bk;
}